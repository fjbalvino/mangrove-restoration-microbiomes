"""Describe community components and gene associations without new hypothesis tests."""
from pathlib import Path
import argparse,datetime,json,hashlib,itertools
import numpy as np
import pandas as pd
import rule_logic as logic

ap=argparse.ArgumentParser();ap.add_argument('--analysis',type=Path,required=True);ap.add_argument('--ora',type=Path,required=True);ap.add_argument('--definitions',type=Path,required=True);ap.add_argument('--cf-data',type=Path,required=True);ap.add_argument('--reference-genes',type=Path,required=True);ap.add_argument('--old-reference',type=Path,required=True);ap.add_argument('--out-root',type=Path,required=True);a=ap.parse_args()
O=a.out_root/('22_environment_blocks_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S'));O.mkdir(parents=True);(a.out_root/'LATEST.txt').write_text(str(O.resolve()))
assert logic.self_tests()==44
defs=json.load(open(a.definitions));entries=defs['entries'];expr={k:v['definition'] for k,v in entries.items()}
ora=pd.read_csv(a.ora/'primary_enriched_modules.tsv',sep='\t');mods=list(ora.module.unique());trees={k:logic.parse(expr[k],expr) for k in mods};keys=sorted(set.union(*(logic.kos(t,True) for t in trees.values())))
meta=pd.read_csv(a.cf_data/'metadata_003b.csv');meta['profile_id']=meta.sample_id.str.rsplit('_',n=1).str[0];assert len(meta)==51 and meta.profile_id.nunique()==17;assert meta.groupby('profile_id').depth_cm.apply(lambda x:sorted(x)==[5,20,40]).all()
meta.to_csv(O/'samples.tsv',sep='\t',index=False)
patterns=json.load(open(a.cf_data/'03_unit_annotations/annotation_patterns.json'))['KO'];ref=[l.strip() for l in a.reference_genes.read_text().splitlines() if l.strip()];assert len(ref)==76
inventory=[];sample_rules=[];sample_steps=[];quant=[]
blocks={k:logic.split_top(expr[k],'space') for k in mods}
for variant,stem,statuses in [('complete','KO_complete',{'complete_consensus'}),('both_mates','KO_complete_both',{'complete_consensus'}),('supported','KO_supported',{'complete_consensus','partial_consensus'})]:
 def read(stem):
  d=pd.read_csv(a.cf_data/'matrices'/stem,sep='\t',index_col=0);assert d.index.is_unique and set(d.columns)==set(meta.sample_id);assert np.isfinite(d.to_numpy()).all() and (d.to_numpy()>=0).all();return d.loc[:,meta.sample_id]
 frac=read(stem+'_counts.tsv');full=read((stem.replace('_both','_full_terms_both') if '_both' in stem else stem+'_full_terms')+'_counts.tsv');rpk=read(stem+'_RPK.tsv');fullrpk=read((stem.replace('_both','_full_terms_both') if '_both' in stem else stem+'_full_terms')+'_RPK.tsv')
 assert frac.index.equals(full.index) and frac.index.equals(rpk.index)
 den=rpk.loc[ref].median(axis=0).to_numpy();assert (den>0).all()
 sizes={k:set() for k in keys}
 for p in patterns:
  if p['status'] in statuses and len(p['labels'])>1:
   for k in set(p['labels'])&set(keys):sizes[k].add(len(p['labels']))
 broad={};strict={}
 for k in keys:
  if k not in frac.index:
   inventory.append(dict(variant=variant,KO=k,available=False));continue
  f=frac.loc[k].to_numpy();u=full.loc[k].to_numpy();tol=1e-6+1e-9*np.maximum(f,u);assert (u>=f-tol).all();d=np.maximum(0,u-f)
  if sizes[k]:lo=np.maximum(0,u-d/(1-1/min(sizes[k])));hi=np.maximum(0,u-d/(1-1/max(sizes[k])))
  else:assert np.allclose(u,f);lo=f.copy();hi=f.copy()
  lo[lo<=tol]=0;hi[hi<=tol]=0;assert (hi>=lo-tol).all()
  broad[k]=(f>tol).astype('int8');strict[k]=np.where(lo>tol,1,np.where(hi<=tol,0,2)).astype('int8')
  inventory.append(dict(variant=variant,KO=k,available=True,fractional_positive=int((f>tol).sum()),single_positive=int((lo>tol).sum()),single_unknown=int(((lo<=tol)&(hi>tol)).sum())))
  # Abundance is retained at KO level; no sum of required subunits is a pathway abundance.
  rf=rpk.loc[k].to_numpy();ru=fullrpk.loc[k].to_numpy();rd=np.maximum(0,ru-rf)
  rl=np.maximum(0,ru-rd/(1-1/min(sizes[k]))) if sizes[k] else rf
  rh=np.maximum(0,ru-rd/(1-1/max(sizes[k]))) if sizes[k] else rf
  for j,s in enumerate(meta.sample_id):quant.append(dict(variant=variant,KO=k,sample_id=s,fractional_RPK=rf[j],fractional_reference_ratio=rf[j]/den[j],single_RPK_lower=rl[j],single_RPK_upper=rh[j],single_reference_lower=rl[j]/den[j],single_reference_upper=rh[j]/den[j]))
 for mid,t in trees.items():
  for evidence,values in [('fractional',broad),('single_KO_bounds',strict)]:
   vv=logic.evaluate(t,values,51);bv=np.stack([logic.evaluate(logic.parse(b,expr),values,51) for b in blocks[mid]])
   for j,s in enumerate(meta.sample_id):sample_rules.append(dict(module=mid,variant=variant,evidence=evidence,sample_id=s,state=['not_met','met','unknown'][vv[j]],steps=len(bv),met_steps=int((bv[:,j]==1).sum()),unknown_steps=int((bv[:,j]==2).sum())))
   for i,b in enumerate(blocks[mid]):
    for j,s in enumerate(meta.sample_id):sample_steps.append(dict(module=mid,step=i+1,expression=b,variant=variant,evidence=evidence,sample_id=s,state=['not_met','met','unknown'][bv[i,j]]))
pd.DataFrame(inventory).to_csv(O/'KO_inventory.tsv',sep='\t',index=False)
q=pd.DataFrame(quant).merge(meta,on='sample_id',validate='many_to_one');q.to_csv(O/'KO_quantification.tsv.gz',sep='\t',index=False)
q.groupby(['variant','KO','restoration4','depth_cm']).agg(n=('sample_id','size'),median=('fractional_reference_ratio','median'),minimum=('fractional_reference_ratio','min'),maximum=('fractional_reference_ratio','max')).to_csv(O/'descriptive_KO_stage_depth.tsv',sep='\t')
rr=pd.DataFrame(sample_rules).merge(meta,on='sample_id',validate='many_to_one');rr.to_csv(O/'rule_sample_evidence.tsv.gz',sep='\t',index=False)
ss=pd.DataFrame(sample_steps);ss.to_csv(O/'step_sample_evidence.tsv.gz',sep='\t',index=False)
rr.groupby(['module','variant','evidence','restoration4','depth_cm','state']).size().rename('samples').to_csv(O/'rule_stage_depth.tsv',sep='\t')
genes=pd.read_csv(a.analysis/'prepared/selected_genes.tsv',sep='\t',keep_default_na=False);assert len(genes)==100000 and genes.gene_numeric_id.is_unique
gk=genes[['gene_numeric_id','historical_cluster_id','KO_set']].assign(KO=lambda d:d.KO_set.str.split(','));gk=gk.explode('KO');gk=gk[gk.KO.isin(keys)].drop_duplicates();gk['single_KO']=~gk.KO_set.str.contains(',')
axes=list(ora.axis.unique());member=[];assoc=[];components=[];checks=[]
old=pd.read_csv(a.old_reference/'kegg_ko_module_links.tsv',sep='\t',header=None,names=['KO','module']);old.KO=old.KO.str.replace('ko:','',regex=False);old.module=old.module.str.replace('md:','',regex=False)
for mid in mods:
 oldset=set(old.loc[old.module.eq(mid),'KO']);newset=logic.kos(trees[mid],True)
 checks.append(dict(module=mid,old_link_KOs=len(oldset),definition_KOs=len(newset),only_old_links=','.join(sorted(oldset-newset)),only_new_definition=','.join(sorted(newset-oldset))))
for axis in axes:
 common=pd.concat([d[d.family.eq('common')] for d in pd.read_csv(a.analysis/'models'/f'{axis}_all_results.tsv.gz',sep='\t',chunksize=100000)],ignore_index=True)
 assert len(common)==100000 and common.gene_numeric_id.is_unique
 joined=gk.merge(common,on='gene_numeric_id',validate='many_to_one');joined.to_csv(O/f'component_gene_estimates_{axis}.tsv.gz',sep='\t',index=False)
 for row in ora[ora.axis.eq(axis)].itertuples():
  mid=row.module;sign=1 if row.direction=='positive' else -1;allko=logic.kos(trees[mid],True);sub=joined[joined.KO.isin(allko)].copy();sub['module']=mid;sub['target_direction']=row.direction
  sub['supported_same_sign']=sub.usable & ~sub.singular & sub.q.le(.1) & (sub.estimate_sd*sign>0);member.append(sub)
  for k in sorted(allko):
   x=sub[sub.KO.eq(k)];valid=x[x.usable & ~x.singular];single=valid[valid.single_KO]
   components.append(dict(module=mid,axis=axis,direction=row.direction,KO=k,tested_genes=x.gene_numeric_id.nunique(),usable_genes=valid.gene_numeric_id.nunique(),same_sign_supported=int(valid.supported_same_sign.sum()),opposite_sign_supported=int((valid.q.le(.1)&(valid.estimate_sd*sign<0)).sum()),single_usable_genes=single.gene_numeric_id.nunique(),single_same_sign_supported=int(single.supported_same_sign.sum()),median_gene_beta=valid.estimate_sd.median(),min_gene_beta=valid.estimate_sd.min(),max_gene_beta=valid.estimate_sd.max()))
  for cutoff in [.1,.05]:
   for evidence in ['all_assignments','single_KO']:
    selected=sub[sub.usable & ~sub.singular & sub.q.le(cutoff)&(sub.estimate_sd*sign>0)]
    if evidence=='single_KO':selected=selected[selected.single_KO]
    values={k:np.array([int(k in set(selected.KO))]) for k in allko}
    for i,b in enumerate(blocks[mid]):
     state=int(logic.evaluate(logic.parse(b,expr),values,1)[0]);assoc.append(dict(module=mid,axis=axis,direction=row.direction,cutoff=cutoff,evidence=evidence,step=i+1,expression=b,associated_step=state==1,associated_gene_count=selected[selected.KO.isin(logic.kos(logic.parse(b,expr),True))].gene_numeric_id.nunique()))
 print('Association overlay',axis,flush=True)
pd.concat(member).to_csv(O/'all_module_gene_estimates.tsv.gz',sep='\t',index=False)
pd.DataFrame(components).to_csv(O/'component_associations.tsv',sep='\t',index=False)
aa=pd.DataFrame(assoc);aa.to_csv(O/'step_associations.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(O/'definition_vs_frozen_links.tsv',sep='\t',index=False)
summary=[]
for row in ora.itertuples():
 mid=row.module;z=rr[(rr.module==mid)&(rr.variant=='complete')];ass=aa[(aa.module==mid)&(aa.axis==row.axis)&(aa.direction==row.direction)&(aa.cutoff==.1)]
 rec=dict(module=mid,name=row.name,axis=row.axis,direction=row.direction,steps=len(blocks[mid]),q_ora=row.q_ora_modules,q_single=row.q_ora_single_KO,single_ora_eligible=pd.notna(row.q_ora_single_KO),fold_enrichment=row.fold_enrichment)
 for e in ['fractional','single_KO_bounds']:
  zz=z[z.evidence==e];rec[e+'_met_samples']=int(zz.state.eq('met').sum());rec[e+'_unknown_samples']=int(zz.state.eq('unknown').sum())
 for e in ['all_assignments','single_KO']:rec[e+'_associated_steps']=int(ass[ass.evidence==e].associated_step.sum())
 summary.append(rec)
pd.DataFrame(summary).to_csv(O/'block_contrast_summary.tsv',sep='\t',index=False)
pd.DataFrame([entries[k] for k in mods]).to_csv(O/'module_rules.tsv',sep='\t',index=False)
assert len(rr)==10*51*3*2 and len(ora)==11
for var in ['complete','both_mates','supported']:
 x=rr[rr.variant.eq(var)].pivot(index=['module','sample_id'],columns='evidence',values='state');assert (x.loc[x.single_KO_bounds.eq('met'),'fractional']=='met').all()
paths=[a.definitions,a.reference_genes,a.analysis/'prepared/selected_genes.tsv',a.ora/'primary_enriched_modules.tsv',a.cf_data/'03_unit_annotations/annotation_patterns.json']+list((a.cf_data/'matrices').glob('KO_*.tsv'))
(O/'validation.json').write_text(json.dumps({'status':'PASS','samples':51,'profiles':17,'modules':10,'selected_contrasts':11,'rule_sample_rows':len(rr),'logic_tests':44,'new_hypothesis_tests':0,'genes_reannotated':False,'selection':'Post-result selection of 10 modules from 11 original ORA contrasts','definition_retrieved_utc':defs['retrieved_utc'],'interpretation':'Community component evidence, not within-genome completeness or activity; association coverage is descriptive','hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}},indent=2))
print(pd.DataFrame(summary).to_string(index=False));print(O)
