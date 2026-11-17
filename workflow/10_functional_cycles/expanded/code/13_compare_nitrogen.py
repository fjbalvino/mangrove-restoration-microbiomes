from pathlib import Path
import json,hashlib,platform
import numpy as np,pandas as pd
from scipy import stats
from stats_cns import Design,adjust,profile_matrix
P=Path(__file__).resolve().parent;D=P.parent/'panel/data';O=P/'13_compare_nitrogen_results';O.mkdir(exist_ok=True)
ref=(P/'uscg_76_kegg_min_2313.lst').read_text().split();assert len(ref)==76 and hashlib.sha256((P/'uscg_76_kegg_min_2313.lst').read_bytes()).hexdigest()=='b34870d49b2d4793554e45274a2ebfdac763f69f286c6425e2f39a607aea7d88'
m=pd.read_csv(D/'metadata_003b.csv');m['profile_id']=m.sample_id.str.rsplit('_',n=1).str[0];m=m.sort_values(['profile_id','depth_cm']);p=m.groupby('profile_id',sort=True).first().reset_index();st=['Degraded','Early restoration','Intermediate restoration','Advanced restoration']
X=np.column_stack([np.ones(len(p)),*(p.locality.eq(l) for l in ['Cozumel','Tuxpan']),*(p.restoration4.eq(s) for s in st)]).astype(float);des=Design(X);assert des.df==10
keys=sorted(pd.read_csv(P.parent/'nitrogen_pilot/KO_sample_evidence.tsv',sep='\t').KO.unique());pat=json.load(open(D/'03_unit_annotations/annotation_patterns.json'))['KO']
def read(n):return pd.read_csv(D/'matrices'/n,sep='\t',index_col=0).loc[:,m.sample_id]
def bounds(r,f,ks,statuses):
 lo={};hi={}
 for k in ks:
  sizes=sorted({len(z['labels']) for z in pat if z['status'] in statuses and k in z['labels'] and len(z['labels'])>1});a=r.loc[k];b=f.loc[k];d=np.maximum(0,b-a)
  lo[k]=np.maximum(0,b-d/(1-1/min(sizes))) if sizes else a
  hi[k]=np.maximum(0,b-d/(1-1/max(sizes))) if sizes else a
 return pd.DataFrame(lo).T,pd.DataFrame(hi).T
r=read('KO_complete_RPK.tsv');f=read('KO_complete_full_terms_RPK.tsv');available=[k for k in keys if k in r.index];lo,hi=bounds(r,f,sorted(set(available+ref)),{'complete_consensus'})
exact=[k for k in available if np.allclose(lo.loc[k],hi.loc[k],rtol=1e-8,atol=1e-7)]
num=r.loc[available].copy();num.loc[exact]=lo.loc[exact];eligible=[k for k in available if (num.loc[k]>1e-7).all()]
registry=pd.DataFrame([dict(KO=k,available=k in available,quantification='single_KO' if k in exact else 'fractional_provisional',eligible=k in eligible,positive_samples=int((num.loc[k]>1e-7).sum()) if k in available else 0) for k in keys]);registry.to_csv(O/'eligibility.tsv',sep='\t',index=False)
# Freeze target and family before outcome tests; abundance-only pilot.
plan=dict(scope='exploratory nitrogen map; no new pathway-specificity validation',primary='mean log2 KO/reference over three horizons per profile ~ locality + categorical restoration',secondary='paired log2 20/5 and 40/5 per profile ~ locality + restoration; restoration omnibus tests whether depth changes differ by stage',family='Holm jointly across K eligible KOs x 3 tests',within_stage_depth='descriptive paired medians and individual profiles; no small-n sign tests',pairwise='adjusted point estimates vs Conserved; no pairwise significance or uncalibrated pairwise CI',gate='available and positive in all 51 primary samples; zeros/missing descriptive only; no pseudocount',eligible=eligible,n_tests=len(eligible)*3,draws=99999,seed=2026111701,units='17 profiles, not 51 independent samples')
(O/'plan.json').write_text(json.dumps(plan,indent=2))
strict=[k for k in ref if np.allclose(lo.loc[k],hi.loc[k],rtol=1e-8,atol=1e-7) and (lo.loc[k]>0).all()];assert len(strict)==56
variants={'primary':num.div(r.loc[ref].median()),'fractional':r.loc[available].div(r.loc[ref].median()),'reference56':num.div(lo.loc[strict].median())}
rb=read('KO_complete_both_RPK.tsv');fb=read('KO_complete_full_terms_both_RPK.tsv');lb,hb=bounds(rb,fb,available,{'complete_consensus'});nb=rb.loc[available].copy();nb.loc[exact]=lb.loc[exact];variants['both_mates']=nb.div(rb.loc[ref].median());rs=read('KO_supported_RPK.tsv');variants['supported_fractional']=rs.loc[available].div(rs.loc[ref].median())
signs=np.random.default_rng(plan['seed']).choice(np.array([-1.,1.]),size=(99999,17));tests=[];effects=[];profile_rows=[];depth_rows=[];diag=[]
for vn,v in variants.items():
 v.to_csv(O/(vn+'_normalized.tsv'),sep='\t');results=[]
 for k in eligible:
  if not (v.loc[k]>1e-7).all():
   for target in ['restoration','depth20_by_stage','depth40_by_stage']:results.append(dict(variant=vn,KO=k,target=target,p=1.,status='not_estimable_zero'))
   continue
  arr=np.log2(profile_matrix(v.loc[[k]],m,p.profile_id)[0]);means=arr.mean(axis=1)
  for target,y in [('restoration',means),('depth20_by_stage',arr[:,1]-arr[:,0]),('depth40_by_stage',arr[:,2]-arr[:,0])]:
   out=des.wild_test(y,signs);results.append(dict(variant=vn,KO=k,target=target,status='estimated',**out));beta,res,cov,stat=des.fit(y)
   for j,s in enumerate(st):effects.append(dict(variant=vn,KO=k,target=target,stage=s,log2_ratio=beta[0,3+j],ratio=2**beta[0,3+j]))
   if vn=='primary':
    residual=res[0];mse=np.sum(residual**2)/des.df;cook=residual**2/(des.p*mse)*des.h/(1-des.h)**2
    diag.append(dict(KO=k,target=target,max_cook=float(cook.max()),most_influential=p.profile_id.iloc[cook.argmax()],max_leverage=float(des.h.max()),residual_shapiro_p=stats.shapiro(residual).pvalue))
  if vn=='primary':
   for i,row in p.iterrows():
    profile_rows.append(dict(KO=k,profile_id=row.profile_id,stage=row.restoration4,locality=row.locality,mean_log2=means[i],geometric_representation=2**means[i]))
    for j,d in [(1,20),(2,40)]:depth_rows.append(dict(KO=k,profile_id=row.profile_id,stage=row.restoration4,locality=row.locality,depth=d,log2_change=arr[i,j]-arr[i,0],ratio=2**(arr[i,j]-arr[i,0])))
 for z,q in zip(results,adjust([z['p'] for z in results])):z['p_holm']=q
 tests.extend(results);print(vn,len(eligible),min(z['p_holm'] for z in results),flush=True)
 pd.DataFrame(tests).to_csv(O/'tests.tsv',sep='\t',index=False)
pd.DataFrame(effects).to_csv(O/'effects.tsv',sep='\t',index=False);pd.DataFrame(profile_rows).to_csv(O/'profile_values.tsv',sep='\t',index=False);pd.DataFrame(depth_rows).to_csv(O/'paired_depth_values.tsv',sep='\t',index=False);pd.DataFrame(diag).to_csv(O/'diagnostics.tsv',sep='\t',index=False)
m.to_csv(O/'samples.tsv',sep='\t',index=False);p.to_csv(O/'profiles.tsv',sep='\t',index=False)
print('DONE',len(eligible),'eligible;',len(exact),'exact available')
