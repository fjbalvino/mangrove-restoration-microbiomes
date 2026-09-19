#!/usr/bin/env python3
# Inputs: convergence808/inputs/; validator also needs a completed run.
# Outputs: results root/808_functional_convergence_<timestamp>/ or reference-run audit.
# Algorithmic provenance: CLR/Aitchison; Freedman-Lane profile permutations;
# locality-stratified profile bootstrap; Holm/BH families in PROTOCOL.md.
# Source SHA-256: 067c2466a71adb9faf324785a7ca21f45dded16a6552595eebde1eaf77d486e5
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
"""Conditional functional convergence; see PROTOCOL.md for the fixed analysis contract."""
import argparse,hashlib,json,platform,shutil,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import null_space
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

STAGES=['Degraded','Early restoration','Intermediate restoration','Advanced restoration']
LOCS=['Carmen','Cozumel','Tuxpan'];DEPTH=[5,20,40]
BPERM=9999;BBOOT=4999

def design(meta,interaction=False,within_depth=False):
 cols=[np.ones(len(meta))];names=['Intercept']
 for loc in sorted(meta.locality.unique())[1:]:cols.append((meta.locality==loc).astype(float));names.append('locality:'+loc)
 if not within_depth:
  for depth in DEPTH[1:]:cols.append((meta.depth_cm==depth).astype(float));names.append('depth:'+str(depth))
 for stage in STAGES[1:]:cols.append((meta.restoration4==stage).astype(float));names.append(stage)
 if interaction:
  for stage in STAGES[1:]:
   for depth in DEPTH[1:]:cols.append(((meta.restoration4==stage)&(meta.depth_cm==depth)).astype(float));names.append(stage+':'+str(depth))
 return np.array(cols,dtype=float).T,names

def correction(p,method='holm'):
 p=np.asarray(p);order=np.argsort(p);n=len(p);out=np.empty(n)
 if method=='holm':vals=np.maximum.accumulate(p[order]*(n-np.arange(n)))
 else:vals=np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1]
 out[order]=np.minimum(1,vals);return out

def fl_test(y,X,C,perms):
 """F test of C beta=0 via its constrained reduced design and synchronized FL draws."""
 rank=np.linalg.matrix_rank(X);assert rank==X.shape[1]
 q=np.linalg.matrix_rank(C);Z=X@null_space(C)
 Q=np.linalg.qr(X,mode='reduced')[0];Q0=np.linalg.qr(Z,mode='reduced')[0]
 fit0=Q0@(Q0.T@y);res=y-fit0;yp=fit0[:,None]+res[perms].T
 def stat(v):
  full=Q@(Q.T@v);red=Q0@(Q0.T@v)
  rss=((v-full)**2).sum(axis=0);ss=((full-red)**2).sum(axis=0)
  return (ss/q)/(rss/(len(y)-rank))
 observed=float(stat(y));null=stat(yp);p=(1+np.count_nonzero(null>=observed-1e-10*max(1,abs(observed))))/(len(perms)+1)
 return observed,p,null,len(y)-rank,q

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--bundle',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--out-root',type=Path);a=ap.parse_args()
 base=a.bundle.resolve();inp=base/'inputs';root=a.out_root or base/'results';out=root/('808_functional_convergence_'+datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%fZ'))
 for f in ['tables','figures','logs','audit','draws']:(out/f).mkdir(parents=True,exist_ok=False)
 started=time.monotonic();checks=[]
 def log(s):print(datetime.now(timezone.utc).strftime('%H:%M:%S')+' | '+s,flush=True)
 def ck(name,ok):
  checks.append({'check':name,'status':'PASS' if ok else 'FAIL'})
  if not ok:raise ValueError(name)
 def save(rows,name):pd.DataFrame(rows).to_csv(out/'tables'/name,sep='\t',index=False,float_format='%.17g')
 shutil.copy2(__file__,out/'logs/808_code_executed.py');shutil.copy2(base/'PROTOCOL.md',out/'PROTOCOL.md')
 log('Auditing metadata and reference profiles')
 m=pd.read_csv(inp/'metadata.csv');hist=pd.read_csv(inp/'historical_metadata.csv')
 m['profile']=m.sample_id.str.replace(r'_[123]$','',regex=True);m=m.sort_values(['locality','profile','depth_cm']).reset_index(drop=True)
 ck('51 unique samples',len(m)==51 and m.sample_id.nunique()==51)
 ck('17 complete profiles',m.profile.nunique()==17 and m.groupby('profile').depth_cm.apply(lambda x:sorted(x)==DEPTH).all())
 ck('historical labels agree',m.set_index('sample_id').restoration4.equals(hist.set_index('sample_id').loc[m.sample_id].restoration4))
 ck('profile labels constant',m.groupby('profile')[['locality','restoration4']].nunique().max().max()==1)
 refs=m[m.restoration4=='Conserved'];refmap={(r.locality,int(r.depth_cm)):r.sample_id for r in refs.itertuples()}
 ck('one conserved reference per locality and depth',len(refs)==9 and len(refmap)==9)
 ck('reference identity',set(refs.profile)=={'CC5','CZ3','TX1'})
 prof=m.drop_duplicates('profile');save(prof[['profile','locality','restoration4','name_site','year_res']],'808_profile_design.tsv')
 save(pd.crosstab(prof.locality,prof.restoration4).reset_index(),'808_stage_replication.tsv')
 matrices={};selected={};distance_frames=[]
 def distances(layer,pc):
  X=matrices[layer];keep=selected[layer];L=np.log(X[:,keep]+pc);Z=L-L.mean(axis=1,keepdims=True)
  ids={s:i for i,s in enumerate(m.sample_id)};ri=np.array([ids[refmap[(r.locality,int(r.depth_cm))]] for r in m.itertuples()])
  d=np.sqrt(((Z-Z[ri])**2).sum(axis=1));ck(layer+' pc '+str(pc)+' reference zero',np.all(d[m.restoration4=='Conserved']==0))
  return d
 for layer in ['KEGG_ko','PFAM']:
  x=pd.read_csv(inp/f'082A_{layer}_abundance_samples_x_functions.tsv.gz',sep='\t',index_col=0)
  ck(layer+' sample IDs match',set(x.index)==set(m.sample_id) and x.index.is_unique and x.columns.is_unique)
  X=x.loc[m.sample_id].to_numpy(float);ck(layer+' nonnegative finite counts',np.isfinite(X).all() and (X>=0).all())
  matrices[layer]=X;prev=(X>0).sum(axis=0);selected[layer]=prev>=26
  save(pd.DataFrame({'function_id':x.columns,'prevalence_n':prev,'included':prev>=26}),'808_'+layer+'_feature_universe.tsv')
 for layer,pc,scenario in [('KEGG_ko',1,'KO_primary'),('KEGG_ko',.5,'KO_pc0.5'),('KEGG_ko',2,'KO_pc2'),('PFAM',1,'PFAM_sensitivity')]:
  d=distances(layer,pc);frame=m[['sample_id','profile','locality','depth_cm','restoration4']].copy();frame['scenario']=scenario;frame['n_features']=selected[layer].sum();frame['reference_sample_id']=[refmap[(r.locality,int(r.depth_cm))] for r in m.itertuples()];frame['distance']=d;frame['included_in_inference']=m.restoration4!='Conserved';distance_frames.append(frame)
 all_d=pd.concat(distance_frames,ignore_index=True);save(all_d,'808_distances_all_scenarios.tsv')
 data=distance_frames[0].query('included_in_inference').reset_index(drop=True);save(data,'808_primary_analytical_data.tsv')
 ck('42 nonreference samples in 14 profiles',len(data)==42 and data.profile.nunique()==14)
 pdx=data.drop_duplicates('profile').reset_index(drop=True);nprof=len(pdx);profile_rows=[np.where(data.profile==v)[0] for v in pdx.profile]
 for r in profile_rows:ck('complete inference profile '+data.profile.iloc[r[0]],list(data.depth_cm.iloc[r])==DEPTH)
 rng=np.random.default_rng(8081709);permprof=np.tile(np.arange(nprof),(BPERM,1));brng=np.random.default_rng(8081710);bootprof=np.tile(np.arange(nprof),(BBOOT,1))
 for loc in LOCS:
  ix=np.where(pdx.locality==loc)[0]
  for j in range(BPERM):permprof[j,ix]=rng.permutation(ix)
  bootprof[:,ix]=brng.choice(ix,size=(BBOOT,len(ix)),replace=True)
 rowblock=np.array(profile_rows);perms=rowblock[permprof].reshape(BPERM,-1);boots=rowblock[bootprof].reshape(BBOOT,-1)
 ck('permutations preserve locality and depth',np.all(data.locality.to_numpy()[perms]==data.locality.to_numpy()) and np.all(data.depth_cm.to_numpy()[perms]==data.depth_cm.to_numpy()))
 np.savez_compressed(out/'draws/808_resampling.npz',sample_ids=data.sample_id.to_numpy(dtype=str),profile_ids=pdx.profile.to_numpy(dtype=str),permutation_indices_base0=perms,bootstrap_indices_base0=boots,profile_permutation_indices_base0=permprof,profile_bootstrap_indices_base0=bootprof)
 y=data.distance.to_numpy();X,names=design(data);XI,inames=design(data,interaction=True);effects=[];tests=[];nulls={}
 def test(name,xx,cols):
  C=np.eye(xx.shape[1])[cols];F,p,null,df,q=fl_test(y,xx,C,perms);nulls[name]=null;tests.append(dict(test=name,F=F,df_numerator=q,df_denominator=df,p_raw=p,n_permutations=BPERM));log(name+' P='+str(p))
 log('Primary five-test family: 9,999 profile-preserving permutations')
 stagecols=[names.index(s) for s in STAGES[1:]];test('Stage_omnibus',X,stagecols)
 for s in STAGES[1:]:test(s+'_vs_Degraded',X,[names.index(s)])
 test('Stage_x_depth',XI,list(range(X.shape[1],XI.shape[1])))
 adj=correction([t['p_raw'] for t in tests]);
 for t,p in zip(tests,adj):t['p_holm_five_tests']=p
 save(tests,'808_primary_tests.tsv')
 log('Bootstrap confidence intervals: 4,999 stratified profile draws')
 beta=np.linalg.lstsq(X,y,rcond=None)[0];bb=np.full((BBOOT,3),np.nan)
 for j,idx in enumerate(boots):
  if np.linalg.matrix_rank(X[idx])==X.shape[1]:bb[j]=np.linalg.lstsq(X[idx],y[idx],rcond=None)[0][stagecols]
 for k,s in enumerate(STAGES[1:]):
  valid=bb[:,k][np.isfinite(bb[:,k])];lo,hi=np.quantile(valid,[.025,.975]);t=tests[k+1]
  effects.append(dict(stage=s,contrast='stage minus Degraded',adjusted_distance_difference=beta[stagecols[k]],ci_low=lo,ci_high=hi,n_bootstrap_valid=len(valid),n_bootstrap_attempted=BBOOT,p_raw=t['p_raw'],p_holm_five_tests=t['p_holm_five_tests'],direction='closer' if beta[stagecols[k]]<0 else 'farther'))
 save(effects,'808_primary_stage_contrasts.tsv');depthrows=[];depthboots=[]
 log('Exploratory depth-specific contrasts')
 for depth in DEPTH:
  ids=np.where(data.depth_cm==depth)[0];md=data.iloc[ids];yd=y[ids];xd,nd=design(md,within_depth=True);sc=[nd.index(s) for s in STAGES[1:]];bd=np.linalg.lstsq(xd,yd,rcond=None)[0];bds=np.full((BBOOT,3),np.nan)
  for j,idx in enumerate(bootprof):
   if np.linalg.matrix_rank(xd[idx])==xd.shape[1]:bds[j]=np.linalg.lstsq(xd[idx],yd[idx],rcond=None)[0][sc]
  for k,s in enumerate(STAGES[1:]):
   F,p,nul,df,q=fl_test(yd,xd,np.eye(xd.shape[1])[[sc[k]]],permprof);nulls[str(depth)+'_'+s]=nul;valid=bds[:,k][np.isfinite(bds[:,k])];lo,hi=np.quantile(valid,[.025,.975]);depthrows.append(dict(depth_cm=depth,stage=s,adjusted_distance_difference=bd[sc[k]],ci_low=lo,ci_high=hi,p_raw=p,n_profiles=14,n_permutations=BPERM,n_bootstrap_valid=len(valid),n_bootstrap_attempted=BBOOT))
  depthboots.append(bds)
 for r,p in zip(depthrows,correction([r['p_raw'] for r in depthrows],'bh')):r['q_BH_nine_contrasts']=p
 save(depthrows,'808_depth_contrasts_exploratory.tsv');np.savez_compressed(out/'draws/808_null_and_bootstrap.npz',primary_bootstrap_beta=bb,depth_bootstrap_beta=np.array(depthboots),**nulls)
 log('Descriptive preprocessing and locality-exclusion sensitivities')
 sens=[]
 for frame in distance_frames:
  sub=frame.query('included_in_inference').reset_index(drop=True);xx,nn=design(sub);b=np.linalg.lstsq(xx,sub.distance,rcond=None)[0]
  for s in STAGES[1:]:sens.append(dict(scenario=sub.scenario.iloc[0],omitted_locality='none',stage=s,n_profiles=sub.profile.nunique(),n_features=int(sub.n_features.iloc[0]),adjusted_distance_difference=b[nn.index(s)]))
 for loc in LOCS:
  sub=data[data.locality!=loc];xx,nn=design(sub);ck('full rank omit '+loc,np.linalg.matrix_rank(xx)==xx.shape[1]);b=np.linalg.lstsq(xx,sub.distance,rcond=None)[0]
  for s in STAGES[1:]:sens.append(dict(scenario='KO_primary',omitted_locality=loc,stage=s,n_profiles=sub.profile.nunique(),n_features=int(sub.n_features.iloc[0]),adjusted_distance_difference=b[nn.index(s)]))
 save(sens,'808_descriptive_sensitivity.tsv')
 save(data.groupby(['locality','depth_cm','restoration4'],sort=False).agg(n_profiles=('profile','nunique'),mean_distance=('distance','mean'),min_distance=('distance','min'),max_distance=('distance','max')).reset_index(),'808_observed_distances_by_stratum.tsv')
 ck('all tests within 0 and 1',all(0<=t['p_raw']<=1 and 0<=t['p_holm_five_tests']<=1 for t in tests))
 ck('bootstrap coverage above 80 percent',np.isfinite(bb).all(axis=1).mean()>.8)
 log('Drawing observed distances and adjusted contrasts')
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
 fig,axs=plt.subplots(1,3,figsize=(10,4),sharey=True)
 colors=['#0072B2','#D55E00','#009E73']
 for ax,depth in zip(axs,DEPTH):
  for j,loc in enumerate(LOCS):
   dd=data[(data.depth_cm==depth)&(data.locality==loc)]
   xpos=[STAGES.index(v)+(j-1)*.16 for v in dd.restoration4]
   ax.scatter(xpos,dd.distance,color=colors[j],label=loc,s=40,alpha=.85)
  ax.set_xticks(range(4),['Degraded','Early','Intermediate','Advanced'],rotation=25,ha='right');ax.set_title(str(depth)+' cm');ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.2)
 axs[0].set_ylabel('Distance to conserved reference (Aitchison)');axs[-1].legend(frameon=False,fontsize=8);fig.tight_layout()
 for ext in ['pdf','svg','png']:fig.savefig(out/'figures'/('808_distances_by_depth.'+ext),dpi=300,facecolor='white')
 plt.close(fig);fig,axs=plt.subplots(2,2,figsize=(9,7),sharex=True)
 for ax,label,rows in [(axs[0,0],'Common stage effects',effects)]+[(a,str(d)+' cm',[r for r in depthrows if r['depth_cm']==d]) for a,d in zip([axs[0,1],axs[1,0],axs[1,1]],DEPTH)]:
  for i,r in enumerate(rows):ax.plot([r['ci_low'],r['ci_high']],[2-i,2-i],color='#0072B2');ax.plot(r['adjusted_distance_difference'],2-i,'o',color='#0072B2')
  ax.axvline(0,color='gray',linestyle='--');ax.set_yticks([2,1,0],['Early − degraded','Intermediate − degraded','Advanced − degraded']);ax.set_title(label);ax.set_xlabel('Adjusted distance difference');ax.set_ylim(-.5,2.5)
 fig.tight_layout()
 for ext in ['pdf','svg','png']:fig.savefig(out/'figures'/('808_stage_contrasts.'+ext),dpi=300,facecolor='white')
 plt.close(fig);save(checks,'808_acceptance_checks.tsv')
 provenance={'status':'PASS','scope':'conditional reference-based functional convergence','n_samples_total':51,'n_profiles_total':17,'n_samples_inference':42,'n_profiles_inference':14,'n_conserved_reference_profiles':3,'n_KO_primary':int(selected['KEGG_ko'].sum()),'n_PFAM_sensitivity':int(selected['PFAM'].sum()),'n_permutations':BPERM,'n_bootstrap_attempted':BBOOT,'n_bootstrap_valid':int(np.isfinite(bb).all(axis=1).sum()),'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inp.iterdir()},'protocol_sha256':hashlib.sha256((base/'PROTOCOL.md').read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-started,'reference_uncertainty_included':False,'R_executed':False}
 (out/'808_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n');(root/'LATEST_808.txt').write_text(str(out)+'\n');log('PASS; output '+str(out))
if __name__=='__main__':main()
