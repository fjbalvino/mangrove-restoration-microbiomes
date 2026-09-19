#!/usr/bin/env python3
# Inputs: convergence808/inputs/; validator also needs a completed run.
# Outputs: results root/808_functional_convergence_<timestamp>/ or reference-run audit.
# Algorithmic provenance: CLR/Aitchison; Freedman-Lane profile permutations;
# locality-stratified profile bootstrap; Holm/BH families in PROTOCOL.md.
# Source SHA-256: 24a7b4ee9c6a98e6f1c4ac483f5f44c90bd3fe8fbed624da41b3915ea3df31ca
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
"""Independent numeric audit using log-ratio distances and Wald-form F statistics."""
from pathlib import Path
import argparse,json,hashlib,math
import numpy as np
import pandas as pd
ap=argparse.ArgumentParser();ap.add_argument('--bundle',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--run',type=Path);a=ap.parse_args();b=a.bundle.resolve();out=a.run or Path((b/'results/LATEST_808.txt').read_text().strip());rows=[]
def ck(name,ok,err=None):
 rows.append({'check':name,'status':'PASS' if ok else 'FAIL','max_error':err})
 if not ok:raise AssertionError((name,err))
def eq(name,x,y,tol=1e-8):
 e=float(np.max(np.abs(np.asarray(x)-np.asarray(y))));ck(name,e<tol,e)
prov=json.loads((out/'808_provenance.json').read_text())
for name,h in prov['input_sha256'].items():ck('hash '+name,hashlib.sha256((b/'inputs'/name).read_bytes()).hexdigest()==h)
ck('executed code hash',hashlib.sha256((out/'logs/808_code_executed.py').read_bytes()).hexdigest()==prov['code_sha256'])
ck('protocol hash',hashlib.sha256((out/'PROTOCOL.md').read_bytes()).hexdigest()==prov['protocol_sha256'])
raw=pd.read_csv(b/'inputs/082A_KEGG_ko_abundance_samples_x_functions.tsv.gz',sep='\t',index_col=0);keep=(raw>0).sum()>=26
data=pd.read_csv(out/'tables/808_primary_analytical_data.tsv',sep='\t');full=pd.read_csv(out/'tables/808_distances_all_scenarios.tsv',sep='\t');full=full[full.scenario=='KO_primary'];ind=[]
for r in full.itertuples():
 v=np.log((raw.loc[r.sample_id,keep].to_numpy()+1)/(raw.loc[r.reference_sample_id,keep].to_numpy()+1));mean=math.fsum(v)/len(v);ind.append(math.sqrt(math.fsum((float(w)-mean)**2 for w in v)))
eq('51 independent log-ratio distances',ind,full.distance)
ck('no conserved rows enter inference',not (data.restoration4=='Conserved').any())
stages=['Early restoration','Intermediate restoration','Advanced restoration']
X=np.column_stack([np.ones(len(data)),pd.get_dummies(data.locality,dtype=float)[['Cozumel','Tuxpan']],pd.get_dummies(data.depth_cm,dtype=float)[[20,40]],pd.get_dummies(data.restoration4,dtype=float)[stages]])
extras=np.column_stack([((data.restoration4==s)&(data.depth_cm==d)).astype(float) for s in stages for d in [20,40]]);XI=np.column_stack([X,extras]);y=data.distance.to_numpy()
resamp=np.load(out/'draws/808_resampling.npz');saved=np.load(out/'draws/808_null_and_bootstrap.npz');idx=resamp['permutation_indices_base0'];obs=pd.read_csv(out/'tables/808_primary_tests.tsv',sep='\t');contr=pd.read_csv(out/'tables/808_primary_stage_contrasts.tsv',sep='\t')
beta=np.linalg.solve(X.T@X,X.T@y);eq('stage coefficients via normal equations',beta[-3:],contr.adjusted_distance_difference)
for row in obs.itertuples():
 xx=XI if row.test=='Stage_x_depth' else X
 ids=list(range(8,14)) if row.test=='Stage_x_depth' else list(range(5,8)) if row.test=='Stage_omnibus' else [5+stages.index(row.test.replace('_vs_Degraded',''))]
 C=np.eye(xx.shape[1])[ids];V=np.linalg.inv(xx.T@xx);M=xx@V;beta0=V@xx.T@y;W=np.linalg.inv(C@V@C.T);restricted=beta0-V@C.T@W@(C@beta0);fit=xx@restricted;Y=fit[:,None]+(y-fit)[idx].T
 def statistic(Y):
  bb=V@xx.T@Y;cb=C@bb;rss=((Y-xx@bb)**2).sum(axis=0)
  num=np.sum(cb*(W@cb),axis=0)/len(ids)
  return num/(rss/(len(y)-xx.shape[1]))
 eq(row.test+' observed Wald F',statistic(y),row.F)
 n=statistic(Y);eq(row.test+' all 9999 Wald null F',n,saved[row.test],1e-7)
 p=(1+np.count_nonzero(n>=row.F-1e-10*max(1,abs(row.F))))/10000;eq(row.test+' P',p,row.p_raw,1e-12)
order=np.argsort(obs.p_raw.to_numpy());h=np.empty(5);h[order]=np.minimum(1,np.maximum.accumulate(obs.p_raw.to_numpy()[order]*np.arange(5,0,-1)));eq('five-test Holm',h,obs.p_holm_five_tests)
boot=resamp['bootstrap_indices_base0'];valid=[]
for i in range(len(boot)):
 ii=boot[i];x=X[ii]
 if np.linalg.matrix_rank(x)==8:valid.append(np.linalg.solve(x.T@x,x.T@y[ii])[-3:])
valid=np.array(valid);eq('valid bootstrap count',len(valid),prov['n_bootstrap_valid']);ci=np.quantile(valid,[.025,.975],axis=0);eq('all primary bootstrap intervals',ci,np.array([contr.ci_low,contr.ci_high]))
depth=pd.read_csv(out/'tables/808_depth_contrasts_exploratory.tsv',sep='\t');pp=depth.p_raw.to_numpy();ix=np.argsort(pp);q=np.empty(9);q[ix]=np.minimum(1,np.minimum.accumulate((pp[ix]*9/np.arange(1,10))[::-1])[::-1]);eq('nine-test BH',q,depth.q_BH_nine_contrasts)
ck('all permutation rows are bijections',np.all(np.sort(idx,axis=1)==np.arange(42)))
ck('permutation locality preserved',np.all(data.locality.to_numpy()[idx]==data.locality.to_numpy()))
ck('permutation depth preserved',np.all(data.depth_cm.to_numpy()[idx]==data.depth_cm.to_numpy()))
ck('all three depths stay within one source profile',np.all(data.profile.to_numpy()[idx].reshape(-1,14,3)==data.profile.to_numpy()[idx].reshape(-1,14,3)[:,:,:1]))
pd.DataFrame(rows).to_csv(out/'audit/808_independent_checks.tsv',sep='\t',index=False);print('PASS',len(rows),'independent checks')
