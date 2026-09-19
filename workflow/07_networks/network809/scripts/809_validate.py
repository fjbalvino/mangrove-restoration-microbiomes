#!/usr/bin/env python3
# Inputs: network809/inputs/{clr_200_KO.tsv.gz,metadata.csv}.
# Diagnostics additionally consume the completed 809 run (objects, draws and tables).
# Outputs: configured results root/809_residualized_functional_networks_*/; see PROTOCOL.md.
# Provenance: fixed-lambda ridge + LIONESS (Kuijjer et al., 2019, doi:10.1016/j.isci.2019.03.021).
# Curation: OpenAI Codex (OpenAI, 2026); documentation only; algorithm body preserved.
# Module filenames remain unchanged to preserve Python imports and historical source-hash checks.
"""Independent network audit via QR residuals and sample-space Woodbury precision."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
from scipy.linalg import solve,expm
from scipy.special import logsumexp
from scipy.spatial.distance import pdist,squareform
ap=argparse.ArgumentParser();ap.add_argument('--bundle',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--run',type=Path);a0=ap.parse_args();base=a0.bundle.resolve();run=a0.run or Path((base/'results/LATEST_809.txt').read_text().strip());rows=[]
def ck(name,ok,error=None):
 rows.append(dict(check=name,status='PASS' if ok else 'FAIL',max_error=error))
 if not ok:raise ValueError((name,error))
def eq(name,x,y,tol=1e-8):
 e=float(np.nanmax(abs(np.asarray(x)-np.asarray(y))));ck(name,e<tol,e)
def read(n):return pd.read_csv(run/'tables'/n,sep='\t',float_precision='round_trip')
SC=['unadjusted','space_depth','space_depth_environment'];IND=['HI','MHI_local'];ENDS=['mean_abs_edge_weight','natural_connectivity_abs_weighted'];meta=read('809_metadata_aligned.tsv');ids=meta.sample_id.to_numpy();source=pd.read_csv(base/'inputs/clr_200_KO.tsv.gz',sep='\t',index_col=0);z0=np.load(run/'objects/809_unadjusted_networks.npz');nodes=z0['nodes'];x=source.loc[ids,nodes].to_numpy();a=z0['node_a_index'];b=z0['node_b_index'];n=200
D=np.column_stack([np.ones(51),(meta.depth_cm==20).astype(float),(meta.depth_cm==40).astype(float),pd.get_dummies(meta.locality,dtype=float)[['Cozumel','Tuxpan']]]);plans=np.load(run/'draws/809_schedules.npz');P=plans['permutations_base0'];tests=read('809_index_associations.tsv');globaltab=read('809_global_configuration.tsv');mechtab=read('809_magnitude_redistribution.tsv')
def ridge_woodbury(x,C):
 Q=np.linalg.qr(C,mode='reduced')[0];r=x-Q@(Q.T@x);Z=(r-r.mean(0))/r.std(0,ddof=1);aa=.1+1e-8;U=np.sqrt(.9/(len(x)-1))*Z;V=np.eye(n)/aa-U.T@solve(np.eye(len(x))+U@U.T/aa,U,assume_a='pos')/(aa*aa);pc=-V/np.sqrt(np.outer(V.diagonal(),V.diagonal()));return pc[a,b]
def select(v):return np.lexsort((nodes[b],nodes[a],-abs(v)))[:995]
def met(v):
 ii=select(v);A=np.zeros((n,n));A[a[ii],b[ii]]=abs(v[ii]);A+=A.T
 return [abs(v[ii]).mean(),logsumexp(np.linalg.eigvalsh(A))-np.log(n)]
def independent_lion(x,C):
 w=ridge_woodbury(x,C);W=[]
 for i in range(len(x)):
  ii=np.arange(len(x))!=i;W.append(len(x)*w-(len(x)-1)*ridge_woodbury(x[ii],C[ii]))
 return np.array(W)
def holm(p):
 ix=np.argsort(p);v=np.empty(len(p));v[ix]=np.minimum(1,np.maximum.accumulate(np.array(p)[ix]*np.arange(len(p),0,-1)));return v
for sid in SC:
 z=np.load(run/'objects'/('809_'+sid+'_networks.npz'));C=z['nuisance_design'];W=z['weights'];Y=z['sample_metrics'];checkw=independent_lion(x,C);eq(sid+' all 1014900 LIONESS weights',checkw,W,1e-9)
 eq(sid+' fitted CLR decomposition',z['residual_CLR']+C@z['nuisance_coefficients'],x);eq(sid+' residual orthogonality',C.T@z['residual_CLR'],0)
 for i in [0,25,50]:
  ii=select(W[i]);A=np.zeros((n,n));A[a[ii],b[ii]]=abs(W[i,ii]);A+=A.T;eq(sid+' natural connectivity matrix exponential sample '+str(i),np.log(np.trace(expm(A))/n),Y[i,1]);eq(sid+' mean weight sample '+str(i),abs(W[i,ii]).mean(),Y[i,0])
 st=read('809_'+sid+'_edge_stability.tsv.gz');det=np.zeros(len(a));pos=det.copy();neg=det.copy()
 for ii in plans['network_bootstrap_base0']:
  w=ridge_woodbury(x[ii],C[ii]);ss=select(w);det[ss]+=1;pos[ss]+=w[ss]>0;neg[ss]+=w[ss]<0
 eq(sid+' all bootstrap selection frequencies',det/200,st.selection_frequency,1e-12);cons=np.divide(np.maximum(pos,neg),det,out=np.zeros(len(a)),where=det>0);eq(sid+' sign consistency',cons,st.sign_consistency,1e-12);ck(sid+' consensus edge identity',np.array_equal((det>=100)&(cons>=.8),z['stable_mask']))
 tt=tests[tests.construction==sid].reset_index(drop=True);draw=np.load(run/'draws'/('809_'+sid+'_conditional_index.npz'));null=draw['null_t'];ystd=(Y-Y.mean(0))/Y.std(0,ddof=1)
 for j,p in enumerate(IND):
  raw=meta[p].to_numpy();xx=(raw-raw.mean())/raw.std(ddof=1);X=np.column_stack([D,xx]);V=np.linalg.inv(X.T@X);bb=V@X.T@ystd;eq(sid+' full OLS beta '+p,bb[-1],tt.loc[tt.predictor==p,'standardized_beta'])
  cross=np.empty((len(P),6,2));cross[:,:5,:]=D.T@ystd;cross[:,-1,:]=xx[P]@ystd;betas=np.einsum('ij,bjk->bik',V,cross);rss=(ystd**2).sum(0)-np.einsum('bij,bij->bj',betas,cross);T=betas[:,-1,:]/np.sqrt(rss/45*V[-1,-1]);eq(sid+' full OLS 99999 null t '+p,T,null[:,j*2:j*2+2],1e-8)
 obs=tt.statistic.to_numpy();pp=(1+(abs(null)>=abs(obs)).sum(0))/(len(P)+1);eq(sid+' index P',pp,tt.p_raw,1e-12);eq(sid+' within-family Holm',holm(pp),tt.p_holm_within_construction,1e-12);eq(sid+' within-family maxT',(1+(abs(null).max(1)[:,None]>=abs(obs)).sum(0))/(len(P)+1),tt.p_maxT_within_construction,1e-12)
 for j in range(4):eq(sid+' conditional bootstrap CI '+str(j),np.quantile(draw['bootstrap_beta'][:,j],[.025,.975]),tt.loc[j,['ci_low','ci_high']].to_numpy(float))
 # Global pseudo-F and R2 through full projection matrices, independent of residualized-index implementation.
 hi=meta.HI.to_numpy();XF=np.column_stack([D,hi]);H=XF@np.linalg.inv(XF.T@XF)@XF.T;H0=D@np.linalg.inv(D.T@D)@D.T;J=np.eye(51)-np.ones((51,51))/51
 pres=np.zeros_like(W,dtype=bool)
 for i in range(51):pres[i,select(W[i])]=True
 norm=np.linalg.norm(W,axis=1);vectors=[abs(W),W,np.log(norm)[:,None],abs(W)/norm[:,None]];gs=[J@(v@v.T)@J for v in vectors];jd=squareform(pdist(pres,'jaccard'));gs.insert(2,-.5*J@jd@J)
 tab=pd.concat([globaltab[globaltab.construction==sid],mechtab[mechtab.construction==sid]],ignore_index=True)
 for i,G in enumerate(gs):
  ss=np.trace((H-H0)@G);res=np.trace((np.eye(51)-H)@G);red=np.trace((np.eye(51)-H0)@G);eq(sid+' global F '+str(i),45*ss/res,tab.statistic.iloc[i]);eq(sid+' global R2 '+str(i),ss/red,tab.partial_R2.iloc[i])
 # Three complete bootstrap reconstructions verify case multiplicities and repeated nuisance fitting.
 rf=np.load(run/'draws/809_reconstruction_bootstrap.npz')['standardized_beta'];si=SC.index(sid)
 for bi in [0,123,498]:
  ii=plans['refit_bootstrap_base0'][bi]
  if not np.isfinite(rf[bi,si]).all():continue
  ww=independent_lion(x[ii],C[ii]);yy=np.array([met(v) for v in ww]);ys=(yy-yy.mean(0))/yy.std(0,ddof=1);bet=[]
  for p in IND:
   xv=meta[p].to_numpy()[ii];xv=(xv-xv.mean())/xv.std(ddof=1);xx=np.column_stack([D[ii],xv]);bet.extend(np.linalg.solve(xx.T@xx,xx.T@ys)[-1])
  eq(sid+' refitted bootstrap '+str(bi),bet,rf[bi,si],1e-8)
new=tests[tests.construction!='unadjusted'];N=np.column_stack([np.load(run/'draws'/('809_'+s+'_conditional_index.npz'))['null_t'] for s in SC[1:]]);obs=new.statistic.to_numpy();pp=(1+(abs(N)>=abs(obs)).sum(0))/(len(P)+1);eq('eight-test Holm',holm(pp),new.p_holm_eight_new_tests,1e-12);eq('eight-test maxT',(1+(abs(N).max(1)[:,None]>=abs(obs)).sum(0))/(len(P)+1),new.p_maxT_eight_new_tests,1e-12)
for file,suffix in [('809_global_configuration.tsv','six_new_configuration_tests'),('809_magnitude_redistribution.tsv','four_new_mechanism_tests')]:
 t=read(file);t=t[t.construction!='unadjusted'];eq(file+' pooled Holm',holm(t.p_raw.to_numpy()),t['p_holm_'+suffix],1e-12)
ck('all 99999 profile permutations respect locality and depth',np.all(meta.locality.to_numpy()[P]==meta.locality.to_numpy()) and np.all(meta.depth_cm.to_numpy()[P]==meta.depth_cm.to_numpy()))
for name in ['conditional_bootstrap_base0','network_bootstrap_base0']:
 ii=plans[name];ck(name+' each sampled profile retains depths',all(np.all((meta.depth_cm.to_numpy()[ii]==depth).sum(1)==17) for depth in [5,20,40]))
prov=json.loads((run/'809_provenance.json').read_text())
for name,h in prov['input_sha256'].items():ck('input hash '+name,hashlib.sha256((base/'inputs'/name).read_bytes()).hexdigest()==h)
ck('executed code hash',hashlib.sha256((run/'logs/809_residual_networks.py').read_bytes()).hexdigest()==prov['code_sha256']);ck('protocol hash',hashlib.sha256((run/'PROTOCOL.md').read_bytes()).hexdigest()==prov['protocol_sha256'])
pd.DataFrame(rows).to_csv(run/'audit/809_independent_checks.tsv',sep='\t',index=False);print('PASS',len(rows),'independent checks')
