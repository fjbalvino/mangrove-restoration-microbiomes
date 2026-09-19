#!/usr/bin/env python3
# Inputs: network809/inputs/{clr_200_KO.tsv.gz,metadata.csv}.
# Diagnostics additionally consume the completed 809 run (objects, draws and tables).
# Outputs: configured results root/809_residualized_functional_networks_*/; see PROTOCOL.md.
# Provenance: fixed-lambda ridge + LIONESS (Kuijjer et al., 2019, doi:10.1016/j.isci.2019.03.021).
# Curation: OpenAI Codex (OpenAI, 2026); minimal-input loader and deterministic schedule regeneration; numerical methods unchanged.
# Module filenames remain unchanged to preserve Python imports and historical source-hash checks.
"""Matched residualized network sensitivity; analysis contract in PROTOCOL.md."""
import argparse,json,hashlib,shutil,time,platform
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import scipy
from scipy.linalg import cho_factor,cho_solve,eigvalsh,null_space
from scipy.special import logsumexp
from scipy.spatial.distance import pdist,squareform
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import network_statistics as stats

SC=['unadjusted','space_depth','space_depth_environment']
ENDS=['mean_abs_edge_weight','natural_connectivity_abs_weighted'];IND=['HI','MHI_local']
ENV=['vegetation_landscape_PC1','water_inundation_PC1','physicochemical_PC1','nutrients_redox_PC1']
LAMBDA=.1;EPS=1e-8;NEDGE=995;BNET=200;BREFIT=499

def ridge_pipeline(x,c):
 coef=np.linalg.lstsq(c,x,rcond=None)[0];r=x-c@coef;sd=r.std(0,ddof=1);bad=sd<=np.sqrt(np.finfo(float).eps)
 z=r-r.mean(0);z[:,~bad]/=sd[~bad];z[:,bad]=0
 S=z.T@z/(len(x)-1);P=cho_solve(cho_factor(.9*S+(.1+EPS)*np.eye(x.shape[1]),lower=True),np.eye(x.shape[1]))
 w=-P/np.sqrt(np.outer(P.diagonal(),P.diagonal()));np.fill_diagonal(w,0)
 return (w+w.T)/2,int(bad.sum()),r,coef

def selected(v,a,b,nodes):return np.lexsort((nodes[b],nodes[a],-abs(v)))[:NEDGE]

def metric(v,sel,a,b,n):
 if len(sel)==0:return np.array([np.nan,np.nan])
 A=np.zeros((n,n));A[a[sel],b[sel]]=abs(v[sel]);A+=A.T
 return np.array([abs(v[sel]).mean(),logsumexp(eigvalsh(A,check_finite=False))-np.log(n)])

def lioness(x,c,a,b,nodes,source_ids=None):
 n=len(x);M,z,r,coef=ridge_pipeline(x,c);w=M[a,b];W=np.empty((n,len(a)));Y=np.empty((n,2));rankfull=np.linalg.matrix_rank(c);minrank=rankfull;bad=z
 cache={}
 for i in range(n):
  key=i if source_ids is None else int(source_ids[i])
  if key in cache:W[i],Y[i]=cache[key];continue
  keep=np.arange(n)!=i;mi,zi,_,_=ridge_pipeline(x[keep],c[keep]);bad=max(bad,zi);minrank=min(minrank,np.linalg.matrix_rank(c[keep]));wi=n*w-(n-1)*mi[a,b];yi=metric(wi,selected(wi,a,b,nodes),a,b,len(nodes));W[i]=wi;Y[i]=yi;cache[key]=(wi,yi)
 return W,Y,w,r,coef,dict(full_rank=int(rankfull),minimum_LOO_rank=int(minrank),maximum_zero_variance_columns=int(bad))

def adjust_family(rows,null,scope):
 obs=np.array([r['statistic'] for r in rows]);p=stats.pvalue(obs,null);h=stats.adjust_p(p,'holm');mt=stats.pvalue(obs,np.abs(null).max(1)[:,None])
 for j,r in enumerate(rows):r['p_raw']=float(p[j]);r['p_holm_'+scope]=float(h[j]);r['p_maxT_'+scope]=float(mt[j])

def index_inference(m,Y,D,perm,boots,predictors=IND):
 X=m[predictors].to_numpy();bb,rr,rho,t=stats.adjusted(X,Y,D);null=stats.perm_index(X,Y,D,perm);boot=np.full((len(boots),len(predictors)*Y.shape[1]),np.nan)
 for k,ii in enumerate(boots):
  if np.linalg.matrix_rank(D[ii])==D.shape[1]:boot[k]=stats.adjusted(X[ii],Y[ii],D[ii])[0].ravel()
 rows=[]
 for j,p in enumerate(predictors):
  for k,e in enumerate(ENDS):
   v=boot[:,j*2+k];v=v[np.isfinite(v)];lo,hi=np.quantile(v,[.025,.975]);rows.append(dict(predictor=p,endpoint=e,standardized_beta=float(bb[j,k]),ci_low=float(lo),ci_high=float(hi),partial_r2=float(rho[j,k]**2),statistic=float(t[j,k]),n_permutations=len(perm),n_bootstrap_valid=len(v),interval_scope='networks fixed'))
 adjust_family(rows,null,'within_construction');return rows,null,boot

def kernel(v):
 if v.ndim==1:v=v[:,None]
 z=v-v.mean(0);return z@z.T

def gower(d):
 z=d*d;return -.5*(z-z.mean(0)-z.mean(1)[:,None]+z.mean())

def global_tests(Gs,names,m,D,perm):
 q=np.linalg.qr(D,mode='reduced')[0];x=m.HI.to_numpy();xr=x-q@(q.T@x);u=xr/np.linalg.norm(xr);xp=x[perm].T;xp-=q@(q.T@xp);xp/=np.linalg.norm(xp,axis=0);df=len(x)-D.shape[1]-1;rows=[];null=[]
 for name,g in zip(names,Gs):
  red=np.trace(g)-np.trace(q.T@g@q);inc=float(u@g@u);F=df*inc/(red-inc);pn=np.einsum('ij,ij->j',xp,g@xp);nf=df*pn/(red-pn);rows.append(dict(response=name,partial_R2=inc/red,total_R2=inc/np.trace(g),statistic=F,n_permutations=len(perm)));null.append(nf)
 null=np.array(null).T;adjust_family(rows,null,'within_construction');return rows,null

def interaction_test(m,y,D,perm):
 x=m.HI.to_numpy();x=(x-x.mean())/x.std(ddof=1);Z=np.column_stack([D,x]);X=np.column_stack([Z,x*(m.depth_cm==20),x*(m.depth_cm==40)]);q0=np.linalg.qr(Z,mode='reduced')[0];q=np.linalg.qr(X,mode='reduced')[0];fit=q0@(q0.T@y);res=y-fit;yp=fit[:,None]+res[perm].T;df=len(y)-X.shape[1]
 def stat(v):return np.sum((q@(q.T@v)-q0@(q0.T@v))**2,axis=0)/2/(np.sum((v-q@(q.T@v))**2,axis=0)/df)
 obs=float(stat(y));nul=stat(yp);return obs,float((1+np.count_nonzero(nul>=obs-1e-10*max(1,obs)))/(len(perm)+1)),nul

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--bundle',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--out-root',type=Path);args=ap.parse_args();base=args.bundle.resolve();inp=base/'inputs';root=args.out_root or base/'results';out=root/('809_residualized_functional_networks_'+datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%fZ'))
 for f in ['tables','objects','draws','figures','audit','logs']:(out/f).mkdir(parents=True,exist_ok=False)
 start=time.monotonic();checks=[]
 def log(s):print(datetime.now(timezone.utc).strftime('%H:%M:%S')+' | '+s,flush=True)
 def ck(name,ok):
  checks.append(dict(check=name,status='PASS' if ok else 'FAIL'))
  if not ok:raise ValueError(name)
 def save(rows,name):pd.DataFrame(rows).to_csv(out/'tables'/name,sep='\t',index=False,float_format='%.17g')
 for p in (base/'scripts').glob('*.py'):shutil.copy2(p,out/'logs'/p.name)
 shutil.copy2(base/'PROTOCOL.md',out/'PROTOCOL.md')
 xdf=pd.read_csv(inp/'clr_200_KO.tsv.gz',sep='\t',index_col=0);ids=xdf.index.to_numpy(dtype=str);nodes=xdf.columns.to_numpy(dtype=str);b,a=np.tril_indices(len(nodes),-1);m=pd.read_csv(inp/'metadata.csv',dtype={'lat_block':str}).set_index('sample_id').loc[ids].reset_index();m['profile']=m.sample_id.str.replace(r'_[123]$','',regex=True)
 xdf=pd.read_csv(inp/'clr_200_KO.tsv.gz',sep='\t',index_col=0).loc[ids,nodes];x=xdf.to_numpy();D=stats.design(m);E=m[ENV].to_numpy();E=(E-E.mean(0))/E.std(0,ddof=1);CE=np.column_stack([D,E]);Cs=[np.ones((51,1)),D,CE]
 ck('input dimensions and IDs',x.shape==(51,200) and len(set(ids))==51 and len(set(nodes))==200);ck('complete profile design',m.profile.nunique()==17 and m.groupby('profile').depth_cm.apply(lambda s:set(s)=={5,20,40}).all());ck('design ranks',np.linalg.matrix_rank(D)==5 and np.linalg.matrix_rank(CE)==9)
 save(m,'809_metadata_aligned.tsv');pd.DataFrame(CE,index=ids,columns=['Intercept','Depth20','Depth40','Cozumel','Tuxpan']+ENV).to_csv(out/'tables/809_environment_design.tsv',sep='\t',index_label='sample_id',float_format='%.17g')
 cov=[]
 for sid,c in zip(SC,Cs):
  for pred in IND:
   v=m[pred].to_numpy();r=v-c@np.linalg.lstsq(c,v,rcond=None)[0];cov.append(dict(construction=sid,predictor=pred,residual_variance_fraction=np.var(r)/np.var(v),design_rank=np.linalg.matrix_rank(c),condition_number=np.linalg.cond(c)))
 save(cov,'809_covariate_diagnostics.tsv');pd.DataFrame(np.corrcoef(E.T),index=ENV,columns=ENV).to_csv(out/'tables/809_environment_correlations.tsv',sep='\t')
 perm=stats.schedule(m,99999,1032,'permutation',True);boots=stats.schedule(m,4999,426,'bootstrap',True);netboots=stats.schedule(m,BNET,809200,'bootstrap',True);small=perm[:4999]
 ck('restricted index permutations',np.all(D[perm]==D[None,:,:]));np.savez_compressed(out/'draws/809_schedules.npz',sample_ids=ids,permutations_base0=perm,conditional_bootstrap_base0=boots,network_bootstrap_base0=netboots,refit_bootstrap_base0=boots[:BREFIT])
 allrows=[];indexnull=[];allstable=[];stabnull=[];globrows=[];globnull=[];mechrows=[];mechnull=[];introws=[];depthrows=[];arrays=[];masks=[];summ=[];allmetrics=[];envrows=[]
 for si,(sid,c) in enumerate(zip(SC,Cs)):
  log(sid+': refitting full and 51 leave-one-out networks')
  W,Y,w,r,coef,diag=lioness(x,c,a,b,nodes);ck(sid+' residual orthogonality',np.max(abs(c.T@r))<1e-8);ck(sid+' finite networks and metrics',np.isfinite(W).all() and np.isfinite(Y).all());ck(sid+' full-rank leave-one-out',diag['minimum_LOO_rank']==c.shape[1]);ck(sid+' no degenerate functions',diag['maximum_zero_variance_columns']==0)
  arrays.append((W,Y));sel0=selected(w,a,b,nodes);det=np.zeros(len(a));pos=det.copy();neg=det.copy();netdiag=[]
  for j,ii in enumerate(netboots):
   M,z,_,_=ridge_pipeline(x[ii],c[ii]);v=M[a,b];sel=selected(v,a,b,nodes);det[sel]+=1;pos[sel]+=v[sel]>0;neg[sel]+=v[sel]<0;netdiag.append(dict(construction=sid,bootstrap_id=j+1,rank=np.linalg.matrix_rank(c[ii]),zero_variance_functions=z))
  consistency=np.divide(np.maximum(pos,neg),det,out=np.zeros(len(a)),where=det>0);mask=(det/BNET>=.5)&(consistency>=.8);masks.append(mask);stable_idx=np.flatnonzero(mask)
  save(dict(node_a=nodes[a],node_b=nodes[b],full_weight=w,selection_frequency=det/BNET,sign_consistency=consistency,stable_edge=mask),'809_'+sid+'_edge_stability.tsv.gz');save(netdiag,'809_'+sid+'_network_bootstrap_diagnostics.tsv')
  np.savez_compressed(out/'objects'/('809_'+sid+'_networks.npz'),weights=W,full_weights=w,sample_metrics=Y,residual_CLR=r,nuisance_coefficients=coef,nuisance_design=c,stable_mask=mask,nodes=nodes,sample_ids=ids,node_a_index=a,node_b_index=b)
  for i,sample in enumerate(ids):allmetrics.append(dict(construction=sid,sample_id=sample,mean_abs_edge_weight=Y[i,0],natural_connectivity_abs_weighted=Y[i,1]))
  summ.append(dict(construction=sid,n_nodes=200,n_sample_edges=NEDGE,n_stable_edges=int(mask.sum()),ridge_lambda=.1,**diag));log(sid+': stable edges '+str(mask.sum())+'; conditional HI/MHI inference')
  rows,nul,boot=index_inference(m,Y,D,perm,boots)
  for rr in rows:rr['construction']=sid
  allrows+=rows;indexnull.append(nul);np.savez_compressed(out/'draws'/('809_'+sid+'_conditional_index.npz'),null_t=nul,bootstrap_beta=boot)
  if mask.sum()>0:
   ys=np.array([metric(v,stable_idx,a,b,200) for v in W]);rows,nul,bb=index_inference(m,ys,D,small,boots,['HI'])
   for rr in rows:rr['construction']=sid;rr['n_stable_edges']=int(mask.sum())
   allstable+=rows;stabnull.append(nul);np.savez_compressed(out/'objects'/('809_'+sid+'_stable_metrics.npz'),sample_ids=ids,metrics=ys)
   F,p,nf=interaction_test(m,ys[:,0],D,small);introws.append(dict(construction=sid,interaction='HI_x_depth',F=F,p_raw=p,n_permutations=4999))
   for depth in [5,20,40]:
    di=np.flatnonzero(m.depth_cm==depth);md=m.iloc[di].reset_index(drop=True);dd=np.column_stack([np.ones(len(di)),pd.get_dummies(md.locality,drop_first=True).to_numpy(float)]);sp=stats.schedule(md,4999,809300+depth,'permutation',True);sb=stats.schedule(md,4999,809400+depth,'bootstrap',True);yy=ys[di];rr,nn,bbb=index_inference(md,yy,dd,sp,sb,['HI']);r0=rr[0];r0['construction']=sid;r0['depth_cm']=depth;depthrows.append(r0)
  else:
   stabnull.append(np.full((4999,2),np.nan));introws.append(dict(construction=sid,interaction='HI_x_depth',F=np.nan,p_raw=np.nan,n_permutations=0))
  for ei,name in enumerate(ENV):
   for k,end in enumerate(ENDS):envrows.append(dict(construction=sid,environment=name,endpoint=end,spearman_rho=spearmanr(E[:,ei],Y[:,k]).statistic,inference='descriptive only; environment used upstream'))
  presence=np.zeros_like(W,dtype=bool)
  for i,v in enumerate(W):presence[i,selected(v,a,b,nodes)]=True
  gs=[kernel(abs(W)),kernel(W),gower(np.sqrt(squareform(pdist(presence,metric='jaccard'))))];gnames=['absolute_weight_euclidean','signed_weight_euclidean','retained_edge_sqrt_jaccard'];rr,nn=global_tests(gs,gnames,m,D,small)
  for v in rr:v['construction']=sid
  globrows+=rr;globnull.append(nn)
  norm=np.linalg.norm(W,axis=1);U=abs(W)/norm[:,None];rr,nn=global_tests([kernel(np.log(norm)),kernel(U)],['log_L2_magnitude','L2_normalized_pattern'],m,D,perm);beta=stats.adjusted(m[['HI']].to_numpy(),np.log(norm)[:,None],D)[0][0,0]
  for k,v in enumerate(rr):v['construction']=sid;v['standardized_beta']=beta if k==0 else np.nan
  mechrows+=rr;mechnull.append(nn)
  np.savez_compressed(out/'draws'/('809_'+sid+'_global.npz'),configuration_null_F=globnull[-1],mechanism_null_F=mechnull[-1])
 for rows,nuls,per,scope in [(allrows,indexnull,4,'eight_new_tests'),(allstable,stabnull,2,'four_new_stable_tests'),(globrows,globnull,3,'six_new_configuration_tests'),(mechrows,mechnull,2,'four_new_mechanism_tests')]:
  new=[r for r in rows if r['construction']!='unadjusted']
  if len(new)==per*2 and np.isfinite(np.column_stack(nuls[1:])).all():adjust_family(new,np.column_stack(nuls[1:]),scope)
 valid=[r for r in introws if r['construction']!='unadjusted' and np.isfinite(r['p_raw'])]
 if valid:
  for r,p in zip(valid,stats.adjust_p([v['p_raw'] for v in valid],'holm')):r['p_holm_two_new_interactions']=float(p)
 valid=[r for r in depthrows if r['construction']!='unadjusted']
 if valid:
  for r,p in zip(valid,stats.adjust_p([v['p_raw'] for v in valid],'bh')):r['q_BH_six_new_slopes']=float(p)
 save(allrows,'809_index_associations.tsv');save(allstable,'809_stable_consensus_HI.tsv');save(introws,'809_HI_depth_interactions.tsv');save(depthrows,'809_depth_HI_exploratory.tsv');save(globrows,'809_global_configuration.tsv');save(mechrows,'809_magnitude_redistribution.tsv');save(envrows,'809_environment_descriptive.tsv');save(allmetrics,'809_sample_network_metrics.tsv');save(summ,'809_network_summary.tsv')
 log('499 paired full reconstructions per construction; features and indices fixed')
 refit=np.full((BREFIT,3,4),np.nan);diagrows=[]
 for bi,ii in enumerate(boots[:BREFIT]):
  for si,(sid,c) in enumerate(zip(SC,Cs)):
   rank=np.linalg.matrix_rank(c[ii])
   if rank<c.shape[1]:diagrows.append(dict(bootstrap_id=bi+1,construction=sid,status='rank deficient',rank=rank));continue
   ww,yy,_,_,_,dg=lioness(x[ii],c[ii],a,b,nodes,ii)
   good=dg['minimum_LOO_rank']==c.shape[1] and dg['maximum_zero_variance_columns']==0 and np.isfinite(yy).all()
   if good:refit[bi,si]=stats.adjusted(m[IND].to_numpy()[ii],yy,D[ii])[0].ravel()
   diagrows.append(dict(bootstrap_id=bi+1,construction=sid,status='valid' if good else 'LOO degenerate',rank=rank,**dg))
  if (bi+1)%25==0 or bi+1==BREFIT:
   np.save(out/'draws/809_refit_bootstrap_checkpoint.npy',refit);log('Paired reconstruction '+str(bi+1)+'/'+str(BREFIT))
 rows=[]
 for si,sid in enumerate(SC):
  for j,pred in enumerate(IND):
   for k,end in enumerate(ENDS):
    v=refit[:,si,j*2+k];valid=np.isfinite(v);lo,hi=np.quantile(v[valid],[.025,.975]);r=dict(construction=sid,predictor=pred,endpoint=end,standardized_beta=allrows[si*4+j*2+k]['standardized_beta'],refit_ci_low=lo,refit_ci_high=hi,n_refit_valid=int(valid.sum()),n_refit_attempted=BREFIT,scope='residualization ridge LIONESS and sample edge selection reconstructed; features and indices fixed')
    if si>0:
     delta=v-refit[:,0,j*2+k];delta=delta[np.isfinite(delta)];dl,dh=np.quantile(delta,[.025,.975]);r.update(beta_change_vs_baseline=r['standardized_beta']-allrows[j*2+k]['standardized_beta'],paired_change_ci_low=dl,paired_change_ci_high=dh,n_paired_valid=len(delta))
    rows.append(r)
 save(rows,'809_reconstruction_bootstrap_intervals.tsv');save(diagrows,'809_reconstruction_bootstrap_diagnostics.tsv');np.savez_compressed(out/'draws/809_reconstruction_bootstrap.npz',standardized_beta=refit,construction=np.array(SC),predictors=np.array(IND),endpoints=np.array(ENDS));ck('at least 90 percent valid reconstructions each construction',np.all(np.isfinite(refit).all(2).mean(0)>.9))
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
 fig,axs=plt.subplots(2,2,figsize=(10,7),sharex=True);colors=['#777777','#0072B2','#D55E00'];labels=['Original CLR','Locality + depth','Locality + depth + environment']
 for ax,(pred,end) in zip(axs.flat,[(p,e) for p in IND for e in ENDS]):
  for si,sid in enumerate(SC):
   r=next(r for r in rows if r['construction']==sid and r['predictor']==pred and r['endpoint']==end);v=next(r for r in allrows if r['construction']==sid and r['predictor']==pred and r['endpoint']==end);ax.plot([r['refit_ci_low'],r['refit_ci_high']],[2-si,2-si],lw=1.3,color=colors[si]);ax.plot([v['ci_low'],v['ci_high']],[2-si+.09,2-si+.09],lw=3,color=colors[si],alpha=.4);ax.plot(r['standardized_beta'],2-si,'o',color=colors[si])
  ax.set_yticks([2,1,0],labels);ax.axvline(0,color='gray',ls='--');ax.set_ylim(-.4,2.4);ax.set_title(('MHI' if pred=='MHI_local' else pred)+' · '+('Mean absolute weight' if end==ENDS[0] else 'Natural connectivity'));ax.set_xlabel('Standardized slope')
 fig.tight_layout()
 for ext in ['pdf','svg','png']:fig.savefig(out/'figures'/('809_residualization_comparison.'+ext),dpi=300,facecolor='white')
 plt.close(fig);save(checks,'809_acceptance_checks.tsv')
 prov=dict(status='PASS',python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,elapsed_seconds=time.monotonic()-start,n_samples=51,n_profiles=17,n_KO=200,conditional_bootstrap=4999,reconstruction_bootstrap=BREFIT,stability_bootstrap=BNET,permutations_primary=99999,input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inp.iterdir()},protocol_sha256=hashlib.sha256((base/'PROTOCOL.md').read_bytes()).hexdigest(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),R_executed=False)
 (out/'809_provenance.json').write_text(json.dumps(prov,indent=2)+'\n');(root/'LATEST_809.txt').write_text(str(out)+'\n');log('PASS; output '+str(out))
if __name__=='__main__':main()
