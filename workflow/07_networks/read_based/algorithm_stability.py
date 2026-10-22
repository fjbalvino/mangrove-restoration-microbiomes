#!/usr/bin/env python3
"""Paired algorithm stability benchmark; no HI/MHI association tests or method tuning."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
import pandas as pd
import sklearn
from scipy.linalg import cho_factor,cho_solve
from scipy.stats import spearmanr
from sklearn.covariance import ledoit_wolf,oas
import network_analysis as core

METHODS=['ridge_010','ledoit_wolf','oas']

def estimate(x,c,method,a,b):
    assert np.linalg.matrix_rank(c)==c.shape[1]
    resid=x-c@np.linalg.lstsq(c,x,rcond=None)[0];sd=resid.std(0,ddof=1)
    assert (sd>np.sqrt(np.finfo(float).eps)).all()
    z=(resid-resid.mean(0))/sd
    lam=.1 if method=='ridge_010' else (ledoit_wolf(z,assume_centered=True)[1] if method=='ledoit_wolf' else oas(z,assume_centered=True)[1])
    S=z.T@z/(len(z)-1);cov=(1-lam)*S+(lam+core.EPS)*np.eye(x.shape[1])
    precision=cho_solve(cho_factor(cov,lower=True),np.eye(x.shape[1]))
    pc=-precision/np.sqrt(np.outer(precision.diagonal(),precision.diagonal()))
    assert np.isfinite(pc).all() and 0<=lam<=1
    return pc[a,b],float(lam)

def lioness(x,c,method,a,b,nodes):
    n=len(x);full,lam=estimate(x,c,method,a,b);W=np.empty((n,len(a)));lams=[]
    for i in range(n):
        keep=np.arange(n)!=i;w,l=estimate(x[keep],c[keep],method,a,b)
        W[i]=n*full-(n-1)*w;lams.append(l)
    Y=np.array([core.metric(w,core.selected(w,a,b,nodes),a,b,len(nodes)) for w in W])
    return full,W,Y,lam,lams

def mask(w,a,b,nodes):
    out=np.zeros(len(w),bool);out[core.selected(w,a,b,nodes)]=True;return out

def compare(v,w,a,b,nodes):
    p=mask(v,a,b,nodes);q=mask(w,a,b,nodes);shared=p&q;nshared=shared.sum()
    return dict(edge_jaccard=float(nshared/(p|q).sum()),
        chance_adjusted_overlap=float((nshared/core.NEDGE-core.NEDGE/len(v))/(1-core.NEDGE/len(v))),
        weight_correlation=float(np.corrcoef(v,w)[0,1]),
        shared_edge_sign_agreement=float((np.sign(v[shared])==np.sign(w[shared])).mean()) if nshared else np.nan)

def compare_samples(W,V,Y,Z,a,b,nodes):
    rows=[compare(v,w,a,b,nodes) for v,w in zip(W,V)]
    return dict(sample_edge_jaccard_median=float(np.median([r['edge_jaccard'] for r in rows])),
        sample_edge_jaccard_min=float(min(r['edge_jaccard'] for r in rows)),
        sample_weight_correlation=float(np.corrcoef(W.ravel(),V.ravel())[0,1]),
        mean_weight_rank_correlation=float(spearmanr(Y[:,0],Z[:,0]).statistic),
        natural_connectivity_rank_correlation=float(spearmanr(Y[:,1],Z[:,1]).statistic))

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--bundle',type=Path,required=True);ap.add_argument('--run',type=Path,required=True)
    args=ap.parse_args();r=args.run;out=r/'algorithm_comparison';out.mkdir(exist_ok=False);(out/'objects').mkdir()
    started=time.monotonic();meta=pd.read_csv(r/'tables/809_metadata_aligned.tsv',sep='\t');ids=meta.sample_id
    frame=pd.read_csv(args.bundle/'inputs/clr_200_KO.tsv.gz',sep='\t',index_col=0).loc[ids]
    nodes=frame.columns.to_numpy(str);x=frame.to_numpy();ppm=pd.read_csv(args.bundle/'inputs/clr_200_KO_per_million.tsv.gz',sep='\t',index_col=0).loc[ids,nodes].to_numpy()
    reference={s:np.load(r/'objects'/f'809_{s}_networks.npz') for s in core.SC}
    a=reference[core.SC[0]]['node_a_index'];b=reference[core.SC[0]]['node_b_index']
    boots=np.load(r/'draws/809_schedules.npz')['network_bootstrap_base0']
    bootrows=[];deletion=[];normalization=[];summaries=[];checks=[]
    for method in METHODS:
        for sid in core.SC:
            c=reference[sid]['nuisance_design'];full,W,Y,lam,lams=lioness(x,c,method,a,b,nodes)
            if method=='ridge_010':
                error=float(np.max(abs(W-reference[sid]['weights'])));assert error<1e-9
                checks.append(dict(check=sid+' ridge reproduction',max_error=error,status='PASS'))
            det=np.zeros(len(a));pos=det.copy();neg=det.copy();blams=[]
            for j,ii in enumerate(boots):
                w,l=estimate(x[ii],c[ii],method,a,b);sel=mask(w,a,b,nodes)
                det+=sel;pos+=sel&(w>0);neg+=sel&(w<0);blams.append(l)
                bootrows.append(dict(method=method,construction=sid,bootstrap_id=j+1,shrinkage=l,**compare(full,w,a,b,nodes)))
            consistency=np.divide(np.maximum(pos,neg),det,out=np.zeros(len(a)),where=det>0)
            stable=(det/len(boots)>=.5)&(consistency>=.8)
            if method=='ridge_010':assert np.array_equal(stable,reference[sid]['stable_mask'])
            summaries.append(dict(method=method,construction=sid,n_KO=len(nodes),n_edges=core.NEDGE,
                full_shrinkage=lam,LOO_shrinkage_min=min(lams),LOO_shrinkage_max=max(lams),
                bootstrap_shrinkage_min=min(blams),bootstrap_shrinkage_median=float(np.median(blams)),bootstrap_shrinkage_max=max(blams),n_stable_edges=int(stable.sum())))
            np.savez_compressed(out/'objects'/f'{method}_{sid}.npz',full_weights=full,weights=W,metrics=Y,
                selection_frequency=det/len(boots),sign_consistency=consistency,stable_mask=stable,nodes=nodes,sample_ids=ids.to_numpy(str))
            print(method,sid,'200 aggregate bootstraps complete',flush=True)
            for profile in meta.profile.unique():
                ii=np.flatnonzero(meta.profile.ne(profile));f,w,y,l,ll=lioness(x[ii],c[ii],method,a,b,nodes)
                deletion.append(dict(method=method,construction=sid,omitted_profile=profile,
                    full_shrinkage=l,**compare(full,f,a,b,nodes),**compare_samples(W[ii],w,Y[ii],y,a,b,nodes)))
            f,w,y,l,ll=lioness(ppm,c,method,a,b,nodes)
            normalization.append(dict(method=method,construction=sid,full_shrinkage=l,
                **compare(full,f,a,b,nodes),**compare_samples(W,w,Y,y,a,b,nodes)))
            print(method,sid,'17 complete-profile deletions and normalization check complete',flush=True)
    B=pd.DataFrame(bootrows);D=pd.DataFrame(deletion);N=pd.DataFrame(normalization);S=pd.DataFrame(summaries)
    bs=B.groupby(['method','construction']).agg(bootstrap_jaccard_median=('edge_jaccard','median'),
        bootstrap_jaccard_q025=('edge_jaccard',lambda x:x.quantile(.025)),bootstrap_jaccard_q975=('edge_jaccard',lambda x:x.quantile(.975)),
        bootstrap_adjusted_overlap_median=('chance_adjusted_overlap','median'),bootstrap_weight_correlation_median=('weight_correlation','median'),
        bootstrap_sign_agreement_median=('shared_edge_sign_agreement','median')).reset_index()
    ds=D.groupby(['method','construction']).agg(deletion_aggregate_jaccard_median=('edge_jaccard','median'),
        deletion_sample_jaccard_median=('sample_edge_jaccard_median','median'),deletion_sample_jaccard_min=('sample_edge_jaccard_median','min'),
        deletion_mean_rank_median=('mean_weight_rank_correlation','median'),deletion_mean_rank_min=('mean_weight_rank_correlation','min'),
        deletion_natural_rank_median=('natural_connectivity_rank_correlation','median')).reset_index()
    ns=N[['method','construction','sample_edge_jaccard_median','edge_jaccard','mean_weight_rank_correlation']].rename(columns={'sample_edge_jaccard_median':'normalization_sample_jaccard','edge_jaccard':'normalization_aggregate_jaccard','mean_weight_rank_correlation':'normalization_mean_rank'})
    S=S.merge(bs).merge(ds).merge(ns);decisions=[]
    for method in METHODS[1:]:
        for sid in core.SC[1:]:
            new=S[(S.method==method)&(S.construction==sid)].iloc[0];old=S[(S.method==METHODS[0])&(S.construction==sid)].iloc[0]
            delta={k:float(new[k]-old[k]) for k in ['bootstrap_jaccard_median','deletion_sample_jaccard_median','normalization_sample_jaccard','deletion_mean_rank_median']}
            passes=delta['bootstrap_jaccard_median']>=.02 and all(delta[k]>=-.02 for k in delta if k!='bootstrap_jaccard_median')
            decisions.append(dict(method=method,construction=sid,passes_practical_criterion=passes,**{'change_'+k:v for k,v in delta.items()}))
    for name,t in [('bootstrap',B),('profile_deletion',D),('normalization',N),('summary',S),('decision',pd.DataFrame(decisions)),('checks',pd.DataFrame(checks))]:
        t.to_csv(out/f'811_algorithm_{name}.tsv',sep='\t',index=False,float_format='%.17g')
    report=dict(status='PASS',methods=METHODS,n_KO=200,n_profiles=17,n_bootstrap=200,n_deletions=17,
        selection_or_evaluation_uses_HI=False,sklearn=sklearn.__version__,elapsed_seconds=time.monotonic()-started,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),protocol_sha256=hashlib.sha256(Path(__file__).with_name('ALGORITHM_BENCHMARK.md').read_bytes()).hexdigest(),
        consistently_more_stable=[method for method in METHODS[1:] if all(z['passes_practical_criterion'] for z in decisions if z['method']==method)])
    (out/'811_provenance.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
