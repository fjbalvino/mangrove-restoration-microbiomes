#!/usr/bin/env python3
"""Prespecified one-factor sensitivities; never select parameters by HI significance."""
from pathlib import Path
import argparse,json,hashlib,time
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist,squareform
import network_analysis as core
import network_statistics as stats

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--bundle',type=Path,required=True);ap.add_argument('--run',type=Path,required=True)
    args=ap.parse_args();out=args.run;inp=args.bundle/'inputs'
    m=pd.read_csv(out/'tables/809_metadata_aligned.tsv',sep='\t',dtype={'lat_block':str});ids=m.sample_id
    D=stats.design(m);E=m[core.ENV].to_numpy();E=(E-E.mean(0))/E.std(0,ddof=1)
    CE=np.c_[D,E];Cs=[np.ones((51,1)),D,CE]
    p=np.load(out/'draws/809_schedules.npz')['permutations_base0'][:4999]
    variants=[('ridge_005',200,.05,.05,False,False),('ridge_020',200,.20,.05,False,False),
              ('density_0025',200,.1,.025,False,False),('density_010',200,.1,.1,False,False),
              ('features_100',100,.1,.05,False,False),('features_300',300,.1,.05,False,False),
              ('four_environment_axes',200,.1,.05,False,True),('RPK_per_million',200,.1,.05,True,False)]
    primary={s:np.load(out/'objects'/f'809_{s}_networks.npz') for s in core.SC}
    primary_ids={s:primary[s]['nodes'] for s in core.SC}
    rows=[];nulls=[];conf=[];confnull=[];mechanism=[];mechnull=[];diagnostics=[]
    for variant,n,lam,density,ppm,four in variants:
        f=inp/(f'clr_{n}_KO_per_million.tsv.gz' if ppm else f'clr_{n}_KO.tsv.gz')
        frame=pd.read_csv(f,sep='\t',index_col=0).loc[ids];x=frame.to_numpy();nodes=frame.columns.to_numpy(str)
        b,a=np.tril_indices(n,-1);core.LAMBDA=lam;core.NEDGE=int(np.floor(len(a)*density+.5))
        constructions=[2] if four else [0,1,2]
        for si in constructions:
            sid=core.SC[si];c=Cs[si]
            if four:c=np.c_[D,np.delete(E,core.ENV.index('moisture_stress_PC1'),axis=1)]
            W,Y,w,resid,coef,dg=core.lioness(x,c,a,b,nodes)
            assert np.isfinite(W).all() and np.isfinite(Y).all() and dg['minimum_LOO_rank']==c.shape[1] and dg['maximum_zero_variance_columns']==0
            equal_nodes=np.array_equal(nodes,primary_ids[sid]);original=primary[sid]['weights']
            agreement=float(np.corrcoef(W.ravel(),original.ravel())[0,1]) if equal_nodes else np.nan
            identical=bool(equal_nodes and np.allclose(W,original,rtol=0,atol=1e-12))
            z=(resid-resid.mean(0))/resid.std(0,ddof=1);eig=np.linalg.eigvalsh(z.T@z/(len(z)-1))
            condition=float(((1-lam)*eig.max()+lam+core.EPS)/((1-lam)*eig.min()+lam+core.EPS))
            diagnostics.append(dict(variant=variant,construction=sid,n_KO=n,ridge_lambda=lam,density=density,
                n_retained_edges=core.NEDGE,regularized_condition_number=condition,
                same_nodes_as_primary=equal_nodes,weight_correlation_to_primary=agreement,
                weights_identical_to_primary=identical,**dg))
            beta,_,rho,t=stats.adjusted(m[core.IND].to_numpy(),Y,D)
            local=[]
            for j,pred in enumerate(core.IND):
                for k,endpoint in enumerate(core.ENDS):
                    local.append(dict(variant=variant,construction=sid,predictor=pred,endpoint=endpoint,
                        standardized_beta=beta[j,k],partial_R2=rho[j,k]**2,statistic=t[j,k],n_permutations=len(p)))
            null=stats.perm_index(m[core.IND].to_numpy(),Y,D,p)
            core.adjust_family(local,null,'within_variant_construction')
            if si>0:rows.extend(local);nulls.append(null)
            # Density changes retained-edge presence but leaves complete edge weights unchanged.
            presence=np.zeros_like(W,dtype=bool)
            for j,v in enumerate(W):presence[j,core.selected(v,a,b,nodes)]=True
            gs=[core.kernel(abs(W)),core.kernel(W),core.gower(np.sqrt(squareform(pdist(presence,'jaccard'))))]
            cr,cn=core.global_tests(gs,['absolute_weight_euclidean','signed_weight_euclidean','retained_edge_sqrt_jaccard'],m,D,p)
            norm=np.linalg.norm(W,axis=1)
            mr,mn=core.global_tests([core.kernel(np.log(norm)),core.kernel(abs(W)/norm[:,None])],['log_L2_magnitude','L2_normalized_pattern'],m,D,p)
            if si>0:
                for row in cr+mr:row.update(variant=variant,construction=sid)
                conf.extend(cr);confnull.append(cn);mechanism.extend(mr);mechnull.append(mn)
        print(variant,'complete',flush=True)
    for records,null,scope in [(rows,nulls,'all_sensitivity_index_tests'),(conf,confnull,'all_sensitivity_configuration_tests'),(mechanism,mechnull,'all_sensitivity_mechanism_tests')]:
        core.adjust_family(records,np.column_stack(null),scope)
    for data,name in [(rows,'index_associations'),(conf,'configuration'),(mechanism,'magnitude_pattern'),(diagnostics,'diagnostics')]:
        pd.DataFrame(data).to_csv(out/'tables'/f'810_parameter_{name}.tsv',sep='\t',index=False,float_format='%.17g')
    provenance={'status':'PASS','variants':[v[0] for v in variants],'synchronized_permutations':4999,
        'index_test_family':len(rows),'configuration_test_family':len(conf),'mechanism_test_family':len(mechanism),
        'parameters_selected_using_results':False,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'Exploratory one-factor sensitivity; no full bootstrap reconstructions of sensitivity variants'}
    (out/'810_sensitivity_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps(provenance,indent=2),flush=True)

if __name__=='__main__':main()
