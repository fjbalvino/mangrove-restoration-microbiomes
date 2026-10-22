#!/usr/bin/env python3
"""Independent QR and direct sklearn covariance audit of alternative LIONESS weights."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from sklearn.covariance import ledoit_wolf,oas

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--bundle',type=Path,required=True);args=ap.parse_args();r=args.run;rows=[]
    meta=pd.read_csv(r/'tables/809_metadata_aligned.tsv',sep='\t');frame=pd.read_csv(args.bundle/'inputs/clr_200_KO.tsv.gz',sep='\t',index_col=0).loc[meta.sample_id];x=frame.to_numpy()
    for method,estimator in [('ledoit_wolf',ledoit_wolf),('oas',oas)]:
        for scenario in ['unadjusted','space_depth','space_depth_environment']:
            original=np.load(r/'objects'/f'809_{scenario}_networks.npz');C=original['nuisance_design'];a=original['node_a_index'];b=original['node_b_index']
            result=np.load(r/'algorithm_comparison/objects'/f'{method}_{scenario}.npz')
            def fit(ii):
                Q=np.linalg.qr(C[ii],mode='reduced')[0];res=x[ii]-Q@(Q.T@x[ii]);z=(res-res.mean(0))/res.std(0,ddof=1)
                cov,shrinkage=estimator(z,assume_centered=True);cov*=len(ii)/(len(ii)-1);cov+=1e-8*np.eye(200)
                P=np.linalg.inv(cov);return (-P/np.sqrt(np.outer(P.diagonal(),P.diagonal())))[a,b]
            full=fit(np.arange(51));W=np.array([51*full-50*fit(np.flatnonzero(np.arange(51)!=i)) for i in range(51)])
            error=float(np.max(abs(W-result['weights'])));assert error<1e-9
            rows.append(dict(method=method,construction=scenario,check='all 1014900 weights via QR and sklearn covariance',max_error=error,status='PASS'))
    d=pd.read_csv(r/'algorithm_comparison/811_algorithm_decision.tsv',sep='\t')
    for _,row in d.iterrows():
        good=row.change_bootstrap_jaccard_median>=.02 and all(row[c]>=-.02 for c in ['change_deletion_sample_jaccard_median','change_normalization_sample_jaccard','change_deletion_mean_rank_median'])
        assert good==row.passes_practical_criterion
    rows.append(dict(check='prespecified stability decision thresholds',status='PASS'))
    pd.DataFrame(rows).to_csv(r/'algorithm_comparison/811_independent_checks.tsv',sep='\t',index=False)
    print('PASS',len(rows),'checks including all 6,089,400 alternative LIONESS weights')

if __name__=='__main__':main()
