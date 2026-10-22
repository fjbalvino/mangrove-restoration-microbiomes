#!/usr/bin/env python3
"""Numerical feasibility screen only; not SPIEC-EASI or sampling-stability inference."""
from pathlib import Path
import argparse,json,time,warnings
import numpy as np
import pandas as pd
from sklearn.covariance import graphical_lasso

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',type=Path,required=True);args=ap.parse_args();rows=[];conditioning=[]
    for scenario in ['unadjusted','space_depth','space_depth_environment']:
        x=np.load(args.run/'objects'/f'809_{scenario}_networks.npz')['residual_CLR'];z=(x-x.mean(0))/x.std(0,ddof=1);S=z.T@z/(len(z)-1)
        correlations=np.corrcoef(x.T)[np.triu_indices(len(S),1)]
        conditioning.append(dict(construction=scenario,residual_rank=int(np.linalg.matrix_rank(x)),
            pairs_abs_r_over_0999=int((abs(correlations)>.999).sum()),pairs_abs_r_over_099=int((abs(correlations)>.99).sum()),
            maximum_absolute_correlation=float(abs(correlations).max())))
        for alpha in [.5,.2,.1,.05]:
            for control in ['default','strict']:
                start=time.monotonic();row=dict(construction=scenario,alpha=alpha,control=control,status='failed')
                with warnings.catch_warnings(record=True) as seen:
                    warnings.simplefilter('always')
                    try:
                        cov,precision,costs,niter=graphical_lasso(S,alpha,max_iter=200 if control=='default' else 1000,
                            tol=1e-4 if control=='default' else 1e-5,enet_tol=1e-4 if control=='default' else 1e-7,
                            eps=np.finfo(float).eps if control=='default' else 1e-6,return_costs=True,return_n_iter=True)
                        gap=float(costs[-1][1]);nedge=int((abs(precision[np.triu_indices(len(S),1)])>1e-8).sum())
                        eig=np.linalg.eigvalsh(precision);valid=np.isfinite(precision).all() and eig.min()>0 and abs(gap)<(1e-4 if control=='default' else 1e-5)
                        row.update(status='converged' if valid else 'not converged',dual_gap=gap,iterations=niter,n_edges=nedge,precision_condition=float(eig.max()/eig.min()))
                    except (FloatingPointError,ValueError,np.linalg.LinAlgError) as exc:row['error']=str(exc)
                row.update(warning_count=len(seen),warnings=' | '.join(sorted(set(str(w.message) for w in seen)))[:1000],elapsed_seconds=time.monotonic()-start)
                rows.append(row);pd.DataFrame(rows).to_csv(args.run/'tables/811_graphical_lasso_screen.tsv',sep='\t',index=False)
                print(scenario,alpha,control,row['status'],round(row['elapsed_seconds'],2),flush=True)
    pd.DataFrame(conditioning).to_csv(args.run/'tables/811_input_conditioning.tsv',sep='\t',index=False)

if __name__=='__main__':main()
