"""Independent high-precision calibration of the replacement interval.

Kept separate so the rejected bootstrap-t interval remains documented.
No ecological outcomes are used.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from calibrar_y_verificar import designs

def main():
    _,p,ds=designs();rng=np.random.default_rng(51239);out=[]
    for name,D in ds.items():
        for case in ['normal','heavy_t3','skew','stage_variance','locality_variance','single_profile_variance']:
            if case=='heavy_t3':y=rng.standard_t(3,(100000,17))/np.sqrt(3)
            elif case=='skew':y=(np.exp(rng.normal(size=(100000,17)))-np.exp(.5))/np.sqrt(np.exp(2)-np.exp(1))
            else:y=rng.normal(size=(100000,17))
            scales=np.ones(17)
            if case=='stage_variance':scales[p.restoration4.eq('Degraded')]=3
            if case=='locality_variance':scales[p.locality.eq('Tuxpan')]=3
            if case=='single_profile_variance':scales[np.argmax(D.h)]=10
            beta,_,cov,_=D.fit(y*scales)
            tt=beta[:,3:]/np.sqrt(np.diagonal(cov,axis1=1,axis2=2));crit=stats.t.ppf(.975,D.df)
            cover=(abs(tt)<=crit).mean(axis=0)
            out.append(dict(design=name,case=case,mincoverage=cover.min(),maxcoverage=cover.max()))
    dest=Path(__file__).resolve().parents[1]/'results/diagnostico_intervalos_HC3_t.tsv'
    pd.DataFrame(out).to_csv(dest,sep='\t',index=False)
    print(pd.DataFrame(out).to_string(index=False))

if __name__=='__main__':main()
