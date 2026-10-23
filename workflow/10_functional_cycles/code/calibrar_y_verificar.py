"""Check inference and profile grouping before viewing empirical associations."""
from pathlib import Path
import json,time
import numpy as np
import pandas as pd
from scipy import stats
from stats_cns import Design,depth_test,profile_matrix,adjust
ROOT=Path(__file__).resolve().parents[1]

def designs():
    m=pd.read_csv(ROOT/'inputs/metadata_003b.csv');m['profile_id']=m.sample_id.str.rsplit('_',n=1).str[0]
    p=m.groupby('profile_id',sort=True).first().reset_index()
    assert len(p)==17 and p.lat_block.nunique()==17
    loc=np.column_stack([np.ones(17),p.locality.eq('Cozumel'),p.locality.eq('Tuxpan')]).astype(float)
    stages=['Degraded','Early restoration','Intermediate restoration','Advanced restoration']
    X=np.column_stack([loc,*[p.restoration4.eq(s) for s in stages]]).astype(float)
    ds={'restoration':Design(X)}
    for v in ['redox_shallow_mv','n_no3_umol_l','n_nh4plus_umol_l','s_2_umol']:
        z=m.groupby('profile_id',sort=True)[v].mean().to_numpy()
        if v!='redox_shallow_mv':assert (z>0).all();z=np.log2(z)
        z=(z-z.mean())/z.std(ddof=1);ds[v]=Design(np.column_stack([loc,z]))
    return m,p,ds

def main(nsim=2000,B=999):
    start=time.time();out=ROOT/'results';m,p,ds=designs();n=len(p);rng=np.random.default_rng(2026102202)
    signs=rng.choice([-1.,1.],size=(B,n));checks=[]
    # Expansion to the original depths must repeat the same profile weight.
    ix={pid:i for i,pid in enumerate(p.profile_id)}
    sample_signs=signs[:,[ix[pid] for pid in m.profile_id]]
    for pid in p.profile_id:assert (np.ptp(sample_signs[:,m.profile_id.eq(pid)],axis=1)==0).all()
    checks.append(dict(check='one_sign_per_profile_shared_by_three_depths',status='PASS',n_profiles=17,n_draws=B))
    # Duplicated profile-depth records must be rejected, not create replication.
    fake=pd.DataFrame(np.ones((1,len(m))),index=['test'],columns=m.sample_id)
    assert profile_matrix(fake,m,p.profile_id).shape==(1,17,3)
    try:profile_matrix(fake,pd.concat([m,m.iloc[[0]]]),p.profile_id)
    except AssertionError:checks.append(dict(check='duplicate_profile_depth_rejected',status='PASS'))
    else:raise AssertionError('Pseudoreplication guard failed')
    # Independently compute HC3 covariance and the explicit null fitted bootstrap.
    for name,D in ds.items():
        y=rng.normal(size=n);beta,res,cov,stat=D.fit(y)
        ia=np.linalg.inv(D.X.T@D.X)
        direct=ia@D.X.T@np.diag((res[0]/(1-D.h))**2)@D.X@ia
        assert np.allclose(cov[0],direct[3:,3:],atol=1e-10)
        nullfit=y@D.H0;rr=(y-nullfit)/np.sqrt(1-D.h0)
        explicit=D.fit(nullfit+signs[:25]*rr)[3];shortcut=D.fit(signs[:25]*rr)[3]
        assert np.allclose(explicit,shortcut,atol=1e-9,rtol=1e-10)
        shift=D.X[:,:3]@rng.normal(size=3)
        assert D.wild_test(y,signs)['p']==D.wild_test(y+shift,signs)['p']
    checks.append(dict(check='HC3_independent_formula_and_null_fitted_invariance',status='PASS',n_designs=len(ds)))
    assert np.allclose(adjust([.01,.04,.03]),[.03,.06,.06])
    assert np.allclose(adjust([.01,.04,.03],'BH'),[.03,.04,.04])
    checks.append(dict(check='multiplicity_known_examples',status='PASS'))
    pd.DataFrame(checks).to_csv(out/'verificacion_pseudoreplicacion.tsv',sep='\t',index=False)
    cases=['normal','heavy_t3','skew','stage_variance','locality_variance','single_profile_variance']
    records=[]
    for name,D in ds.items():
        for case in cases:
            if case=='heavy_t3':Y=rng.standard_t(3,size=(nsim,n))/np.sqrt(3)
            elif case=='skew':Y=(np.exp(rng.normal(size=(nsim,n)))-np.exp(.5))/np.sqrt(np.exp(2)-np.exp(1))
            else:Y=rng.normal(size=(nsim,n))
            scale=np.ones(n)
            if case=='stage_variance':scale[p.restoration4.eq('Degraded')]=3
            if case=='locality_variance':scale[p.locality.eq('Tuxpan')]=3
            if case=='single_profile_variance':scale[np.argmax(D.h)]=10
            Y*=scale
            beta,res,cov,T=D.fit(Y);se=np.sqrt(np.diagonal(cov,axis1=1,axis2=2))
            pv=[];covers=[]
            for j in range(0,nsim,20):
                y=Y[j:j+20];nrun=len(y)
                rr=(y@D.R0)/np.sqrt(1-D.h0)
                bb=(rr[:,None,:]*signs[None,:,:]).reshape(nrun*B,n)
                tw=D.fit(bb)[3].reshape(nrun,B)
                pv.extend((1+(tw>=T[j:j+nrun,None]-1e-12).sum(axis=1))/(B+1))
                # Pointwise unrestricted studentized intervals, calibrated separately.
                u=(y@D.R)/np.sqrt(1-D.h)
                yy=(u[:,None,:]*signs[None,:,:]).reshape(nrun*B,n)
                b,_,v,_=D.fit(yy)
                tb=(b[:,3:]/np.sqrt(np.diagonal(v,axis1=1,axis2=2))).reshape(nrun,B,D.q)
                qs=np.quantile(tb,[.025,.975],axis=1)
                lo=beta[j:j+nrun,3:]-qs[1]*se[j:j+nrun]
                hi=beta[j:j+nrun,3:]-qs[0]*se[j:j+nrun]
                covers.extend((lo<=0)&(hi>=0))
            rate=float((np.asarray(pv)<=.05).mean());cover=np.asarray(covers).mean(axis=0)
            records.append(dict(design=name,scenario=case,simulations=nsim,draws=B,type_I_error=rate,pointwise_CI_coverage_min=float(cover.min()),pointwise_CI_coverage_max=float(cover.max()),max_leverage=float(D.h.max())))
            pd.DataFrame(records).to_csv(out/'calibracion_diseno_y_intervalos.tsv',sep='\t',index=False)
        print(name,'completed',round(time.time()-start,1),'seconds',flush=True)
    (out/'calibracion_estado.json').write_text(json.dumps(dict(n_profiles=n,n_simulations_per_scenario=nsim,n_designs=len(ds),n_scenarios=len(cases),draws=B,ecological_response_used=False),indent=2))
    print(pd.DataFrame(records).to_string(index=False),flush=True)

if __name__=='__main__':main()
