"""Exact Gaussian random-intercept REML with analytic Satterthwaite derivatives.

Balanced profiles have three observations. Restricted likelihood is evaluated
in an orthonormal residual space; scalar variance-ratio optimization is batched.
This evaluates the same REML likelihood and observed Hessian used by lmerTest.
No asymptotic normal substitution is used. Validation against R is required.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg, stats
import json,time

AXES=['MHI_local','vegetation_landscape_PC1','water_inundation_PC1','moisture_stress_PC1','physicochemical_PC1','nutrients_redox_PC1']
ROOT=Path(__file__).resolve().parent

def design(m,axis,interaction):
    x=m[axis].to_numpy(float);x=(x-x.mean())/x.std(ddof=1)
    d20=m.depth_cm.eq(20).to_numpy(float);d40=m.depth_cm.eq(40).to_numpy(float)
    X=np.column_stack([np.ones(len(m)),x,d20,d40,m.locality.eq('Cozumel').to_numpy(float),m.locality.eq('Tuxpan').to_numpy(float)])
    if interaction:X=np.column_stack([X,x*d20,x*d40])
    assert np.linalg.matrix_rank(X)==X.shape[1]
    z=pd.get_dummies(m.profile_id,dtype=float).to_numpy();assert np.all(z.sum(0)==3)
    return X,z@z.T/3

def fit(Y,X,B,L):
    # Y: response standardized by its across-sample SD; rows are genes.
    n,p=X.shape;df=n-p;f=len(Y)
    K=linalg.null_space(X.T);ev,E=linalg.eigh(K.T@(3*B)@K)
    ev[np.abs(ev)<1e-10]=0;ev[np.abs(ev-3)<1e-10]=3
    unique=[];which=[]
    for e in ev:
        if not unique or abs(e-unique[-1])>1e-9:unique.append(e)
        which.append(len(unique)-1)
    U=Y@(K@E)
    vals=np.array(unique);mult=np.bincount(which).astype(float)
    ss=np.stack([np.sum(U[:,np.array(which)==j]**2,axis=1) for j in range(len(vals))],axis=1)
    def objective(s):
        t=np.expm1(s)/3
        den=1+t[:,None]*vals
        Q=np.sum(ss/den,axis=1)
        return np.sum(np.log(den)*mult,axis=1)+df*np.log(Q/df)
    grid=np.linspace(0,18,121)
    fv=np.stack([objective(np.full(f,s)) for s in grid],axis=1)
    ix=np.argmin(fv,axis=1)
    assert not np.any(ix==len(grid)-1),'Variance ratio upper bound; expand explicitly'
    low=grid[np.maximum(ix-1,0)];high=grid[np.minimum(ix+1,len(grid)-1)]
    phi=(np.sqrt(5)-1)/2
    c=high-phi*(high-low);d=low+phi*(high-low);fc=objective(c);fd=objective(d)
    for _ in range(65):
        left=fc<fd
        high=np.where(left,d,high);low=np.where(left,low,c)
        c=high-phi*(high-low);d=low+phi*(high-low);fc=objective(c);fd=objective(d)
    s=(low+high)/2
    boundary=objective(np.zeros(f))<=objective(s)+1e-11
    s[boundary]=0;t=np.expm1(s)/3
    # Refine the score root; objective-only minimization loses precision near
    # zero variance. Parameterize curvature by variance ratio, not its square
    # root, to avoid a degenerate Hessian as the fitted intercept tends to zero.
    for _ in range(5):
        den=1+t[:,None]*vals
        Q=np.sum(ss/den,axis=1);u=np.sum(ss*vals/den**2,axis=1);v2=np.sum(ss*vals**2/den**3,axis=1)
        score_t=np.sum(mult*vals/den,axis=1)-df*u/Q
        curvature=-np.sum(mult*vals**2/den**2,axis=1)+df*(2*v2/Q-(u/Q)**2)
        step=np.divide(score_t,curvature,out=np.zeros_like(t),where=curvature>0)
        t=np.where(boundary,0,np.maximum(0,t-step))
    theta=np.sqrt(t);s=np.log1p(3*t)
    den=1+t[:,None]*vals
    Q=np.sum(ss/den,axis=1);sigma2=Q/df;sigma=np.sqrt(sigma2)
    q1=np.sum(-ss*vals/den**2,axis=1)
    q2=np.sum(ss*2*vals**2/den**3,axis=1)
    ld2=np.sum(-mult*vals**2/den**2,axis=1)
    H=np.zeros((f,2,2));H[:,0,0]=ld2+q2/sigma2;H[:,0,1]=H[:,1,0]=-2*q1/sigma**3;H[:,1,1]=4*df/sigma2
    evals,evecs=np.linalg.eigh(H);ok_h=(evals>1e-8).all(1)
    vinv=np.where(evals>1e-8,1/np.maximum(evals,1e-300),0)
    A=2*np.einsum('fij,fj,fkj->fik',evecs,vinv,evecs)
    singular=theta<1e-4
    A[singular]=0;A[singular,1,1]=sigma2[singular]/(2*df);ok_h[singular]=True
    chol=linalg.cholesky(X.T@X,lower=True);invl=linalg.solve_triangular(chol,np.eye(p),lower=True)
    r,V=linalg.eigh(invl@X.T@B@X@invl.T)
    r[np.abs(r)<1e-10]=0;r[np.abs(r-1)<1e-10]=1
    S=invl.T@V;T=X@S
    alpha=3*t/(1+3*t);idg=1/(1-alpha[:,None]*r)
    score=(Y@T-alpha[:,None]*(Y@(B@T)))*idg
    beta=score@S.T
    W=L@S
    estimate=beta@L.T
    base=np.einsum('lj,fj,mj->flm',W,idg,W)
    cov=sigma2[:,None,None]*base
    dt=3/(1+3*t)**2
    J0=sigma2[:,None,None]*np.einsum('lj,fj,mj->flm',W,dt[:,None]*r*idg**2,W)
    J1=2*sigma[:,None,None]*base
    v=np.diagonal(cov,axis1=1,axis2=2);g=np.stack([np.diagonal(J0,axis1=1,axis2=2),np.diagonal(J1,axis1=1,axis2=2)],axis=2)
    var_v=np.einsum('fli,fij,flj->fl',g,A,g)
    nu=2*v*v/var_v;se=np.sqrt(v);tt=estimate/se
    pv=2*stats.t.sf(np.abs(tt),nu);qcrit=stats.t.ppf(.975,nu)
    usable=(~singular)&ok_h&np.isfinite(nu).all(1)&(nu>0).all(1)
    ret=dict(estimate=estimate,se=se,df=nu,p=pv,ci_low=estimate-qcrit*se,ci_high=estimate+qcrit*se,theta=theta,sigma=sigma,icc=t/(1+t),singular=singular,hessian_ok=ok_h,usable=usable,beta=beta,restricted_deviance=objective(s))
    if len(L)==3:
        D=np.array([[-1,1,0],[-1,0,1]])
        C=np.einsum('ij,fjk,lk->fil',D,cov,D);b=estimate@D.T
        vv,PP=np.linalg.eigh(C)
        jt=[np.einsum('ij,fjk,lk->fil',D,J,D) for J in [J0,J1]]
        gg=np.stack([np.einsum('fji,fjk,fki->fi',PP,J,PP) for J in jt],axis=2)
        nnu=2*vv*vv/np.einsum('fli,fij,flj->fl',gg,A,gg)
        weighted=np.einsum('fji,fj->fi',PP,b)
        F=np.sum(weighted**2/vv,axis=1)/2
        EE=np.sum(nnu/(nnu-2),axis=1)
        ddf=np.where((nnu<=2).any(1),2,np.where(np.abs(nnu[:,1]-nnu[:,0])<1e-8,nnu.mean(1),2*EE/(EE-2)))
        ret.update(interaction_F=F,interaction_df=ddf,interaction_p=stats.f.sf(F,2,ddf))
    return ret

def run():
    out=ROOT/'models';out.mkdir(exist_ok=True)
    m=pd.read_csv(ROOT/'prepared/metadata.csv');raw=np.load(ROOT/'prepared/selected_CLR.npy');sd=raw.std(axis=1,ddof=1);Y=(raw-raw.mean(axis=1,keepdims=True))/sd[:,None]
    np.save(out/'response_sd.npy',sd)
    for axis in AXES:
        for kind in ['common','interaction']:
            dest=out/(axis+'_'+kind+'.npz')
            if dest.exists():continue
            X,B=design(m,axis,kind=='interaction')
            if kind=='common':L=np.zeros((1,6));L[0,1]=1
            else:
                L=np.zeros((3,8));L[:,1]=1;L[1,6]=1;L[2,7]=1
            pieces=[];start=time.time()
            for begin in range(0,len(Y),10000):
                pieces.append(fit(Y[begin:begin+10000],X,B,L))
            result={k:np.concatenate([r[k] for r in pieces]) for k in pieces[0]}
            np.savez_compressed(dest,**result)
            print(axis,kind,'seconds',round(time.time()-start,1),'usable',int(result['usable'].sum()),'singular',int(result['singular'].sum()),flush=True)
    (out/'STATUS.txt').write_text('FITTED_VALIDATION_PENDING\n')

if __name__=='__main__':run()
