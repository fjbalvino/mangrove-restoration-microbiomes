"""Profile-level inference. Bootstrap draws never become biological replicates."""
import numpy as np
from scipy import stats

def adjust(p,method='holm'):
    p=np.asarray(p,float);n=len(p);order=np.argsort(p);v=p[order]
    if method=='holm':v=np.maximum.accumulate(v*(n-np.arange(n)))
    elif method=='BH':v=np.minimum.accumulate((v*n/(np.arange(n)+1))[::-1])[::-1]
    else:raise ValueError(method)
    out=np.empty(n);out[order]=np.minimum(v,1);return out

class Design:
    def __init__(self,X,nuisance=3):
        self.X=np.asarray(X,float);self.n,self.p=self.X.shape;self.nuisance=nuisance
        assert np.linalg.matrix_rank(self.X)==self.p
        self.df=self.n-self.p;self.q=self.p-nuisance
        assert self.df>0 and self.q>0
        self.A=np.linalg.solve(self.X.T@self.X,self.X.T)
        self.H=self.X@self.A;self.R=np.eye(self.n)-self.H;self.h=np.diag(self.H)
        X0=self.X[:,:nuisance];A0=np.linalg.solve(X0.T@X0,X0.T)
        self.H0=X0@A0;self.R0=np.eye(self.n)-self.H0;self.h0=np.diag(self.H0)
        self.C=self.A[nuisance:]
        assert self.h.max()<1 and self.h0.max()<1

    def fit(self,Y):
        Y=np.atleast_2d(Y);beta=Y@self.A.T;res=Y@self.R
        C=self.C;u=res/(1-self.h)
        cov=np.einsum('ki,mi,li->mkl',C,u*u,C,optimize=True)
        b=beta[:,self.nuisance:]
        if self.q==1:stat=b[:,0]**2/cov[:,0,0]
        else:stat=np.einsum('mi,mi->m',b,np.linalg.solve(cov,b[...,None])[...,0])/self.q
        return beta,res,cov,np.maximum(stat,0)

    def wild_test(self,y,signs,batch=4096):
        y=np.asarray(y);obs=self.fit(y)[3][0]
        assert signs.shape[1]==self.n
        residual=(y@self.R0)/np.sqrt(1-self.h0)
        exceed=0
        for start in range(0,len(signs),batch):
            yy=signs[start:start+batch]*residual
            exceed+=int((self.fit(yy)[3]>=obs-1e-12).sum())
        return dict(statistic=float(obs),p=(exceed+1)/(len(signs)+1),exceedances=exceed,bootstrap_draws=len(signs),n_profiles=self.n,residual_df=self.df)

    def wild_ci(self,y,signs,batch=4096):
        beta,res,cov,_=self.fit(y);b=beta[0,self.nuisance:];se=np.sqrt(np.diag(cov[0]))
        u=res[0]/np.sqrt(1-self.h);ts=[]
        for start in range(0,len(signs),batch):
            bb,_,vv,_=self.fit(signs[start:start+batch]*u)
            ts.append(bb[:,self.nuisance:]/np.sqrt(np.diagonal(vv,axis1=1,axis2=2)))
        quant=np.quantile(np.concatenate(ts),[.025,.975],axis=0)
        return b,se,b-quant[1]*se,b-quant[0]*se

    def permutation_test(self,y,perms):
        Q=np.linalg.qr(self.X)[0];y=np.asarray(y)
        null=np.sum((y@self.H0)**2);total=y@y
        full=np.sum((y@Q)**2);obs=max(0,(full-null)/self.q)/((total-full)/self.df)
        yp=y[perms];fullb=np.sum((yp@Q)**2,axis=1)
        fb=np.maximum(0,(fullb-null)/self.q)/np.maximum((total-fullb)/self.df,1e-15)
        return (int((fb>=obs-1e-12).sum())+1)/(len(perms)+1)

def depth_test(delta):
    x=np.asarray(delta,float);n=len(x);zero=np.isclose(x,0,rtol=0,atol=1e-12);nz=x[~zero]
    p=stats.binomtest(int((nz>0).sum()),len(nz),.5).pvalue if len(nz) else 1.
    k=max([k for k in range(1,n//2+1) if 2*stats.binom.cdf(k-1,n,.5)<=.05],default=0)
    sx=np.sort(x);lo=sx[k-1] if k else -np.inf;hi=sx[n-k] if k else np.inf
    mean=x.mean();se=x.std(ddof=1)/np.sqrt(n);crit=stats.t.ppf(.975,n-1)
    return dict(n_profiles=n,n_nonzero=len(nz),n_positive=int((nz>0).sum()),n_zero=int(zero.sum()),median_log2=float(np.median(x)),ci_low=float(lo),ci_high=float(hi),interval_coverage=float(1-2*stats.binom.cdf(k-1,n,.5)),p=float(p),mean_log2=float(mean),mean_ci_low=float(mean-crit*se),mean_ci_high=float(mean+crit*se))

def profile_matrix(values,meta,profiles):
    """Validate exactly one observation per profile-depth and preserve pairing."""
    out=[]
    assert meta.sample_id.is_unique and not meta.duplicated(['profile_id','depth_cm']).any()
    for pid in profiles:
        block=meta.loc[meta.profile_id.eq(pid)].sort_values('depth_cm')
        assert block.depth_cm.tolist()==[5,20,40]
        out.append(values.loc[:,block.sample_id].to_numpy())
    # markers x profiles x depth
    return np.stack(out,axis=1)
