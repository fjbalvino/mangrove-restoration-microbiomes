#!/usr/bin/env python3
"""Numerical helpers retained from legacy803; original inference functions unchanged."""
import argparse,hashlib,json,platform,shutil,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata
PRIMARY='KEGG_ko__clr_variance_N200'
END=['natural_connectivity_abs_weighted','mean_abs_edge_weight']
IND=['MHI_local','HI']
ENV=['vegetation_landscape_pc1','water_inundation_pc1','physicochemical_pc1','nutrients_redox_pc1']
B=4999;HIGH=99999

def digest(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def design(m):
    return np.column_stack([np.ones(len(m)),(m.depth_cm.to_numpy()==20).astype(float),(m.depth_cm.to_numpy()==40).astype(float),pd.get_dummies(m.locality,drop_first=True).to_numpy(float)])
def residual(a,d):
    q=np.linalg.qr(d,mode='reduced')[0]
    return a-q@(q.T@a)
def adjusted(X,Y,D):
    xr=residual(X,D);yr=residual(Y,D)
    cross=xr.T@yr;ssx=(xr*xr).sum(0);ssy=(yr*yr).sum(0)
    raw=cross/ssx[:,None]
    beta=raw*X.std(0,ddof=1)[:,None]/Y.std(0,ddof=1)[None,:]
    rho=cross/np.sqrt(ssx[:,None]*ssy[None,:]);df=len(X)-np.linalg.matrix_rank(D)-1
    t=rho*np.sqrt(df/(1-rho**2))
    return beta,raw*X.std(0,ddof=1)[:,None],rho,t

def correlation(X,Y,method):
    if method=='spearman':X=rankdata(X,axis=0);Y=rankdata(Y,axis=0)
    x=X-X.mean(0);y=Y-Y.mean(0)
    return (x.T@y)/np.sqrt((x*x).sum(0)[:,None]*(y*y).sum(0)[None,:])

def adjust_p(p,method):
    p=np.asarray(p,float);order=np.argsort(p);n=len(p);sp=p[order]
    a=np.maximum.accumulate(sp*np.arange(n,0,-1)) if method=='holm' else np.minimum.accumulate((sp*n/np.arange(1,n+1))[::-1])[::-1]
    out=np.empty(n);out[order]=np.minimum(a,1);return out

def pvalue(obs,null):return (1+(np.abs(null)>=np.abs(obs)).sum(0))/(len(null)+1)

def schedule(meta,n,seed,kind,within):
    blocks=sorted(meta.lat_block.unique());br=[np.flatnonzero(meta.lat_block.to_numpy()==b) for b in blocks]
    br=np.array([v[np.argsort(meta.depth_cm.iloc[v])] for v in br])
    loc=np.array([meta.locality.iloc[v[0]] for v in br]);groups=[np.flatnonzero(loc==v) for v in sorted(set(loc))] if within else [np.arange(len(blocks))]
    rng=np.random.default_rng(seed)
    if kind=='permutation':
        idx=np.tile(np.arange(len(meta)),(n,1))
        for g in groups:
            source=rng.permuted(np.broadcast_to(g,(n,len(g))).copy(),axis=1)
            idx[:,br[g].ravel()]=br[source].reshape(n,-1)
    else:
        draws=np.concatenate([rng.choice(g,size=(n,len(g)),replace=True) for g in groups],axis=1)
        idx=br[draws].reshape(n,-1)
    return idx.astype(np.uint8)

def perm_index(X,Y,D,idx):
    q=np.linalg.qr(D,mode='reduced')[0];yr=residual(Y,D);ssy=(yr*yr).sum(0);df=len(X)-np.linalg.matrix_rank(D)-1
    stats=np.empty((len(idx),X.shape[1]*Y.shape[1]))
    for start in range(0,len(idx),5000):
        ii=idx[start:start+5000]
        for j in range(X.shape[1]):
            xp=X[ii,j].T;xr=xp-q@(q.T@xp)
            r=(yr.T@xr).T/np.sqrt((xr*xr).sum(0)[:,None]*ssy[None,:])
            stats[start:start+len(ii),j*Y.shape[1]:(j+1)*Y.shape[1]]=r*np.sqrt(df/(1-r*r))
    return stats

def perm_env(X,Y,idx,method):
    if method=='spearman':X=rankdata(X,axis=0);Y=rankdata(Y,axis=0)
    yc=Y-Y.mean(0);ssy=(yc*yc).sum(0);out=np.empty((len(idx),X.shape[1]*Y.shape[1]))
    for j in range(X.shape[1]):
        xp=X[idx,j].T;xp-=xp.mean(0)
        r=(yc.T@xp).T/np.sqrt((xp*xp).sum(0)[:,None]*ssy[None,:])
        out[:,j*Y.shape[1]:(j+1)*Y.shape[1]]=r
    return out
