#!/usr/bin/env python3
"""Analysis 611. Exact weighted ordinal trend with locality blocks and LOLO.

Python implementation v1.0.0. No random permutation sampling, asymptotic P,
feature selection by core, or use of a held-out locality for training direction.
"""
from pathlib import Path
from itertools import permutations, combinations
from functools import lru_cache
import argparse, gzip, hashlib, io, json, math, platform, sys, traceback
import numpy as np
import pandas as pd

VERSION='1.0.0'
STAGES=['Degraded','Early restoration','Intermediate restoration','Advanced restoration','Conserved']
LOCS=['Carmen','Cozumel','Tuxpan']; DEPTHS=[5,20,40]
VARIANTS=[('main_CLR',5e-7,True),('CLR_0.1ppm',1e-7,True),('CLR_1ppm',1e-6,True),('relative_abundance',0.,True),('without_conserved',5e-7,False)]

def write_table(df,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    raw=df.to_csv(None,sep='\t',index=False,na_rep='NA',float_format='%.15g').encode('utf-8')
    data=gzip.compress(raw,mtime=0) if path.suffix=='.gz' else raw
    tmp=path.with_name(path.name+'.tmp');tmp.write_bytes(data)
    assert tmp.read_bytes()==data;tmp.replace(path)

def adjusted(p,method='BH'):
    p=np.asarray(p,float);ok=np.isfinite(p);v=np.where(ok,p,1.)
    order=np.argsort(v,kind='stable');m=len(v);q=np.minimum.accumulate((v[order]*m/np.arange(1,m+1))[::-1])[::-1]
    if method=='BY':q*=np.sum(1/np.arange(1,m+1))
    result=np.empty(m);result[order]=np.minimum(q,1);result[~ok]=np.nan
    return result

def signed(v): return np.sign(np.asarray(v)).astype(int)

@lru_cache(None)
def design(counts):
    k=len(counts);n=sum(counts);scale=k*(k-1)
    assert all(v>0 for v in counts)
    assert max(counts)<=2 and sum(v==2 for v in counts)<=1,'Unsupported block counts; revise exact lattice'
    labels=tuple(g for g,c in enumerate(counts) for _ in range(c))
    perms=np.array(sorted(set(permutations(labels))),dtype=int)
    ii,jj=np.triu_indices(n,1)
    weights=np.zeros((len(perms),len(ii)),dtype=np.int16)
    for t,(i,j) in enumerate(zip(ii,jj)):
        a=perms[:,i];b=perms[:,j];den=np.array(counts)[a]*np.array(counts)[b]
        weights[:,t]=np.where(a==b,0,np.sign(b-a)*2//den)
    return scale,ii,jj,weights

def tie_key(y):
    s=np.sort(y);breaks=np.r_[0,np.flatnonzero(np.diff(s)!=0)+1,len(s)]
    return tuple(np.diff(breaks).astype(int))

@lru_cache(None)
def null_hist(counts,ties):
    scale,ii,jj,w=design(counts)
    canonical=np.repeat(np.arange(len(ties)),ties)
    score=w @ np.sign(canonical[jj]-canonical[ii]).astype(np.int16)
    hist=np.bincount(score+scale,minlength=2*scale+1).astype(np.int64)
    assert hist.sum()==len(w)
    return hist

@lru_cache(None)
def exact_lookup(descriptors):
    hist=np.array([1],dtype=np.int64);bound=0
    for counts,ties in descriptors:
        hist=np.convolve(hist,null_hist(counts,ties));bound+=len(counts)*(len(counts)-1)
    values=np.arange(-bound,bound+1);total=int(hist.sum())
    expected=math.prod(math.factorial(sum(c))//math.prod(math.factorial(x) for x in c) for c,t in descriptors)
    assert total==expected
    tail=np.array([hist[np.abs(values)>=a].sum()/total for a in range(bound+1)])
    return tail,total

def local_evidence(z,stage):
    k=int(stage.max())+1;counts=tuple(int(np.sum(stage==g)) for g in range(k))
    scale,ii,jj,w=design(counts)
    a=stage[ii];b=stage[jj];den=np.array(counts)[a]*np.array(counts)[b]
    weights=np.where(a==b,0,np.sign(b-a)*2//den).astype(np.int16)
    score=weights @ np.sign(z[jj]-z[ii]).astype(np.int16)
    descriptors=[(counts,tie_key(z[:,f])) for f in range(z.shape[1])]
    return score.astype(int),descriptors,scale

def exact_p(scores,descriptors,indices):
    qsum=sum(scores[j] for j in indices)
    p=np.empty(len(qsum));nperm=np.empty(len(qsum),int)
    for f,obs in enumerate(qsum):
        key=tuple(descriptors[j][f] for j in indices)
        tail,n=exact_lookup(key);p[f]=tail[abs(int(obs))];nperm[f]=n
    return p,nperm

def pava(y):
    # Equal weights for stage means. Used only to summarize fitted shape.
    values=[];weights=[]
    for val in y:
        values.append(float(val));weights.append(1)
        while len(values)>1 and values[-2]>values[-1]:
            ww=weights[-2]+weights[-1];vv=(values[-2]*weights[-2]+values[-1]*weights[-1])/ww
            values[-2:]=[vv];weights[-2:]=[ww]
    return np.repeat(values,weights)

def transform(x,delta):
    p=x/x.sum(1,keepdims=True)
    if delta==0:return p,p
    z=np.log2(p+delta);z-=z.mean(1,keepdims=True)
    return p,z

def analyze_slice(meta,x,features,ontology,depth,variant,delta,with_conserved):
    p,z=transform(x,delta)
    # Deterministic numerical equality for cross-language rank comparisons.
    # Proportions use 15 decimals; CLR uses 10 decimals (not detection thresholds).
    zr=np.round(z,15 if delta==0 else 10)
    use=meta.depth_cm.eq(depth).to_numpy()
    if not with_conserved:use &= meta.restoration4.ne('Conserved').to_numpy()
    m=meta.loc[use].reset_index(drop=True);xx=x[use];pp=p[use];zz=z[use];rankz=zr[use]
    stage=np.array([STAGES.index(v) for v in m.restoration4]);k=5 if with_conserved else 4;scale=k*(k-1)
    positive=xx>0;npos=positive.sum(0);eligible=(npos>=3)&(np.ptp(rankz,axis=0)>0)
    scores=[];descriptors=[];local_ok=[];local_pos=[]
    means=[];means_rel=[]
    for loc in LOCS:
        ix=m.locality.eq(loc).to_numpy()
        qq,dd,ss=local_evidence(rankz[ix],stage[ix]);assert ss==scale
        scores.append(qq);descriptors.append(dd);local_pos.append(positive[ix].sum(0))
        local_ok.append((positive[ix].sum(0)>=3)&(np.ptp(rankz[ix],axis=0)>0))
        means.append(np.array([zz[ix & (stage==g)].mean(0) for g in range(k)]))
        means_rel.append(np.array([pp[ix & (stage==g)].mean(0) for g in range(k)]))
    pv,nperm=exact_p(scores,descriptors,[0,1,2]);pv[~eligible]=np.nan
    S=sum(scores)/(3*scale);direction=signed(S)
    mean=np.mean(means,axis=0);mean_rel=np.mean(means_rel,axis=0)
    adj=direction[None,:]*np.diff(mean,axis=0)
    monotone=(adj.min(0)>=-1e-10)&(direction!=0)&eligible
    main=pd.DataFrame(dict(ontology=ontology,feature_id=features,depth_cm=depth,variant=variant,n_profiles=len(m),n_positive=npos,zero_fraction=1-npos/len(m),eligible=eligible,S=S,direction=direction,p_exact=pv,n_assignments=nperm,monotone_observed=monotone,max_adjacent_reversal=np.maximum(0,-adj.min(0)),endpoint_difference=mean[-1]-mean[0],endpoint_relative_difference=mean_rel[-1]-mean_rel[0]))
    for j,loc in enumerate(LOCS):
        main['S_'+loc]=scores[j]/scale;main['n_positive_'+loc]=local_pos[j];main['eligible_'+loc]=local_ok[j]
    if variant!='main_CLR':return main,None,None
    ref=meta.restoration4.eq('Conserved').to_numpy();det=p>=.0001
    main['core_reference_9of9']=det[ref].all(0)
    main['core_depth_3of3']=det[ref & meta.depth_cm.eq(depth).to_numpy()].all(0)
    localdir=np.array([signed(q) for q in scores]);main['local_direction_concordance_n']=np.sum((localdir==direction)&(direction!=0)&np.array(local_ok),axis=0)
    main['opposite_local_direction']=np.any(localdir*direction<0,axis=0)
    fit=np.array([pava(mean[:,f]*direction[f])*direction[f] if direction[f] else np.repeat(mean[:,f].mean(),k) for f in range(len(features))]).T
    main['isotonic_rmse']=np.sqrt(np.mean((mean-fit)**2,axis=0))
    stages=[]
    for g in range(k):
        stages.append(pd.DataFrame(dict(ontology=ontology,feature_id=features,depth_cm=depth,stage_order=g,restoration4=STAGES[g],mean_CLR=mean[g],isotonic_CLR=fit[g],mean_relative=mean_rel[g])))
        main['mean_CLR_'+str(g)]=mean[g]
        if g:main['adjacent_difference_'+str(g)]=mean[g]-mean[g-1]
    folds=[]
    for omitted,loc in enumerate(LOCS):
        ix=m.locality.ne(loc).to_numpy();train=[j for j in range(3) if j!=omitted]
        ptrain,ntrainperm=exact_p(scores,descriptors,train)
        Strain=sum(scores[j] for j in train)/(2*scale)
        etrain=(positive[ix].sum(0)>=3)&(np.ptp(rankz[ix],axis=0)>0);ptrain[~etrain]=np.nan
        train_dir=signed(Strain);held=scores[omitted]/scale
        reftrain=ref & meta.locality.ne(loc).to_numpy()
        folds.append(pd.DataFrame(dict(ontology=ontology,feature_id=features,depth_cm=depth,omitted_locality=loc,n_training=int(ix.sum()),n_heldout=int((~ix).sum()),n_positive_training=positive[ix].sum(0),n_positive_heldout=local_pos[omitted],eligible_training=etrain,eligible_heldout=local_ok[omitted],S_training=Strain,S_heldout=held,training_direction=train_dir,p_training=ptrain,n_training_assignments=ntrainperm,training_matches_full=(train_dir==direction)&(direction!=0)&etrain,heldout_matches_training=(signed(held)==train_dir)&(train_dir!=0)&local_ok[omitted]&etrain,core_training_6of6=det[reftrain].all(0),core_training_depth_2of2=det[reftrain & meta.depth_cm.eq(depth).to_numpy()].all(0))))
    return main,pd.concat(folds,ignore_index=True),pd.concat(stages,ignore_index=True)

def audit_inputs(inputs):
    m=pd.read_csv(inputs/'metadata.tsv',sep='\t')
    assert not m.sample_id.duplicated().any() and len(m)==51
    assert m.profile_id.nunique()==17
    assert set(m.depth_cm)==set(DEPTHS) and set(m.locality)==set(LOCS)
    assert m.groupby('profile_id').depth_cm.apply(lambda v:sorted(v)==DEPTHS).all()
    assert (m.groupby('profile_id').restoration4.nunique()==1).all()
    assert (m.groupby('profile_id').locality.nunique()==1).all()
    for dep in DEPTHS:
        tab=pd.crosstab(m[m.depth_cm==dep].locality,m[m.depth_cm==dep].restoration4).reindex(index=LOCS,columns=STAGES)
        assert tab.values.tolist()==[[1,1,2,1,1],[1,1,2,1,1],[1,1,1,1,1]]
    return m

def run(inputs,out):
    inputs=Path(inputs);out=Path(out)
    if out.exists():raise ValueError('Output already exists; use a new directory')
    out.mkdir(parents=True);(out/'STATUS.txt').write_text('RUNNING\n')
    try:
        meta=audit_inputs(inputs);write_table(meta,out/'metadata_aligned.tsv')
        provenance=[];all_main=[];all_sens=[];all_folds=[];all_stages=[]
        for ontology in ['KO','PFAM']:
            file=inputs/(ontology+'.tsv.gz');df=pd.read_csv(file,sep='\t',index_col=0)
            assert not df.index.duplicated().any() and not df.columns.duplicated().any()
            assert set(df.index)==set(meta.sample_id)
            df=df.loc[meta.sample_id];x=df.to_numpy(float);features=df.columns.to_numpy()
            assert np.isfinite(x).all() and (x>=0).all() and (x.sum(1)>0).all()
            provenance.append(dict(file=file.name,sha256=hashlib.sha256(file.read_bytes()).hexdigest(),n_samples=len(df),n_features=len(features)))
            ont_m=[];ont_s=[];ont_f=[];ont_g=[]
            for variant,delta,inc in VARIANTS:
                for dep in DEPTHS:
                    print(ontology,variant,dep,flush=True)
                    tab,fold,stage=analyze_slice(meta,x,features,ontology,dep,variant,delta,inc)
                    if variant=='main_CLR':ont_m.append(tab);ont_f.append(fold);ont_g.append(stage)
                    else:ont_s.append(tab)
            main=pd.concat(ont_m,ignore_index=True);sens=pd.concat(ont_s,ignore_index=True);folds=pd.concat(ont_f,ignore_index=True)
            main['q_BH']=adjusted(main.p_exact);main['q_BY']=adjusted(main.p_exact,'BY')
            main['significant_BH']=main.q_BH.le(.05);main['significant_BY']=main.q_BY.le(.05)
            for v,ix in sens.groupby('variant',sort=False).groups.items():
                sens.loc[ix,'q_BH']=adjusted(sens.loc[ix,'p_exact']);sens.loc[ix,'q_BY']=adjusted(sens.loc[ix,'p_exact'],'BY')
            for loc,ix in folds.groupby('omitted_locality',sort=False).groups.items():
                folds.loc[ix,'q_training_BH']=adjusted(folds.loc[ix,'p_training']);folds.loc[ix,'q_training_BY']=adjusted(folds.loc[ix,'p_training'],'BY')
            folds['selected_training']=folds.q_training_BH.le(.05)
            folds['selected_and_transferred']=folds.selected_training & folds.heldout_matches_training
            keys=['ontology','feature_id','depth_cm']
            agg=folds.groupby(keys,sort=False).agg(lolo_sign_matches_n=('training_matches_full','sum'),heldout_matches_training_n=('heldout_matches_training','sum'),training_selected_n=('selected_training','sum'),training_selected_transferred_n=('selected_and_transferred','sum'),min_S_training=('S_training','min'),max_S_training=('S_training','max')).reset_index()
            main=main.merge(agg,on=keys,validate='one_to_one',sort=False)
            main['direction_stable']=main.lolo_sign_matches_n.eq(3)&main.heldout_matches_training_n.eq(3)&main.local_direction_concordance_n.eq(3)
            main['stable_monotone_candidate']=main.significant_BH&main.monotone_observed&main.direction_stable
            for var in ['CLR_0.1ppm','CLR_1ppm','relative_abundance','without_conserved']:
                sub=sens[sens.variant==var][keys+['S','p_exact','q_BH','eligible','monotone_observed']].rename(columns={c:var+'__'+c for c in ['S','p_exact','q_BH','eligible','monotone_observed']})
                main=main.merge(sub,on=keys,validate='one_to_one',sort=False)
            main['zero_sign_stable']=signed(main['CLR_0.1ppm__S']).__eq__(main.direction)&signed(main['CLR_1ppm__S']).__eq__(main.direction)&main['CLR_0.1ppm__eligible']&main['CLR_1ppm__eligible']&main.eligible
            main['zero_selection_stable']=main['CLR_0.1ppm__q_BH'].le(.05).eq(main.significant_BH)&main['CLR_1ppm__q_BH'].le(.05).eq(main.significant_BH)
            all_main.append(main);all_sens.append(sens);all_folds.append(folds);all_stages.append(pd.concat(ont_g,ignore_index=True))
        main=pd.concat(all_main,ignore_index=True);sens=pd.concat(all_sens,ignore_index=True);folds=pd.concat(all_folds,ignore_index=True);stages=pd.concat(all_stages,ignore_index=True)
        write_table(main,out/'611_main_results.tsv.gz');write_table(sens,out/'611_sensitivity_results.tsv.gz');write_table(folds,out/'611_LOLO_results.tsv.gz');write_table(stages,out/'611_stage_means.tsv.gz')
        selected=main[main.significant_BH].sort_values(['ontology','depth_cm','q_BH','feature_id'])
        write_table(selected,out/'611_selected_functions.tsv')
        summary=main.groupby(['ontology','depth_cm'],sort=False).agg(n_features=('feature_id','size'),n_eligible=('eligible','sum'),n_BH=('significant_BH','sum'),n_BY=('significant_BY','sum'),n_monotone_BH=('monotone_observed',lambda x:int((x & main.loc[x.index,'significant_BH']).sum())),n_stable_monotone=('stable_monotone_candidate','sum'),n_core_reference=('core_reference_9of9','sum'),min_p=('p_exact','min'),min_q=('q_BH','min')).reset_index()
        write_table(summary,out/'611_summary.tsv');write_table(pd.DataFrame(provenance),out/'input_manifest.tsv')
        params=dict(version=VERSION,engine='Python',engine_version=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,core_relative_floor=.0001,core_reference_n=9,min_positive=3,rank_round_CLR=10,rank_round_relative=15,delta_main=5e-7,delta_sensitivity=[1e-7,1e-6],family='ontology x all functions x 3 depths; untestable P=1 internally',direction_test='two-sided absolute S',exact_null='local label distributions convolved with integer counts',profile_independence='Assumed conditional on canonical profile IDs; intervention-site replication not independently verified',cache_distributions=exact_lookup.cache_info().currsize)
        (out/'parameters.json').write_text(json.dumps(params,indent=2));(out/'STATUS.txt').write_text('COMPLETE\n');print(summary.to_string(index=False),flush=True)
    except Exception:
        (out/'STATUS.txt').write_text('FAILED\n'+traceback.format_exc());raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',required=True);p.add_argument('--out',required=True)
    args=p.parse_args();run(args.inputs,args.out)
