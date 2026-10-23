"""Execute the locked CNS panel with profiles as the inferential unit."""
from pathlib import Path
import json,time,hashlib
import numpy as np
import pandas as pd
from scipy import stats
from stats_cns import Design,adjust,depth_test,profile_matrix
from calibrar_y_verificar import designs
from output_io import save_npz

ROOT=Path(__file__).resolve().parents[1]
GENES={'K02588':'nifH','K00368':'nirK','K15864':'nirS','K00376':'nosZ','K03385':'nrfA','K04561':'norB','K11180':'dsrA','K00394':'aprA','K17224':'soxB','K01601':'rbcL','K14138':'acsB'}
CYCLE={k:('C' if k in ['K01601','K14138'] else 'S' if k in ['K11180','K00394','K17224'] else 'N') for k in GENES}
STAGES=['Degraded','Early restoration','Intermediate restoration','Advanced restoration']

def bounds(r,f,patterns,keys):
    lo={};hi={}
    for k in keys:
        sizes=sorted({len(p['labels']) for p in patterns if k in p['labels'] and len(p['labels'])>1})
        delta=np.maximum(0,f.loc[k]-r.loc[k])
        if sizes:
            lo[k]=np.maximum(0,f.loc[k]-delta/(1-1/min(sizes)))
            hi[k]=np.maximum(0,f.loc[k]-delta/(1-1/max(sizes)))
        else:
            assert np.allclose(r.loc[k],f.loc[k],rtol=1e-8,atol=1e-7)
            lo[k]=hi[k]=r.loc[k]
    return pd.DataFrame(lo).T,pd.DataFrame(hi).T

def linear_bounds(weights,low,high):
    return np.maximum(weights,0)@low+np.minimum(weights,0)@high,np.maximum(weights,0)@high+np.minimum(weights,0)@low

def prepare():
    inp=ROOT/'inputs';out=ROOT/'results';spec=json.loads((ROOT/'analysis_specification.json').read_text())
    m,p,design=designs();m=m.sort_values(['profile_id','depth_cm']).reset_index(drop=True)
    assert m.locality.notna().all() and m.restoration4.notna().all()
    assert m.groupby('profile_id').locality.nunique().eq(1).all() and m.groupby('profile_id').restoration4.nunique().eq(1).all()
    assert m.groupby('profile_id').latitude.nunique().eq(1).all() and m.groupby('profile_id').longitude.nunique().eq(1).all()
    samples=m.sample_id.tolist();us=(inp/'uscg_76_kegg_min_2313.lst').read_text().split();panel=list(GENES)
    exact=spec['post_audit_quantification_groups']['exact_single_KO_positive_primary_candidates']
    mats={q.stem:pd.read_csv(q,sep='\t',index_col=0).loc[:,samples] for q in inp.glob('KO_*.tsv')}
    patterns=json.loads((inp/'annotation_patterns_selected.json').read_text());complete=[a for a in patterns if a['status']=='complete_consensus']
    keys=sorted(set(us+panel));r=mats['KO_complete_RPK'];f=mats['KO_complete_full_terms_RPK']
    lo,hi=bounds(r,f,complete,keys);den=r.loc[us].median();dlo=lo.loc[us].median();dhi=hi.loc[us].median()
    assert (dlo>0).all() and np.allclose(lo.loc[exact],hi.loc[exact],rtol=1e-8,atol=1e-7)
    num=r.loc[panel].copy();num.loc[exact]=lo.loc[exact]
    strict=[k for k in us if np.allclose(lo.loc[k],hi.loc[k],rtol=1e-8,atol=1e-7) and (lo.loc[k]>0).all()]
    assert len(strict)==56
    variants={'primary_76':num.div(den,axis=1),'fractional_76':r.loc[panel].div(den,axis=1),'strict_reference_56':num.div(lo.loc[strict].median(),axis=1)}
    rb=mats['KO_complete_both_RPK'];fb=mats['KO_complete_full_terms_both_RPK'];lb,hb=bounds(rb,fb,complete,keys)
    assert np.allclose(lb.loc[exact],hb.loc[exact],rtol=1e-8,atol=1e-7)
    nb=rb.loc[panel].copy();nb.loc[exact]=lb.loc[exact]
    variants['both_mates_76']=nb.div(rb.loc[us].median(),axis=1)
    rs=mats['KO_supported_RPK'];variants['supported_fractional_76']=rs.loc[panel].div(rs.loc[us].median(),axis=1)
    nl=lo.loc[panel].div(dhi,axis=1);nh=hi.loc[panel].div(dlo,axis=1)
    assert (nl>0).all().all() and (nh>=nl-1e-12).all().all()
    nl.to_csv(out/'single_KO_normalized_lower.tsv',sep='\t');nh.to_csv(out/'single_KO_normalized_upper.tsv',sep='\t')
    for name,values in variants.items():
        assert np.isfinite(values).all().all() and (values>0).all().all()
        values.to_csv(out/(name+'_normalized.tsv'),sep='\t')
    registry=pd.DataFrame([dict(KO=k,gene=g,cycle=CYCLE[k],quantification='resolved_single_KO' if k in exact else 'provisional_fractional',sequence_specificity='not_validated') for k,g in GENES.items()])
    registry.to_csv(out/'marker_registry.tsv',sep='\t',index=False);m.to_csv(out/'sample_manifest.tsv',sep='\t',index=False);p.to_csv(out/'profile_manifest.tsv',sep='\t',index=False)
    profile_env=[]
    for v,D in design.items():
        if v=='restoration':continue
        raw=m.groupby('profile_id',sort=True)[v].mean()
        for i,pid in enumerate(p.profile_id):profile_env.append(dict(profile_id=pid,variable=v,mean_raw=raw.loc[pid],standardized_predictor=D.X[i,-1],transformation='identity' if v=='redox_shallow_mv' else 'log2_of_profile_mean'))
    pd.DataFrame(profile_env).to_csv(out/'profile_environment.tsv',sep='\t',index=False)
    allmarkers=spec['candidate_primary_markers']+spec['support_markers']+spec['descriptive_or_pending_markers']
    detection=[]
    for k in allmarkers:
        cc=mats['KO_complete_counts'].loc[k] if k in r.index else pd.Series(0,index=samples)
        detection.append(dict(KO=k,positive_samples=int((cc>0).sum()),n_samples=51,total_fractional_fragments=float(cc.sum()),median_fractional_fragments=float(cc.median()),tested=k in panel))
    pd.DataFrame(detection).to_csv(out/'marker_detection.tsv',sep='\t',index=False)
    return spec,m,p,design,registry,variants,nl,nh

def main():
    start=time.time();out=ROOT/'results';spec,m,p,design,registry,variants,nl,nh=prepare();D=design['restoration'];panel=list(GENES)
    B=spec['bootstrap_draws'];rng=np.random.default_rng(spec['bootstrap_seed']);signs=rng.choice(np.array([-1.,1.]),size=(B,17))
    save_npz(out/'bootstrap_profile_signs.npz',signs=signs.astype('int8'),profile_id=p.profile_id.to_numpy(str))
    blocks=[np.flatnonzero(p.locality.eq(l)) for l in sorted(p.locality.unique())];perms=np.tile(np.arange(17),(19999,1))
    for row in perms:
        for block in blocks:row[block]=rng.permutation(block)
    assert (p.locality.to_numpy()[perms]==p.locality.to_numpy()).all()
    save_npz(out/'permutation_profiles.npz',indices=perms,profile_id=p.profile_id.to_numpy(str))
    # Explicit proof at sample resolution; there are 17 independent weights, not 51.
    ix={pid:i for i,pid in enumerate(p.profile_id)};sample_weights=signs[:1000,[ix[x] for x in m.profile_id]]
    for pid in p.profile_id:assert np.all(np.ptp(sample_weights[:,m.profile_id.eq(pid)],axis=1)==0)
    alltests=[];effects=[];diagnostics=[];profiles_long=[];environment=[];leaveouts=[];bounds_records=[]
    larr=np.log2(profile_matrix(nl,m,p.profile_id));harr=np.log2(profile_matrix(nh,m,p.profile_id));plow=larr.mean(axis=2);phigh=harr.mean(axis=2)
    for variant,frame in variants.items():
        arr=np.log2(profile_matrix(frame,m,p.profile_id));means=arr.mean(axis=2)
        tests=[]
        for j,k in enumerate(panel):
            result=D.wild_test(means[j],signs)
            tests.append(dict(variant=variant,KO=k,gene=GENES[k],test='restoration_omnibus',**result))
            for dd,depth in [(1,20),(2,40)]:
                z=depth_test(arr[j,:,dd]-arr[j,:,0]);tests.append(dict(variant=variant,KO=k,gene=GENES[k],test=f'depth_{depth}_vs_5',**z))
        assert len(tests)==33
        adjusted=adjust([x['p'] for x in tests])
        for a,b in zip(tests,adjusted):a['p_holm_33']=b
        alltests.extend(tests)
        for j,k in enumerate(panel):
            beta,res,cov,_=D.fit(means[j]);se=np.sqrt(np.diag(cov[0]))
            for q,stage in enumerate(STAGES):
                effects.append(dict(variant=variant,KO=k,gene=GENES[k],contrast=stage+' vs Conserved',estimate_log2=beta[0,3+q],HC3_se=se[q],n_profiles=17))
            for variable in spec['environment_candidates'][k]:
                E=design[variable];er=E.wild_test(means[j],signs);eb,eres,ecov,_=E.fit(means[j])
                ese=float(np.sqrt(ecov[0,0,0]));est=float(eb[0,-1]);crit=stats.t.ppf(.975,E.df)
                environment.append(dict(variant=variant,KO=k,gene=GENES[k],variable=variable,estimate_log2_per_SD=est,HC3_se=ese,ci_low=est-crit*ese,ci_high=est+crit*ese,ci_method='pointwise_HC3_t13',p_HC3_t_sensitivity=float(2*stats.t.sf(abs(est/ese),E.df)),**er))
            if variant!='primary_76':continue
            mse=float(np.sum(res**2)/D.df);cook=res[0]**2/(D.p*mse)*D.h/(1-D.h)**2
            diagnostics.append(dict(model='restoration',KO=k,max_cooks_distance=float(cook.max()),most_influential_profile=p.profile_id.iloc[cook.argmax()],residual_shapiro_p=stats.shapiro(res[0]).pvalue,max_leverage=D.h.max(),residual_df=D.df))
            for pi,pid in enumerate(p.profile_id):
                profiles_long.append(dict(KO=k,gene=GENES[k],cycle=CYCLE[k],profile_id=pid,locality=p.locality.iloc[pi],restoration=p.restoration4.iloc[pi],mean_log2=means[j,pi],depth_5_log2=arr[j,pi,0],depth_20_log2=arr[j,pi,1],depth_40_log2=arr[j,pi,2]))
            for variable in spec['environment_candidates'][k]:
                E=design[variable];eb,er,ec,_=E.fit(means[j]);mse=float(np.sum(er**2)/E.df);cook=er[0]**2/(E.p*mse)*E.h/(1-E.h)**2
                diagnostics.append(dict(model=variable,KO=k,max_cooks_distance=float(cook.max()),most_influential_profile=p.profile_id.iloc[cook.argmax()],residual_shapiro_p=stats.shapiro(er[0]).pvalue,max_leverage=E.h.max(),residual_df=E.df))
                low,high=linear_bounds(E.A[-1],plow[j],phigh[j]);bounds_records.append(dict(KO=k,target=variable,allocation_estimate_low=low,allocation_estimate_high=high,scope='between_profile_slope_not_confidence_interval'))
            for dd,depth in [(1,20),(2,40)]:
                low=larr[j,:,dd]-harr[j,:,0];high=harr[j,:,dd]-larr[j,:,0]
                bounds_records.append(dict(KO=k,target=f'depth_{depth}_vs_5',allocation_estimate_low=np.median(low),allocation_estimate_high=np.median(high),n_certain_positive=int((low>0).sum()),n_certain_negative=int((high<0).sum()),scope='median_paired_change_not_confidence_interval'))
            for omitted in list(p.profile_id)+list(sorted(p.locality.unique())):
                keep=~p.profile_id.eq(omitted) if omitted in set(p.profile_id) else ~p.locality.eq(omitted)
                for dd,depth in [(1,20),(2,40)]:
                    delta=arr[j,keep,dd]-arr[j,keep,0]
                    leaveouts.append(dict(KO=k,omitted=omitted,kind='profile' if omitted in set(p.profile_id) else 'locality',target=f'depth_{depth}_vs_5',estimate_log2=float(np.median(delta)),n_profiles=int(keep.sum())))
                for variable in spec['environment_candidates'][k]:
                    X=design[variable].X[keep];zero=np.all(np.abs(X)<1e-14,axis=0);X=X[:,~zero]
                    if np.linalg.matrix_rank(X)<X.shape[1]:
                        # With Carmen omitted, intercept is the sum of both locality dummies.
                        X=np.delete(X,1,axis=1)
                    b=np.linalg.lstsq(X,means[j,keep],rcond=None)[0]
                    leaveouts.append(dict(KO=k,omitted=omitted,kind='profile' if omitted in set(p.profile_id) else 'locality',target=variable,estimate_log2=b[-1],n_profiles=int(keep.sum())))
        pd.DataFrame(alltests).to_csv(out/'all_primary_and_sensitivity_tests.tsv',sep='\t',index=False)
        print(variant,'finished',round(time.time()-start,1),'seconds',flush=True)
    env=pd.DataFrame(environment)
    for variant,ixs in env.groupby('variant').groups.items():
        assert len(ixs)==19;env.loc[ixs,'q_BH_19']=adjust(env.loc[ixs,'p'],'BH')
        env.loc[ixs,'q_BH_HC3_sensitivity_19']=adjust(env.loc[ixs,'p_HC3_t_sensitivity'],'BH')
    env.to_csv(out/'environment_tests.tsv',sep='\t',index=False)
    pd.DataFrame(effects).to_csv(out/'restoration_contrasts.tsv',sep='\t',index=False)
    pd.DataFrame(diagnostics).to_csv(out/'model_diagnostics.tsv',sep='\t',index=False)
    pd.DataFrame(profiles_long).to_csv(out/'profile_marker_values.tsv',sep='\t',index=False)
    pd.DataFrame(leaveouts).to_csv(out/'leave_one_out_effects.tsv',sep='\t',index=False)
    pd.DataFrame(bounds_records).to_csv(out/'allocation_effect_bounds.tsv',sep='\t',index=False)
    tests=pd.DataFrame(alltests);tests[tests.variant.eq('primary_76')].to_csv(out/'primary_tests.tsv',sep='\t',index=False)
    # Random-intercept REML cross-check; this estimates mean depth changes, not medians.
    from fit_random_intercept import fit
    arr=np.log2(profile_matrix(variants['primary_76'],m,p.profile_id));Y=arr.reshape(11,51)
    depth=np.tile(np.array([[0,0],[1,0],[0,1]]),(17,1));X=np.column_stack([np.repeat(D.X,3,axis=0),depth]);Z=np.repeat(np.eye(17),3,axis=0)
    L=np.eye(9)[3:];sd=Y.std(axis=1,ddof=1);fitres=fit((Y-Y.mean(axis=1,keepdims=True))/sd[:,None],X,Z@Z.T/3,L)
    mixed=[]
    for j,k in enumerate(panel):
        assert np.allclose(fitres['estimate'][j,:4]*sd[j],D.fit(arr.mean(axis=2)[j])[0][0,3:],rtol=1e-6,atol=1e-7)
        for q,contrast in enumerate([s+' vs Conserved' for s in STAGES]+['mean_depth_20_vs_5','mean_depth_40_vs_5']):
            mixed.append(dict(KO=k,contrast=contrast,estimate_log2=fitres['estimate'][j,q]*sd[j],ci_low=fitres['ci_low'][j,q]*sd[j],ci_high=fitres['ci_high'][j,q]*sd[j],p_sensitivity=fitres['p'][j,q],df_satterthwaite=fitres['df'][j,q],ICC=fitres['icc'][j],singular=bool(fitres['singular'][j]),usable=bool(fitres['usable'][j])))
    pd.DataFrame(mixed).to_csv(out/'mixed_model_crosscheck.tsv',sep='\t',index=False)
    sums=dict(n_samples=51,n_profiles=17,n_localities=3,n_primary_tests=33,n_environment_tests=19,bootstrap_draws=B,independent_weights_per_draw=17,depths_kept_together=True,locality_permutation_draws=len(perms),model_df_stage=10,model_df_environment=13,landa_used=False,main_significant_tests=tests[tests.variant.eq('primary_76')&tests.p_holm_33.lt(.05)][['KO','test','p','p_holm_33']].to_dict('records'),exploratory_environment_q_below_005=env[env.variant.eq('primary_76')&env.q_BH_19.lt(.05)][['KO','variable','p','q_BH_19']].to_dict('records'),runtime_seconds=round(time.time()-start,2))
    (out/'run_summary.json').write_text(json.dumps(sums,indent=2));print(json.dumps(sums,indent=2),flush=True)

if __name__=='__main__':main()
