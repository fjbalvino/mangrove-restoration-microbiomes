#!/usr/bin/env python3
"""Validate archived Figure 1 against unchanged counts and current metadata.

The original 301 permutation P values are retained, not simulated again.
Deterministic B/C quantities and the historical D fractions are independently
recomputed. D is then recomputed with the explicitly supplied current ECI.
No taxa are removed, relabelled, or reselected.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
import warnings
import numpy as np
import pandas as pd
import rdata
from scipy import stats

STAGES = ['Degraded','Early restoration','Intermediate restoration',
          'Advanced restoration','Conserved']
COMMUNITIES = ['01_full_species_community_all51',
               '03_raw_nonunimodal_candidates_1031_available']

def atomic(path, payload):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as f:
        temporary = f.name; f.write(payload); f.flush(); os.fsync(f.fileno())
    os.replace(temporary, path)
    assert hashlib.sha256(path.read_bytes()).digest() == hashlib.sha256(payload).digest()

def csv(path, frame):
    atomic(path, frame.to_csv(index=False, float_format='%.15g').encode())

def find(root, name):
    matches = list(Path(root).rglob(name))
    if len(matches) != 1:
        raise ValueError(f'Expected exactly one {name}; found {len(matches)}')
    return matches[0]

def dummy(series):
    return pd.get_dummies(series, drop_first=True, dtype=float).to_numpy()

def design(md, fields):
    arrays = [np.ones((len(md), 1))]
    arrays += [md[[s]].to_numpy(float) if s == 'ECI_PC1' else dummy(md[s]) for s in fields]
    return np.column_stack(arrays)

def fitted_projector(x):
    return x @ np.linalg.pinv(x)

def model_stats(y, full, reduced):
    hf, hr = fitted_projector(full), fitted_projector(reduced)
    gram = y @ y.T
    ss = float(np.sum((hf-hr) * gram))
    res = float(np.sum((np.eye(len(y))-hf) * gram))
    total = float(np.sum((y-y.mean(axis=0))**2))
    df1 = np.linalg.matrix_rank(full)-np.linalg.matrix_rank(reduced)
    df2 = len(y)-np.linalg.matrix_rank(full)
    return dict(R2=ss/total, pseudo_F=(ss/df1)/(res/df2), df_num=int(df1),df_den=int(df2))

def adjusted_r2(y, x):
    resid = y - fitted_projector(x) @ y
    ssr = np.sum(resid**2); sst = np.sum((y-y.mean(axis=0))**2)
    return 1-(ssr/sst)*(len(y)-1)/(len(y)-np.linalg.matrix_rank(x))

def unique_fractions(y, md):
    all_fields = ['restoration4','ECI_PC1','locality']
    total_adj = adjusted_r2(y, design(md, all_fields))
    return {f: total_adj-adjusted_r2(y,design(md,[x for x in all_fields if x != f])) for f in all_fields}

def clr(counts):
    z = np.log1p(counts)
    return z-z.mean(axis=1, keepdims=True)

def geometric_median(a):
    center = a.mean(axis=0)
    for _ in range(20000):
        dist = np.linalg.norm(a-center,axis=1)
        if (dist < 1e-13).any():
            center += 1e-10
            continue
        weights = 1/dist
        new = (a*weights[:,None]).sum(axis=0)/weights.sum()
        if np.linalg.norm(new-center) < 1e-11:
            return new
        center = new
    raise ValueError('Spatial median did not converge')

def tax_frame(parsed):
    a = parsed.object.attributes
    while a is not None and a.value is not None:
        tag = rdata.conversion.convert(a.tag) if a.tag is not None else None
        if tag == 'tax_table':
            attrs = rdata.conversion.convert(a.value[0].attributes)
            values = rdata.conversion.convert(a.value[0])
            return pd.DataFrame(np.array(values).reshape(attrs['dim'],order='F'),
                                index=attrs['dimnames'][0],columns=attrs['dimnames'][1])
        a = a.value[1]
    raise ValueError('No taxonomy slot')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--phyloseq',required=True)
    ap.add_argument('--metadata',required=True)
    ap.add_argument('--archive',required=True)
    ap.add_argument('--raw-metrics',required=True)
    ap.add_argument('--membership',required=True)
    ap.add_argument('--output',required=True)
    args=ap.parse_args(); out=Path(args.output); checks=[]
    used={}
    def read(name):
        p=find(args.archive,name);used[name]=p
        return pd.read_csv(p)
    def equal(label,a,b,atol=1e-7,rtol=1e-8):
        a,b=np.asarray(a,float),np.asarray(b,float)
        ok=bool(np.allclose(a,b,atol=atol,rtol=rtol,equal_nan=True))
        checks.append(dict(check=label,passed=ok,max_abs_difference=float(np.nanmax(np.abs(a-b)))))
        if not ok: raise ValueError(f'Validation failed: {label}; {checks[-1]}')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        parsed=rdata.parser.parse_file(args.phyloseq)
        ps=rdata.conversion.convert(parsed)
        tax=tax_frame(parsed)
    counts=pd.DataFrame(ps.otu_table.values.T,index=ps.otu_table.coords['dim_1'].values,
                        columns=ps.otu_table.coords['dim_0'].values)
    assert counts.shape == (51,12761) and counts.index.is_unique and counts.columns.is_unique
    assert np.isfinite(counts.to_numpy()).all() and (counts.to_numpy()>=0).all()
    assert np.equal(counts.to_numpy(),np.round(counts.to_numpy())).all()
    original=read('301_input_md5.csv').set_index('input')
    md5=hashlib.md5(Path(args.phyloseq).read_bytes()).hexdigest()
    assert md5 == original.loc['ps_rds','md5']
    for role,path in [('raw_metrics_csv',args.raw_metrics),('membership_csv',args.membership)]:
        assert hashlib.md5(Path(path).read_bytes()).hexdigest()==original.loc[role,'md5'],role
        used[role]=Path(path)
    meta=pd.read_csv(args.metadata).set_index('sample_id',drop=False)
    assert meta.index.is_unique and set(meta.index)==set(counts.index)
    meta=meta.loc[counts.index].copy()
    meta['profile_id']=meta['locality']+'::'+meta['lat_block'].astype(str)
    assert not meta[['restoration4','depth_cm','locality','lat_block','ECI_PC1']].isna().any().any()
    assert set(meta.restoration4)==set(STAGES)
    profiles=meta.groupby('profile_id',sort=True)
    assert len(profiles)==17
    assert all(set(g.depth_cm)=={5,20,40} and len(g)==3 and g.restoration4.nunique()==1 and g.locality.nunique()==1 for _,g in profiles)
    pcoa=read('301_PCoA_source_data.csv')
    pcoa=pcoa[pcoa.community.isin(COMMUNITIES)&pcoa.metric.eq('Aitchison_primary')].copy()
    for community,g in pcoa.groupby('community'):
        assert set(g.sample_id)==set(meta.index)
        g=g.set_index('sample_id').loc[meta.index]
        for col in ['locality','restoration4','depth_cm','profile_id']:
            assert np.array_equal(g[col].astype(str).to_numpy(),meta[col].astype(str).to_numpy()),col
    current_eci=meta.ECI_PC1.to_numpy(float)
    old_eci=read('302_auditoria_ECI_PC1.csv').set_index('sample_id').loc[meta.index,'ECI_PC1']
    csv(out/'environment_versions.csv',pd.DataFrame({'sample_id':meta.index,'historical_ECI_PC1':old_eci.values,'current_ECI_PC1':current_eci}))
    raw=read('202_counts_raw1031_muestras_por_taxa.csv').set_index('sample_id').loc[meta.index]
    assert raw.shape==(51,1031)
    screen=pd.read_csv(args.raw_metrics,dtype={'Original_Taxon_ID':str})
    selected=screen.loc[screen.observed_joint_raw.astype(str).str.lower().eq('true'),'Original_Taxon_ID']
    assert selected.is_unique and set(selected)==set(raw.columns)
    equal('Archived candidate counts equal unchanged phyloseq subset',counts.loc[:,raw.columns],raw,atol=0,rtol=0)
    abundance=counts.to_numpy(); probability=abundance/abundance.sum(axis=1,keepdims=True)
    logs=np.zeros_like(probability);np.log(probability,out=logs,where=probability>0)
    hill=pd.DataFrame({'sample_id':meta.index,'profile_id':meta.profile_id,
                       'q0_richness':(abundance>0).sum(axis=1),
                       'q1_exp_shannon':np.exp(-np.sum(probability*logs,axis=1)),
                       'q2_inverse_simpson':1/np.sum(probability**2,axis=1)})
    metrics=['q0_richness','q1_exp_shannon','q2_inverse_simpson']
    hp=hill.groupby('profile_id')[metrics].mean()
    old_hp=read('301_diversidad_Hill_completa_por_perfil.csv').set_index('profile_id').loc[hp.index]
    equal('All 51 profile-mean Hill values',hp[metrics],old_hp[metrics])
    pm=meta.groupby('profile_id').first().loc[hp.index]
    ht=read('301_pruebas_Hill_comunidad_completa.csv')
    ht=ht[ht.analysis_role.eq('primary_observed_Hill_profile_mean')].set_index('response')
    em=read('301_Hill_estimaciones_ajustadas_IC95_por_etapa.csv')
    x=design(pm,['locality','restoration4']);xr=design(pm,['locality'])
    for metric in metrics:
        y=hp[[metric]].to_numpy();fit=model_stats(y,x,xr)
        equal(f'{metric}: stage F',fit['pseudo_F'],ht.loc[metric,'observed_F'])
        coef=np.linalg.lstsq(x,y,rcond=None)[0].ravel()
        variance=float(np.sum((y-x@coef[:,None])**2)/(len(y)-np.linalg.matrix_rank(x)))
        covariance=variance*np.linalg.inv(x.T@x)
        # Reuse the fitted design's category levels when forming each EMM.
        loc_levels=sorted(pm.locality.unique());stage_levels=sorted(pm.restoration4.unique())
        for stage in STAGES:
            rows=[]
            for loc in loc_levels:
                rows.append([1]+[int(loc==v) for v in loc_levels[1:]]+[int(stage==v) for v in stage_levels[1:]])
            contrast=np.mean(rows,axis=0);estimate=contrast@coef
            se=np.sqrt(contrast@covariance@contrast);half=stats.t.ppf(.975,10)*se
            old=em[(em.response==metric)&(em.restoration4==stage)].iloc[0]
            equal(f'{metric} {stage}: adjusted mean and CI',[estimate,estimate-half,estimate+half],old[['estimate','conf_low','conf_high']])
    result=read('301_resultados_restauracion_orden_manuscrito.csv')
    result=result[result.community.isin(COMMUNITIES)&result.metric.eq('Aitchison_primary')&result.model_role.eq('primary_adjusted_locality_depth')]
    dispersion=read('301_dispersion_restauracion.csv')
    dispersion=dispersion[dispersion.community.isin(COMMUNITIES)&dispersion.metric.eq('Aitchison_primary')]
    clrs={COMMUNITIES[0]:clr(abundance),COMMUNITIES[1]:clr(raw.to_numpy())}
    disp_checks=[]
    for community,y in clrs.items():
        fit=model_stats(y,design(meta,['locality','depth_cm','restoration4']),design(meta,['locality','depth_cm']))
        old=result[result.community==community].iloc[0]
        equal(f'{community}: marginal R2 and F',[fit['R2'],fit['pseudo_F']],old[['R2','pseudo_F']])
        yc=y-y.mean(axis=0);eig,vectors=np.linalg.eigh(yc@yc.T)
        idx=np.argsort(eig)[::-1];eig=eig[idx];vectors=vectors[:,idx]
        points=vectors[:,:2]*np.sqrt(eig[:2])
        old_pcoa=pcoa[pcoa.community==community].set_index('sample_id').loc[meta.index]
        for k,col in enumerate(['PCoA1','PCoA2']):
            sign=np.sign(np.dot(points[:,k],old_pcoa[col]))
            equal(f'{community}: {col} coordinates',points[:,k]*sign,old_pcoa[col])
        equal(f'{community}: PCoA explained fractions',eig[:2]/eig[eig>0].sum(),old_pcoa.iloc[0][['PCoA1_variance_positive','PCoA2_variance_positive']])
        coordinates=vectors[:,eig>1e-8]*np.sqrt(eig[eig>1e-8])
        old_disp=dispersion[dispersion.community==community].iloc[0]
        for center_kind in ['median','centroid']:
            distances=np.zeros(len(meta))
            for stage in STAGES:
                ii=np.flatnonzero(meta.restoration4.to_numpy()==stage);a=coordinates[ii]
                center=geometric_median(a) if center_kind=='median' else a.mean(axis=0)
                distances[ii]=np.linalg.norm(a-center,axis=1)*np.sqrt(len(ii)/(len(ii)-1))
            f=model_stats(distances[:,None],design(meta,['restoration4']),design(meta,[]))['pseudo_F']
            disp_checks.append(dict(community=community,center=center_kind,recomputed_F=f,archived_F=old_disp.pseudo_F))
            # vegan and this independent Weiszfeld solver use different stopping
            # tolerances for the spatial median; verify F to four decimal places.
            if center_kind=='median':equal(f'{community}: spatial-median dispersion F',f,old_disp.pseudo_F,atol=1e-4)
    vp=read('302_particion_varianza_por_profundidad.csv')
    names={'restoration4':'restoration_pure','ECI_PC1':'ECI_PC1_pure','locality':'locality_pure'}
    rows=[]
    for depth in [5,20,40]:
        ii=np.flatnonzero(meta.depth_cm.to_numpy()==depth);md=meta.iloc[ii].copy();y=clrs[COMMUNITIES[0]][ii]
        for version in ['historical','current']:
            md['ECI_PC1']=old_eci.iloc[ii].to_numpy() if version=='historical' else current_eci[ii]
            fractions=unique_fractions(y,md)
            for field,value in fractions.items():
                if version=='historical':
                    old=vp[(vp.community=='full_species_community')&(vp.depth_cm==depth)&(vp.predictor_set==names[field])]
                    equal(f'{depth} cm {field}: historical variance fraction',value,old['Adj.R.square'].iloc[0])
                rows.append(dict(depth_cm=depth,version=version,component=field,unique_adjusted_R2=value,percentage_total_variation=100*value))
    csv(out/'variance_partition_comparison.csv',pd.DataFrame(rows))
    csv(out/'dispersion_center_audit.csv',pd.DataFrame(disp_checks))
    csv(out/'hill_profiles.csv',old_hp.reset_index())
    csv(out/'hill_adjusted_means.csv',em)
    csv(out/'hill_tests.csv',ht.reset_index())
    csv(out/'pcoa_scores.csv',pcoa)
    csv(out/'composition_tests.csv',result)
    csv(out/'dispersion_tests.csv',dispersion)
    csv(out/'metadata_audited.csv',meta.reset_index(drop=True))
    csv(out/'profile_counts.csv',pm.groupby(['locality','restoration4']).size().rename('n_profiles').reset_index())
    taxa=[]
    for universe,ids in [('complete',counts.columns),('candidate',raw.columns)]:
        taxgroup=tax.loc[ids,'kingdom'].fillna('Unassigned').value_counts()
        for kingdom,n in taxgroup.items():taxa.append(dict(universe=universe,kingdom=kingdom,n_taxa=int(n)))
    csv(out/'existing_taxonomy_inventory.csv',pd.DataFrame(taxa))
    csv(out/'numeric_validation.csv',pd.DataFrame(checks))
    provenance=[]
    for name,p in dict(phyloseq=Path(args.phyloseq),canonical_metadata=Path(args.metadata),**used).items():
        b=p.read_bytes();provenance.append(dict(role=name,path=str(p.resolve()),size_bytes=len(b),sha256=hashlib.sha256(b).hexdigest()))
    csv(out/'input_provenance.csv',pd.DataFrame(provenance))
    summary=dict(n_samples=51,n_profiles=17,n_full_taxa=12761,n_candidate_taxa=1031,
                 validation_checks=len(checks),all_checks_passed=True,
                 metadata_stage_locality_depth_profile_concordance=True,
                 historical_current_eci_correlation=float(np.corrcoef(old_eci,current_eci)[0,1]),
                 permutation_p_values='Retained from archived R 301 analysis; not regenerated',
                 depth_partition='Recomputed with current metadata_003b ECI; historical fractions independently reproduced',
                 taxonomic_selection='Unchanged; existing tax_table includes viral assignments',
                 depth_significance='No sequential-test asterisk displayed on unique adjusted fractions')
    atomic(out/'audit_summary.json',(json.dumps(summary,indent=2)+'\n').encode())
    print(json.dumps(summary,indent=2));print(pd.DataFrame(rows).to_string(index=False))

if __name__=='__main__':main()
