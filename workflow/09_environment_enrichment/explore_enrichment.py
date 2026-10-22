"""Exploratory ORA of Figure 4a common-slope discoveries; no model refitting.

python3 explore_enrichment.py --source project/Figura4 --output results
Primary correction: BH across both ontologies, five axes and both directions.
The gene-q=0.05 sensitivity is corrected as a separate family.
"""
from pathlib import Path
import argparse, json, hashlib, platform
import numpy as np
import pandas as pd
import scipy
from scipy.stats import hypergeom, fisher_exact

AXES = ['vegetation_landscape_PC1', 'water_inundation_PC1',
        'moisture_stress_PC1', 'physicochemical_PC1', 'nutrients_redox_PC1']
VARIABLES = ['shdi','ndvi_p50','ndvi_sd','ndwi_p50','ndwi_sd','mndwi_p50',
             'mndwi_sd','water_prop_mndwi','ndmi_p50','ndmi_sd','temperature_c',
             'salinity_ups','ph','redox_shallow_mv','n_no2_umol_l','n_no3_umol_l',
             'n_nh3_mg_l_campo','n_nh4plus_umol_l','p_po4_3_umol_l','s_2_umol']

def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p, kind='stable')
    q = np.ones(len(p))
    if len(p):
        q[order] = np.minimum(1, np.minimum.accumulate(
            (p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1])
    return q

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--single-ko',action='store_true',help='KO sensitivity: exclude multi-KO clusters from the KO universe and hits; PFAM unchanged')
    args = ap.parse_args()
    source, out = args.source, args.output
    out.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(source/'prepared/selected_genes.tsv', sep='\t', keep_default_na=False)
    assert len(genes) == 100000 and genes.gene_numeric_id.is_unique
    annotations = {}
    for ont in ['KO','PFAM']:
        ann = genes[['gene_numeric_id',ont+'_set']].rename(columns={ont+'_set':'term'})
        if ont=='KO' and args.single_ko:
            ann=ann.loc[~ann.term.str.contains(',')].copy()
        ann = ann.assign(term=ann.term.str.split(',')).explode('term')
        ann['term'] = ann.term.str.strip()
        ann = ann.loc[~ann.term.isin(['','-','NA'])].drop_duplicates()
        annotations[ont] = ann
    records, universes, discoveries = [], [], []
    for axis in AXES:
        # Read all chunks and retain only common-slope rows, rather than assuming ordering.
        common = pd.concat([c.loc[c.family.eq('common')] for c in pd.read_csv(
            source/'models'/f'{axis}_all_results.tsv.gz', sep='\t', chunksize=100000)],ignore_index=True)
        assert len(common)==100000 and common.gene_numeric_id.is_unique
        assert set(common.gene_numeric_id)==set(genes.gene_numeric_id)
        assert common.usable.dtype == bool and common.singular.dtype == bool
        assert not (common.usable & common.singular).any()
        usable = common.loc[common.usable & ~common.singular].copy()
        selected = usable.loc[usable.q.le(.10)].copy()
        selected['direction'] = np.where(selected.estimate_sd.gt(0),'positive','negative')
        selected['also_gene_q005'] = selected.q.le(.05)
        selected = selected.merge(genes, on='gene_numeric_id',validate='one_to_one')
        assert selected.feature_id.equals(selected.historical_cluster_id)
        discoveries.append(selected)
        for ont, annotation in annotations.items():
            ann = annotation.loc[annotation.gene_numeric_id.isin(usable.gene_numeric_id)]
            background = pd.Index(ann.gene_numeric_id.unique())
            N = len(background)
            counts = ann.groupby('term').gene_numeric_id.nunique()
            counts = counts[counts.ge(10)]  # independent of outcome selection
            for cutoff in [.10,.05]:
                for direction, sign in [('negative',-1),('positive',1)]:
                    hit = usable.loc[usable.q.le(cutoff) & usable.estimate_sd.mul(sign).gt(0)]
                    n = hit.gene_numeric_id.isin(background).sum()
                    a = ann.loc[ann.gene_numeric_id.isin(hit.gene_numeric_id)].groupby('term').gene_numeric_id.nunique().reindex(counts.index,fill_value=0)
                    K = counts
                    b, c, d = n-a, K-a, N-K-n+a
                    assert min(a.min(),b.min(),c.min(),d.min())>=0
                    p = hypergeom.sf(a.to_numpy()-1,N,K.to_numpy(),n) if n else np.ones(len(K))
                    records.append(pd.DataFrame({
                        'gene_q_cutoff':cutoff,'axis':axis,'ontology':ont,'direction':direction,
                        'term':K.index,'background_genes':N,'selected_genes':n,
                        'background_term_genes':K.to_numpy(),'selected_term_genes':a.to_numpy(),
                        'expected_term_genes':n*K.to_numpy()/N,
                        'fold_enrichment':(a.to_numpy()/n)/(K.to_numpy()/N) if n else np.nan,
                        'odds_ratio_haldane':((a+.5)*(d+.5)/((b+.5)*(c+.5))).to_numpy(),
                        'p_fisher_greater':p}))
                    universes.append({'axis':axis,'ontology':ont,'direction':direction,
                        'gene_q_cutoff':cutoff,'tested_genes':len(common),'usable_genes':len(usable),
                        'ontology_annotated_usable_genes':N,'selected_all_genes':len(hit),
                        'selected_annotated_genes':n,'tested_terms':len(K)})
        print(axis, 'usable',len(usable),'selected q.1',len(selected),flush=True)
    results = pd.concat(records,ignore_index=True)
    for _, ix in results.groupby('gene_q_cutoff').groups.items():
        results.loc[ix,'q_ora_global'] = bh(results.loc[ix,'p_fisher_greater'])
    for _, ix in results.groupby(['gene_q_cutoff','ontology']).groups.items():
        results.loc[ix,'q_ora_within_ontology'] = bh(results.loc[ix,'p_fisher_greater'])
    results['enriched_exploratory'] = results.q_ora_global.le(.1) & results.selected_term_genes.ge(5)
    # Meaningful implementation check against an independent exact-test function.
    for row in results.sample(min(120,len(results)),random_state=4701).itertuples():
        a=row.selected_term_genes; b=row.selected_genes-a; c=row.background_term_genes-a
        d=row.background_genes-row.background_term_genes-b
        expected=fisher_exact([[a,b],[c,d]],alternative='greater').pvalue
        assert np.isclose(expected,row.p_fisher_greater,rtol=1e-9,atol=1e-14)
    results.to_csv(out/'all_ORA_tests.tsv.gz',sep='\t',index=False)
    sig=results.loc[results.enriched_exploratory].sort_values(['gene_q_cutoff','axis','ontology','q_ora_global'])
    sig.to_csv(out/'enriched_terms.tsv',sep='\t',index=False)
    primary=sig.loc[sig.gene_q_cutoff.eq(.1)].copy()
    sensitivity=sig.loc[sig.gene_q_cutoff.eq(.05),['axis','ontology','direction','term']].assign(retained_gene_q005=True)
    primary=primary.merge(sensitivity,how='left',on=['axis','ontology','direction','term'],validate='one_to_one')
    primary['retained_gene_q005']=primary.retained_gene_q005.eq(True)
    primary.to_csv(out/'primary_enriched_terms.tsv',sep='\t',index=False)
    pd.DataFrame(universes).to_csv(out/'tested_universes.tsv',sep='\t',index=False)
    discoveries=pd.concat(discoveries,ignore_index=True)
    discoveries.to_csv(out/'associated_genes_Figure4A.tsv.gz',sep='\t',index=False)
    counts=primary.groupby(['axis','ontology','direction']).agg(enriched_terms=('term','size'),
        retained_gene_q005=('retained_gene_q005','sum')).reset_index()
    counts.to_csv(out/'enrichment_counts.tsv',sep='\t',index=False)
    # Descriptive correlations only. They are not original PCA loadings or tests.
    meta=pd.read_csv(source/'prepared/metadata.csv')
    corr=[]
    for axis in AXES:
        frame=meta.groupby('profile_id')[VARIABLES+AXES].mean() if axis in AXES[:3] else meta
        for v in VARIABLES:
            complete=frame[[v,axis]].dropna()
            corr.append({'axis':axis,'variable':v,'pearson_r':complete[v].corr(complete[axis]),
                         'observations':len(complete),'unit':'profile_mean' if axis in AXES[:3] else 'sample',
                         'interpretation':'descriptive correlation, not original PCA loading'})
    pd.DataFrame(corr).to_csv(out/'axis_orientation_descriptive.tsv',sep='\t',index=False)
    meta[AXES].corr().to_csv(out/'between_axis_correlations.tsv',sep='\t')
    params={'analysis':'Exploratory ORA, Figure4A common slopes only',
        'single_KO_only_for_KO_tests':args.single_ko,
        'axes':AXES,'gene_q_primary':.1,'gene_q_sensitivity':.05,'minimum_background_term_genes':10,
        'minimum_selected_term_genes_after_BH':5,'ora_q':.1,
        'multiplicity_primary':'BH across KO and PFAM, five axes and two signs, separately for each gene-q cutoff',
        'tests_per_cutoff':results.groupby('gene_q_cutoff').size().to_dict(),
        'MHI':'omitted: no common-slope primary discoveries at gene q<=.10',
        'duplicate_rule':'one distinct historical gene-cluster ID per term, multiannotated genes may belong to multiple terms',
        'limitations':['ORA treats distinct gene clusters as independent units; gene coabundance may make nominal FDR optimistic.',
            'Selection is limited to the 100000 most variable genes.',
            'Axis models fitted separately; association is not a unique causal environmental effect.',
            'DNA abundance is relative functional potential, not expression, activity or process rate.'],
        'validation':'120 sampled Fisher exact tests agree with hypergeometric; ID uniqueness and hit/background containment checked',
        'software':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__}}
    (out/'parameters.json').write_text(json.dumps(params,indent=2,ensure_ascii=False)+'\n')
    print(counts.to_string(index=False),flush=True)

if __name__=='__main__':
    main()
