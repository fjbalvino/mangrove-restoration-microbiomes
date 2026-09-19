#!/usr/bin/env python3
# ============================================================
# 07_00_select_features.py
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Historical 801 selection implementation; requires its original input layout. Current 809 bundle freezes the resulting 200 KO.
# Inputs (source expressions; complete list in docs/contracts/07_00_select_features.json):
#   return hashlib.sha256(p.read_bytes()).hexdigest()
#   summary=pd.read_csv(source/'082A_run_summary.tsv',sep='\t').iloc[0]
#   hist=pd.read_csv(source/'082A_metadata_aligned_51samples.csv',dtype={'lat_block':str})
#   meta=pd.read_csv(base/'inputs/metadata/003b_metadata_integrada_canon_51.csv',dtype={'lat_block':str})
#   d=pd.read_csv(f,sep='\t',index_col=0)
# Outputs (source expressions; complete list in contract):
#   def save(rows,name,sub='tables'):
#   df.to_csv(out/sub/name,sep='\t',index=False,float_format='%.17g')
#   save(meta,'801_metadata_aligned.tsv')
#   save(ranking,f'801_{layer}_{method}_ranking.tsv')
#   mat.to_csv(dest,sep='\t',index_label='sample_id',float_format='%.17g')
#   save(layer_rows,'801a_functional_matrix_summary.tsv','audit')
#   pd.DataFrame(checks).rename(columns={'pass_':'pass'}).to_csv(out/'audit/801a_acceptance_checks.tsv',sep='\t',index=False)
#   save(pd.concat(selections),'801_selected_functions.tsv')
#   save(manifest,'801_network_input_scenario_manifest.tsv')
#   detail=pd.DataFrame(loo); save(detail,'801_leave_one_profile_out_detail.tsv')
# Algorithmic provenance:
# Technical feature selection by prevalence and variance or abundance; 200-KO primary set.
#   Kuijjer et al. (2019), doi:10.1016/j.isci.2019.03.021. Fixed lambda is a study choice, not analytic shrinkage estimation. No Louvain algorithm is used.
# Source SHA-256: 7a8d6b4c9214db9c989872b8582a021e9e78febb67a2b5f0cabc8a21d835477b
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================

"""Computational selection/LOPO counterpart of 801; no networks or RDS produced.
See README for historical provenance and differences from the R reporting layer.
"""
import argparse, hashlib, json, math, platform, time, shutil
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd

PRIMARY = 'KEGG_ko__clr_variance_N200'
METHODS = ['clr_variance', 'mean_relative_abundance']
CAPS = [100, 200, 300]

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def clr(x):
    logged = np.log(x + 1.0)
    return logged - logged.mean(axis=1, keepdims=True)

def candidates(frame, layer):
    """Relative abundance denominator is the complete layer, before prevalence."""
    x = frame.to_numpy(float)
    prev = (x > 0).sum(axis=0)
    eligible = prev >= math.ceil(.5 * len(frame))
    z = clr(x[:, eligible])
    return pd.DataFrame(dict(
        annotation_layer=layer, feature_id=frame.columns[eligible],
        prevalence_n=prev[eligible], prevalence_prop=prev[eligible]/len(frame),
        total_abundance=x.sum(axis=0)[eligible],
        mean_relative_abundance=(x/x.sum(axis=1, keepdims=True)).mean(axis=0)[eligible],
        clr_variance=z.var(axis=0, ddof=1)))

def rank(stats, method):
    return stats.sort_values([method, 'prevalence_n','total_abundance','feature_id'],
                           ascending=[False,False,False,True], kind='stable').reset_index(drop=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bundle', type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument('--out-root', type=Path)
    args=ap.parse_args(); base=args.bundle.resolve()
    root=args.out_root or base/'results'
    out=root/('801_seleccionar_features_red_funcional_python_'+datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%fZ'))
    for sub in ['tables','matrices','audit','logs']: (out/sub).mkdir(parents=True,exist_ok=False)
    (out/"logs/code_executed").mkdir()
    for code in (base/"scripts").glob("*"):
        if code.is_file(): shutil.copy2(code, out/"logs/code_executed"/code.name)
    started=time.monotonic(); checks=[]
    def check(name, ok, critical=True):
        checks.append(dict(check=name,pass_=bool(ok),critical=critical))
        if critical and not ok: raise ValueError(name)
    def save(rows,name,sub='tables'):
        df=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows)
        df.to_csv(out/sub/name,sep='\t',index=False,float_format='%.17g')
    source=base/'inputs/082A/tables'
    summary=pd.read_csv(source/'082A_run_summary.tsv',sep='\t').iloc[0]
    hist=pd.read_csv(source/'082A_metadata_aligned_51samples.csv',dtype={'lat_block':str})
    meta=pd.read_csv(base/'inputs/metadata/003b_metadata_integrada_canon_51.csv',dtype={'lat_block':str})
    check('51 unique metadata samples',len(meta)==51 and meta.sample_id.is_unique)
    check('same historical and updated sample universe',set(meta.sample_id)==set(hist.sample_id))
    # Preserve historical row order. Health/environment fields never enter ranking.
    meta=meta.set_index('sample_id').loc[hist.sample_id].reset_index()
    for col in ['depth_cm','locality','lat_block']:
        if col in ['depth_cm','lat_block']:
            ok=np.allclose(pd.to_numeric(meta[col]),pd.to_numeric(hist[col]),rtol=0,atol=1e-9)
        else: ok=np.array_equal(meta[col].astype(str),hist[col].astype(str))
        check('historical design concordance: '+col,ok)
    groups=meta.groupby('lat_block',sort=False)
    check('17 complete locality-nested profiles',len(groups)==17 and all(
        len(g)==3 and set(g.depth_cm)=={5,20,40} and g.locality.nunique()==1 for _,g in groups))
    check('historical summary PASS',summary.overall_status=='PASS')
    for col in ['indices_used_for_aggregation','functional_cap_applied','transformation_applied','network_inference_performed']:
        check('historical declaration: '+col,str(summary[col]).lower()=='false')
    save(meta,'801_metadata_aligned.tsv')
    frames={}; layer_rows=[]; selected={}; selections=[]; manifest=[]
    for layer in ['KEGG_ko','PFAM']:
        f=source/f'082A_{layer}_abundance_samples_x_functions.tsv.gz'
        d=pd.read_csv(f,sep='\t',index_col=0)
        check(layer+' sample and feature identifiers',d.index.is_unique and d.columns.is_unique and set(d.index)==set(meta.sample_id))
        d=d.loc[meta.sample_id]; x=d.to_numpy(float)
        check(layer+' nonnegative finite complete counts',np.isfinite(x).all() and (x>=0).all() and (x.sum(axis=1)>0).all())
        frames[layer]=d; st=candidates(d,layer)
        check(layer+' totals match historical summary',d.shape[1]==int(summary['n_complete_'+layer]) and len(st)==int(summary['n_eligible_'+layer]))
        layer_rows.append(dict(annotation_layer=layer,n_samples=len(d),n_complete_features=d.shape[1],n_features_prevalence_ge_threshold=len(st),min_prevalence_n=26))
        for method in METHODS:
            ranking=rank(st,method); ranking['global_rank']=np.arange(1,len(ranking)+1)
            save(ranking,f'801_{layer}_{method}_ranking.tsv')
            for cap in CAPS:
                sid=f'{layer}__{method}_N{cap}'; chosen=ranking.head(cap).copy()
                chosen['scenario_id']=sid; chosen['selection_method']=method
                selected[sid]=chosen.feature_id.tolist(); selections.append(chosen)
                c=d.loc[:,selected[sid]]; z=clr(c.to_numpy(float))
                check(sid+' correct size/prevalence',len(chosen)==cap and (chosen.prevalence_n>=26).all())
                check(sid+' CLR centered',np.abs(z.mean(axis=1)).max()<1e-10)
                # Ratio reconstruction detects loss of the original count ratios.
                check(sid+' CLR ratio identity',np.allclose(z[:,1:]-z[:,[0]],np.log((c.to_numpy()[:,1:]+1)/(c.to_numpy()[:,[0]]+1)),atol=1e-12,rtol=1e-12))
                for name,mat in [('counts',c),('clr',pd.DataFrame(z,index=c.index,columns=c.columns))]:
                    dest=out/'matrices'/f'801_{sid}_{name}.tsv.gz'
                    mat.to_csv(dest,sep='\t',index_label='sample_id',float_format='%.17g')
                manifest.append(dict(scenario_id=sid,annotation_layer=layer,selection_method=method,n_functions=cap,n_samples=51,is_primary=sid==PRIMARY,
                    counts_file=f'matrices/801_{sid}_counts.tsv.gz',clr_file=f'matrices/801_{sid}_clr.tsv.gz'))
        print(layer,'complete',d.shape[1],'eligible',len(st),flush=True)
    save(layer_rows,'801a_functional_matrix_summary.tsv','audit')
    # This is a new input audit, not the unavailable historical acceptance table.
    pd.DataFrame(checks).rename(columns={'pass_':'pass'}).to_csv(out/'audit/801a_acceptance_checks.tsv',sep='\t',index=False)
    save(pd.concat(selections),'801_selected_functions.tsv')
    save(manifest,'801_network_input_scenario_manifest.tsv')
    loo=[]; loo_sets={sid:[] for sid in selected}; frequencies=[]
    for number,(block,g) in enumerate(groups,1):
        train=meta.loc[meta.lat_block!=block,'sample_id']
        for layer,d in frames.items():
            st=candidates(d.loc[train],layer)
            check(f'{block}/{layer} supports largest cap',len(st)>=max(CAPS))
            for method in METHODS:
                ranked=rank(st,method)
                for cap in CAPS:
                    sid=f'{layer}__{method}_N{cap}'; ids=ranked.head(cap).feature_id.tolist()
                    a=set(ids); b=set(selected[sid]); overlap=len(a&b)
                    loo_sets[sid].append(a)
                    loo.append(dict(scenario_id=sid,omitted_lat_block=block,omitted_locality=g.locality.iloc[0],n_train_samples=len(train),
                        min_prevalence_n=24,n_eligible_features=len(st),overlap_with_full_n=overlap,
                        retention_of_full=overlap/cap,jaccard_with_full=overlap/len(a|b),novel_vs_full_n=len(a-b)))
                    for pos,f in enumerate(ids,1): frequencies.append(dict(scenario_id=sid,omitted_lat_block=block,rank=pos,feature_id=f))
        print(f'LOPO {number}/17: omitted {block}',flush=True)
    detail=pd.DataFrame(loo); save(detail,'801_leave_one_profile_out_detail.tsv')
    save(frequencies,'801_leave_one_profile_out_selected_functions.tsv.gz')
    rows=[]; freqrows=[]
    for sid,g in detail.groupby('scenario_id',sort=False):
        v=g.jaccard_with_full
        rows.append(dict(scenario_id=sid,n_profiles_omitted=len(g),min_jaccard=v.min(),median_jaccard=v.median(),max_jaccard=v.max(),median_retention_of_full=g.retention_of_full.median()))
        for f in sorted(set(selected[sid]).union(*loo_sets[sid])):
            n=sum(f in ss for ss in loo_sets[sid])
            freqrows.append(dict(scenario_id=sid,feature_id=f,selected_in_full=f in selected[sid],n_leave_one_profile_out_selected=n,leave_one_profile_out_selection_frequency=n/17))
    stability=pd.DataFrame(rows); save(stability,'801_leave_one_profile_out_summary.tsv')
    save(freqrows,'801_function_selection_stability.tsv')
    for layer in frames:
        for method in METHODS:
            check(f'{layer}/{method} nested caps',set(selected[f'{layer}__{method}_N100'])<=set(selected[f'{layer}__{method}_N200'])<=set(selected[f'{layer}__{method}_N300']))
        sid=f'{layer}__clr_variance_N200'
        check(sid+' median Jaccard >=0.70',float(stability.set_index('scenario_id').loc[sid,'median_jaccard'])>=.70,False)
    check('12 scenarios and 204 profile-scenario results',len(selected)==12 and len(detail)==204)
    overlap=[]
    for i,(sa,a) in enumerate(selected.items()):
        for sb,b in list(selected.items())[i+1:]:
            if sa.split('__')[0]==sb.split('__')[0]:
                overlap.append(dict(scenario_a=sa,scenario_b=sb,intersection=len(set(a)&set(b)),jaccard=len(set(a)&set(b))/len(set(a)|set(b))))
    save(overlap,'801_selection_pairwise_overlap.tsv')
    save(pd.DataFrame(checks).rename(columns={'pass_':'pass'}),'801_acceptance_checks.tsv')
    status='PASS' if all(c['pass_'] for c in checks) else 'PASS_WITH_WARNINGS'
    provenance=dict(status=status,execution_language='Python',R_executed=False,python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,
        elapsed_seconds=time.monotonic()-started,primary=PRIMARY,prevalence=.5,pseudocount=1,
        ranking_clr_universe='prevalence eligible functions within each full/training set',network_clr_universe='selected functions in each scenario',
        networks_inferred=False,original_aggregation_repeated=False,source='historical 082A recovered from 605a',
        input_sha256={str(f.relative_to(base)):digest(f) for f in sorted((base/'inputs').rglob('*')) if f.is_file()},
        code_sha256={str(f.relative_to(base)):digest(f) for f in sorted((base/'scripts').glob('*')) if f.is_file()})
    (out/'801_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    save([dict(overall_status=status,n_samples=51,n_profiles=17,n_scenarios=12,n_loo_results=204,n_checks=len(checks),n_failed=sum(not c['pass_'] for c in checks))],'801_run_summary.tsv')
    (root/'LATEST_801_python.txt').write_text(str(out.resolve())+'\n')
    print(stability.to_string(index=False),flush=True)
    print('STATUS',status,'OUTPUT',out.resolve(),flush=True)
if __name__=='__main__': main()
