"""Post hoc audit of the comparison with previous analyses; no new hypothesis tests."""
from pathlib import Path
import json, hashlib, argparse
import numpy as np
import pandas as pd
from marker_annotations import ANNOTATIONS, ORDER

ROOT=Path(__file__).resolve().parents[1];I=ROOT/'inputs';R=ROOT/'results'

def main(verify_editorial_baseline=False):
    previous=pd.read_csv(I/'previous_marker_restoration_tests.tsv',sep='\t')
    net=pd.read_csv(I/'previous_network_200_KOs.tsv',sep='\t')
    meta=pd.read_csv(R/'sample_manifest.tsv',sep='\t')
    samples=meta.sample_id.tolist()
    frac=pd.read_csv(I/'KO_complete_RPK.tsv',sep='\t',index_col=0).loc[:,samples]
    uscgs=(I/'uscg_76_kegg_min_2313.lst').read_text().split()
    reference=frac.loc[uscgs].median()
    total=pd.read_csv(I/'full_KO_RPK_totals.tsv',sep='\t',index_col=0).total_KO_RPK.loc[samples]
    primary=pd.read_csv(R/'primary_76_normalized.tsv',sep='\t',index_col=0).loc[:,samples]
    exact_numerator=primary.mul(reference,axis=1)
    variants={'primary_reference76':primary,
              'same_numerator_totalKO':exact_numerator.div(total,axis=1),
              'fractional_reference76':frac.loc[primary.index].div(reference,axis=1),
              'fractional_totalKO':frac.loc[primary.index].div(total,axis=1)}
    records=[]
    for variant,values in variants.items():
        for gene in ORDER:
            k=ANNOTATIONS[gene]['KO']
            wide=pd.DataFrame({'profile':meta.profile_id,'depth':meta.depth_cm,
                               'value':np.log2(values.loc[k,samples].values)}).pivot(index='profile',columns='depth',values='value')
            for depth in [20,40]:
                delta=wide[depth]-wide[5]
                records.append(dict(gene=gene,KO=k,cycle=ANNOTATIONS[gene]['cycle'],variant=variant,
                    contrast=f'{depth}_vs_5',n_profiles=len(delta),n_lower=int((delta<0).sum()),
                    median_log2=float(delta.median()),median_ratio=float(np.exp2(delta.median())),
                    scope='post_hoc_descriptive_denominator_audit_no_new_significance_test'))
    pd.DataFrame(records).to_csv(R/'normalization_comparison.tsv',sep='\t',index=False)
    membership=pd.DataFrame([dict(gene=g,KO=ANNOTATIONS[g]['KO'],cycle=ANNOTATIONS[g]['cycle'],
        included_in_previous_network=ANNOTATIONS[g]['KO'] in set(net.KO),
        previous_ordinal_tests=int(previous.feature_id.eq(ANNOTATIONS[g]['KO']).sum()),
        minimum_previous_restoration_q=float(previous.loc[previous.feature_id.eq(ANNOTATIONS[g]['KO']),'q_BH'].min())) for g in ORDER])
    membership.to_csv(R/'previous_analysis_marker_comparison.tsv',sep='\t',index=False)
    if verify_editorial_baseline:
        hashes=json.loads((R/'numeric_hashes_before_editorial_revision.json').read_text())
        assert all(hashlib.sha256((R/name).read_bytes()).hexdigest()==h for name,h in hashes.items())
    summary=dict(primary_numerical_outputs_unchanged=True if verify_editorial_baseline else None,previous_network_KOs=len(net),
        panel_markers_in_previous_network=int(membership.included_in_previous_network.sum()),
        previous_targeted_restoration_tests=len(previous),
        previous_targeted_restoration_BH_below_010=int((previous.q_BH<.1).sum()),
        analysis_status='Targeted exploratory follow-up; not independent validation or original-study preregistration',
        adjustment_scope='Holm covers the declared 33-test marker family, not all analyses attempted in the project',
        normalization_diagnostic='The total-KO denominator changes the estimand to composition; this diagnostic reports all 11 markers without new P values')
    (R/'comparison_with_previous_analyses.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify-editorial-baseline',action='store_true',help='Check byte identity against the pre-edit files, only for this editorial revision; omit after a scientific rerun.')
    args=parser.parse_args();main(args.verify_editorial_baseline)
