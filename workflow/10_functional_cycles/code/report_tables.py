"""Generate the two essential supplementary source tables, retaining every marker."""
from pathlib import Path
import pandas as pd
from marker_annotations import export_guide, ORDER

ROOT = Path(__file__).resolve().parents[1]


def main():
    export_guide(ROOT)
    result = ROOT/'results'
    guide = pd.read_csv(result/'marker_function_guide.tsv', sep='\t')
    tests = pd.read_csv(result/'primary_tests.tsv', sep='\t')
    depth = pd.read_csv(result/'depth_effects_for_reporting.tsv', sep='\t')
    table = guide[['gene', 'KO', 'cycle', 'label', 'allocation_status', 'limitation', 'source_url']].copy()
    restoration = tests.loc[tests.test.eq('restoration_omnibus'), ['KO', 'statistic', 'p', 'p_holm_33']]
    table = table.merge(restoration.rename(columns={c: 'restoration_'+c for c in ['statistic', 'p', 'p_holm_33']}), on='KO', validate='one_to_one')
    for cm in [20, 40]:
        names = ['median_ratio', 'ratio_ci_low', 'ratio_ci_high', 'p', 'p_holm_33', 'n_profiles', 'n_positive', 'n_zero']
        subset = depth.loc[depth.test.eq(f'depth_{cm}_vs_5'), ['KO', *names]]
        table = table.merge(subset.rename(columns={c: f'depth_{cm}_vs_5_'+c for c in names}), on='KO', validate='one_to_one')
    table = table.set_index('gene').loc[ORDER].reset_index()
    assert len(table) == 11
    table.to_csv(result/'Supplementary_Table16_CNS.tsv', sep='\t', index=False)
    environment = pd.read_csv(result/'environment_tests.tsv', sep='\t')
    environment = environment.loc[environment.variant.eq('primary_76')].copy()
    environment = environment.merge(guide[['KO', 'cycle', 'allocation_status']], on='KO', validate='many_to_one')
    assert len(environment) == 19
    environment.to_csv(result/'Supplementary_Table17_CNS.tsv', sep='\t', index=False)
    print('Wrote all 11 marker rows and all 19 environment rows; no significance filtering.')


if __name__ == '__main__':
    main()
