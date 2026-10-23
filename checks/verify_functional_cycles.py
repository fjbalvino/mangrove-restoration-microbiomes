#!/usr/bin/env python3
"""Independent profile, multiplicity, regression and source-table checks for CNS."""
from pathlib import Path
import argparse
import json
import math
import sys
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy import stats
from PIL import Image

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO/'workflow/10_functional_cycles/code'))
from stats_cns import Design, depth_test, profile_matrix
from check_functional_cycles_inputs import validate_cohort


def rejects(call):
    try:
        call()
    except (ValueError, AssertionError):
        return
    raise AssertionError('Invalid input was accepted')


def self_test():
    rng = np.random.default_rng(442017)
    X = np.column_stack([np.ones(17), np.arange(17)<6, np.arange(17)>=12, rng.normal(size=17)])
    design = Design(X)
    y = rng.normal(size=17)
    beta, residual, covariance, _ = design.fit(y)
    inverse = np.linalg.inv(X.T@X)
    hc3 = inverse@X.T@np.diag((residual[0]/(1-np.diag(X@inverse@X.T)))**2)@X@inverse
    np.testing.assert_allclose(covariance[0], hc3[3:, 3:], atol=1e-12)
    signs = rng.choice([-1., 1.], (99, 17))
    nullfit = X[:, :3]@np.linalg.lstsq(X[:, :3], y, rcond=None)[0]
    rescaled = (y-nullfit)/np.sqrt(1-design.h0)
    explicit = design.fit(nullfit+rescaled*signs)[3]
    direct_p = (1+np.count_nonzero(explicit>=design.fit(y)[3][0]-1e-12))/100
    assert design.wild_test(y, signs)['p'] == direct_p
    assert design.df == 13 and design.n == 17
    zero = depth_test(np.zeros(17))
    assert zero['p'] == 1 and zero['n_nonzero'] == 0 and zero['n_zero'] == 17
    one_zero = depth_test(np.r_[0., np.ones(16)])
    assert one_zero['n_nonzero'] == 16 and one_zero['p'] == 2/2**16
    rows = []
    stages = ['Degraded', 'Early restoration', 'Intermediate restoration', 'Advanced restoration', 'Conserved']
    for p in range(17):
        for depth in [5, 20, 40]:
            rows.append(dict(sample_id=f'profile{p}_{depth}', lat_block=f'profile{p}', depth_cm=depth,
                             locality=['Carmen','Cozumel','Tuxpan'][min(p//6,2)],
                             restoration4=stages[p%5], latitude=p, longitude=p,
                             redox_shallow_mv=0., n_no3_umol_l=1., n_nh4plus_umol_l=1., s_2_umol=1.))
    meta = validate_cohort(pd.DataFrame(rows))
    values = pd.DataFrame(rng.uniform(1,2,(2,51)), index=['a','b'], columns=meta.sample_id)
    profiles = sorted(meta.profile_id.unique())
    base = profile_matrix(values, meta, profiles)
    shuffled = meta.sample(frac=1, random_state=4)
    np.testing.assert_array_equal(base, profile_matrix(values[shuffled.sample_id], shuffled, profiles))
    bad = meta.copy(); bad.loc[1,'depth_cm'] = 5
    rejects(lambda: validate_cohort(bad))
    rejects(lambda: profile_matrix(values, bad, profiles))
    rejects(lambda: validate_cohort(pd.concat([meta, meta.iloc[[0]]])))
    bad = meta.copy(); bad.loc[0,'lat_block'] = 'profile1'
    rejects(lambda: validate_cohort(bad))
    assert np.isclose(depth_test(np.arange(17.))['interval_coverage'], 1-2*sum(math.comb(17,k) for k in range(5))/2**17)
    return ['HC3 matches independent matrix formula', 'restricted bootstrap matches explicit null fits',
            'zeros and ties handled without pseudocounts', 'pairing invariant to sample order',
            'duplicate samples/depths and inconsistent profile identities rejected']


def verify_run(run, reference=None):
    run = Path(run); result = run/'results'
    tests = pd.read_csv(result/'all_primary_and_sensitivity_tests.tsv',sep='\t')
    env = pd.read_csv(result/'environment_tests.tsv',sep='\t')
    assert len(tests) == 165 and len(env) == 95
    assert tests.variant.nunique() == env.variant.nunique() == 5
    for _, frame in tests.groupby('variant'):
        assert len(frame) == 33 and len(frame[['KO','test']].drop_duplicates()) == 33
        p = frame.p.to_numpy(); order = np.argsort(p); adjusted = np.empty(33)
        adjusted[order] = [min(1., max((33-j)*p[order[j]] for j in range(i+1))) for i in range(33)]
        np.testing.assert_allclose(frame.p_holm_33, adjusted, rtol=1e-10, atol=1e-12)
        assert frame.loc[frame.test.eq('restoration_omnibus'),'residual_df'].eq(10).all()
    for _, frame in env.groupby('variant'):
        assert len(frame)==19 and frame.residual_df.eq(13).all()
        np.testing.assert_allclose(frame.q_BH_19, stats.false_discovery_control(frame.p.to_numpy()), rtol=1e-10, atol=1e-12)
    meta = validate_cohort(pd.read_csv(run/'inputs/metadata_003b.csv'))
    samples = meta.sample_id.tolist()
    for variant in tests.variant.unique():
        normalized = pd.read_csv(result/(variant+'_normalized.tsv'),sep='\t',index_col=0)
        assert normalized.shape == (11,51) and (normalized>0).all().all()
        for ko in normalized.index:
            wide = pd.DataFrame({'profile':meta.profile_id, 'depth':meta.depth_cm,
                                 'value':np.log2(normalized.loc[ko,samples].to_numpy())}).pivot(index='profile',columns='depth',values='value')
            for depth in [20,40]:
                delta = (wide[depth]-wide[5]).to_numpy()
                row = tests.loc[tests.variant.eq(variant)&tests.KO.eq(ko)&tests.test.eq(f'depth_{depth}_vs_5')].iloc[0]
                nonzero = delta[~np.isclose(delta,0,rtol=0,atol=1e-12)]
                n = len(nonzero); k = min((nonzero>0).sum(),(nonzero<0).sum())
                exact = min(1.,2*sum(math.comb(n,j) for j in range(k+1))/2**n) if n else 1.
                np.testing.assert_allclose([row.median_log2,row.ci_low,row.ci_high,row.p],
                                           [np.median(delta),np.sort(delta)[4],np.sort(delta)[12],exact], rtol=1e-10,atol=1e-12)
    archive = np.load(result/'bootstrap_profile_signs.npz')
    weights, profiles = archive['signs'], archive['profile_id']
    assert weights.shape == (99999,17) and set(np.unique(weights)) == {-1,1}
    index = {pid:i for i,pid in enumerate(profiles)}
    expanded = weights[:1000,[index[pid] for pid in meta.profile_id]]
    for pid in profiles:
        assert np.ptp(expanded[:,meta.profile_id.eq(pid)],axis=1).max() == 0
    guide = pd.read_csv(result/'marker_function_guide.tsv',sep='\t')
    assert guide.cycle.value_counts().to_dict() == {'Nitrogen':6,'Sulfur':3,'Carbon':2}
    assert guide.allocation_status.eq('provisional_fractional').sum() == 5
    table16 = pd.read_csv(result/'Supplementary_Table16_CNS.tsv',sep='\t')
    table17 = pd.read_csv(result/'Supplementary_Table17_CNS.tsv',sep='\t')
    assert len(table16)==11 and table16.KO.is_unique and len(table17)==19
    expected_env = env.loc[env.variant.eq('primary_76')]
    pd.testing.assert_frame_equal(table17[expected_env.columns].reset_index(drop=True),expected_env.reset_index(drop=True),check_dtype=False)
    for depth in [20,40]:
        expected = tests.loc[tests.variant.eq('primary_76')&tests.test.eq(f'depth_{depth}_vs_5')].set_index('KO').loc[table16.KO]
        np.testing.assert_allclose(table16[f'depth_{depth}_vs_5_median_ratio'],np.exp2(expected.median_log2),rtol=1e-10,atol=1e-12)
        np.testing.assert_allclose(table16[f'depth_{depth}_vs_5_p_holm_33'],expected.p_holm_33,rtol=1e-10,atol=1e-12)
    pd.testing.assert_frame_equal(table16[['KO','cycle']].sort_values('KO').reset_index(drop=True),guide[['KO','cycle']].sort_values('KO').reset_index(drop=True))
    matched = []
    if reference:
        reference = Path(reference)
        for path in sorted(result.glob('*.tsv')):
            other = reference/path.name
            if not other.exists():
                if path.name not in {'Supplementary_Table16_CNS.tsv','Supplementary_Table17_CNS.tsv'}:
                    raise ValueError('Reference lacks generated table: '+path.name)
                continue
            pd.testing.assert_frame_equal(pd.read_csv(path,sep='\t'),pd.read_csv(other,sep='\t'),
                                          check_dtype=False,check_exact=False,rtol=1e-9,atol=1e-11)
            matched.append(path.name)
        if len(matched)<20:
            raise ValueError('Insufficient reference tables for a complete numerical comparison')
    figures = []
    execution = json.loads((run/'execution_manifest.json').read_text())
    if execution['figures_requested']:
        for stem in ['Figura_CNS_profundidad','Figura_CNS_sensibilidad']:
            png = run/'figures'/(stem+'.png')
            with Image.open(png) as image:
                image.verify()
            with Image.open(png) as image:
                image.load()
                assert min(image.size)>3000
            pdf = run/'figures'/(stem+'.pdf')
            content = pdf.read_bytes()
            assert content.startswith(b'%PDF-') and content.rstrip().endswith(b'%%EOF')
            svg = run/'figures'/(stem+'.svg')
            assert ET.parse(svg).getroot().tag.endswith('svg')
            figures += [png.name,pdf.name,svg.name]
    report = {'status':'PASS','profiles':17,'samples':51,'primary_tests_per_variant':33,
              'environment_tests_per_variant':19,'quantification_variants':5,
              'independently_reconstructed_depth_tests':110,'reference_tables_compared':matched,
              'complete_figure_files_verified':figures,
              'scope':'All generated scientific tables compared when a reference is supplied; simulations and raw-read reconstruction are separate.'}
    (result/'repository_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--run-root',type=Path)
    parser.add_argument('--reference-results',type=Path)
    args=parser.parse_args()
    if not __debug__:
        parser.error('Assertions must remain enabled.')
    if not args.self_test and not args.run_root:
        parser.error('Supply --self-test or --run-root.')
    if args.reference_results and not args.run_root:
        parser.error('--reference-results requires --run-root.')
    if args.self_test:
        print(json.dumps({'status':'PASS','checks':self_test()},indent=2))
    if args.run_root:
        print(json.dumps(verify_run(args.run_root,args.reference_results),indent=2))


if __name__=='__main__':
    main()
