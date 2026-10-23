"""Post-run independent checks and compact derived reporting tables."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy import stats
from stats_cns import depth_test,adjust,profile_matrix
from calibrar_y_verificar import designs
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results'

def main():
    m,p,designs_=designs();m=m.sort_values(['profile_id','depth_cm']).reset_index(drop=True)
    vals=pd.read_csv(R/'primary_76_normalized.tsv',sep='\t',index_col=0)
    a=np.log2(profile_matrix(vals,m,p.profile_id));t=pd.read_csv(R/'primary_tests.tsv',sep='\t')
    assert len(t)==33 and len(t[['KO','test']].drop_duplicates())==33
    checks=[]
    for j,k in enumerate(vals.index):
        for dd,depth in [(1,20),(2,40)]:
            row=t[t.KO.eq(k)&t.test.eq(f'depth_{depth}_vs_5')].iloc[0]
            # Recreate pairing directly by sample IDs rather than relying on reshape.
            wide=pd.DataFrame({'profile':m.profile_id,'depth':m.depth_cm,'value':np.log2(vals.loc[k,m.sample_id].values)}).pivot(index='profile',columns='depth',values='value')
            d=wide[depth]-wide[5]
            assert np.isclose(row.median_log2,np.median(d),atol=1e-12)
            assert np.isclose(row.p,stats.binomtest(int((d>0).sum()),len(d),.5).pvalue,rtol=1e-10,atol=1e-15)
    assert np.allclose(t.p_holm_33,adjust(t.p),rtol=1e-10,atol=1e-12)
    checks.append(dict(check='independent_pairing_and_all_22_depth_tests',status='PASS'))
    checks.append(dict(check='Holm_family_all_33_tests',status='PASS'))
    signs=np.load(R/'bootstrap_profile_signs.npz')['signs'];assert signs.shape==(99999,17)
    assert set(np.unique(signs))=={-1,1}
    checks.append(dict(check='99999_draws_each_with_17_profile_weights',status='PASS'))
    perms=np.load(R/'permutation_profiles.npz')['indices'];assert perms.shape==(19999,17)
    assert (p.locality.to_numpy()[perms]==p.locality.to_numpy()).all()
    D=designs_['restoration'];pp=[]
    for j,k in enumerate(vals.index):pp.append(dict(KO=k,p_permutation_sensitivity=D.permutation_test(a[j].mean(axis=1),perms),n_profiles=17,n_permutations=19999))
    pd.DataFrame(pp).to_csv(R/'restoration_permutation_sensitivity.tsv',sep='\t',index=False)
    # All signed synthetic outcomes share one biological profile unit; no B-sized df.
    assert (t[t.test.eq('restoration_omnibus')].residual_df==10).all()
    checks.append(dict(check='model_df_uses_17_profiles_not_draw_count',status='PASS'))
    tests=pd.read_csv(R/'all_primary_and_sensitivity_tests.tsv',sep='\t');env=pd.read_csv(R/'environment_tests.tsv',sep='\t')
    for v,g in env.groupby('variant'):
        assert len(g)==19 and np.allclose(g.q_BH_19,adjust(g.p,'BH'))
    checks.append(dict(check='BH_family_all_19_prespecified_environment_pairs',status='PASS'))
    sens=tests[tests.test.eq('depth_40_vs_5')].groupby(['KO','gene']).agg(variants_tested=('variant','nunique'),variants_Holm_below_005=('p_holm_33',lambda x:int((x<.05).sum())),median_log2_min=('median_log2','min'),median_log2_max=('median_log2','max')).reset_index()
    loo=pd.read_csv(R/'leave_one_out_effects.tsv',sep='\t');lo=loo[loo.target.eq('depth_40_vs_5')].groupby('KO').estimate_log2.agg(loo_median_min='min',loo_median_max='max').reset_index();sens=sens.merge(lo,on='KO')
    sens.to_csv(R/'sensitivity_summary_40_vs_5.tsv',sep='\t',index=False)
    dep=t[t.test.str.startswith('depth')].copy()
    dep['median_ratio']=np.exp2(dep.median_log2);dep['ratio_ci_low']=np.exp2(dep.ci_low);dep['ratio_ci_high']=np.exp2(dep.ci_high);dep['percent_change']=100*(dep.median_ratio-1)
    dep.to_csv(R/'depth_effects_for_reporting.tsv',sep='\t',index=False)
    # Actual changes in the input matrices, not just stochastic resampling differences.
    changes=[]
    for name in ['primary_76','fractional_76','strict_reference_56','both_mates_76','supported_fractional_76']:
        f=R/(name+'_normalized.tsv')
        z=pd.read_csv(f,sep='\t',index_col=0).loc[vals.index,vals.columns]
        for k in vals.index:changes.append(dict(variant=f.stem.replace('_normalized',''),KO=k,max_absolute_relative_change=float(np.max(np.abs(z.loc[k]/vals.loc[k]-1)))))
    pd.DataFrame(changes).to_csv(R/'sensitivity_input_changes.tsv',sep='\t',index=False)
    pd.DataFrame(checks).to_csv(R/'verification_results.tsv',sep='\t',index=False)
    print(pd.DataFrame(checks).to_string(index=False));print(sens.to_string(index=False))

if __name__=='__main__':main()
