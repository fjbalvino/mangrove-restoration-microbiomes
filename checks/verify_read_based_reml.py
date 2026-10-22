#!/usr/bin/env python3
"""Refit the archived validation subset and compare with external lmerTest results.

This validates 714 model fits (1,785 contrasts), not a new full-catalogue run.
The archived 20 boundary-sensitive R refits replace only their df and P values.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();root=args.analysis_root;validation=root/'validation'
    spec=importlib.util.spec_from_file_location('read_based_reml',ROOT/'workflow/06_gene_associations/read_based/fit_random_intercept.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    metadata=pd.read_csv(root/'prepared/metadata.csv')
    responses=pd.read_csv(validation/'responses.csv').set_index('feature_index')
    selected=pd.read_csv(validation/'selected_validation_models.csv')
    reference=pd.read_csv(validation/'lmerTest_results.csv')
    stable=pd.read_csv(validation/'R_stable_derivatives.csv')
    keys=['axis','feature_index','model','contrast']
    if len(selected)!=357 or len(reference)!=1785 or selected.duplicated().any():
        raise AssertionError('Unexpected validation subset')
    if not np.allclose(responses.std(axis=1,ddof=1),1,rtol=1e-12):
        raise AssertionError('Validation responses must be standardized')
    records=[]
    for axis,selection in selected.groupby('axis',sort=False):
        indices=selection.feature_index.to_numpy()
        y=responses.loc[indices].to_numpy(float)
        for kind in ['common','interaction']:
            X,B=module.design(metadata,axis,kind=='interaction')
            L=np.zeros((1,6) if kind=='common' else (3,8));L[:,1]=1
            if kind=='interaction':L[1,6]=1;L[2,7]=1
            result=module.fit(y,X,B,L)
            for i,feature in enumerate(indices):
                for j,contrast in enumerate(['common'] if kind=='common' else ['slope5','slope20','slope40']):
                    row=dict(axis=axis,feature_index=feature,model=kind,contrast=contrast,
                             singular=bool(result['singular'][i]))
                    row.update({key:result[key][i,j] for key in ['estimate','se','df','p']})
                    records.append(row)
                if kind=='interaction':
                    records.append(dict(axis=axis,feature_index=feature,model=kind,contrast='joint',
                        singular=bool(result['singular'][i]),df=result['interaction_df'][i],p=result['interaction_p'][i]))
    current=pd.DataFrame(records)
    expected=reference.merge(stable[keys+['df','p']],on=keys,how='left',suffixes=('','_stable'),validate='one_to_one')
    expected['boundary_refit']=expected.df_stable.notna()
    for key in ['df','p']:expected[key]=expected[key+'_stable'].fillna(expected[key])
    joined=current.merge(expected,on=keys,suffixes=('_python','_R'),validate='one_to_one')
    if len(joined)!=1785:raise AssertionError('Incomplete contrast comparison')
    errors={key:float((joined[key+'_python']-joined[key+'_R']).abs().max()) for key in ['estimate','se','df','p']}
    boundary=joined.loc[joined.boundary_refit]
    stable_errors={key:float((boundary[key+'_python']-boundary[key+'_R']).abs().max()) for key in ['df','p']}
    mismatch=int(joined.singular_python.ne(joined.singular_R).sum())
    passed=(errors['estimate']<1e-7 and errors['se']<1e-7 and errors['df']<.001 and errors['p']<1e-6
            and stable_errors['df']<1e-5 and stable_errors['p']<1e-7 and mismatch==0)
    report={'status':'PASS' if passed else 'REVIEW','model_fits':714,'contrasts':len(joined),
        'singularity_classification_mismatches':mismatch,'max_absolute_errors':errors,
        'stable_R_max_absolute_errors':stable_errors,
        'scope':'Fresh Python subset fits compared with archived independent R validation; full models not refitted'}
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'reml_subset_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if not passed:raise SystemExit(1)

if __name__=='__main__':main()
