#!/usr/bin/env python3
"""LOLO for every prespecified sensitivity, using the same discovery families."""
from pathlib import Path
import argparse,importlib.util
import pandas as pd

spec=importlib.util.spec_from_file_location('engine611',Path(__file__).with_name('611_monotonic.py'))
eng=importlib.util.module_from_spec(spec);spec.loader.exec_module(eng)

def run(inputs,results):
    inputs=Path(inputs);results=Path(results);meta=eng.audit_inputs(inputs);out=[]
    for ont in ['KO','PFAM']:
        df=pd.read_csv(inputs/(ont+'.tsv.gz'),sep='\t',index_col=0).loc[meta.sample_id]
        for variant,delta,inc in eng.VARIANTS[1:]:
            rows=[]
            for depth in eng.DEPTHS:
                print('SENSITIVITY LOLO',ont,variant,depth,flush=True)
                # The main mode returns detailed folds; delta and stage selection
                # define the sensitivity. No inferential selection occurs here.
                _,fold,_=eng.analyze_slice(meta,df.to_numpy(float),df.columns.to_numpy(),ont,depth,'main_CLR',delta,inc)
                fold['variant']=variant;rows.append(fold)
            tab=pd.concat(rows,ignore_index=True)
            for loc,ix in tab.groupby('omitted_locality').groups.items():
                tab.loc[ix,'q_training_BH']=eng.adjusted(tab.loc[ix,'p_training'])
                tab.loc[ix,'q_training_BY']=eng.adjusted(tab.loc[ix,'p_training'],'BY')
            tab['selected_training']=tab.q_training_BH.le(.05)
            tab['selected_and_transferred']=tab.selected_training&tab.heldout_matches_training
            out.append(tab)
    eng.write_table(pd.concat(out,ignore_index=True),results/'611_sensitivity_LOLO.tsv.gz')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',required=True);p.add_argument('--results',required=True)
    a=p.parse_args();run(a.inputs,a.results)
