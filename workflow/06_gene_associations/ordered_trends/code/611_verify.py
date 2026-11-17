#!/usr/bin/env python3
"""Meaningful exact-test checks and cross-language verification of all outputs."""
from pathlib import Path
from itertools import permutations
from collections import Counter
import argparse,importlib.util,json,math
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent.parent
spec=importlib.util.spec_from_file_location('engine611',Path(__file__).with_name('611_monotonic.py'))
eng=importlib.util.module_from_spec(spec);spec.loader.exec_module(eng)

def independent_score(y,labels):
    groups=sorted(set(labels));terms=[]
    for a in groups:
        for b in groups:
            if a>=b:continue
            terms.append(np.mean([np.sign(y[j]-y[i]) for i in range(len(y)) for j in range(len(y)) if labels[i]==a and labels[j]==b]))
    return float(np.mean(terms))

def unit():
    evidence=[]
    for labels,values in [([0,1,2,2,3,4],[0,0,1,1,1,2]),([0,1,2,3,4],[1,4,2,0,3]),([0,1,2,2,3],[0,0,1,2,3]),([0,1,2,3],[1,1,1,1])]:
        labels=np.array(labels);y=np.array(values,float);k=max(labels)+1;scale=k*(k-1)
        q,d,s=eng.local_evidence(y[:,None],labels)
        assert np.isclose(q[0]/scale,independent_score(y,labels))
        brute=Counter(int(round(scale*independent_score(y,p))) for p in set(permutations(labels)))
        h=eng.null_hist(*d[0]);got={i-scale:int(n) for i,n in enumerate(h) if n}
        assert dict(brute)==got
        q2,d2,_=eng.local_evidence((y+1000)[:,None],labels);assert q2[0]==q[0]
        q3,_,_=eng.local_evidence((-y)[:,None],labels);assert q3[0]==-q[0]
    evidence.append('Independent group-pair score matches integer implementation, with ties, constants and 4/5 stages')
    evidence.append('Exact local distributions match brute label enumeration; locality shifts preserve score and value reversal reverses its sign')
    a=((1,1,2,1,1),(1,1,1,1,1,1));b=((1,1,1,1,1),(1,1,1,1,1))
    tail,n=eng.exact_lookup((a,a,b));assert n==15552000 and tail[60]==2/n
    h1=eng.null_hist(*a);h2=eng.null_hist(*b)
    enum=Counter()
    for i,ni in enumerate(h1):
        for j,nj in enumerate(h2):enum[i+j-40]+=int(ni*nj)
    p,n=eng.exact_lookup((a,b))
    for v in range(41):assert np.isclose(p[v],sum(c for q,c in enum.items() if abs(q)>=v)/43200)
    evidence.append('Convolved exact tails match independent two-block enumeration; perfect three-block trend has P=2/15552000')
    assert np.allclose(eng.adjusted([.001,.02,np.nan]),[.003,.03,np.nan],equal_nan=True)
    y=[3.,1.,2.,5.,4.];iso=eng.pava(y);assert np.allclose(iso,[2,2,2,4.5,4.5])
    evidence.append('Planned FDR denominator retains untestable functions; PAVA permits pooled plateaus')
    return evidence

def parity(py,rr):
    files={
      '611_main_results.tsv.gz':['ontology','feature_id','depth_cm'],
      '611_sensitivity_results.tsv.gz':['ontology','feature_id','depth_cm','variant'],
      '611_LOLO_results.tsv.gz':['ontology','feature_id','depth_cm','omitted_locality'],
      '611_sensitivity_LOLO.tsv.gz':['ontology','feature_id','depth_cm','omitted_locality','variant'],
      '611_stage_means.tsv.gz':['ontology','feature_id','depth_cm','stage_order'],
      '611_summary.tsv':['ontology','depth_cm']}
    result={}
    for name,keys in files.items():
        a=pd.read_csv(py/name,sep='\t').sort_values(keys).reset_index(drop=True)
        b=pd.read_csv(rr/name,sep='\t').sort_values(keys).reset_index(drop=True)
        assert set(a.columns)==set(b.columns),(name,'columns')
        assert len(a)==len(b);maxdiff=0
        for col in a.columns:
            x=a[col];y=b[col]
            if pd.api.types.is_numeric_dtype(x) and not pd.api.types.is_bool_dtype(x):
                assert np.allclose(x,y,rtol=2e-11,atol=2e-11,equal_nan=True),(name,col,np.nanmax(abs(x-y)))
                d=np.nanmax(abs(x-y)) if np.isfinite(x).any() else 0;maxdiff=max(maxdiff,float(d))
            else:assert x.fillna('NA').astype(str).equals(y.fillna('NA').astype(str)),(name,col)
        result[name]=dict(rows=len(a),maximum_absolute_numeric_difference=maxdiff,status='PASS')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--unit-only',action='store_true');a=p.parse_args()
    report={'status':'PASS','exact_test_checks':unit()}
    if not a.unit_only:report['Python_R']=parity(ROOT/'results_python',ROOT/'results_R')
    dest=ROOT/'verification'/('unit_tests.json' if a.unit_only else 'Python_R_validation.json');dest.parent.mkdir(exist_ok=True)
    dest.write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
