#!/usr/bin/env python3
"""Check declared multiplicity families, published totals and optional archive agreement."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

def require(condition,message):
    if not condition:raise AssertionError(message)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,required=True)
    parser.add_argument('--reference-results',type=Path,help='Optional external verified result directory')
    args=parser.parse_args();root=args.results
    terms=pd.read_csv(root/'primary_enriched_terms_named.tsv',sep='\t')
    modules=pd.read_csv(root/'primary_enriched_modules.tsv',sep='\t')
    all_terms=pd.read_csv(root/'all_ORA_tests.tsv.gz',sep='\t')
    all_modules=pd.read_csv(root/'all_module_ORA_tests.tsv.gz',sep='\t')
    require(all_terms.groupby('gene_q_cutoff').size().to_dict()=={.05:42770,.1:42770},'KO/PFAM family size changed')
    require(all_modules.groupby('gene_q_cutoff').size().to_dict()=={.05:2606,.1:2606},'Module family size changed')
    require(len(terms)==102 and len(terms[['ontology','term']].drop_duplicates())==96,'Primary term counts changed')
    require(terms.ontology.value_counts().to_dict()=={'PFAM':87,'KO':15},'Ontology counts changed')
    require(int(terms.retained_gene_q005.sum())==45,'Gene-q sensitivity changed')
    require(len(modules)==11 and modules.module.nunique()==10,'Candidate-module counts changed')
    require(int(modules.retained_gene_q005.sum())==1,'Module gene-q sensitivity changed')
    require(int(terms.loc[terms.ontology.eq('KO'),'retained_single_KO'].eq(True).sum())==11,'Single-KO term sensitivity changed')
    require(not modules.retained_single_KO.eq(True).any(),'Single-KO module result changed')
    require(terms.q_ora_global.le(.1).all() and terms.selected_term_genes.ge(5).all(),'Invalid primary term filter')
    require(modules.q_ora_modules.le(.1).all() and modules.selected_term_genes.ge(5).all()
            and modules.selected_distinct_KO.ge(2).all(),'Invalid primary module filter')
    require(terms.groupby('axis').size().to_dict()=={
        'water_inundation_PC1':11,'moisture_stress_PC1':44,
        'physicochemical_PC1':5,'nutrients_redox_PC1':42},'Axis contrast counts changed')
    compared=[]
    if args.reference_results:
        for sub in ['', 'single_KO_functions', 'single_KO_modules']:
            names=['all_module_ORA_tests.tsv.gz'] if sub=='single_KO_modules' else ['all_ORA_tests.tsv.gz']
            if not sub:names+=['all_module_ORA_tests.tsv.gz','primary_enriched_terms_named.tsv','primary_enriched_modules.tsv']
            for name in names:
                a=pd.read_csv(root/sub/name,sep='\t')
                b=pd.read_csv(args.reference_results/sub/name,sep='\t')
                keys=['axis','direction','gene_q_cutoff']+(['ontology','term'] if 'term' in a.columns else ['module'])
                a=a.sort_values(keys).reset_index(drop=True).sort_index(axis=1)
                b=b.sort_values(keys).reset_index(drop=True).sort_index(axis=1)
                pd.testing.assert_frame_equal(a,b,check_dtype=False,rtol=1e-12,atol=1e-14)
                compared.append(str(Path(sub)/name))
    report={'status':'PASS','function_tests_per_cutoff':42770,'module_tests_per_cutoff':2606,
        'function_contrasts':len(terms),'candidate_module_contrasts':len(modules),
        'external_result_tables_compared':compared,
        'scope':'Enrichment output only; no new full-catalogue reconstruction or full model refit is certified'}
    (root/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
