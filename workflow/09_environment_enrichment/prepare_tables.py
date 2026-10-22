"""Join reference labels, annotation sensitivities and auditable gene memberships."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--results',type=Path,required=True)
parser.add_argument('--reference',type=Path,required=True)
args=parser.parse_args();R=args.results.resolve();reference=args.reference.resolve()
t=pd.read_csv(R/'primary_enriched_terms.tsv',sep='\t')
ko_names=pd.read_csv(reference/'kegg_ko_names.tsv',sep='\t',header=None,names=['term','name'])
t=t.merge(ko_names,on='term',how='left',validate='many_to_one')
t['name']=t.name.fillna(t.term)
pfam_names=pd.read_csv(reference/'Pfam-A.clans.tsv.gz',sep='\t',header=None,
    names=['accession','clan','clan_name','term','description'])
t=t.merge(pfam_names[['term','accession','description']],on='term',how='left',validate='many_to_one')
t['name']=t.description.fillna(t.name)
alltests=pd.read_csv(R/'all_ORA_tests.tsv.gz',sep='\t')
sens=alltests.loc[alltests.gene_q_cutoff.eq(.05),['axis','ontology','direction','term','selected_term_genes','q_ora_global']].rename(columns={'selected_term_genes':'selected_term_genes_q005','q_ora_global':'q_ora_gene_q005'})
t=t.merge(sens,on=['axis','ontology','direction','term'],validate='one_to_one')
single=pd.read_csv(R/'single_KO_functions/all_ORA_tests.tsv.gz',sep='\t')
single=single.loc[single.gene_q_cutoff.eq(.1)&single.ontology.eq('KO'),
    ['axis','ontology','direction','term','q_ora_global','enriched_exploratory']].rename(columns={
        'q_ora_global':'q_ora_single_KO','enriched_exploratory':'retained_single_KO'})
t=t.merge(single,on=['axis','ontology','direction','term'],how='left',validate='one_to_one')
t.to_csv(R/'primary_enriched_terms_named.tsv',sep='\t',index=False)
m=pd.read_csv(R/'primary_enriched_modules.tsv',sep='\t').drop(columns=[
    'selected_term_genes_q005','selected_distinct_KO_q005','q_ora_gene_q005'],errors='ignore')
mall=pd.read_csv(R/'all_module_ORA_tests.tsv.gz',sep='\t')
ms=mall.loc[mall.gene_q_cutoff.eq(.05),['axis','direction','module','selected_term_genes','selected_distinct_KO','q_ora_modules']].rename(columns={'selected_term_genes':'selected_term_genes_q005','selected_distinct_KO':'selected_distinct_KO_q005','q_ora_modules':'q_ora_gene_q005'})
m=m.merge(ms,on=['axis','direction','module'],validate='one_to_one')
singlem=pd.read_csv(R/'single_KO_modules/all_module_ORA_tests.tsv.gz',sep='\t')
singlem=singlem.loc[singlem.gene_q_cutoff.eq(.1),['axis','direction','module','q_ora_modules','enriched_exploratory']].rename(columns={
    'q_ora_modules':'q_ora_single_KO','enriched_exploratory':'retained_single_KO'})
m=m.drop(columns=['q_ora_single_KO','retained_single_KO'],errors='ignore').merge(singlem,on=['axis','direction','module'],how='left',validate='one_to_one')
m.to_csv(R/'primary_enriched_modules.tsv',sep='\t',index=False)
d=pd.read_csv(R/'associated_genes_Figure4A.tsv.gz',sep='\t')
members=[]
for ont in ['KO','PFAM']:
    a=d[['axis','direction','gene_numeric_id','feature_id','estimate_sd','q',ont+'_set']].copy()
    a['term']=a[ont+'_set'].str.split(',');a=a.explode('term').drop(columns=ont+'_set').drop_duplicates()
    a=a.merge(t.loc[t.ontology.eq(ont),['axis','direction','term']],on=['axis','direction','term'],validate='many_to_one')
    a['ontology']=ont;members.append(a)
membership=pd.concat(members,ignore_index=True)
membership.to_csv(R/'genes_supporting_enriched_terms.tsv.gz',sep='\t',index=False)
assert membership.groupby(['axis','ontology','direction','term']).gene_numeric_id.nunique().sort_index().equals(t.set_index(['axis','ontology','direction','term']).selected_term_genes.sort_index())
links=pd.read_csv(reference/'kegg_ko_module_links.tsv',sep='\t',header=None,names=['KO','module'])
links.KO=links.KO.str.replace('ko:','',regex=False);links.module=links.module.str.replace('md:','',regex=False)
ma=d[['axis','direction','gene_numeric_id','feature_id','estimate_sd','q','KO_set']].copy()
ma['KO']=ma.KO_set.str.split(',');ma=ma.explode('KO').merge(links,on='KO')
ma=ma.merge(m[['axis','direction','module']],on=['axis','direction','module'],validate='many_to_one').drop_duplicates()
ma.to_csv(R/'genes_supporting_enriched_modules.tsv.gz',sep='\t',index=False)
assert ma.groupby(['axis','direction','module']).gene_numeric_id.nunique().sort_index().equals(m.set_index(['axis','direction','module']).selected_term_genes.sort_index())

# Gene-level annotation ambiguity is exposed, not resolved by inventing new assignments.
qc=d.assign(multi_KO=d.KO_set.str.contains(',',na=False)).groupby('axis').agg(associated_genes=('gene_numeric_id','size'),multi_KO_genes=('multi_KO','sum'))
qc.to_csv(R/'annotation_ambiguity.tsv',sep='\t')

audit={'source_archives':{'Matrices_genicas_Figura4.zip':'5580e569dc7e0f1546a6c7c83ff0af3a0889efc9227ec5ae215da71c96b96d83','Analisis_Figura4.zip':'dd6eceafdd2349b59fc89eba8f6a6b6e73388c061bce0e446315fa3c9208af2f'},
'primary_function_axis_sign_contrasts':len(t),'distinct_primary_terms':int(t[['ontology','term']].drop_duplicates().shape[0]),
'primary_retained_at_gene_q005':int(t.retained_gene_q005.sum()),'module_contrasts':len(m),
'module_retained_at_gene_q005':int(m.retained_gene_q005.sum()),
'primary_KO_contrasts_retained_single_KO':int(t.loc[t.ontology.eq('KO'),'retained_single_KO'].eq(True).sum()),
'module_contrasts_retained_single_KO':int(m.retained_single_KO.eq(True).sum()),
'associated_gene_axis_rows':len(d),'membership_count_checks':'PASS',
'scope':'Exploratory gene-environment functional enrichment; modules are potential candidates'}
(R/'audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(audit,indent=2,ensure_ascii=False))
