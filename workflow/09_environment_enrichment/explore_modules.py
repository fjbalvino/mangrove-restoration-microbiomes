"""Secondary KEGG-module ORA; separate multiplicity family from KO/PFAM."""
from pathlib import Path
import argparse,json
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from explore_enrichment import AXES,bh

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source',type=Path,required=True)
parser.add_argument('--reference',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--single-ko',action='store_true',help='Exclude multi-KO clusters from both selected and background genes')
args=parser.parse_args()
source=args.source.resolve();reference=args.reference.resolve();out=args.output.resolve()
out.mkdir(parents=True,exist_ok=True)
genes=pd.read_csv(source/'prepared/selected_genes.tsv',sep='\t',keep_default_na=False)
if args.single_ko:
    genes=genes.loc[~genes.KO_set.str.contains(',')].copy()
gk=genes[['gene_numeric_id','KO_set']].assign(KO=lambda d:d.KO_set.str.split(',')).explode('KO')
links=pd.read_csv(reference/'kegg_ko_module_links.tsv',sep='\t',header=None,names=['KO','module'])
links['KO']=links.KO.str.replace('ko:','',regex=False)
links['module']=links.module.str.replace('md:','',regex=False)
gm=gk[['gene_numeric_id','KO']].merge(links,on='KO').drop_duplicates()
names=pd.read_csv(reference/'kegg_module_names.tsv',sep='\t',header=None,names=['module','name'])
rows=[]
for axis in AXES:
    common=pd.concat([x.loc[x.family.eq('common')] for x in pd.read_csv(source/'models'/f'{axis}_all_results.tsv.gz',sep='\t',chunksize=100000)],ignore_index=True)
    common=common.loc[common.usable & ~common.singular]
    b=gm.loc[gm.gene_numeric_id.isin(common.gene_numeric_id)]
    ann=b[['gene_numeric_id','module']].drop_duplicates()
    N=ann.gene_numeric_id.nunique()
    K=ann.groupby('module').gene_numeric_id.nunique()
    kb=b.groupby('module').KO.nunique()
    K=K[K.ge(10) & kb.ge(2)]
    for cut in [.1,.05]:
        for direction,sign in [('negative',-1),('positive',1)]:
            ids=common.loc[common.q.le(cut)&common.estimate_sd.mul(sign).gt(0),'gene_numeric_id']
            n=ann.loc[ann.gene_numeric_id.isin(ids),'gene_numeric_id'].nunique()
            selected=b.loc[b.gene_numeric_id.isin(ids)]
            a=selected.groupby('module').gene_numeric_id.nunique().reindex(K.index,fill_value=0)
            ka=selected.groupby('module').KO.nunique().reindex(K.index,fill_value=0)
            koids=selected.groupby('module').KO.agg(lambda s:','.join(sorted(s.unique()))).reindex(K.index,fill_value='')
            rows.append(pd.DataFrame({'axis':axis,'direction':direction,'gene_q_cutoff':cut,'module':K.index,
                'background_genes':N,'selected_genes':n,'background_term_genes':K.values,'selected_term_genes':a.values,
                'background_distinct_KO':kb.reindex(K.index).values,'selected_distinct_KO':ka.values,'selected_KO_list':koids.values,
                'fold_enrichment':(a.values/n)/(K.values/N) if n else np.nan,
                'p_fisher_greater':hypergeom.sf(a.values-1,N,K.values,n) if n else np.ones(len(K))}))
r=pd.concat(rows,ignore_index=True).merge(names,on='module',how='left',validate='many_to_one')
for _,ix in r.groupby('gene_q_cutoff').groups.items():r.loc[ix,'q_ora_modules']=bh(r.loc[ix,'p_fisher_greater'])
r['enriched_exploratory']=r.q_ora_modules.le(.1)&r.selected_term_genes.ge(5)&r.selected_distinct_KO.ge(2)
r.to_csv(out/'all_module_ORA_tests.tsv.gz',sep='\t',index=False)
s=r.loc[r.enriched_exploratory].sort_values(['gene_q_cutoff','axis','q_ora_modules'])
s.to_csv(out/'enriched_modules.tsv',sep='\t',index=False)
primary=s.loc[s.gene_q_cutoff.eq(.1)].copy()
sens=s.loc[s.gene_q_cutoff.eq(.05),['axis','direction','module']].assign(retained_gene_q005=True)
primary=primary.merge(sens,how='left',on=['axis','direction','module'],validate='one_to_one')
primary['retained_gene_q005']=primary.retained_gene_q005.eq(True)
primary.to_csv(out/'primary_enriched_modules.tsv',sep='\t',index=False)
print(primary.to_string(index=False),flush=True)
print('Module tests per cutoff:',r.groupby('gene_q_cutoff').size().to_dict(),flush=True)
(out/'module_parameters.json').write_text(json.dumps({'single_KO_only':args.single_ko,
    'tests_per_cutoff':r.groupby('gene_q_cutoff').size().to_dict(),
    'background_min_genes':10,'background_min_KO':2,'selected_min_genes':5,'selected_min_KO':2,
    'ORA_BH_q':0.1,'gene_selection_q':[0.1,0.05],
    'multiplicity':'all five axes and both signs, separately for each gene-q cutoff'},indent=2)+'\n')
