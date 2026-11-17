from pathlib import Path
import re,json,datetime
import pandas as pd
P=Path(__file__).resolve().parent
maps={'map01200':'C','map00710':'C','map00720':'C','map00680':'CH4_C1','map00910':'N','map00920':'S','map00440':'P'}
links=pd.read_csv(P.parent/'panel/kegg_KO_pathway_links.tsv',sep='\t',names=['KO','map']);links.KO=links.KO.str.replace('ko:','',regex=False);links['map']=links['map'].str.replace('path:','',regex=False)
rows=[dict(KO=z.KO,cycle=maps[z.map],source=z.map) for z in links.itertuples() if z.map in maps]
process=pd.read_csv(P.parent/'panel/curated_processes.tsv',sep='\t').set_index('id');rules=pd.read_csv(P.parent/'panel/operational_rules.tsv',sep='\t')
for z in rules.itertuples():
 cycle=process.loc[z.process_id,'cycle']
 if cycle in ['C','CH4_C1','N','S','P']:
  rows.extend(dict(KO=k,cycle=cycle,source=z.rule_id) for k in sorted(set(re.findall(r'K\d{5}',z.expression))))
reg=pd.DataFrame(rows).drop_duplicates();reg.to_csv(P/'scope.tsv',sep='\t',index=False)
keys=sorted(reg.KO.unique());stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S');out=P/f'15_compare_expanded_{stamp}';out.mkdir();(P/'LATEST.txt').write_text(str(out));(out/'scope.json').write_text(json.dumps(dict(maps=maps,unique_KOs=len(keys),selection='All KOs on seven ecologically selected maps plus all operational curated rules for C,CH4_C1,N,S,P. Fixed before testing expanded outcomes; shared KOs counted once. Existing nitrogen analysis exploratory.'),indent=2))
src=(P.parent/'nitrogen_analysis/13_compare_nitrogen.py').read_text()
src=src.replace("O=P/'13_compare_nitrogen_results'",'O=out')
src=src.replace("keys=sorted(pd.read_csv(P.parent/'nitrogen_pilot/KO_sample_evidence.tsv',sep='\\t').KO.unique());",'keys=keys;')
src=src.replace("exploratory nitrogen map; no new pathway-specificity validation","exploratory expanded five-cycle map and curated-rule KO union; no new pathway-specificity validation")
# The same random sign vectors and methods preserve comparison with the nitrogen pilot.
src=src.replace("out=des.wild_test(y,signs);results.append", "testout=des.wild_test(y,signs);results.append").replace('**out));beta','**testout));beta')
src=src.replace("for k in eligible:\n", "for ki,k in enumerate(eligible):\n  if ki%25==0: print(vn,ki,'/',len(eligible),flush=True)\n")
exec(compile(src,str(P/'15_compare_expanded.py'),'exec'))
