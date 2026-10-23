"""Readable cycle/function labels; annotations describe potential, not rates."""
from pathlib import Path
import pandas as pd
ORDER=['rbcL','acsB','nifH','nirK','nirS','norB','nosZ','nrfA','dsrA','aprA','soxB']
PROVISIONAL={'nifH','nirK','nosZ','soxB','acsB'}
DATA=[
 ('rbcL','K01601','Carbon','CO₂ incorporation\n(RuBisCO)','RuBisCO catalyses CO2 incorporation. This KO also occurs outside a complete Calvin cycle; taxonomic origin and enzyme form were not resolved.'),
 ('acsB','K14138','Carbon','Acetyl-CoA synthesis','Associated with the Wood–Ljungdahl pathway. The marker alone establishes neither a complete pathway nor its net direction.'),
 ('nifH','K02588','Nitrogen','Nitrogen fixation','Nitrogenase iron-protein candidate; assignment does not establish a complete nitrogenase or exclude related sequences.'),
 ('nirK','K00368','Nitrogen','Nitrite → nitric oxide','NO-forming nitrite reduction, a step associated with denitrification and some other nitrogen pathways; not proof of complete denitrification.'),
 ('nirS','K15864','Nitrogen','Nitrite / hydroxylamine\nreduction','KEGG includes both NO-forming nitrite reductase and hydroxylamine reductase in K15864. The KO-level signal is not exclusive to denitrification.'),
 ('norB','K04561','Nitrogen','Nitric oxide → N₂O','Nitric oxide reductase subunit; a step associated with denitrification, not a measurement of N2O production.'),
 ('nosZ','K00376','Nitrogen','N₂O → nitrogen gas','Nitrous oxide reduction; the marker alone does not establish that its carrier performs all other denitrification steps.'),
 ('nrfA','K03385','Nitrogen','Nitrite → ammonium','Nitrite reduction to ammonium, associated with dissimilatory nitrate reduction to ammonium (DNRA).'),
 ('dsrA','K11180','Sulfur','Sulfite transformations','Dissimilatory sulfite reductase alpha subunit. Dsr systems occur in sulfur reduction and oxidation; direction is unresolved from this KO alone.'),
 ('aprA','K00394','Sulfur','Activated sulfate\ntransformations','Adenylylsulfate (APS) reductase subunit A, associated with APS/sulfite transformations. Net sulfur reduction versus oxidation is unresolved.'),
 ('soxB','K17224','Sulfur','Thiosulfate oxidation','Sox-system sulfohydrolase associated with thiosulfate oxidation; does not establish complete oxidation to sulfate by itself.')]
ANNOTATIONS={g:dict(gene=g,KO=k,cycle=c,label=f,limitation=n,source_url='https://www.kegg.jp/entry/'+k) for g,k,c,f,n in DATA}
def export_guide(root):
 rows=[]
 for gene in ORDER:
  row=ANNOTATIONS[gene].copy();row['label']=row['label'].replace('\n',' ')
  row['allocation_status']='provisional_fractional' if gene in PROVISIONAL else 'recoverable_single_KO'
  row['sequence_specificity']='not_validated';rows.append(row)
 pd.DataFrame(rows).to_csv(Path(root)/'results/marker_function_guide.tsv',sep='\t',index=False)
def display_rows():
 rows=[]
 for cycle in ['Carbon','Nitrogen','Sulfur']:
  rows.append(('header',cycle))
  rows.extend(('gene',g) for g in ORDER if ANNOTATIONS[g]['cycle']==cycle)
 return rows

