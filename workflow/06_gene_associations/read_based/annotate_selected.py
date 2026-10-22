from pathlib import Path
import numpy as np,pandas as pd,gzip,json
import argparse
parser=argparse.ArgumentParser(description='Attach historical KO/PFAM sets to the selected gene IDs.')
parser.add_argument('--input',type=Path,required=True)
parser.add_argument('--analysis-root',type=Path,required=True)
args=parser.parse_args();R=args.analysis_root.resolve();INP=args.input.resolve()
ids=np.load(R/'prepared/selected_gene_numeric_ids.npy');targets=set(map(int,ids));rows=[];focal=[]
old={'CZ5_2.contigs_935909','CZ1_3.contigs_1220634','CC1_3.contigs_1003408','CC2_3.contigs_980998','CZ3_2.contigs_1844408'}
with gzip.open(INP/'02_annotations/historical_gene_dictionary.tsv.gz','rt') as f:
    header=next(f).rstrip('\n').split('\t')
    for i,line in enumerate(f,1):
        s=line.rstrip('\n').split('\t');n=int(s[0])
        if n in targets:rows.append(s)
        if s[1] in old:focal.append(s)
        if i%5000000==0:print('DICTIONARY',i,'selected',len(rows),flush=True)
assert len(rows)==100000
a=pd.DataFrame(rows,columns=header);a['gene_numeric_id']=a.gene_numeric_id.astype(int)
sel=pd.read_csv(R/'prepared/selection_numeric.tsv',sep='\t').merge(a,on='gene_numeric_id',validate='one_to_one',sort=False)
sel.to_csv(R/'prepared/selected_genes.tsv',sep='\t',index=False)
pd.DataFrame(focal,columns=header).to_csv(R/'prepared/previous_candidates.tsv',sep='\t',index=False)
print('ANNOTATION_PASS',len(rows),'dictionary_genes',i,'focal',focal,flush=True)
