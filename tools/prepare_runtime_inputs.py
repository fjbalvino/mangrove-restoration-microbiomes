#!/usr/bin/env python3
"""Stage exact input copies at legacy bundle paths; do not run analyses."""
import argparse
import hashlib
import shutil
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--group',choices=['all','network','convergence'],default='all')
    group=parser.parse_args().group
    pairs=[]
    if group in {'all','network'}:
        base='workflow/07_networks/network809/inputs/'
        pairs += [('source_data/indices/metadata_003b.csv',base+'metadata.csv'),
                  ('source_data/functional/clr_200_KO.tsv.gz',base+'clr_200_KO.tsv.gz')]
    if group in {'all','convergence'}:
        base='workflow/05_functional_reference/convergence808/inputs/'
        pairs += [('source_data/indices/metadata_003b.csv',base+'metadata.csv'),
                  ('source_data/indices/convergence808_historical_metadata.csv',base+'historical_metadata.csv')]
        for layer in ['KEGG_ko','PFAM']:
            name=f'082A_{layer}_abundance_samples_x_functions.tsv.gz'
            pairs.append(('source_data/functional/tables/'+name,base+name))
    for source,target in pairs:
        src,dst=ROOT/source,ROOT/target
        expected=hashlib.sha256(src.read_bytes()).hexdigest()
        if dst.exists() and hashlib.sha256(dst.read_bytes()).hexdigest()!=expected:
            raise ValueError(f'Existing runtime input differs; inspect it before replacing: {dst}')
        dst.parent.mkdir(parents=True,exist_ok=True)
        if not dst.exists():shutil.copy2(src,dst)
        print(f'READY {target} SHA256={expected}')

if __name__=='__main__':main()
