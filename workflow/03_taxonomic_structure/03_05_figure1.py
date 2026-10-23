#!/usr/bin/env python3
"""Validate original taxonomic inputs and regenerate Figure 1 panels B-D locally."""
import argparse
import csv
import hashlib
from pathlib import Path
import subprocess
import sys

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',required=True)
    parser.add_argument('--output-dir',required=True)
    args=parser.parse_args()
    here=Path(__file__).resolve().parent/'figure1'
    inputs=Path(args.input_dir).resolve();output=Path(args.output_dir).resolve()
    if output.exists():raise SystemExit(f'Refusing to overwrite existing output: {output}')
    with (here/'inputs_required.tsv').open() as f:
        expected=list(csv.DictReader(f,delimiter='\t'))
    errors=[]
    for item in expected:
        p=inputs/item['relative_path']
        if not p.is_file():errors.append(f'Missing {p}');continue
        payload=p.read_bytes()
        if len(payload)!=int(item['size_bytes']) or hashlib.sha256(payload).hexdigest()!=item['sha256']:
            errors.append(f'Input does not match verified version: {p}')
    if errors:raise SystemExit('\n'.join(errors))
    output.mkdir(parents=True)
    print(f'Validated {len(expected)} external files',flush=True)
    command=[sys.executable,str(here/'audit_figure1.py'),
             '--phyloseq',str(inputs/'001_phyloseq_canon_51.rds'),
             '--metadata',str(inputs/'metadata_003b.csv'),
             '--archive',str(inputs/'archived_tables'),
             '--raw-metrics',str(inputs/'201_metricas_bimodalidad_todas_especies.csv'),
             '--membership',str(inputs/'202_membresia_universos_84_12_3.csv'),
             '--output',str(output/'audit')]
    subprocess.run(command,check=True)
    command=[sys.executable,str(here/'render_figure1.py'),'--audit',str(output/'audit'),
             '--output',str(output/'figures')]
    subprocess.run(command,check=True)
    print(f'Completed Figure 1 panels B-D: {output}',flush=True)

if __name__=='__main__':main()

