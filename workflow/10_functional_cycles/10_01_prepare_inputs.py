#!/usr/bin/env python3
"""Extract the reviewed compact CNS input package from the corrected KO export."""
from pathlib import Path
import argparse
import json
import shutil
import sys
import pandas as pd

SECTION = Path(__file__).resolve().parent
REPO = SECTION.parents[1]
sys.path.insert(0,str(REPO/'checks'))
from check_functional_cycles_inputs import digest, verify, MANIFEST


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',required=True,type=Path,
                        help='Corrected export input directory containing matrices/, 03_unit_annotations/ and metadata_003b.csv.')
    parser.add_argument('--reference-list',required=True,type=Path,
                        help='Externally supplied official 76-KO reference list; exact hash is verified.')
    parser.add_argument('--output-dir',required=True,type=Path)
    args=parser.parse_args()
    manifest=json.loads(MANIFEST.read_text())
    if args.output_dir.exists():
        parser.error('Output already exists; choose a fresh directory.')
    ref=next(x for x in manifest['core'] if x['path']=='uscg_76_kegg_min_2313.lst')
    if not args.reference_list.is_file() or digest(args.reference_list)!=ref['sha256']:
        parser.error('Reference list differs from the reviewed 76-KO list.')
    for entry in manifest['full_source']:
        path=args.source/entry['path']
        if not path.is_file() or digest(path)!=entry['sha256']:
            parser.error('Full source is missing or changed: '+entry['path'])
    spec=json.loads((SECTION/'analysis_specification.json').read_text())
    selected=set(args.reference_list.read_text().split()+spec['candidate_primary_markers']+
                 spec['support_markers']+spec['descriptive_or_pending_markers'])
    args.output_dir.mkdir(parents=True,exist_ok=False)
    for entry in manifest['full_source']:
        if 'subset_path' not in entry:
            continue
        frame=pd.read_csv(args.source/entry['path'],sep='\t',index_col=0)
        frame.loc[frame.index.intersection(sorted(selected))].to_csv(
            args.output_dir/entry['subset_path'],sep='\t',float_format='%.15g')
    patterns=json.loads((args.source/'03_unit_annotations/annotation_patterns.json').read_text())['KO']
    kept=[p for p in patterns if set(p['labels'])&selected]
    (args.output_dir/'annotation_patterns_selected.json').write_text(json.dumps(kept,indent=2))
    shutil.copy2(args.source/'metadata_003b.csv',args.output_dir/'metadata_003b.csv')
    shutil.copy2(args.reference_list,args.output_dir/'uscg_76_kegg_min_2313.lst')
    report,_=verify(args.output_dir)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
