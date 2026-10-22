#!/usr/bin/env python3
"""Check the external files used for the current Figure 4A enrichment.

Checksums pin the archived model tables and annotation snapshots. They do not
download databases or certify independently refitted models.
"""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-root',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    args=parser.parse_args()
    roots={'analysis':args.analysis_root,'reference':args.reference}
    manifest=json.loads((ROOT/'docs/ENVIRONMENT_ENRICHMENT_INPUTS.json').read_text())
    failures=[]
    for entry in manifest['files']:
        path=roots[entry['root']]/entry['path']
        if not path.is_file():
            failures.append(f'Missing external input: {path}')
            continue
        with path.open('rb') as stream:
            digest=hashlib.file_digest(stream,'sha256').hexdigest()
        if path.stat().st_size!=entry['bytes'] or digest!=entry['sha256']:
            failures.append(f'Input differs from the verified snapshot: {path}')
    metadata=args.analysis_root/'prepared/metadata.csv'
    if not metadata.is_file():
        failures.append(f'Missing external metadata: {metadata}')
    else:
        m=pd.read_csv(metadata)
        required=['sample_id','profile_id','depth_cm','locality','restoration4']
        if not set(required).issubset(m.columns):
            failures.append('Metadata is missing cohort identifiers')
        else:
            valid=(len(m)==51 and m.sample_id.is_unique and m.profile_id.nunique()==17)
            valid=valid and all(sorted(g.depth_cm)==[5,20,40] and
                g.locality.nunique()==g.restoration4.nunique()==1 for _,g in m.groupby('profile_id'))
            if not valid:failures.append('Metadata does not describe 17 complete three-depth profiles')
    print(json.dumps({'status':'FAIL' if failures else 'PASS','files_checked':len(manifest['files']),
                      'failures':failures},indent=2))
    if failures:raise SystemExit(1)

if __name__=='__main__':main()
