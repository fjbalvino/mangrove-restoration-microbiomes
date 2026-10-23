#!/usr/bin/env python3
"""Verify external CNS inputs without publishing or downloading their contents."""
from pathlib import Path
import argparse
import hashlib
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'docs/FUNCTIONAL_CYCLES_INPUTS.json'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def validate_cohort(meta):
    """Fail explicitly on duplicated samples, depth units or inconsistent profiles."""
    required = {'sample_id', 'lat_block', 'depth_cm', 'locality', 'restoration4',
                'latitude', 'longitude', 'redox_shallow_mv', 'n_no3_umol_l',
                'n_nh4plus_umol_l', 's_2_umol'}
    if not required.issubset(meta.columns):
        raise ValueError('Missing metadata columns: ' + ', '.join(sorted(required-set(meta.columns))))
    if len(meta) != 51 or meta.sample_id.isna().any() or not meta.sample_id.is_unique:
        raise ValueError('Expected exactly 51 unique nonmissing sample IDs')
    meta = meta.copy()
    meta['profile_id'] = meta.sample_id.str.rsplit('_', n=1).str[0]
    if meta.profile_id.nunique() != 17 or meta.lat_block.nunique() != 17:
        raise ValueError('Expected 17 complete profiles')
    if meta.duplicated(['profile_id', 'depth_cm']).any():
        raise ValueError('Duplicate profile-depth observation')
    expected_stages = {'Degraded', 'Early restoration', 'Intermediate restoration',
                       'Advanced restoration', 'Conserved'}
    if set(meta.restoration4) != expected_stages or set(meta.locality) != {'Carmen', 'Cozumel', 'Tuxpan'}:
        raise ValueError('Unexpected locality or restoration labels')
    for _, block in meta.groupby('profile_id'):
        if sorted(block.depth_cm.tolist()) != [5, 20, 40]:
            raise ValueError('Every profile must contain one 5, 20 and 40 cm observation')
        for column in ['lat_block', 'locality', 'restoration4', 'latitude', 'longitude']:
            if block[column].isna().any() or block[column].nunique() != 1:
                raise ValueError('Inconsistent profile identity: ' + column)
    if not meta.groupby('lat_block').profile_id.nunique().eq(1).all():
        raise ValueError('Profile IDs and lat_block are not one-to-one')
    return meta


def verify(input_dir, with_previous_comparison=False):
    input_dir = Path(input_dir).resolve()
    manifest = json.loads(MANIFEST.read_text())
    entries = list(manifest['core'])
    if with_previous_comparison:
        entries += manifest['previous_comparison']
    errors = []
    for item in entries:
        path = input_dir / item['path']
        if not path.is_file():
            errors.append('Missing external input: ' + item['path'])
        elif path.stat().st_size != item['bytes'] or digest(path) != item['sha256']:
            errors.append('External input differs from the reviewed release: ' + item['path'])
    if errors:
        raise ValueError('\n'.join(errors))
    meta = validate_cohort(pd.read_csv(input_dir/'metadata_003b.csv'))
    report = {'status': 'PASS', 'files_verified': len(entries), 'n_samples': len(meta),
              'n_profiles': meta.profile_id.nunique(), 'with_previous_comparison': with_previous_comparison}
    return report, entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--with-previous-comparison', action='store_true')
    args = parser.parse_args()
    try:
        report, _ = verify(args.input_dir, args.with_previous_comparison)
    except ValueError as error:
        parser.exit(1, str(error)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
