#!/usr/bin/env python3
"""Read-only, bounded collection of explicitly inventoried analytical inputs.

Inputs: docs/MISSING_INPUTS.tsv and the core/result files currently listed there.
Outputs: timestamped collection directory, manifest.json and .tar.gz.
No analyses, R deserialisation, deletion or LATEST mutation.
AI-assisted implementation: OpenAI Codex (OpenAI, 2026).
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import tarfile
from datetime import datetime, timezone

def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser()
    ap.add_argument('--group',choices=['core'],default='core')
    ap.add_argument('--manifest',type=Path,default=root/'docs/MISSING_INPUTS.tsv')
    ap.add_argument('--out-root',type=Path,default=Path('/home/fjbalvino/Tipping_points/auditoria_github'))
    args=ap.parse_args()
    records=[r for r in csv.DictReader(args.manifest.open(),delimiter='\t') if r['collection_group']==args.group]
    if not records: raise ValueError('No files in requested collection group')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    dest=args.out_root/f'GITHUB_INPUTS_{args.group}_{stamp}_{os.getpid()}'
    dest.mkdir(parents=True,exist_ok=False)
    (dest/'files').mkdir()
    report=[]
    for i,r in enumerate(records):
        p=Path(r['server_path']);item=dict(r)
        try:
            before=p.stat()
            if not p.is_file(): raise ValueError('Input is not a regular file')
            # Deliberately copy only explicit files, never traverse directories.
            out=dest/'files'/f'{i:02d}_{p.name}'
            digest=hashlib.sha256()
            with p.open('rb') as source,out.open('xb') as target:
                while True:
                    block=source.read(1024*1024)
                    if not block: break
                    target.write(block);digest.update(block)
            after=p.stat()
            if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
                raise ValueError('Source changed during copy')
            item.update(status='COPIED',collected_path=str(out.relative_to(dest)),observed_bytes=out.stat().st_size,collected_sha256=digest.hexdigest(),same_size_as_inventory=(not r['size_bytes'] or before.st_size==int(r['size_bytes'])))
            if r['sha256'] and r['sha256']!=digest.hexdigest():item['status']='HASH_MISMATCH'
        except Exception as error:
            item.update(status='FAILED',error=str(error))
        report.append(item)
        print(f'{i+1}/{len(records)} {p.name}: {item["status"]}',flush=True)
    (dest/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    archive=dest.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz',compresslevel=1) as tar:tar.add(dest,arcname=dest.name)
    digest=hashlib.sha256()
    with archive.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):digest.update(block)
    print(f'ARCHIVE_TO_SHARE={archive}\nSHA256={digest.hexdigest()}',flush=True)
    return 0 if all(r['status']=='COPIED' and r['same_size_as_inventory'] for r in report) else 1

if __name__=='__main__':raise SystemExit(main())
