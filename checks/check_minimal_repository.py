#!/usr/bin/env python3
"""Check shipped input hashes and Python/Bash syntax; no model fitting."""
import ast,csv,hashlib,json,subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def main():
    failures=[];checked=0
    with (ROOT/'docs/INPUTS_MINIMAL.tsv').open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            path=ROOT/row['path'];checked+=1
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:
                failures.append('Input missing or modified: '+row['path'])
    for folder in ['workflow','checks','tools']:
        for path in (ROOT/folder).rglob('*'):
            if not path.is_file():continue
            if path.suffix=='.py':
                checked+=1
                try:ast.parse(path.read_text())
                except SyntaxError as error:failures.append(f'{path.relative_to(ROOT)}: {error}')
            elif path.suffix=='.sh':
                checked+=1
                result=subprocess.run(['bash','-n',str(path)],capture_output=True,text=True)
                if result.returncode:failures.append(result.stderr)
    print(json.dumps({'status':'FAIL' if failures else 'PASS','checks':checked,'failures':failures,'scope':'Shipped inputs and syntax only; missing server inputs and R execution are not validated.'},indent=2))
    if failures:raise SystemExit(1)

if __name__=='__main__':main()
