#!/usr/bin/env python3
"""Run the paired nitrogen-depth pilot from hash-verified external inputs."""
import argparse,datetime,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
P=Path(__file__).resolve().parent/'nitrogen_depth'
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--pilot-dir',type=Path,required=True);ap.add_argument('--out-root',type=Path,default=Path('local_runs/nitrogen_depth'));ap.add_argument('--check-only',action='store_true');a=ap.parse_args();source=a.pilot_dir.resolve();manifest=json.loads((P/'inputs.json').read_text())
 for row in manifest:
  f=source/row['path']
  if not f.is_file():raise SystemExit('Missing '+str(f))
  if hashlib.sha256(f.read_bytes()).hexdigest()!=row['sha256']:raise SystemExit('Hash mismatch '+str(f))
 print('Verified',len(manifest),'input files',flush=True)
 if a.check_only:return
 run=a.out_root.resolve()/('10_02_run_nitrogen_depth_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S'));run.mkdir(parents=True,exist_ok=False)
 for row in manifest:(run/row['path']).symlink_to(source/row['path'])
 for f in (P/'code').glob('*.py'):shutil.copy2(f,run/f.name)
 (run/'launch.json').write_text(json.dumps({'source':str(source),'input_manifest':manifest,'code_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (P/'code').glob('*.py')}},indent=2))
 subprocess.run([sys.executable,str(run/'19_nitrogen_depth.py')],cwd=run,env=dict(os.environ,OPENBLAS_NUM_THREADS='1'),check=True)
 print('Results:',(run/'LATEST.txt').read_text())
if __name__=='__main__':main()
