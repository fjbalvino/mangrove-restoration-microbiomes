#!/usr/bin/env python3
"""Stage verified external inputs and run the preserved expanded-panel producer."""
import argparse,datetime,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
P=Path(__file__).resolve().parent

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--bundle-root',type=Path,required=True,help='reproducibilidad directory extracted from the expanded-panel external archive')
 ap.add_argument('--out-root',type=Path,default=Path('local_runs/expanded_cycles'))
 ap.add_argument('--check-only',action='store_true')
 ap.add_argument('--prepare-only',action='store_true')
 a=ap.parse_args();source=a.bundle_root.resolve();manifest=json.loads((P/'expanded/inputs.json').read_text());fail=[]
 for item in manifest:
  f=source/item['path']
  if not f.is_file():fail.append('Missing: '+item['path'])
  elif hashlib.sha256(f.read_bytes()).hexdigest()!=item['sha256']:fail.append('Modified: '+item['path'])
 if fail:raise SystemExit('\n'.join(fail))
 import pandas as pd
 m=pd.read_csv(source/'panel/data/metadata_003b.csv');m['profile_id']=m.sample_id.str.rsplit('_',n=1).str[0]
 assert len(m)==51 and m.sample_id.is_unique and m.profile_id.nunique()==17
 for _,g in m.groupby('profile_id'):
  assert sorted(g.depth_cm)==[5,20,40] and g.restoration4.nunique()==1 and g.locality.nunique()==1
 print(json.dumps({'input_files_verified':len(manifest),'samples':51,'profiles':17}),flush=True)
 if a.check_only:return
 stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S')
 run=a.out_root.resolve()/('10_01_run_expanded_'+stamp);run.mkdir(parents=True,exist_ok=False)
 for item in manifest:
  rel=Path(item['path']);dest=run/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.symlink_to(source/rel)
 (run/'expanded_analysis').mkdir(exist_ok=True);(run/'nitrogen_analysis').mkdir(exist_ok=True)
 for f in (P/'expanded/code').glob('*.py'):
  dst=run/('nitrogen_analysis' if f.name=='13_compare_nitrogen.py' else 'expanded_analysis')/f.name
  shutil.copy2(f,dst)
 shutil.copy2(P/'expanded/manuscript_edits.json',run/'expanded_analysis/manuscript_edits.json')
 (run/'launch.json').write_text(json.dumps({'bundle_root':str(source),'code_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (P/'expanded/code').glob('*.py')}},indent=2))
 env=dict(os.environ,OPENBLAS_NUM_THREADS='1')
 print('Run:',run,flush=True)
 if a.prepare_only:return
 for name in ['15_compare_expanded.py','17_present_expanded.py','16_build_atlas.py','18_report.py']:
  subprocess.run([sys.executable,str(run/'expanded_analysis'/name)],cwd=run,env=env,check=True)
 print('Results:',(run/'expanded_analysis/LATEST.txt').read_text(),flush=True)
if __name__=='__main__':main()
