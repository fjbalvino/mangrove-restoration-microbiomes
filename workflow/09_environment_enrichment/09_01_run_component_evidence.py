"""Reproduce ORA, audit archived gene families, then evaluate component evidence."""
from pathlib import Path
import argparse,subprocess,datetime,json,hashlib
ap=argparse.ArgumentParser()
for n in ['analysis','reference','cf-data','definitions','reference-genes','out-root']:ap.add_argument('--'+n,type=Path,required=True)
ap.add_argument('--reference-results',type=Path)
a=ap.parse_args();here=Path(__file__).resolve().parent;repo=here.parents[1]
out=a.out_root/('09_01_run_component_evidence_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S'));out.mkdir(parents=True,exist_ok=False)
commands=[['bash',str(here/'09_00_run_enrichment.sh'),str(a.analysis),str(a.reference),str(out/'09_enrichment')],['python3',str(here/'23_verify_gene_families.py'),'--analysis',str(a.analysis),'--output',str(out/'23_gene_validation')],['python3',str(repo/'checks/verify_read_based_reml.py'),'--analysis-root',str(a.analysis),'--output',str(out/'reml_validation')]]
if a.reference_results:commands.append(['python3',str(repo/'checks/verify_environment_enrichment.py'),'--results',str(out/'09_enrichment'),'--reference-results',str(a.reference_results)])
commands.append(['python3',str(here/'22_environment_blocks.py'),'--analysis',str(a.analysis),'--ora',str(out/'09_enrichment'),'--definitions',str(a.definitions),'--cf-data',str(a.cf_data),'--reference-genes',str(a.reference_genes),'--old-reference',str(a.reference),'--out-root',str(out)])
with (out/'execution.log').open('w') as log:
 for cmd in commands:
  print('Running',cmd[1],flush=True);subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT)
 block=Path((out/'LATEST.txt').read_text().strip());cmd=['python3',str(here/'24_plot_environment_blocks.py'),'--ora',str(out/'09_enrichment'),'--blocks',str(block),'--out-root',str(out)];commands.append(cmd);subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT)
(out/'execution.json').write_text(json.dumps({'status':'COMPLETE','commands':commands,'scope':'Enrichment rerun and component overlay; all archived BH families audited; 714 model fits rechecked, not a full catalogue/model refit','definitions_sha256':hashlib.sha256(a.definitions.read_bytes()).hexdigest()},indent=2))
print(out)
