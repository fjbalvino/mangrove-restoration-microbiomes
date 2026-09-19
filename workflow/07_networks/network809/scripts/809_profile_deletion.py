#!/usr/bin/env python3
# Inputs: network809/inputs/{clr_200_KO.tsv.gz,metadata.csv}.
# Diagnostics additionally consume the completed 809 run (objects, draws and tables).
# Outputs: configured results root/809_residualized_functional_networks_*/; see PROTOCOL.md.
# Provenance: fixed-lambda ridge + LIONESS (Kuijjer et al., 2019, doi:10.1016/j.isci.2019.03.021).
# Curation: OpenAI Codex (OpenAI, 2026); documentation only; algorithm body preserved.
# Module filenames remain unchanged to preserve Python imports and historical source-hash checks.
"""Post-hoc diagnostic explicitly described in AMENDMENT_01.md; no hypothesis tests."""
import argparse,time
from pathlib import Path
import pandas as pd,numpy as np
import network_statistics as stats
from importlib import import_module
core=import_module('809_residual_networks')
ap=argparse.ArgumentParser();ap.add_argument('--bundle',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--run',type=Path,required=True);args=ap.parse_args();b=args.bundle.resolve();out=args.run;m=pd.read_csv(out/'tables/809_metadata_aligned.tsv',sep='\t');x0=pd.read_csv(b/'inputs/clr_200_KO.tsv.gz',sep='\t',index_col=0);rows=[]
for sid in core.SC:
 z=np.load(out/'objects'/('809_'+sid+'_networks.npz'));x=x0.loc[m.sample_id,z['nodes']].to_numpy();c=z['nuisance_design'];a=z['node_a_index'];bb=z['node_b_index'];Y=z['sample_metrics'];nodes=z['nodes']
 for profile in m.profile.unique():
  ii=np.flatnonzero(m.profile!=profile);mm=m.iloc[ii].reset_index(drop=True);D=stats.design(mm);W,newY,_,_,_,diag=core.lioness(x[ii],c[ii],a,bb,nodes);rebuilt=stats.adjusted(mm[core.IND].to_numpy(),newY,D)[0];frozen=stats.adjusted(mm[core.IND].to_numpy(),Y[ii],D)[0]
  for j,p in enumerate(core.IND):
   for k,e in enumerate(core.ENDS):rows.append(dict(construction=sid,omitted_profile=profile,predictor=p,endpoint=e,beta_network_reconstructed=rebuilt[j,k],beta_network_fixed=frozen[j,k],n_samples=48,n_profiles=16,diagnostic='post hoc; no P values',**diag))
 print(sid,'17 deletion reconstructions complete',flush=True)
d=pd.DataFrame(rows);d.to_csv(out/'tables/809_profile_deletion_diagnostic.tsv',sep='\t',index=False,float_format='%.17g');d.groupby(['construction','predictor','endpoint']).agg(rebuilt_beta_min=('beta_network_reconstructed','min'),rebuilt_beta_median=('beta_network_reconstructed','median'),rebuilt_beta_max=('beta_network_reconstructed','max'),frozen_beta_min=('beta_network_fixed','min'),frozen_beta_max=('beta_network_fixed','max')).reset_index().to_csv(out/'tables/809_profile_deletion_summary.tsv',sep='\t',index=False,float_format='%.17g')
