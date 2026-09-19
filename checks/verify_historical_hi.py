#!/usr/bin/env python3
"""Reconstruct archived HI without redefining its calibration.
Inputs: source_data/indices/{historical_CLR,historical_metadata,historical_HI,metadata_003b}.csv
Output: validation/HI_reconstruction_audit.json
Method: categorical nuisance regression by QR; centered, unscaled PCA by SVD.
AI-assisted curation/implementation: OpenAI Codex (OpenAI, 2026).
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
src = root / 'source_data/indices'
meta = pd.read_csv(src/'historical_metadata.csv').set_index('sample_id')
x = pd.read_csv(src/'historical_CLR.csv', index_col=0).T.loc[meta.index]
archived = pd.read_csv(src/'historical_HI.csv').set_index('sample_id')
current = pd.read_csv(src/'metadata_003b.csv').set_index('sample_id')
d = np.column_stack([np.ones(len(meta)), pd.get_dummies(meta.locality, drop_first=True).to_numpy(float), pd.get_dummies(meta.depth_cm, drop_first=True).to_numpy(float)])
q = np.linalg.qr(d, mode='reduced')[0]
r = x.to_numpy() - q @ (q.T @ x.to_numpy())
u, s, vt = np.linalg.svd(r-r.mean(0), full_matrices=False)
score = u[:, 0] * s[0]
flip = -1 if np.median(score[meta.collapsed_stage=='Preserved']) < np.median(score[meta.collapsed_stage=='Degraded']) else 1
score *= flip
observed = archived.loc[meta.index, 'RBMA1'].to_numpy()
errs = {'QR_residual_SVD_to_archived_HI':float(np.max(abs(score-observed))), 'archived_HI_to_current_metadata':float(np.max(abs(archived.loc[current.index,'RBMA1']-current.HI)))}
res = {'status':'PASS' if max(errs.values())<1e-9 else 'FAIL','samples':len(meta),'taxa':x.shape[1],'design_rank':int(np.linalg.matrix_rank(d)), 'depth_coding':'categorical','PCA_scale':False,'index_standardized_after_PCA':False,'PC1_variance_fraction':float(s[0]**2/(s*s).sum()),'max_absolute_errors':errs,'historical_reference_stage':'Preserved','historical_taxon_set_is_current_1031_screen':False}
out=root/'validation/HI_reconstruction_audit.json'
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(res,indent=2)+'\n')
print(out.read_text())
assert res['status']=='PASS'
