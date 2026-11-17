"""Recompute BH families from archived P values, not a full model refit."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
ap=argparse.ArgumentParser();ap.add_argument('--analysis',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
axes=['MHI_local','vegetation_landscape_PC1','water_inundation_PC1','moisture_stress_PC1','physicochemical_PC1','nutrients_redox_PC1'];data={};counts=[];genes=set();maxerr=0.
def bh(p):
 order=np.argsort(p,kind='stable');q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
for axis in axes:
 d=pd.read_csv(a.analysis/'models'/f'{axis}_all_results.tsv.gz',sep='\t',usecols=['gene_numeric_id','family','depth_cm','p','q','usable','singular','estimate_sd'])
 assert len(d)==500000
 for fam in ['common','interaction','slopes']:
  x=d[d.family==fam];data[axis,fam]=(np.where(x.usable,x.p,1),x.q.to_numpy(),x.usable.to_numpy())
  for depth,y in x.groupby('depth_cm'):
   h=y[y.q.le(.1)];counts.append(dict(axis=axis,family=fam,depth_cm=int(depth),n_q010=len(h),n_q005=int(y.q.le(.05).sum())))
   if axis!='MHI_local':genes.update(h.gene_numeric_id)
 print(axis,flush=True)
for group in [axes[:1],axes[1:]]:
 for fam in ['common','interaction','slopes']:
  p=np.concatenate([data[x,fam][0] for x in group]);q=bh(p);offset=0
  for axis in group:
   pp,old,u=data[axis,fam];new=q[offset:offset+len(pp)];assert np.isnan(old[~u]).all();err=np.max(np.abs(new[u]-old[u]));maxerr=max(maxerr,float(err));assert err<1e-12;offset+=len(pp)
c=pd.DataFrame(counts);env=c[c.axis!='MHI_local'];tot=env.groupby('family').n_q010.sum().to_dict();assert tot=={'common':42222,'interaction':1234,'slopes':63942};assert len(genes)==48580;assert c[c.axis=='MHI_local'].n_q010.sum()==0
a.output.mkdir(parents=True,exist_ok=True);c.to_csv(a.output/'gene_counts_rechecked.tsv',sep='\t',index=False)
r={'status':'PASS','environment_counts':tot,'distinct_environment_genes':len(genes),'MHI_discoveries':0,'max_BH_error':maxerr,'scope':'All archived P values and BH families checked; no new full model fit'};(a.output/'gene_family_verification.json').write_text(json.dumps(r,indent=2));print(r)
