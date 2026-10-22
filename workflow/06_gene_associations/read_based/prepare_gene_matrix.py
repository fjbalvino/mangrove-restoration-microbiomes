"""Reconstruct read-derived gene inputs; filter independently of predictors.
Primary response: log(RPK+1), centred across the full retained catalogue per sample.
Filtering: nonzero unique fragment counts in >=5 samples, total >=10 fragments.
Selection: 100000 largest sample CLR variances, with numeric gene ID tie-break.
"""
from pathlib import Path
import json,time,gzip
import numpy as np
import pandas as pd

import argparse,hashlib
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input',type=Path,required=True,help='Reassembled read-derived gene export directory')
parser.add_argument('--analysis-root',type=Path,required=True)
args=parser.parse_args()
ROOT=args.analysis_root.resolve();ROOT.mkdir(parents=True,exist_ok=True)
INP=args.input.resolve();OUT=ROOT/'prepared';CACHE=ROOT/'sample_cache'
metadata_hash=hashlib.sha256((INP/'analysis_metadata/metadata_003b.csv').read_bytes()).hexdigest()
if metadata_hash!='d9fb0bef5be343d76cdf0bd4ef5b4c7059386e9102e3a670dd67b34a1b3aaf1b':
    raise SystemExit('Metadata differs from the verified 51-sample input; review lineage before analysis')
OUT.mkdir(exist_ok=True);CACHE.mkdir(exist_ok=True)
m=pd.read_csv(INP/'analysis_metadata/metadata_003b.csv')
m['profile_id']=m.sample_id.str.replace(r'_[123]$','',regex=True)
assert len(m)==51 and m.sample_id.is_unique and m.profile_id.nunique()==17
for _,g in m.groupby('profile_id'):
    assert sorted(g.depth_cm)==[5,20,40] and g.locality.nunique()==g.restoration4.nunique()==1
m.to_csv(OUT/'metadata.csv',index=False)
cap=30000000
prev=np.zeros(cap,dtype=np.uint8);total=np.zeros(cap,dtype=np.uint64)
qc=[]
for j,sid in enumerate(m.sample_id):
    p=CACHE/(sid+'.npz')
    if not p.exists():
        a=pd.read_csv(INP/'samples'/sid/'historical_gene_counts.tsv.gz',sep='\t',dtype={'gene_numeric_id':np.uint32,'unique_total':np.uint32,'RPK':np.float64,'both_RPK':np.float64},usecols=['gene_numeric_id','unique_total','RPK','both_RPK'])
        assert a.gene_numeric_id.is_unique and np.isfinite(a.to_numpy()).all() and (a.unique_total>0).all() and (a.RPK>0).all()
        np.savez(p,ids=a.gene_numeric_id.to_numpy(),counts=a.unique_total.to_numpy(),rpk=a.RPK.to_numpy(),both_rpk=a.both_RPK.to_numpy())
        del a
    with np.load(p) as a:
        ids=a['ids'];counts=a['counts'];rpk=a['rpk']
        if ids.max()>=len(prev):
            cap=int(ids.max())+1000000;prev=np.pad(prev,(0,cap-len(prev)));total=np.pad(total,(0,cap-len(total)))
        prev[ids]+=1;total[ids]+=counts
        stats=json.loads((INP/'samples'/sid/'COMPLETE.json').read_text())['stats']
        assert int(counts.sum(dtype=np.uint64))==stats['historical_gene_assigned_fragments'] and len(ids)==stats['historical_gene_nonzero']
        qc.append({'sample_id':sid,'nonzero_genes':len(ids),'assigned_fragments':int(counts.sum(dtype=np.uint64)),'all_unique_fragments':stats['unique_fragments'],'gene_assigned_fraction':int(counts.sum(dtype=np.uint64))/stats['unique_fragments'],'RPK_sum':float(rpk.sum())})
    print('INPUT',j+1,51,sid,len(ids),flush=True)
keep=(prev>=5)&(total>=10);keep_ids=np.flatnonzero(keep).astype(np.uint32);nkeep=len(keep_ids)
np.save(OUT/'filtered_gene_numeric_ids.npy',keep_ids)
np.save(OUT/'prevalence.npy',prev[keep]);np.save(OUT/'total_unique_fragments.npy',total[keep])
print('FILTERED',nkeep,'DETECTED_UNION',int((prev>0).sum()),flush=True)
sumlog=np.zeros(len(prev));sumlog2=np.zeros(len(prev));cross=np.zeros(len(prev));centers=[]
for j,sid in enumerate(m.sample_id):
    with np.load(CACHE/(sid+'.npz')) as a:
        ids=a['ids'];mask=keep[ids];ids=ids[mask];v=np.log1p(a['rpk'][mask]);c=float(v.sum()/nkeep)
        centers.append(c);sumlog[ids]+=v;sumlog2[ids]+=v*v;cross[ids]+=v*c
        qc[j].update(filtered_detected_genes=len(ids),filtered_fragments=int(a['counts'][mask].sum(dtype=np.uint64)),filtered_RPK=float(a['rpk'][mask].sum()),CLR_center=c)
    print('CENTER',j+1,51,sid,flush=True)
centers=np.array(centers)
mean=(sumlog[keep]-centers.sum())/51
variance=(sumlog2[keep]-2*cross[keep]+(centers**2).sum()-51*mean**2)/50
assert np.isfinite(variance).all() and variance.min()>-1e-12
candidate=np.argpartition(variance,-100000)[-100000:]
candidate=candidate[np.lexsort((keep_ids[candidate],-variance[candidate]))]
sel_ids=keep_ids[candidate];np.save(OUT/'selected_gene_numeric_ids.npy',sel_ids)
index=np.full(len(prev),-1,dtype=np.int32);index[sel_ids]=np.arange(100000)
raw=np.zeros((100000,51));counts_out=np.zeros((100000,51),dtype=np.uint32);both=np.zeros((100000,51))
for j,sid in enumerate(m.sample_id):
    with np.load(CACHE/(sid+'.npz')) as a:
        pos=index[a['ids']];mask=pos>=0
        raw[pos[mask],j]=a['rpk'][mask];counts_out[pos[mask],j]=a['counts'][mask];both[pos[mask],j]=a['both_rpk'][mask]
clr=np.log1p(raw)-centers[None,:]
assert np.max(np.abs(clr.var(axis=1,ddof=1)-variance[candidate]))<1e-10
np.save(OUT/'selected_RPK.npy',raw);np.save(OUT/'selected_unique_counts.npy',counts_out);np.save(OUT/'selected_both_RPK.npy',both);np.save(OUT/'selected_CLR.npy',clr)
selection=pd.DataFrame({'rank':np.arange(1,100001),'gene_numeric_id':sel_ids,'variance_clr':variance[candidate],'sd_clr':np.sqrt(variance[candidate]),'prevalence':prev[sel_ids],'total_fragments':total[sel_ids]})
selection.to_csv(OUT/'selection_numeric.tsv',sep='\t',index=False)
pd.DataFrame(qc).to_csv(OUT/'sample_qc.tsv',sep='\t',index=False)
summary={'status':'PASS','samples':51,'profiles':17,'detected_gene_union':int((prev>0).sum()),'filtered_genes':nkeep,'selected_genes':100000,'min_prevalence':5,'min_total_unique_fragments':10,'response':'RPK','CLR_pseudocount':1,'CLR_center':'full_filtered_catalogue','selection':'global_CLR_variance','predictors_used_for_selection':False,'stage_used_for_selection':False,'metadata_sha256':'d9fb0bef5be343d76cdf0bd4ef5b4c7059386e9102e3a670dd67b34a1b3aaf1b'}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary),flush=True)
