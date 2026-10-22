#!/usr/bin/env python3
"""Prepare corrected read-derived KO network inputs without using any index in selection."""
from pathlib import Path
import argparse,hashlib,json,math
import numpy as np
import pandas as pd

EXPECTED={
 'KO_complete_RPK.tsv':'40ffa464ac84e97955c39f268e461dfe81b237feb29d9b3e302c5762564ed408',
 'KO_complete_counts.tsv':'5eb8b36e3f4c31454289061277380eda465e4d09f382006499e6848a53b19b0b'}
META_SHA='d9fb0bef5be343d76cdf0bd4ef5b4c7059386e9102e3a670dd67b34a1b3aaf1b'

def clr(x):
    log=np.log1p(x);return log-log.mean(axis=1,keepdims=True)

def ranking(rpk,counts):
    prevalence=(counts>0).sum(0);total=counts.sum(0)
    keep=(prevalence>=math.ceil(len(counts)/2))&(total>=10)&(np.log1p(rpk).var(0,ddof=1)>0)
    z=clr(rpk.loc[:,keep].to_numpy())
    return pd.DataFrame({'KO':rpk.columns[keep],'prevalence':prevalence[keep].to_numpy(),
        'total_fragments':total[keep].to_numpy(),'total_RPK':rpk.loc[:,keep].sum(0).to_numpy(),
        'CLR_variance':z.var(0,ddof=1)}).sort_values(['CLR_variance','prevalence','total_RPK','KO'],ascending=[False,False,False,True],kind='stable').reset_index(drop=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrices',type=Path,required=True);p.add_argument('--metadata',type=Path,required=True)
    p.add_argument('--bundle',type=Path,required=True);args=p.parse_args()
    if args.bundle.exists():raise SystemExit('Choose a fresh runtime bundle directory')
    hashes={}
    for name,expected in EXPECTED.items():
        path=args.matrices/name;digest=hashlib.sha256(path.read_bytes()).hexdigest()
        if digest!=expected:raise ValueError(f'Unexpected input hash: {path}')
        hashes[name]=digest
    if hashlib.sha256(args.metadata.read_bytes()).hexdigest()!=META_SHA:raise ValueError('Metadata calibration or cohort changed')
    m=pd.read_csv(args.metadata,dtype={'lat_block':str});m['profile']=m.sample_id.str.replace(r'_[123]$','',regex=True)
    assert len(m)==51 and m.sample_id.is_unique and m.profile.nunique()==17
    assert all(len(g)==3 and sorted(g.depth_cm)==[5,20,40] and g.locality.nunique()==1 and g.lat_block.nunique()==1 for _,g in m.groupby('profile'))
    r=pd.read_csv(args.matrices/'KO_complete_RPK.tsv',sep='\t',index_col=0).T
    c=pd.read_csv(args.matrices/'KO_complete_counts.tsv',sep='\t',index_col=0).T
    assert r.index.is_unique and r.columns.is_unique and set(r.index)==set(m.sample_id)
    assert r.index.equals(c.index) and r.columns.equals(c.columns)
    r=r.loc[m.sample_id];c=c.loc[m.sample_id]
    assert np.isfinite(r).all().all() and (r>=0).all().all() and np.isfinite(c).all().all() and (c>=0).all().all()
    assert np.array_equal(r.to_numpy()>0,c.to_numpy()>0)
    rank=ranking(r,c);assert len(rank)==7802
    out=args.bundle;inp=out/'inputs';inp.mkdir(parents=True);selection=out/'selection';selection.mkdir()
    rank.to_csv(selection/'KO_ranking.tsv',sep='\t',index=False)
    m.to_csv(inp/'metadata.csv',index=False)
    for n in [100,200,300]:
        ids=rank.KO.head(n);z=clr(r.loc[:,ids].to_numpy())
        assert np.max(np.abs(z.mean(1)))<1e-12
        pd.DataFrame(z,index=m.sample_id,columns=ids).to_csv(inp/f'clr_{n}_KO.tsv.gz',sep='\t',index_label='sample_id',float_format='%.17g')
    ppm=r.div(r.sum(1),axis=0)*1e6
    z=clr(ppm.loc[:,rank.KO.head(200)].to_numpy())
    pd.DataFrame(z,index=m.sample_id,columns=rank.KO.head(200)).to_csv(inp/'clr_200_KO_per_million.tsv.gz',sep='\t',index_label='sample_id',float_format='%.17g')
    detail=[];chosen=set(rank.KO.head(200))
    for profile in m.profile.unique():
        ii=m.loc[m.profile.ne(profile),'sample_id'];candidate=ranking(r.loc[ii],c.loc[ii]);new=set(candidate.KO.head(200));overlap=len(new&chosen)
        detail.append({'omitted_profile':profile,'eligible_KO':len(candidate),'overlap_with_primary':overlap,'retention':overlap/200,'jaccard':overlap/len(new|chosen)})
    pd.DataFrame(detail).to_csv(selection/'profile_deletion_selection.tsv',sep='\t',index=False)
    report={'status':'PASS','n_samples':51,'n_profiles':17,'complete_KO':len(r.columns),'eligible_KO':len(rank),'selected_KO':200,
       'primary_RPK':True,'pseudocount':1,'ranking_CLR_center':'all eligible KOs','network_CLR_center':'selected KOs',
       'input_sha256':dict(hashes,metadata=META_SHA),'selection_uses_predictors':False,
       'minimum_selection_retention':min(d['retention'] for d in detail),'median_selection_retention':float(np.median([d['retention'] for d in detail]))}
    (selection/'preparation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
