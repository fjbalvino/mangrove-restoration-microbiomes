"""Descriptive concordance of supporting subunits; no extra hypothesis tests."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from ejecutar_cns import bounds
from calibrar_y_verificar import designs
from stats_cns import profile_matrix
ROOT=Path(__file__).resolve().parents[1]

def main():
    r=pd.read_csv(ROOT/'inputs/KO_complete_RPK.tsv',sep='\t',index_col=0)
    f=pd.read_csv(ROOT/'inputs/KO_complete_full_terms_RPK.tsv',sep='\t',index_col=0)
    us=(ROOT/'inputs/uscg_76_kegg_min_2313.lst').read_text().split();den=r.loc[us].median()
    patterns=[p for p in json.loads((ROOT/'inputs/annotation_patterns_selected.json').read_text()) if p['status']=='complete_consensus']
    mapping=[('K02588','K02586','nifD'),('K02588','K02591','nifK'),('K11180','K11181','dsrB'),('K00394','K00395','aprB'),('K01601','K00855','prk'),('K04561','K02305','norC')]
    keys=[k for _,k,_ in mapping];lo,hi=bounds(r,f,patterns,keys)
    m,p,_=designs();v=pd.read_csv(ROOT/'results/primary_76_normalized.tsv',sep='\t',index_col=0)
    main=np.log2(profile_matrix(v,m,p.profile_id));results=[];normalized={}
    for target,k,gene in mapping:
        exact=np.allclose(lo.loc[k],hi.loc[k],rtol=1e-8,atol=1e-7)
        x=(lo.loc[k] if exact else r.loc[k])/den;assert (x>0).all();normalized[k]=x
        a=np.log2(profile_matrix(pd.DataFrame([x],index=[k]),m,p.profile_id)[0]);j=v.index.get_loc(target)
        delta=a[:,2]-a[:,0]
        results.append(dict(primary_KO=target,support_KO=k,gene=gene,quantification='single_KO' if exact else 'fractional_candidate',median_log2_40_vs_5=float(np.median(delta)),median_ratio_40_vs_5=float(np.exp2(np.median(delta))),negative_profiles=int((delta<0).sum()),n_profiles=17,profile_mean_spearman=float(spearmanr(main[j].mean(axis=1),a.mean(axis=1)).statistic),hypothesis_test_performed=False))
    pd.DataFrame(results).to_csv(ROOT/'results/supporting_subunits.tsv',sep='\t',index=False)
    pd.DataFrame(normalized).T.to_csv(ROOT/'results/supporting_subunits_normalized.tsv',sep='\t')
    print(pd.DataFrame(results).to_string(index=False))

if __name__=='__main__':main()
