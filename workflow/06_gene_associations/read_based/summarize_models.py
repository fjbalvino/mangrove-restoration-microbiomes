from pathlib import Path
import numpy as np,pandas as pd,json
from fit_random_intercept import AXES
import argparse
parser=argparse.ArgumentParser(description='Adjust declared gene-model families and export model tables.')
parser.add_argument('--analysis-root',type=Path,required=True)
args=parser.parse_args();R=args.analysis_root.resolve();O=R/'models';sel=pd.read_csv(R/'prepared/selected_genes.tsv',sep='\t');ids=sel.historical_cluster_id.to_numpy();n=len(ids)
def bh(p):
    p=np.asarray(p);idx=np.argsort(p,kind='stable');q=np.minimum.accumulate((p[idx]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];out=np.empty_like(q);out[idx]=np.minimum(q,1);return out
res={}
for axis in AXES:
    with np.load(O/(axis+'_common.npz')) as a:
        res[axis,'common']={k:a[k].copy() for k in a.files}
    with np.load(O/(axis+'_interaction.npz')) as a:
        res[axis,'interaction']={k:a[k].copy() for k in a.files}
for group in [AXES[:1],AXES[1:]]:
    for family,key,model in [('common','p','common'),('interaction','interaction_p','interaction'),('slopes','p','interaction')]:
        pp=[];ppd=[];sizes=[]
        for axis in group:
            a=res[axis,model];p=a[key];u=a['usable'] if p.ndim==1 else a['usable'][:,None];diag=a['hessian_ok'] if p.ndim==1 else a['hessian_ok'][:,None]
            pp.append(np.where(u,p,1).ravel());ppd.append(np.where(diag,p,1).ravel());sizes.append(p.size)
        q=bh(np.concatenate(pp));qd=bh(np.concatenate(ppd));off=0
        for axis,size in zip(group,sizes):
            a=res[axis,model];p=a[key];u=a['usable'] if p.ndim==1 else a['usable'][:,None]
            a[family+'_q']=np.where(u,q[off:off+size].reshape(p.shape),np.nan)
            a[family+'_diagnostic_q']=qd[off:off+size].reshape(p.shape);off+=size
hits=[];summ=[];diagnostics=[]
for axis in AXES:
    rows=[]
    for family,model in [('common','common'),('interaction','interaction'),('slopes','interaction')]:
        a=res[axis,model];q=a[family+'_q'];qd=a[family+'_diagnostic_q']
        for j,dep in enumerate([5,20,40] if family=='slopes' else [0]):
            idx=j if family=='slopes' else 0; joint=family=='interaction'
            qq=q if q.ndim==1 else q[:,idx];qqd=qd if qd.ndim==1 else qd[:,idx]
            pp=a['interaction_p'] if joint else a['p'][:,idx]
            estimate=np.full(n,np.nan) if joint else a['estimate'][:,idx]
            tab=pd.DataFrame({'gene_numeric_id':sel.gene_numeric_id,'feature_id':ids,'axis':axis,'family':family,'depth_cm':dep,'estimate_sd':estimate,'se_sd':np.full(n,np.nan) if joint else a['se'][:,idx],'df_den':a['interaction_df'] if joint else a['df'][:,idx],'statistic':a['interaction_F'] if joint else a['estimate'][:,idx]/a['se'][:,idx],'ci_low_sd':np.full(n,np.nan) if joint else a['ci_low'][:,idx],'ci_high_sd':np.full(n,np.nan) if joint else a['ci_high'][:,idx],'p':pp,'q':qq,'q_including_singular_diagnostic':qqd,'usable':a['usable'],'singular':a['singular'],'icc':a['icc']})
            rows.append(tab);hits.append(tab[tab.q<=.1])
            summ.append({'axis':axis,'family':family,'depth_cm':dep,'n_tests':n,'n_usable':int(a['usable'].sum()),'n_singular':int(a['singular'].sum()),'n_q_0_10':int((qq<=.1).sum()),'n_q_0_05':int((qq<=.05).sum()),'negative_q_0_10':int(((qq<=.1)&(estimate<0)).sum()),'positive_q_0_10':int(((qq<=.1)&(estimate>0)).sum()),'min_q':float(np.nanmin(qq)),'n_diagnostic_q_0_10':int((qqd<=.1).sum())})
    pd.concat(rows,ignore_index=True).to_csv(O/(axis+'_all_results.tsv.gz'),sep='\t',index=False,compression={'method':'gzip','compresslevel':1})
allhits=pd.concat(hits,ignore_index=True).merge(sel[['gene_numeric_id','KO_set','PFAM_set']],on='gene_numeric_id',validate='many_to_one')
allhits.to_csv(O/'discoveries_q010.tsv.gz',sep='\t',index=False,compression={'method':'gzip','compresslevel':1})
pd.DataFrame(summ).to_csv(O/'discovery_counts.tsv',sep='\t',index=False)
print(pd.DataFrame(summ).to_string(index=False),flush=True)
summary={'status':'FITTED_R_VALIDATION_PENDING','total_discoveries_q010':len(allhits),'distinct_genes_q010':allhits.feature_id.nunique(),'environment_distinct_genes_q010':allhits[allhits.axis!='MHI_local'].feature_id.nunique(),'MHI_distinct_genes_q010':allhits[allhits.axis=='MHI_local'].feature_id.nunique()}
(O/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary,flush=True)
