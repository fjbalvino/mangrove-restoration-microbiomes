from pathlib import Path
import argparse,datetime,json,hashlib
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
ap=argparse.ArgumentParser();ap.add_argument('--expanded-dir',type=Path,required=True);ap.add_argument('--depth-dir',type=Path,required=True);ap.add_argument('--out-root',type=Path,default=Path('figure6_depth'));a=ap.parse_args();S=a.expanded_dir;D=a.depth_dir;O=a.out_root/('21_figure6_depth_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S'));O.mkdir(parents=True,exist_ok=False);a.out_root.mkdir(exist_ok=True);(a.out_root/'LATEST.txt').write_text(str(O.resolve()))
t=pd.read_csv(D/'depth_tests.tsv',sep='\t');p=t[t.variant.eq('primary')];r=pd.read_csv(S/'eligibility.tsv',sep='\t').set_index('KO');e=pd.read_csv(S/'effects.tsv',sep='\t');v=pd.read_csv(D/'primary_normalized.tsv',sep='\t',index_col=0);m=pd.read_csv(D/'samples.tsv',sep='\t')
info=[('K00261','gdhA','Glutamate dehydrogenase, NAD(P)'),('K00262','gdhA','Glutamate dehydrogenase, NADP'),('K00284','gltS','Glutamate synthesis'),('K00371','narH/nxrB','Nitrate/nitrite transformations'),('K00368','nirK','Nitrite reduction'),('K15864','nirS','Nitrite/hydroxylamine'),('K05601','hcp','Hydroxylamine reduction'),('K01429','ureB','Urea hydrolysis, beta subunit'),('K01430','ureA','Urea hydrolysis, gamma subunit'),('K01501','nitrilase','Nitrile hydrolysis'),('K01725','cynS','Cyanate breakdown'),('K15578','nrtC','Nitrate/nitrite transport')]
keys=[z[0] for z in info];assert set(keys)==set(p.loc[(p.depth==40)&(p.p_holm_98<.05),'KO']) and len(keys)==12
selected=json.load(open(S/'display_selection.json'));st=['Degraded','Early restoration','Intermediate restoration','Advanced restoration'];short=['Degraded','Early','Intermediate','Advanced'];q=e[(e.variant=='primary')&(e.target=='restoration')].pivot(index='KO',columns='stage',values='log2_ratio');M=q.reindex(index=[x[1] for x in selected],columns=st).to_numpy();whole=q.reindex(index=r.index[r.eligible],columns=st)
rob=t.pivot(index=['KO','depth'],columns='variant',values='p_holm_1409');stable=set(k for (k,d) in rob[(rob<.05).all(axis=1)].index if d==40)
cmap=LinearSegmentedColormap.from_list('div',['#2166ac','white','#b2182b']);cmap.set_bad('#ddd');plt.rcParams.update({'font.size':8,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(7.5,8.2));ax=fig.add_axes([.20,.62,.41,.33]);bx=fig.add_axes([.72,.62,.13,.33]);cx=fig.add_axes([.46,.10,.46,.40]);ca=fig.add_axes([.89,.66,.02,.25]);im=ax.imshow(M,cmap=cmap,vmin=-2,vmax=2,aspect='auto');ax.set_yticks(range(12),[cy+' | '+name.split(' · ')[0]+(' †' if r.loc[k,'available'] and r.loc[k,'quantification']!='single_KO' else '') for cy,k,name in selected],fontsize=7.5);ax.set_xticks(range(4),['D','E','I','A']);ax.set_title('a  Stage / conserved',loc='left',fontsize=10)
for i in range(12):
 for j in range(4):ax.text(j,i,'—' if np.isnan(M[i,j]) else f'{2**M[i,j]:.2f}',ha='center',va='center',fontsize=7.3)
bx.imshow(whole,cmap=cmap,vmin=-2,vmax=2,aspect='auto',interpolation='nearest');bx.set_yticks([]);bx.set_xticks(range(4),['D','E','I','A']);bx.set_title('b  437 KOs',loc='left',fontsize=10);cb=fig.colorbar(im,cax=ca);cb.set_ticks([-2,-1,0,1,2],labels=['¼','½','1','2','4']);fig.text(.20,.575,'D Degraded · E Early · I Intermediate · A Advanced',fontsize=8)
labels=[f'{gene} | {func}\n{k}'+(' †' if r.loc[k,'quantification']!='single_KO' else '') for k,gene,func in info];yy=np.arange(12)
for depth,offset,col in [(20,-.16,'#1685a4'),(40,.16,'#b44b36')]:
 z=p[p.depth.eq(depth)].set_index('KO').loc[keys]
 for i,(k,row) in enumerate(z.iterrows()):
  cx.plot([row.ratio_low,row.ratio_high],[i+offset]*2,color=col,lw=1)
  marker='D' if row.p_holm_1409<.05 else 'o';cx.scatter(row.ratio,i+offset,s=24,marker=marker,facecolors=col if row.p_holm_98<.05 else 'white',edgecolors=col,zorder=3)
cx.set_xscale('log',base=2);cx.set_xlim(.025,2);cx.set_xticks([.03125,.0625,.125,.25,.5,1,2],['1/32','1/16','1/8','1/4','1/2','1','2']);cx.axvline(1,color='#888',lw=.8);cx.set_ylim(11.6,-.6);cx.set_yticks(yy,labels,fontsize=7)
for tick,k in zip(cx.get_yticklabels(),keys):
 if k in stable:tick.set_fontweight('bold')
cx.set_title('c  Paired nitrogen depth changes',loc='left',fontsize=10,pad=13);cx.set_xlabel('Depth / 5 cm (median paired ratio)',fontsize=9);cx.legend(handles=[Line2D([],[],color='#1685a4',label='20 / 5 cm'),Line2D([],[],color='#b44b36',label='40 / 5 cm')],loc='lower left',bbox_to_anchor=(-.1,1.10),ncol=2,frameon=False,fontsize=8)
fig.text(.05,.025,'Open: unsupported · Filled: Holm98 < 0.05 · Diamond: Holm1409 < 0.05\nBold c labels: Holm1409 support in all five variants · † Provisional quantification',fontsize=8)
fig.savefig(O/'Figure6.png',dpi=350);fig.savefig(O/'Figure6.pdf');plt.close(fig)
# A larger independent panel for easy review.
# Paired trajectories preserve actual depth spacing and show all 17 profiles.
rows=[]
for k in keys:
 for pid,g in m.groupby('profile_id'):
  g=g.sort_values('depth_cm');assert g.depth_cm.tolist()==[5,20,40]
  values=v.loc[k,g.sample_id].to_numpy();rel=values/values[0]
  for ix,z in enumerate(g.itertuples()):rows.append(dict(KO=k,profile_id=pid,sample_id=z.sample_id,depth_cm=z.depth_cm,locality=z.locality,stage=z.restoration4,normalized_abundance=values[ix],ratio_to_5=rel[ix]))
tr=pd.DataFrame(rows);assert len(tr)==612 and np.allclose(tr.loc[tr.depth_cm.eq(5),'ratio_to_5'],1);tr.to_csv(O/'trajectory_data.tsv',sep='\t',index=False)
p[p.KO.isin(keys)].to_csv(O/'panel_c_data.tsv',sep='\t',index=False);pd.DataFrame(info,columns=['KO','gene','function']).to_csv(O/'display_selection.tsv',sep='\t',index=False)
lo=2**np.floor(np.log2(tr.ratio_to_5.min()));hi=2**np.ceil(np.log2(tr.ratio_to_5.max()));colors={'Carmen':'#d95f02','Cozumel':'#1b9e77','Tuxpan':'#7570b3'}
fig,axes=plt.subplots(4,3,figsize=(11,12),sharex=True,sharey=True)
for ax,(k,gene,func) in zip(axes.flat,info):
 a=tr[tr.KO.eq(k)]
 for pid,g in a.groupby('profile_id'):ax.plot(g.depth_cm,g.ratio_to_5,'o-',color=colors[g.locality.iloc[0]],lw=.7,alpha=.60,ms=2.5)
 med=a.groupby('depth_cm').ratio_to_5.median();ax.scatter(med.index,med.values,color='black',marker='D',s=20,zorder=5);ax.axhline(1,color='#777',lw=.6,ls='--');ax.set_yscale('log',base=2);ax.set_ylim(lo,hi);ax.set_yticks([1/128,1/32,1/8,1/2,1,2,8],['1/128','1/32','1/8','1/2','1','2','8']);ax.set_xticks([5,20,40]);ax.set_title(f'{gene} · {k}'+(' †' if r.loc[k,'quantification']!='single_KO' else '')+'\n'+func,fontsize=9);ax.grid(axis='y',alpha=.12)
for ax in axes[-1]:ax.set_xlabel('Depth (cm)')
for ax in axes[:,0]:ax.set_ylabel('Representation relative to 5 cm')
fig.legend(handles=[Line2D([],[],color=c,label=l) for l,c in colors.items()]+[Line2D([],[],color='black',marker='D',linestyle='none',label='Median across profiles')],loc='upper center',ncol=4,bbox_to_anchor=(.5,.945),frameon=False)
fig.suptitle('Nitrogen depth trajectories in 17 profiles\nLines connect paired observations; no linear trend was fitted',fontsize=14,y=.99);fig.tight_layout(rect=[0,.025,1,.91]);fig.text(.03,.012,'All panels share the same scale. Every profile starts at 1 at 5 cm. † Provisional fractional quantification.',fontsize=9);fig.savefig(O/'Supplementary_nitrogen_trajectories.png',dpi=250);fig.savefig(O/'Supplementary_nitrogen_trajectories.pdf');plt.close(fig)
(O/'validation.json').write_text(json.dumps(dict(selected_by='Primary Holm98 < 0.05 for 40/5; post-result selection explicitly disclosed',selected_KOs=12,eligible_panel=49,paired_trajectories=204,trajectory_points=612,baseline_all_one=True,linear_models_fitted=False,stable_combined_keys=sorted(stable),stats_recomputed=False,source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [S/'effects.tsv',D/'depth_tests.tsv',D/'primary_normalized.tsv',D/'samples.tsv']}),indent=2));print(O)
