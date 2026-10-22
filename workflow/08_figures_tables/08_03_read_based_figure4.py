from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D
AXES=['MHI_local','vegetation_landscape_PC1','water_inundation_PC1','moisture_stress_PC1','physicochemical_PC1','nutrients_redox_PC1']
import argparse
parser=argparse.ArgumentParser(description='Render corrected read-based Figure 4 in PNG, PDF and SVG.')
group=parser.add_mutually_exclusive_group(required=True)
group.add_argument('--analysis-root',type=Path,help='Contains prepared/ and models/ from the corrected gene workflow')
group.add_argument('--figure-data',type=Path,help='Verified frozen Figure4 association-count and candidate-statistic tables')
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args();O=args.output.resolve();O.mkdir(parents=True,exist_ok=True)
D=O/'figure_data';D.mkdir(exist_ok=True)
candidates=[('CC1_1.contigs_267171','moisture_stress_PC1','Amino-acid N-acetyltransferase','K00619','NDMI'),('CC1_1.contigs_1063266','nutrients_redox_PC1','Phosphate acyltransferase','K03621','Nutrients / sulfide')]
if args.figure_data:
    summary=pd.read_csv(args.figure_data/'Figure4_association_counts.tsv',sep='\t')
    candidate_data=pd.read_csv(args.figure_data/'Figure4_candidate_statistics.tsv',sep='\t')
else:
    R=args.analysis_root.resolve()
    summary=pd.read_csv(R/'models/discovery_counts.tsv',sep='\t')
    sel=pd.read_csv(R/'prepared/selected_genes.tsv',sep='\t')
    m=pd.read_csv(R/'prepared/metadata.csv')
    count=np.load(R/'prepared/selected_unique_counts.npy')
    parts=[]
    for fid,axis,title,ko,label in candidates:
        for chunk in pd.read_csv(R/'models'/(axis+'_all_results.tsv.gz'),sep='\t',chunksize=100000):
            a=chunk[chunk.feature_id==fid].copy()
            if len(a):parts.append(a)
    candidate_data=pd.concat(parts,ignore_index=True)
    for fid,ax,*_ in candidates:
        i=sel.index[sel.historical_cluster_id==fid][0]
        for dep in [5,20,40]:candidate_data.loc[(candidate_data.feature_id==fid)&candidate_data.depth_cm.eq(dep),'n_detected']=int((count[i,m.depth_cm==dep]>0).sum())
candidate_data.to_csv(D/'Figure4_candidate_statistics.tsv',sep='\t',index=False);summary.to_csv(D/'Figure4_association_counts.tsv',sep='\t',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8.5,'axes.titlesize':9,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'xtick.major.width':.65,'ytick.major.width':.65})
fig=plt.figure(figsize=(185/25.4,172/25.4),facecolor='white')
labels=['MHI','Vegetation / landscape','Surface-water signal','Moisture (NDMI)','Physicochemistry','Nutrients / sulfide']
A=fig.add_axes([.255,.555,.315,.335]);B=fig.add_axes([.65,.555,.315,.335])
y=np.arange(6);common=summary[summary.family.eq('common')].set_index('axis').loc[AXES]
A.barh(y,-common.negative_q_0_10,color='#0072B2',height=.62);A.barh(y,common.positive_q_0_10,color='#D55E00',height=.62)
A.axvline(0,color='#333333',lw=.65);A.set_yticks(y,labels);A.set_ylim(5.6,-.6);A.set_xlim(-23000,10500);A.set_xticks([-20000,-10000,0,10000]);A.xaxis.set_major_formatter(FuncFormatter(lambda v,pos:f'{v/1000:g}'))
A.set_xlabel('Associations (×10³)',labelpad=5);A.tick_params(axis='y',length=0,pad=6);A.grid(axis='x',color='#e6e9ec',lw=.5);A.set_axisbelow(True)
A.legend(handles=[Line2D([0],[0],color='#0072B2',lw=6,label='Negative'),Line2D([0],[0],color='#D55E00',lw=6,label='Positive')],loc='upper center',bbox_to_anchor=(.5,-.21),ncol=2,frameon=False,fontsize=8,handlelength=1,columnspacing=1.1)
for i,(_,row) in enumerate(common.iterrows()):
    if row.n_q_0_10==0:A.text(700,i,'0',va='center',fontsize=8,color='#555555')
M=np.array([[int(summary[(summary.axis==ax)&summary.family.eq('slopes')&summary.depth_cm.eq(d)].n_q_0_10.iloc[0]) for d in [5,20,40]] for ax in AXES])
M05=np.array([[int(summary[(summary.axis==ax)&summary.family.eq('slopes')&summary.depth_cm.eq(d)].n_q_0_05.iloc[0]) for d in [5,20,40]] for ax in AXES])
cmap=LinearSegmentedColormap.from_list('counts',['#ffffff','#deebf7','#9ecae1','#3182bd','#08519c'])
im=B.imshow(M,aspect='auto',cmap=cmap,vmin=0,vmax=20000,extent=[-.5,2.5,5.5,-.5]);B.set_ylim(5.6,-.6);B.set_yticks([]);B.set_xticks([0,1,2],['5 cm','20 cm','40 cm']);B.tick_params(axis='x',length=0,pad=6)
for i in range(6):
 for j in range(3):B.text(j,i,f'{M[i,j]:,}\n({M05[i,j]:,})',ha='center',va='center',fontsize=8,color='white' if M[i,j]>16000 else '#152b31',linespacing=1.25)
for edge in np.arange(-.5,6,.999999):B.axhline(edge,color='#d1d8da',lw=.4)
for edge in [-.5,.5,1.5,2.5]:B.axvline(edge,color='#d1d8da',lw=.4)
B.spines[['left','bottom','top','right']].set_visible(False)
cbax=fig.add_axes([.69,.485,.235,.011]);cb=fig.colorbar(im,cax=cbax,orientation='horizontal',ticks=[0,10000,20000]);cb.ax.tick_params(labelsize=8,length=2,pad=2);cb.ax.xaxis.set_major_formatter(FuncFormatter(lambda v,pos:f'{v/1000:g}k' if v else '0'));cb.outline.set_visible(False)
fig.text(.02,.94,'a',weight='bold',fontsize=11);fig.text(.255,.94,'Common slopes',weight='bold',fontsize=10)
fig.text(.61,.94,'b',weight='bold',fontsize=11);fig.text(.65,.94,'Slopes by depth',weight='bold',fontsize=10)
fig.text(.65,.908,'q ≤ 0.10  (q ≤ 0.05)',fontsize=8)
colors=['#0072B2','#E69F00','#009E73'];axes=[]
for k,(fid,axis,title,ko,label) in enumerate(candidates):
    x0=.145 if k==0 else .64
    ax=fig.add_axes([x0,.14,.265,.20]);axes.append(ax)
    d=candidate_data[(candidate_data.feature_id==fid)&candidate_data.family.eq('slopes')].sort_values('depth_cm');qint=float(candidate_data[(candidate_data.feature_id==fid)&candidate_data.family.eq('interaction')].q.iloc[0])
    for j,(_,v) in enumerate(d.iterrows()):
        ax.errorbar(v.estimate_sd,j,xerr=[[v.estimate_sd-v.ci_low_sd],[v.ci_high_sd-v.estimate_sd]],fmt='o',color=colors[j],mfc=colors[j] if v.q<=.1 else 'white',mew=1,ms=5,elinewidth=1.1,capsize=2.5)
        ax.text(1.04,j,f'{int(v.n_detected)}/17',transform=ax.get_yaxis_transform(),va='center',ha='left',fontsize=8)
    ax.axvline(0,color='#7b858c',lw=.7,ls='--');ax.set_yticks([0,1,2],['5 cm','20 cm','40 cm']);ax.set_ylim(2.55,-.55);ax.set_xlim(-1.15,.8);ax.set_xticks([-1,-.5,0,.5]);ax.set_xlabel('Standardized slope (95% CI)',labelpad=6);ax.tick_params(axis='y',length=0,pad=5)
    ax.text(1.04,1.08,'Detected',transform=ax.transAxes,fontsize=8,ha='left')
    fig.text(.02 if k==0 else .52,.42,'c' if k==0 else 'd',weight='bold',fontsize=11)
    fig.text(.06 if k==0 else .56,.42,title,weight='bold',fontsize=9)
    fig.text(.06 if k==0 else .56,.39,f'{ko} · {label}',fontsize=8.5)
    fig.text(.06 if k==0 else .56,.36,f'Interaction q = {qint:.4f}',fontsize=8)
fig.text(.145,.055,'●  q ≤ 0.10     ○  q > 0.10',fontsize=8)
fig.text(.64,.055,'51 samples · 17 profiles',fontsize=8)
for ext in ['pdf','svg','png']:fig.savefig(O/('Figure4.'+ext),dpi=600 if ext=='png' else 150,facecolor='white')
fig.savefig(O/'Figure4_preview.png',dpi=160,facecolor='white');plt.close(fig)
(D/'Figure4_specification.json').write_text(json.dumps({'width_mm':185,'height_mm':172,'font_min_pt':8,'font_max_pt':11,'candidates':candidates,'primary_q':.1,'secondary_q':.05},indent=2)+'\n')
print('FIGURE_READY',O,flush=True)
