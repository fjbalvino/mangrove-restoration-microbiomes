#!/usr/bin/env python3
"""Assemble manuscript figures from recalculated tables at final print size.
No statistical refitting and no illustrative/artificial data.
"""
from pathlib import Path
import json, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import f as f_dist

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.labelsize':9,'axes.titlesize':10,
 'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8.5,'axes.linewidth':0.6,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
STAGES=['Degraded','Early restoration','Intermediate restoration','Advanced restoration','Conserved']
SHORT=['Degraded','Early','Intermediate','Advanced','Conserved']
COLORS=dict(zip(STAGES,['#C54C43','#D59127','#8870B3','#26867B','#3273A5']))
MARKERS={'Carmen':'o','Cozumel':'s','Tuxpan':'^'}
SIZES={5:16,20:29,40:46}
MM=1/25.4

def clean(ax):
    ax.spines[['top','right']].set_visible(False)
    for sp in ['left','bottom']:ax.spines[sp].set_color('#9BA1A6')
    ax.tick_params(length=2.5,width=.6,pad=2)

def legend(fig,y=.974):
    handles=[Line2D([],[],marker='o',linestyle='',markerfacecolor=COLORS[s],markeredgecolor='white',markersize=6,label=l) for s,l in zip(STAGES,SHORT)]
    fig.legend(handles=handles,ncol=5,loc='upper center',bbox_to_anchor=(.52,y),frameon=False,columnspacing=1.35,handletextpad=.3,borderaxespad=0)

def ellipse(ax,xy,color):
    # Same normal ellipse convention as ggplot2::stat_ellipse(type='norm',level=.68).
    xy=np.asarray(xy);n=len(xy)
    cov=np.cov(xy,rowvar=False);vals,vecs=np.linalg.eigh(cov);ix=vals.argsort()[::-1];vals=vals[ix];vecs=vecs[:,ix]
    radius=np.sqrt(2*f_dist.ppf(.68,2,n-1))
    angle=np.degrees(np.arctan2(vecs[1,0],vecs[0,0]))
    ax.add_patch(Ellipse(xy.mean(0),*(2*radius*np.sqrt(vals)),angle=angle,facecolor=color,edgecolor=color,alpha=.085,lw=.8,zorder=1))
    ax.add_patch(Ellipse(xy.mean(0),*(2*radius*np.sqrt(vals)),angle=angle,facecolor='none',edgecolor=color,alpha=.65,lw=.7,zorder=2))

def verify_fonts(fig,name,out):
    from matplotlib.text import Text
    fig.canvas.draw()
    sizes=[t.get_fontsize() for t in fig.findobj(Text) if t.get_visible() and t.get_text()]
    assert min(sizes)>=8-1e-6 and max(sizes)<=11+1e-6,(name,min(sizes),max(sizes))
    return dict(figure=name,width_mm=round(fig.get_figwidth()/MM,3),height_mm=round(fig.get_figheight()/MM,3),minimum_pt=min(sizes),maximum_pt=max(sizes),text_objects=len(sizes))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out_root',type=Path,default=Path(__file__).resolve().parents[2]/'local_runs/functional_reference')
    a=ap.parse_args()
    r=a.out_root;out=r/'figures';out.mkdir(parents=True,exist_ok=True)
    sc=pd.read_csv(r/'603_partial_dbRDA/tables/603_sample_scores.tsv',sep='\t');sc=sc[sc.pseudocount==1]
    axvar=pd.read_csv(r/'603_partial_dbRDA/tables/603_axis_variance.tsv',sep='\t');axvar=axvar[axvar.pseudocount==1]
    tests=pd.read_csv(r/'603_partial_dbRDA/tables/603_partial_dbRDA_results.tsv',sep='\t').query('pseudocount==1').set_index('ontology')
    delta=pd.read_csv(r/'604_depth_reference_contrasts/tables/604_distances_and_local_contrasts.tsv',sep='\t')
    fontlog=[]
    fig=plt.figure(figsize=(185*MM,230*MM))
    legend(fig,.982)
    loc_h=[Line2D([],[],marker=m,linestyle='',markerfacecolor='white',markeredgecolor='#3C4145',markersize=5,label=l) for l,m in MARKERS.items()]
    fig.legend(handles=loc_h,ncol=3,loc='upper left',bbox_to_anchor=(.11,.95),frameon=False,columnspacing=1.25,handletextpad=.3,borderaxespad=0)
    dep_h=[Line2D([],[],marker='o',linestyle='',markerfacecolor='#ADB4B9',markeredgecolor='#555555',markersize=np.sqrt(SIZES[d]),label=f'{d} cm') for d in SIZES]
    fig.legend(handles=dep_h,ncol=3,loc='upper right',bbox_to_anchor=(.95,.95),frameon=False,columnspacing=.9,handletextpad=.3,borderaxespad=0)
    for col,ontology in enumerate(['KEGG_ko','PFAM']):
        left=.10+.50*col
        ax=fig.add_axes([left,.59,.35,.285]);clean(ax)
        z=sc[sc.ontology==ontology]
        for st in STAGES:ellipse(ax,z.loc[z.restoration4==st,['axis1','axis2']],COLORS[st])
        for row in z.itertuples():
            ax.scatter(row.axis1,row.axis2,s=SIZES[row.depth_cm],marker=MARKERS[row.locality],facecolor=COLORS[row.restoration4],edgecolor='#363B3F',linewidth=.45,zorder=3,alpha=.95)
        ax.axhline(0,color='#DDDFE1',lw=.55,zorder=0);ax.axvline(0,color='#DDDFE1',lw=.55,zorder=0)
        ax.set_aspect('equal',adjustable='datalim')
        ev=axvar[axvar.ontology==ontology].set_index('axis').percent_constrained
        ax.set_xlabel(f'dbRDA1 ({ev["dbRDA1"]:.1f}%)',labelpad=3)
        ax.set_ylabel(f'dbRDA2 ({ev["dbRDA2"]:.1f}%)',labelpad=3)
        t=tests.loc[ontology];layer='KO' if col==0 else 'PFAM'
        fig.text(left-.07,.900,'ab'[col],fontsize=11,fontweight='bold')
        fig.text(left,.900,f'{layer} composition',fontsize=10,fontweight='bold')
        fig.text(left,.881,f'{int(t.n_features):,} features  |  Stage R² = {100*t.r2_unique_fraction_of_total:.1f}%',fontsize=8.5)
        fig.text(left,.535,f'Stage P = {t.p_value:.4f}  |  Dispersion P = {t.prior_dispersion_p:.4f}',fontsize=8)
        fig.text(left,.499,f'{layer} distance contrasts',fontsize=10,fontweight='bold')
        fig.text(left-.07,.499,'cd'[col],fontsize=11,fontweight='bold')
        fig.text(left-.066,.294,'Δ distance to Conserved',rotation=90,fontsize=9,va='center')
        for i,depth in enumerate([5,20,40]):
            for j,loc in enumerate(MARKERS):
                ax=fig.add_axes([left+j*.119,.365-i*.116,.105,.096]);clean(ax)
                ax.axhspan(-50,0,color='#F1F6F8',zorder=0)
                ax.axhline(0,color='#939A9F',lw=.65,ls=(0,(3,2)),zorder=1)
                group=delta[(delta.ontology==ontology)&(delta.locality==loc)&(delta.depth_cm==depth)]
                for row in group.itertuples():
                    ax.scatter(row.x,row.delta_distance,s=25,marker=MARKERS[loc],facecolor=COLORS[row.restoration4],edgecolor='#343A3F',linewidth=.45,zorder=3)
                ax.set_xlim(.65,4.35);ax.set_ylim(-50,38);ax.set_yticks([-40,-20,0,20])
                ax.set_xticks([1,2,3,4]);ax.set_xticklabels(['D','E','I','A'] if i==2 else [])
                if j!=0:ax.set_yticklabels([]);ax.tick_params(axis='y',length=0);ax.spines['left'].set_visible(False)
                if i==0:ax.set_title(loc,fontsize=8.5,pad=5)
                if j==2:ax.text(1.06,.5,f'{depth} cm',transform=ax.transAxes,rotation=90,va='center',fontsize=8.5)
        fig.text(left+.175,.096,'Restoration stage',ha='center',fontsize=9)
    fig.text(.52,.063,'D  Degraded     E  Early     I  Intermediate     A  Advanced',ha='center',fontsize=8)
    fig.text(.52,.040,'Negative Δ: closer to Conserved than the local Degraded baseline',ha='center',fontsize=8.5)
    fig.text(.52,.020,'51 samples · 17 profiles · matched references by locality and depth',ha='center',fontsize=8)
    fontlog.append(verify_fonts(fig,'Figure 3',out))
    for ext in ['pdf','svg','png']:fig.savefig(out/f'Figure_3.{ext}',dpi=400)
    plt.close(fig)

    cov=pd.read_csv(r/'605_reference_core/core_coverage_by_sample.tsv',sep='\t')
    cp=pd.read_csv(r/'605_reference_core/core_distances_by_profile.tsv',sep='\t')
    csum=pd.read_csv(r/'605_reference_core/core_summary.tsv',sep='\t').set_index('layer')
    plt.rcParams.update({'font.size':9,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9})
    with PdfPages(out/'Supplementary_Figure_9.pdf') as pdf:
        fig=plt.figure(figsize=(185*MM,220*MM));legend(fig,.975)
        fig.text(.07,.92,'Conserved reference core',fontsize=11,fontweight='bold')
        for col,layer in enumerate(['KO','PFAM']):
            left=.11+.49*col;z=cov[cov.layer==layer]
            ax=fig.add_axes([left,.56,.34,.28]);clean(ax)
            for k,st in enumerate(STAGES):
                g=z[z.restoration4==st].sort_values(['locality','profile_id','depth_cm'])
                xs=k+np.linspace(-.20,.20,len(g))
                for x,row in zip(xs,g.itertuples()):
                    ax.scatter(x,row.coverage_percent,s=SIZES[row.depth_cm],marker=MARKERS[row.locality],facecolor=COLORS[st],edgecolor='#41474B',linewidth=.4,alpha=.8)
                ax.plot([k-.23,k+.23],[g.coverage_percent.median()]*2,color='black',lw=1)
            ax.set_ylim(85,100.8);ax.set_yticks([85,90,95,100]);ax.set_xlim(-.5,4.5)
            ax.set_xticks(range(5));ax.set_xticklabels(['D','E','I','A','C']);ax.set_xlabel('Restoration stage')
            ax.set_ylabel('Core detected (%)' if col==0 else '')
            if col:ax.set_yticklabels([])
            fig.text(left-.07,.875,'ab'[col],fontsize=11,fontweight='bold');fig.text(left,.875,f'{layer}  |  {int(csum.loc[layer,"core_features"]):,} core features',fontsize=10,fontweight='bold')
            fig.text(left,.855,f'Degraded median: {csum.loc[layer,"degraded_median_percent"]:.2f}%',fontsize=9)
            fig.text(left-.07,.452,'cd'[col],fontsize=11,fontweight='bold');fig.text(left,.452,f'{layer} core proximity',fontsize=10,fontweight='bold')
            for j,loc in enumerate(MARKERS):
                ax=fig.add_axes([left+j*.117,.18,.105,.22]);clean(ax)
                ax.axhline(0,color='#969EA3',lw=.7,ls=(0,(3,2)))
                zz=cp[(cp.layer==layer)&(cp.locality==loc)]
                for k,st in enumerate(STAGES[:4]):
                    g=zz[zz.restoration4==st].sort_values('profile_id');offset=np.linspace(-.12,.12,len(g)) if len(g)>1 else [0]
                    for x,row in zip(k+np.asarray(offset),g.itertuples()):ax.scatter(x,row.distance_reduction,s=28,marker=MARKERS[loc],facecolor=COLORS[st],edgecolor='#343A3F',linewidth=.45)
                ax.set_xticks(range(4));ax.set_xticklabels(['D','E','I','A']);ax.set_xlim(-.5,3.5);ax.set_ylim(-4,6);ax.set_yticks([-4,-2,0,2,4,6]);ax.set_title(loc,fontsize=9)
                if j:ax.set_yticklabels([]);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
            fig.text(left-.065,.29,'Mean distance reduction',rotation=90,va='center',fontsize=9)
        fig.text(.52,.12,'Positive reduction: closer to Conserved than Degraded',ha='center',fontsize=9)
        fig.text(.52,.085,'Core: ≥0.01% relative abundance in all nine Conserved samples',ha='center',fontsize=9)
        fig.text(.52,.055,'Top: sample coverage and stage medians. Bottom: three-depth profile means.',ha='center',fontsize=9)
        fig.text(.52,.030,'D Degraded · E Early · I Intermediate · A Advanced · C Conserved',ha='center',fontsize=9)
        fontlog.append(verify_fonts(fig,'Supplementary Figure 9 page 1',out));pdf.savefig(fig);fig.savefig(out/'Supplementary_Figure_9_page1.png',dpi=300);plt.close(fig)
        intersections=pd.read_csv(r/'605_reference_core/stage_core_exact_intersections.tsv',sep='\t',dtype={'membership_pattern':str})
        sets=pd.read_csv(r/'605_reference_core/stage_core_sizes.tsv',sep='\t')
        fig=plt.figure(figsize=(185*MM,220*MM))
        fig.text(.05,.972,'Overlap of stage cores',fontsize=11,fontweight='bold')
        fig.text(.05,.947,'Exact intersections; ≥0.01% relative abundance in every sample of each stage',fontsize=9)
        for i,layer in enumerate(['KO','PFAM']):
            top=.84-i*.43
            z=intersections[intersections.layer==layer].sort_values(['features','membership_pattern'],ascending=[False,True]).reset_index(drop=True)
            xs=np.arange(len(z));ax=fig.add_axes([.25,top-.12,.70,.17]);clean(ax)
            colors=['#3273A5' if k=='11111' else '#8196A2' for k in z.membership_pattern]
            ax.bar(xs,z.features,color=colors,width=.70);ax.set_yscale('log');ax.set_ylim(.8,2500);ax.set_yticks([1,10,100,1000]);ax.set_yticklabels(['1','10','100','1,000']);ax.set_xlim(-.7,len(z)-.3);ax.set_xticks([]);ax.set_ylabel('Functions (log scale)',fontsize=9)
            for k in range(3):ax.text(k,z.features.iloc[k]*1.18,f'{z.features.iloc[k]:,}',ha='center',fontsize=9)
            fig.text(.05,top+.065,'ef'[i],fontsize=11,fontweight='bold');fig.text(.10,top+.065,layer,fontsize=10,fontweight='bold')
            dm=fig.add_axes([.25,top-.285,.70,.145]);dm.set_xlim(ax.get_xlim());dm.set_ylim(4.6,-.6);dm.axis('off')
            for j,st in enumerate(STAGES):
                if j%2==0:dm.axhspan(j-.5,j+.5,color='#F4F5F6')
                n=int(sets[(sets.layer==layer)&(sets.stage==st)].core_features.iloc[0])
                dm.text(-1.5,j,f'{SHORT[j]}  ({n:,})',ha='right',va='center',fontsize=9)
            for k,pat in enumerate(z.membership_pattern):
                inc=np.array([j for j,v in enumerate(pat.zfill(5)) if v=='1'])
                dm.scatter([k]*5,range(5),s=7,color='#D7DDE1',zorder=2)
                dm.plot([k,k],[inc.min(),inc.max()],color=colors[k],lw=.8,zorder=3)
                dm.scatter([k]*len(inc),inc,s=12,color=colors[k],zorder=4)
        fig.text(.50,.03,'Stage-core sizes are reported in Supplementary Table 14.',ha='center',fontsize=9)
        fontlog.append(verify_fonts(fig,'Supplementary Figure 9 page 2',out));pdf.savefig(fig);fig.savefig(out/'Supplementary_Figure_9_page2.png',dpi=300);plt.close(fig)
    (out/'figure_font_audit.json').write_text(json.dumps(fontlog,indent=2)+'\n')
    print(json.dumps(fontlog,indent=2))

if __name__=='__main__':main()
