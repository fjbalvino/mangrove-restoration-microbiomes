"""Cycle-labelled scientific figures from unchanged numerical results."""
from pathlib import Path
from io import BytesIO
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, NullLocator
from marker_annotations import ORDER, PROVISIONAL, ANNOTATIONS, display_rows, export_guide
from output_io import write_binary

ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results';F=ROOT/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'axes.titlesize':10,'xtick.labelsize':8.5,'ytick.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','mathtext.default':'regular'})
LOC={'Carmen':'#5266A4','Cozumel':'#27806A','Tuxpan':'#C47A32'}
DEP={20:'#C47A32',40:'#246D8C'}
ROWS=display_rows();Y={name:i for i,(kind,name) in enumerate(ROWS) if kind=='gene'}

def save(fig,name):
    for fmt in ['png','pdf','svg']:
        stream=BytesIO()
        fig.savefig(stream,format=fmt,dpi=600 if fmt=='png' else 300,bbox_inches='tight',pad_inches=.08)
        write_binary(F/(name+'.'+fmt),stream.getvalue())
    plt.close(fig)

def labels_axis(ax):
    ax.set_xlim(0,1);ax.set_ylim(len(ROWS)-.5,-.55);ax.axis('off')
    for i,(kind,name) in enumerate(ROWS):
        if kind=='header':
            ax.text(.01,i,name+' cycle',va='center',fontsize=10,fontweight='bold')
        else:
            ax.text(.035,i,name+(' †' if name in PROVISIONAL else ''),va='center',fontstyle='italic',fontsize=9.5)
            ax.text(.29,i,ANNOTATIONS[name]['label'],va='center',fontsize=9,linespacing=1.05)

def main():
    export_guide(ROOT)
    profiles=pd.read_csv(R/'profile_marker_values.tsv',sep='\t')
    tests=pd.read_csv(R/'depth_effects_for_reporting.tsv',sep='\t')
    fig=plt.figure(figsize=(7.3,8.15))
    top=fig.add_gridspec(1,3,left=.105,right=.985,bottom=.690,top=.895,wspace=.40)
    for i,gene in enumerate(['rbcL','nirS','dsrA']):
        ax=fig.add_subplot(top[0,i]);df=profiles[profiles.gene.eq(gene)]
        vals=np.exp2(df[['depth_5_log2','depth_20_log2','depth_40_log2']].to_numpy())
        for (_,row),v in zip(df.iterrows(),vals):
            ax.plot([5,20,40],v,color=LOC[row.locality],alpha=.70,lw=.8,marker='o',ms=2.6)
        ax.set_yscale('log');ax.set_xticks([5,20,40]);ax.set_xlabel('Depth (cm)',labelpad=4)
        ax.set_title(ANNOTATIONS[gene]['cycle']+' cycle\n'+r'$\it{'+gene+'}$',fontsize=10,pad=8)
        candidates=np.array([.001,.003,.01,.03,.1,.3,1.])
        ticks=candidates[(candidates>=vals.min()*.75)&(candidates<=vals.max()*1.25)]
        if len(ticks)>4:ticks=ticks[::2]
        ax.set_yticks(ticks);ax.yaxis.set_major_formatter(FuncFormatter(lambda v,_:f'{v:g}'));ax.yaxis.set_minor_locator(NullLocator())
        ax.grid(axis='y',alpha=.14);ax.set_xlim(2,43)
        if i==0:ax.set_ylabel('Marker/reference ratio\n(log scale)',fontsize=9)
    handles=[Line2D([0],[0],color=c,marker='o',ms=3,lw=.9,label=l) for l,c in LOC.items()]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.55,.995),ncol=3,frameon=False,fontsize=9)
    fig.text(.02,.945,'a',fontsize=13,fontweight='bold')
    fig.text(.02,.618,'b',fontsize=13,fontweight='bold')
    low=fig.add_gridspec(1,2,left=.025,right=.985,bottom=.115,top=.605,width_ratios=[.47,.53],wspace=.035)
    lab=fig.add_subplot(low[0,0]);labels_axis(lab)
    ax=fig.add_subplot(low[0,1]);ax.set_ylim(len(ROWS)-.5,-.55)
    ax.set_xscale('log',base=2);ax.axvline(1,color='#777777',lw=.9,ls='--',zorder=0)
    for i,(kind,name) in enumerate(ROWS):
        if kind=='header':
            ax.axhline(i,color='#DDDDDD',lw=.8,zorder=0);continue
        for depth,offset in [(20,-.17),(40,.17)]:
            row=tests[tests.gene.eq(name)&tests.test.eq(f'depth_{depth}_vs_5')].iloc[0];y=i+offset
            ax.plot([row.ratio_ci_low,row.ratio_ci_high],[y,y],color=DEP[depth],lw=1.5,zorder=2)
            ax.scatter(row.median_ratio,y,s=26,edgecolor=DEP[depth],facecolor='white' if name in PROVISIONAL else DEP[depth],linewidth=.95,zorder=3)
            if row.p_holm_33<.05:ax.annotate('*',(row.median_ratio,y),xytext=(0,3.5),textcoords='offset points',ha='center',va='bottom',color=DEP[depth],fontweight='bold',fontsize=10)
    ax.set_xlim(.055,2.25);ax.set_xticks([.0625,.125,.25,.5,1,2],['1/16','1/8','1/4','1/2','1','2'])
    ax.set_yticks([]);ax.set_xlabel('Ratio relative to 5 cm (1 = no change)',fontsize=9,labelpad=6)
    ax.grid(axis='x',alpha=.14);ax.spines['left'].set_visible(False)
    h=[Line2D([0],[0],color=DEP[x],marker='o',lw=1.3,ms=4,label=f'{x} vs 5 cm') for x in [20,40]]
    ax.legend(handles=h,loc='lower right',bbox_to_anchor=(1,1.015),ncol=2,frameon=False,fontsize=8.5,handlelength=1.3,columnspacing=1.1)
    fig.text(.025,.049,'17 paired profiles  •  95.1% pointwise intervals  •  * Holm P < 0.05 across 33 tests',fontsize=8)
    fig.text(.025,.027,'† Unresolved fractional assignment. Functions indicate annotated potential, not measured rates.',fontsize=8)
    save(fig,'Figura_CNS_profundidad')

    alltests=pd.read_csv(R/'all_primary_and_sensitivity_tests.tsv',sep='\t')
    variants=['primary_76','fractional_76','strict_reference_56','both_mates_76','supported_fractional_76']
    a=alltests[alltests.test.eq('depth_40_vs_5')]
    values=a.pivot(index='gene',columns='variant',values='median_log2').loc[ORDER,variants]
    ps=a.pivot(index='gene',columns='variant',values='p_holm_33').loc[ORDER,variants]
    matrix=np.full((len(ROWS),5),np.nan)
    for gene in ORDER:matrix[Y[gene]]=values.loc[gene]
    cmap=LinearSegmentedColormap.from_list('blue_white_ochre',['#246D8C','#F8F8F5','#C47A32']);cmap.set_bad('white')
    fig=plt.figure(figsize=(7.3,6.55))
    grid=fig.add_gridspec(1,2,left=.025,right=.985,bottom=.235,top=.875,width_ratios=[.47,.53],wspace=.035)
    lab=fig.add_subplot(grid[0,0]);labels_axis(lab)
    ax=fig.add_subplot(grid[0,1]);im=ax.imshow(matrix,cmap=cmap,norm=TwoSlopeNorm(vmin=-2,vcenter=0,vmax=2),aspect='auto')
    ax.set_xticks(range(5),['Primary','Fractional\nlabels','56 reference\ngenes','Both\nmates','Broader\nannotation'],fontsize=8)
    ax.xaxis.tick_top();ax.tick_params(axis='x',pad=5)
    for gene in ORDER:
        for j in range(5):
            v=values.loc[gene].iloc[j]
            ax.text(j,Y[gene],f'{np.exp2(v):.2f}'+('*' if ps.loc[gene].iloc[j]<.05 else ''),ha='center',va='center',fontsize=8.8,color='white' if abs(v)>1.25 else '#222222')
    ax.set_yticks([]);ax.tick_params(length=0);ax.spines[:].set_visible(False)
    for i,(kind,name) in enumerate(ROWS):
        if kind=='header':ax.axhline(i,color='#DDDDDD',lw=.8)
    fig.text(.025,.973,'Sensitivity of depth patterns to quantification',fontsize=11,fontweight='bold')
    fig.text(.025,.94,'Median 40/5 cm ratios across the same 17 profiles',fontsize=9)
    cbax=fig.add_axes([.50,.145,.40,.022]);cbar=fig.colorbar(im,cax=cbax,orientation='horizontal')
    cbar.set_ticks([-2,-1,0,1,2],labels=['1/4','1/2','1','2','4'])
    cbar.set_label('Representation relative to 5 cm (ratio)',fontsize=8.5,labelpad=3)
    fig.text(.025,.058,'* Holm P < 0.05 within each 33-test family.  † Unresolved fractional assignment.',fontsize=8)
    fig.text(.025,.031,'Alternative analyses of the same samples provide sensitivity checks, not independent replication.',fontsize=8)
    save(fig,'Figura_CNS_sensibilidad')
    print('Both cycle/function-labelled figures exported as PNG, PDF and SVG')

if __name__=='__main__':main()

