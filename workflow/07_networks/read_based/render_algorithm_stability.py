#!/usr/bin/env python3
"""Render algorithm stability without displaying or selecting on HI associations."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    r=args.run/'algorithm_comparison';out=args.output;out.mkdir(exist_ok=True,parents=True)
    B=pd.read_csv(r/'811_algorithm_bootstrap.tsv',sep='\t');D=pd.read_csv(r/'811_algorithm_profile_deletion.tsv',sep='\t');N=pd.read_csv(r/'811_algorithm_normalization.tsv',sep='\t')
    methods=['ridge_010','ledoit_wolf','oas'];labels=['Ridge 0.10','Ledoit–Wolf','OAS'];colors=['#788795','#2171b5','#c46a28']
    scenarios=['space_depth','space_depth_environment'];sl=['Locality + depth','+ 5 environmental axes']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,
        'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axs=plt.subplots(2,2,figsize=(7.3,6.2));fig.subplots_adjust(left=.10,right=.98,bottom=.20,top=.87,hspace=.60,wspace=.40)
    for k,(ax,df,col,title,ylabel) in enumerate(zip(axs.ravel(),[B,D,N,D],
        ['edge_jaccard','sample_edge_jaccard_median','sample_edge_jaccard_median','mean_weight_rank_correlation'],
        ['Bootstrap stability','Profile deletion stability','Abundance scaling sensitivity','Sample metric stability'],
        ['Aggregate edge overlap','Sample edge overlap','Sample edge overlap','Mean-weight rank correlation'])):
        for mi,method in enumerate(methods):
            for si,sid in enumerate(scenarios):
                v=df[(df.method==method)&(df.construction==sid)][col].to_numpy();x=si+(mi-1)*.21
                if len(v)>1:
                    low,median,high=np.quantile(v,[.025,.5,.975]);ax.plot([x,x],[low,high],color=colors[mi],lw=1.1)
                    ax.plot([x-.035,x+.035],[low,low],color=colors[mi],lw=.9);ax.plot([x-.035,x+.035],[high,high],color=colors[mi],lw=.9)
                else:median=v[0]
                ax.plot(x,median,'o',color=colors[mi],ms=5)
                ax.annotate(f'{median:.2f}',(x,median),xytext=(0,9 if k!=3 else -14),textcoords='offset points',ha='center',fontsize=6.5)
        ax.set_xticks([0,1],['L + D','L + D + E']);ax.set_xlim(-.43,1.43);ax.set_ylim(0,1.04);ax.set_yticks([0,.25,.5,.75,1]);ax.set_ylabel(ylabel)
        ax.set_title(title,loc='left',fontsize=9,fontweight='bold',pad=11);ax.text(-.22,1.15,'abcd'[k],transform=ax.transAxes,fontweight='bold',fontsize=12)
        ax.grid(axis='y',color='#e5e8eb',lw=.6);ax.set_axisbelow(True)
    fig.suptitle('Network stability across inference algorithms',x=.10,y=.97,ha='left',fontsize=11,fontweight='bold')
    fig.legend(handles=[Line2D([],[],color=c,marker='o',label=l,linestyle='') for c,l in zip(colors,labels)],
        loc='lower center',bbox_to_anchor=(.52,.095),ncol=3,frameon=False)
    fig.text(.10,.035,'Higher values indicate greater repeatability · 200 KOs · 995 selected edges\nPoints: medians; lines: 2.5th–97.5th percentiles of paired resampling diagnostics',fontsize=7.5,linespacing=1.6)
    for ext in ['png','pdf','svg']:fig.savefig(out/f'Network_algorithm_stability.{ext}',dpi=450,facecolor='white')
    fig.savefig(out/'Network_algorithm_stability_preview.png',dpi=150,facecolor='white');plt.close(fig)
    print('Rendered algorithm stability comparison',flush=True)

if __name__=='__main__':main()
