#!/usr/bin/env python3
"""Assemble Figure 2 from audited TITAN2 results; no decorative assets required."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

PANELS = [
    ('MHI_local', 'Environmental index (MHI)',
     'Farther from local\nconserved conditions', 'Closer to local\nconserved conditions'),
    ('moisture_stress_PC1', 'Moisture (NDMI)',
     'Lower NDMI\nmoisture signal', 'Higher NDMI\nmoisture signal'),
    ('nutrients_redox_PC1', 'Nutrients and sulfide',
     'Lower nutrients\nand sulfide', 'Higher nutrients\nand sulfide'),
    ('physicochemical_PC1', 'Physicochemical conditions',
     'Lower salinity\nhigher redox potential', 'Higher salinity\nlower redox potential'),
    ('vegetation_landscape_PC1', 'Vegetation and landscape',
     'Higher median NDVI\nlower heterogeneity', 'Lower median NDVI\ngreater heterogeneity'),
    ('water_inundation_PC1', 'Surface-water signal',
     'Stronger surface-water\nsignal', 'Weaker surface-water\nsignal'),
]
COLOURS = {'z-': '#D59127', 'z+': '#26867B'}
MARKERS = {'z-': 'o', 'z+': '^'}
STRATA = ['Global', '5', '20', '40']

def validate(data):
    required = ['stratum','gradient','direction','n_indicators','cp','q05','q10','q50','q90','q95']
    if not set(required).issubset(data):
        raise ValueError('Missing figure-data columns')
    expected = {(s,g,d) for s in STRATA for g, *_ in PANELS for d in COLOURS}
    actual = list(data[['stratum','gradient','direction']].itertuples(index=False,name=None))
    if len(actual) != 48 or len(set(actual)) != 48 or set(actual) != expected:
        raise ValueError('Figure requires exactly six gradients, four strata and two responses')
    if ((data.n_indicators < 0) | (data.n_indicators % 1 != 0)).any():
        raise ValueError('Indicator counts must be nonnegative integers')
    qcols = ['q05','q10','q50','q90','q95']
    kept = data.loc[data.n_indicators > 0]
    if not np.isfinite(kept[['cp']+qcols].to_numpy()).all():
        raise ValueError('Retained indicators require finite thresholds and quantiles')
    if (np.diff(kept[qcols].to_numpy(),axis=1) < -1e-12).any():
        raise ValueError('Bootstrap quantiles are not ordered')
    if data.loc[data.n_indicators == 0, ['cp']+qcols].notna().any().any():
        raise ValueError('Do not plot filtered thresholds when no indicator was retained')

def build_figure(data):
    validate(data)
    plt.rcParams.update({
        'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,
        'xtick.labelsize':7.5,'ytick.labelsize':8,'axes.linewidth':.55,
        'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
        'axes.unicode_minus':True,'savefig.facecolor':'white',
    })
    fig = plt.figure(figsize=(185/25.4,226/25.4))
    # Each gradient has one common axis for the profile mean and all depths.
    xleft = [.093,.582]; top = [.85,.565,.28]
    width=.325; height=.150
    for i,(g,title,low,high) in enumerate(PANELS):
        col=i%2; row=i//2; bottom=top[row]-height
        ax=fig.add_axes([xleft[col],bottom,width,height])
        tab=data.loc[data.gradient.eq(g)].copy()
        valid=tab.loc[tab.n_indicators.gt(0)]
        lower=min(valid.q05.min(),valid.cp.min()); upper=max(valid.q95.max(),valid.cp.max())
        pad=max((upper-lower)*.065,.015);ax.set_xlim(lower-pad,upper+pad)
        ax.set_ylim(3.5,-.6)
        ax.axhspan(-.47,.47,color='#F2F4F5',zorder=0)
        ax.axhline(.55,color='#D9DEE1',lw=.6,zorder=0)
        if lower-pad < 0 < upper+pad:ax.axvline(0,color='#A8ADB0',lw=.6,ls=(0,(3,3)),zorder=1)
        for y,st in enumerate(STRATA):
            group=tab.loc[tab.stratum.eq(st)]
            counts={r.direction:int(r.n_indicators) for r in group.itertuples()}
            ax.text(1.045,y,f"{counts['z-']} / {counts['z+']}",transform=ax.get_yaxis_transform(),
                    ha='left',va='center',fontsize=7.5,color='#30383C',clip_on=False)
            for d,dy in [('z-',-.12),('z+',.12)]:
                r=group.loc[group.direction.eq(d)].iloc[0]
                if r.n_indicators == 0:continue
                yy=y+dy;c=COLOURS[d]
                ax.plot([r.q05,r.q95],[yy,yy],lw=1.0,color=c,alpha=.40,solid_capstyle='butt',zorder=2)
                ax.plot([r.q10,r.q90],[yy,yy],lw=2.7,color=c,solid_capstyle='butt',zorder=3)
                ax.scatter(r.cp,yy,s=29 if d=='z-' else 36,marker=MARKERS[d],
                           facecolor=c,edgecolor='#263A42',linewidth=.65,zorder=4)
        ax.set_yticks(range(4),['Global','5 cm','20 cm','40 cm'])
        ax.tick_params(axis='y',length=0,pad=5)
        ax.tick_params(axis='x',length=3,width=.55,pad=3)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=3))
        ax.spines[['top','right','left']].set_visible(False)
        ax.spines['bottom'].set_color('#59636A')
        ax.grid(axis='x',color='#E8EBED',linewidth=.45,zorder=0)
        fig.text(xleft[col]-.067,top[row]+.025,chr(97+i),weight='bold',fontsize=11,va='bottom')
        fig.text(xleft[col],top[row]+.025,title,weight='bold',fontsize=9,va='bottom')
        ax.text(1.045,1.025,'n− / n+',transform=ax.transAxes,fontsize=7.5,ha='left',va='bottom')
        ylabels=bottom-.060
        ax.annotate('',xy=(.80,-.267),xytext=(.20,-.267),xycoords='axes fraction',
                    arrowprops=dict(arrowstyle='->',color='#53656F',lw=.65),annotation_clip=False)
        fig.text(xleft[col],ylabels,low,fontsize=7,ha='left',va='top',linespacing=1.15,color='#414D54')
        fig.text(xleft[col]+width,ylabels,high,fontsize=7,ha='right',va='top',linespacing=1.15,color='#414D54')
    handles=[Line2D([0],[0],marker=MARKERS[d],color=COLOURS[d],markersize=5.5,
                    markeredgecolor='#263A42',markeredgewidth=.65,lw=1.8,label=label)
             for d,label in [('z-','Decreasing (z−)'),('z+','Increasing (z+)')]]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.093,.991),
               frameon=False,ncol=2,columnspacing=2.5,handlelength=1.7,fontsize=8)
    fig.text(.093,.933,'Points: community change points     Lines: central 80% (thick) / 90% (thin) bootstrap ranges',
             fontsize=7.4,color='#414D54')
    fig.text(.093,.912,'Global: mean of three depths per profile     n− / n+: retained decreasing / increasing indicators',
             fontsize=7.4,color='#414D54')
    fig.text(.5,.028,'Environmental gradient value (original scale; increasing from left to right)',
             ha='center',fontsize=8.2)
    return fig

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); data=pd.read_csv(a.data,dtype={'stratum':str});fig=build_figure(data)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    for ext in ('png','pdf','svg'):
        fig.savefig(a.output.with_suffix('.'+ext),dpi=600)
    plt.close(fig)

if __name__=='__main__':main()
