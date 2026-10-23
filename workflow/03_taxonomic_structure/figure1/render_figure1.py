#!/usr/bin/env python3
"""Render validated Figure 1 panels B-D from audited numerical tables."""
from __future__ import annotations
import argparse
import io
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
from matplotlib.ticker import FuncFormatter, MaxNLocator
from scipy.stats import f as f_dist
from PIL import Image
from audit_figure1 import atomic, STAGES, COMMUNITIES

# Match the canonical Figure 3 renderer, 05_06_assemble_figures.py,
# git blob aeffa71527a927d89dcbeb4020d4188b8b3287be.
COLORS=dict(zip(STAGES,['#C54C43','#D59127','#8870B3','#26867B','#3273A5']))
ELLIPSE_LEVEL=.68
SHORT=['Degr.','Early','Inter.','Adv.','Cons.']
LOCS={'Carmen':'o','Cozumel':'s','Tuxpan':'^'}
DEPTH_SIZES={5:16,20:29,40:46}
COMPONENTS={'restoration4':('Restoration stage','#43886D'),
            'ECI_PC1':('Environment','#C47A32'),
            'locality':('Locality','#246D8C')}

def ellipse(ax,xy,color):
    """Descriptive normal data ellipse, matching Figure 3 (not a mean CI)."""
    xy=np.asarray(xy,dtype=float)
    assert xy.ndim==2 and xy.shape[1]==2 and len(xy)>2 and np.isfinite(xy).all()
    vals,vecs=np.linalg.eigh(np.cov(xy,rowvar=False))
    order=vals.argsort()[::-1];vals=vals[order];vecs=vecs[:,order]
    assert vals.min()>=-1e-10
    radius=np.sqrt(2*f_dist.ppf(ELLIPSE_LEVEL,2,len(xy)-1))
    angle=np.degrees(np.arctan2(vecs[1,0],vecs[0,0]))
    diameters=2*radius*np.sqrt(np.maximum(vals,0))
    ax.add_patch(Ellipse(xy.mean(0),*diameters,angle=angle,facecolor=color,
                         edgecolor=color,alpha=.085,lw=.8,zorder=1))
    ax.add_patch(Ellipse(xy.mean(0),*diameters,angle=angle,facecolor='none',
                         edgecolor=color,alpha=.65,lw=.7,zorder=2))

def display_fractions(vp):
    """Retain raw estimates; truncate negatives only for the 100% display."""
    shown=vp.copy()
    shown['nonnegative_unique_R2']=shown.unique_adjusted_R2.clip(lower=0)
    totals=shown.groupby('depth_cm').nonnegative_unique_R2.transform('sum')
    assert (totals>0).all()
    shown['relative_share_percent']=100*shown.nonnegative_unique_R2/totals
    assert np.allclose(shown.groupby('depth_cm').relative_share_percent.sum(),100)
    return shown

def export(fig,path,formats=('png','pdf','svg')):
    for ext in formats:
        buffer=io.BytesIO()
        metadata={'Creator':'','Date':None} if ext=='svg' else {'Creator':'','Producer':'','CreationDate':None} if ext=='pdf' else {'Software':'Scientific figure'}
        fig.savefig(buffer,format=ext,dpi=600,facecolor='white',metadata=metadata)
        atomic(path.with_suffix('.'+ext),buffer.getvalue())
    if 'png' in formats:
        with Image.open(path.with_suffix('.png')) as im: im.verify()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--audit',required=True)
    ap.add_argument('--output',required=True)
    args=ap.parse_args();root=Path(args.audit);out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    import json
    assert json.loads((root/'audit_summary.json').read_text())['all_checks_passed']
    hp=pd.read_csv(root/'hill_profiles.csv');em=pd.read_csv(root/'hill_adjusted_means.csv')
    scores=pd.read_csv(root/'pcoa_scores.csv');tests=pd.read_csv(root/'composition_tests.csv')
    vp=pd.read_csv(root/'variance_partition_comparison.csv');vp=vp[vp.version.eq('current')]
    shown=display_fractions(vp)
    atomic(root/'variance_partition_display.csv',shown.to_csv(index=False).encode())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.labelsize':9,
                         'axes.titlesize':10,'xtick.labelsize':8,'ytick.labelsize':8,
                         'axes.spines.top':False,'axes.spines.right':False,
                         'axes.linewidth':.6,'xtick.major.width':.5,'ytick.major.width':.5,
                         'xtick.major.size':2.5,'ytick.major.size':2.5,
                         'svg.fonttype':'none','pdf.fonttype':42})
    width,height=7.3,7.85
    fig=plt.figure(figsize=(width,height))
    def axbox(x,y,w,h):return fig.add_axes([x/width,y/height,w/width,h/height])
    def text(x,y,s,**kw):return fig.text(x/width,y/height,s,**kw)
    def panel(letter,x,y,title):
        text(x,y,letter,weight='bold',size=12,va='baseline')
        text(x+.27,y,title,size=10,weight='bold',va='baseline')
    panel('b',.13,7.59,'Diversity within sediment profiles')
    text(7.10,7.59,'Stage tests: all BH q = 0.968',ha='right',size=8,color='#444444')
    metrics=[('q0_richness','Richness'),('q1_exp_shannon','Exponential Shannon'),
             ('q2_inverse_simpson','Inverse Simpson')]
    for k,(metric,title) in enumerate(metrics):
        ax=axbox(.72+k*2.20,5.98,1.75,1.26)
        ax.set_title(title,pad=8,weight='bold',size=8.5)
        ax.grid(axis='y',color='#E5E8EA',lw=.55,zorder=0)
        for j,stage in enumerate(STAGES):
            rows=hp[hp.restoration4.eq(stage)].sort_values('profile_id')
            offsets=np.linspace(-.16,.16,len(rows))
            ax.scatter(j+offsets,rows[metric],s=20,c=COLORS[stage],edgecolors='#3B4248',linewidths=.5,zorder=3)
            row=em[em.response.eq(metric)&em.restoration4.eq(stage)].iloc[0]
            ax.errorbar(j,row.estimate,yerr=[[row.estimate-row.conf_low],[row.conf_high-row.estimate]],
                        fmt='D',markersize=4.5,mfc=COLORS[stage],mec='#20272C',mew=.75,
                        ecolor='#20272C',elinewidth=1.1,capsize=2,zorder=5)
        ax.set_xlim(-.48,4.48)
        ax.set_xticks(range(5),[f'{s}\n(n={n})' for s,n in zip(SHORT,[3,3,5,3,3])])
        ax.tick_params(axis='x',pad=4,length=0)
        ax.yaxis.set_major_locator(MaxNLocator(4,integer=True))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:,.0f}'))
        ax.margins(y=.13)
        if k==0:ax.set_ylabel('Profile-mean Hill diversity',labelpad=6)
    text(.57,5.48,'Circles: profiles     Diamonds and bars: locality-adjusted means and 95% CIs',size=8,color='#4D545B')
    panel('c',.13,5.12,'Taxonomic composition')
    for k,community in enumerate(COMMUNITIES):
        x=.57+k*3.52
        d=scores[scores.community.eq(community)]
        name='Complete community' if k==0 else 'Candidate bimodal subset'
        count='12,761 taxa' if k==0 else '1,031 taxa · exploratory'
        text(x,4.85,name,weight='bold',size=9)
        text(x,4.67,count,size=7.7,color='#4D545B')
        ax=axbox(x,2.31,2.64,2.23)
        for stage in STAGES:
            ellipse(ax,d.loc[d.restoration4.eq(stage),['PCoA1','PCoA2']],COLORS[stage])
        for stage in STAGES:
            for loc,marker in LOCS.items():
                sub=d[d.restoration4.eq(stage)&d.locality.eq(loc)]
                ax.scatter(sub.PCoA1,sub.PCoA2,s=[DEPTH_SIZES[v] for v in sub.depth_cm],
                           marker=marker,c=COLORS[stage],edgecolors='#303B43',lw=.55,alpha=.95,zorder=4)
        ax.axhline(0,c='#C6CCD0',lw=.65,zorder=1)
        ax.axvline(0,c='#C6CCD0',lw=.65,zorder=1)
        ax.grid(color='#E9ECEE',lw=.45)
        ax.set_aspect('equal',adjustable='datalim')
        ax.margins(.12)
        ax.xaxis.set_major_locator(MaxNLocator(5));ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.set_xlabel(f'PCoA1 ({100*d.PCoA1_variance_positive.iloc[0]:.1f}%)',labelpad=3)
        ax.set_ylabel(f'PCoA2 ({100*d.PCoA2_variance_positive.iloc[0]:.1f}%)',labelpad=3)
        row=tests[tests.community.eq(community)].iloc[0]
        text(x,1.77,f'Stage: R² = {row.R2:.3f}; F = {row.pseudo_F:.2f}; P = {row.p_value:.3f}'+(' (nominal)' if k else ''),size=8)
    stage_handles=[Line2D([],[],marker='o',color='none',mfc=COLORS[s],mec='#4D545B',mew=.5,ms=5,label=s.replace(' restoration','')) for s in STAGES]
    fig.legend(handles=stage_handles,loc='center',bbox_to_anchor=(.52,1.47/height),ncol=5,
               frameon=False,handletextpad=.4,columnspacing=1.2,fontsize=8)
    loc_handles=[Line2D([],[],marker=m,color='none',mfc='#DBE0E3',mec='#4D545B',ms=5,label=l) for l,m in LOCS.items()]
    fig.legend(handles=loc_handles,loc='center left',bbox_to_anchor=(.055,1.19/height),ncol=3,
               frameon=False,handletextpad=.4,columnspacing=1.0,fontsize=8)
    depth_handles=[Line2D([],[],marker='o',color='none',mfc='#DBE0E3',mec='#4D545B',mew=.5,
                          ms=np.sqrt(s),label=f'{d} cm') for d,s in DEPTH_SIZES.items()]
    fig.legend(handles=depth_handles,loc='center right',bbox_to_anchor=(.97,1.19/height),ncol=3,
               frameon=False,handletextpad=.4,columnspacing=1.0,fontsize=8)
    text(.57,.89,'Ellipses show stage distributions (68% normal ellipses)',size=8,color='#4D545B')
    # Append D below B/C in a second canvas, keeping physical type sizes.
    dheight=2.15
    fd=plt.figure(figsize=(width,dheight))
    fd.text(.13/width,1.90/dheight,'d',weight='bold',size=12)
    fd.text(.40/width,1.90/dheight,'Relative shares of unique variation',weight='bold',size=10)
    ax=fd.add_axes([.71/width,.43/dheight,6.22/width,1.15/dheight])
    ax.grid(axis='x',color='#E5E8EA',lw=.5,zorder=0)
    for k,depth in enumerate([5,20,40]):
        left=0.
        for field,(label,color) in COMPONENTS.items():
            value=shown[(shown.depth_cm==depth)&(shown.component==field)].relative_share_percent.iloc[0]
            ax.barh(k,value,left=left,height=.60,color=color,edgecolor='white',lw=.7,zorder=3)
            if value>0:
                ax.text(left+value/2,k,f'{value:.1f}%',va='center',ha='center',size=8,color='white',weight='bold')
            left+=value
        assert np.isclose(left,100)
    ax.set_yticks([0,1,2],['5 cm','20 cm','40 cm']);ax.invert_yaxis()
    ax.set_xlim(0,100);ax.set_xticks([0,25,50,75,100]);ax.set_xlabel('Relative share of unique adjusted variation (%)',labelpad=4)
    ax.tick_params(axis='y',length=0,pad=5);ax.spines['left'].set_visible(False)
    handles=[Line2D([],[],color=c,lw=5,label=l) for l,c in COMPONENTS.values()]
    fd.legend(handles=handles,loc='center',bbox_to_anchor=(.53,1.69/dheight),ncol=3,frameon=False,fontsize=8,columnspacing=2)
    # Crop the unused bottom of the B/C canvas; retain vector objects in SVG/PDF
    # by composing SVGs rather than rasterizing B-D.
    export(fig,out/'Figure1_BC',formats=('svg',));export(fd,out/'Figure1_D',formats=('svg',))
    import xml.etree.ElementTree as ET
    import cairosvg
    bc=ET.fromstring((out/'Figure1_BC.svg').read_bytes())
    dd=ET.fromstring((out/'Figure1_D.svg').read_bytes())
    ns='http://www.w3.org/2000/svg';ET.register_namespace('',ns)
    bc_height=7.0*72
    final_height=bc_height+dheight*72
    svg=ET.Element('{'+ns+'}svg',{'width':f'{width*72}pt','height':f'{final_height}pt','viewBox':f'0 0 {width*72} {final_height}'})
    ET.SubElement(svg,'{'+ns+'}rect',{'width':'100%','height':'100%','fill':'white'})
    for source,offset,prefix in [(bc,0,'bc'),(dd,bc_height,'d')]:
        # matplotlib IDs repeat between standalone figures; namespace them.
        xml=ET.tostring(source,encoding='unicode')
        import re
        ids=re.findall(r'id="([^"]+)"',xml)
        for old in sorted(ids,key=len,reverse=True):
            xml=xml.replace(f'id="{old}"',f'id="{prefix}_{old}"').replace(f'url(#{old})',f'url(#{prefix}_{old})').replace(f'href="#{old}"',f'href="#{prefix}_{old}"')
        source=ET.fromstring(xml)
        group=ET.SubElement(svg,'{'+ns+'}g',{'transform':f'translate(0 {offset})'})
        for child in source:group.append(child)
    svg_bytes=ET.tostring(svg,encoding='utf-8',xml_declaration=True)
    atomic(out/'Figure1_B-D.svg',svg_bytes)
    atomic(out/'Figure1_B-D.pdf',cairosvg.svg2pdf(bytestring=svg_bytes))
    atomic(out/'Figure1_B-D.png',cairosvg.svg2png(bytestring=svg_bytes,output_width=4380))
    # Intermediate standalone BC and D are not user deliverables.
    for stem in ['Figure1_BC','Figure1_D']:
        (out/(stem+'.svg')).unlink()
    path=out/'Figure1_B-D.png'
    with Image.open(path) as im:
        im.load();buffer=io.BytesIO();im.save(buffer,format='PNG',dpi=(600,600))
        atomic(path,buffer.getvalue())
        preview=im.copy();preview.thumbnail((1000,1600));preview.save(out/'Figure1_B-D_preview.png')
        print('Figure1_B-D',im.size)
    plt.close('all')

if __name__=='__main__':main()

