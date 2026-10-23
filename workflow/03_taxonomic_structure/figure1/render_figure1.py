#!/usr/bin/env python3
"""Render validated Figure 1 B-D; optionally retain the supplied raster panel A."""
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
from matplotlib.ticker import FuncFormatter, MaxNLocator
from PIL import Image
from audit_figure1 import atomic, STAGES, COMMUNITIES

COLORS=dict(zip(STAGES,['#E41A1C','#FDAE61','#FFF7A8','#A6DDA0','#2B83BA']))
SHORT=['Degr.','Early','Inter.','Adv.','Cons.']
LOCS={'Carmen':'o','Cozumel':'s','Tuxpan':'^'}
DEPTH_SIZES={5:19,20:33,40:51}
COMPONENTS={'restoration4':('Restoration stage','#43886D'),
            'ECI_PC1':('Environment','#C47A32'),
            'locality':('Locality','#246D8C')}

def export(fig,path):
    for ext in ['png','pdf','svg']:
        buffer=io.BytesIO()
        metadata={'Creator':'','Date':None} if ext=='svg' else {'Creator':'','Producer':'','CreationDate':None} if ext=='pdf' else {'Software':'Scientific figure'}
        fig.savefig(buffer,format=ext,dpi=600,facecolor='white',metadata=metadata)
        atomic(path.with_suffix('.'+ext),buffer.getvalue())
    with Image.open(path.with_suffix('.png')) as im: im.verify()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--audit',required=True)
    ap.add_argument('--output',required=True);ap.add_argument('--source-figure')
    ap.add_argument('--panel-a-bottom',type=int,default=460)
    args=ap.parse_args();root=Path(args.audit);out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    import json
    assert json.loads((root/'audit_summary.json').read_text())['all_checks_passed']
    hp=pd.read_csv(root/'hill_profiles.csv');em=pd.read_csv(root/'hill_adjusted_means.csv')
    scores=pd.read_csv(root/'pcoa_scores.csv');tests=pd.read_csv(root/'composition_tests.csv')
    vp=pd.read_csv(root/'variance_partition_comparison.csv');vp=vp[vp.version.eq('current')]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,
                         'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,
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
    text(.57,5.48,'Circles: profiles     Diamonds and bars: locality-adjusted means and 95% CIs',size=7.4,color='#4D545B')
    panel('c',.13,5.12,'Taxonomic composition')
    for k,community in enumerate(COMMUNITIES):
        x=.57+k*3.52
        d=scores[scores.community.eq(community)]
        name='Complete community' if k==0 else 'Candidate bimodal subset'
        count='12,761 taxa' if k==0 else '1,031 taxa · exploratory'
        text(x,4.85,name,weight='bold',size=9)
        text(x,4.67,count,size=7.7,color='#4D545B')
        ax=axbox(x,2.19,2.64,2.29)
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
        text(x,1.77,f'Stage: R² = {row.R2:.3f}; F = {row.pseudo_F:.2f}; P = {row.p_value:.3f}'+(' (nominal)' if k else ''),size=7.1)
    stage_handles=[Line2D([],[],marker='o',color='none',mfc=COLORS[s],mec='#4D545B',mew=.5,ms=5,label=s.replace(' restoration','')) for s in STAGES]
    fig.legend(handles=stage_handles,loc='center',bbox_to_anchor=(.52,1.47/height),ncol=5,
               frameon=False,handletextpad=.4,columnspacing=1.2,fontsize=7.3)
    loc_handles=[Line2D([],[],marker=m,color='none',mfc='#DBE0E3',mec='#4D545B',ms=5,label=l) for l,m in LOCS.items()]
    fig.legend(handles=loc_handles,loc='center left',bbox_to_anchor=(.055,1.19/height),ncol=3,
               frameon=False,handletextpad=.4,columnspacing=1.0,fontsize=7.1)
    depth_handles=[Line2D([],[],marker='o',color='none',mfc='#DBE0E3',mec='#4D545B',mew=.5,
                          ms=np.sqrt(s),label=f'{d} cm') for d,s in DEPTH_SIZES.items()]
    fig.legend(handles=depth_handles,loc='center right',bbox_to_anchor=(.97,1.19/height),ncol=3,
               frameon=False,handletextpad=.4,columnspacing=1.0,fontsize=7.1)
    # Append D below B/C in a second canvas, keeping physical type sizes.
    dheight=2.15
    fd=plt.figure(figsize=(width,dheight))
    fd.text(.13/width,1.90/dheight,'d',weight='bold',size=12)
    fd.text(.40/width,1.90/dheight,'Unique contributions to whole-community variation',weight='bold',size=10)
    ax=fd.add_axes([.71/width,.43/dheight,6.22/width,1.15/dheight])
    ax.axvline(0,color='#808A92',lw=.8,zorder=2)
    ax.grid(axis='x',color='#E5E8EA',lw=.5,zorder=0)
    for j,(field,(label,color)) in enumerate(COMPONENTS.items()):
        for k,depth in enumerate([5,20,40]):
            value=vp[(vp.depth_cm==depth)&(vp.component==field)].percentage_total_variation.iloc[0]
            yy=k+(j-1)*.24
            ax.barh(yy,value,height=.19,color=color,zorder=3)
            ax.text(value+(.4 if value>=0 else -.35),yy,f'{value:.1f}',va='center',ha='left' if value>=0 else 'right',size=7.1,color=color)
    ax.set_yticks([0,1,2],['5 cm','20 cm','40 cm']);ax.invert_yaxis()
    ax.set_xlim(-4,32);ax.set_xticks([0,10,20,30]);ax.set_xlabel('Unique adjusted variation (% of total)',labelpad=4)
    ax.tick_params(axis='y',length=0,pad=5);ax.spines['left'].set_visible(False)
    handles=[Line2D([],[],color=c,lw=5,label=l) for l,c in COMPONENTS.values()]
    fd.legend(handles=handles,loc='center',bbox_to_anchor=(.53,1.69/dheight),ncol=3,frameon=False,fontsize=7.5,columnspacing=2)
    # Crop the unused bottom of the B/C canvas; retain vector objects in SVG/PDF
    # by composing SVGs rather than rasterizing B-D.
    export(fig,out/'Figure1_BC');export(fd,out/'Figure1_D')
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
        for ext in ['png','pdf','svg']:(out/(stem+'.'+ext)).unlink()
    if args.source_figure:
        import base64
        source=Image.open(args.source_figure).convert('RGB')
        crop=source.crop((0,0,source.width,args.panel_a_bottom))
        a_buf=io.BytesIO();crop.save(a_buf,format='PNG')
        aw=width*72; ah=aw*crop.height/crop.width
        complete=ET.Element('{'+ns+'}svg',{'width':f'{aw}pt','height':f'{ah+final_height}pt','viewBox':f'0 0 {aw} {ah+final_height}'})
        ET.SubElement(complete,'{'+ns+'}rect',{'width':'100%','height':'100%','fill':'white'})
        ET.SubElement(complete,'{'+ns+'}image',{'width':str(aw),'height':str(ah),
                      'href':'data:image/png;base64,'+base64.b64encode(a_buf.getvalue()).decode()})
        group=ET.SubElement(complete,'{'+ns+'}g',{'transform':f'translate(0 {ah})'})
        for child in svg:group.append(child)
        data=ET.tostring(complete,encoding='utf-8',xml_declaration=True)
        atomic(out/'Figure1_complete.svg',data)
        atomic(out/'Figure1_complete.pdf',cairosvg.svg2pdf(bytestring=data))
        atomic(out/'Figure1_complete.png',cairosvg.svg2png(bytestring=data,output_width=4380))
    for name in ['Figure1_B-D','Figure1_complete']:
        path=out/(name+'.png')
        if path.exists():
            with Image.open(path) as im:
                im.load();buffer=io.BytesIO();im.save(buffer,format='PNG',dpi=(600,600))
                atomic(path,buffer.getvalue())
                preview=im.copy();preview.thumbnail((1000,1600));preview.save(out/(name+'_preview.png'))
                print(name,im.size)
    plt.close('all')

if __name__=='__main__':main()
