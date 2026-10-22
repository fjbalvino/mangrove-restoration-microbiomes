from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D

parser = argparse.ArgumentParser(description='Render Figure 5: functional enrichment and potential modules.')
parser.add_argument('--results',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args();SOURCE=args.results.resolve();OUT=args.output.resolve()
OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,
 'svg.fonttype':'none','axes.linewidth':.5,'axes.spines.top':False,'axes.spines.right':False})
blue = LinearSegmentedColormap.from_list('enrichment_blue',['#deebf7','#9ecae1','#4292c6','#08519c','#08306b'])
terms = pd.read_csv(SOURCE/'primary_enriched_terms_named.tsv',sep='\t')
ko = terms[terms.ontology.eq('KO')].copy()
modules = pd.read_csv(SOURCE/'primary_enriched_modules.tsv',sep='\t')
axes = ['vegetation_landscape_PC1','water_inundation_PC1','moisture_stress_PC1','physicochemical_PC1','nutrients_redox_PC1']
axis_labels = ['Vegetation /\nlandscape','Surface-water\nsignal','Moisture\n(NDMI)','Physico-\nchemistry','Nutrients /\nsulfide']
ko_names = {
 'K02456':'GspG secretion protein','K21572':'SusD starch-binding protein',
 'K02027':'Multiple-sugar binding protein','K03088':'ECF sigma factor',
 'K02026':'Multiple-sugar permease','K00854':'Xylulokinase','K01179':'Endoglucanase',
 'K01715':'Enoyl-CoA hydratase','K02035':'Peptide / nickel binding protein',
 'K00231':'Porphyrinogen oxidase','K03427':'Restriction-system methylase',
 'K05889':'Polyvinyl alcohol dehydrogenase','K08738':'Cytochrome c','K00799':'Glutathione S-transferase'}
mod_names = {'M00014':'Uronate pathway','M00001':'Glycolysis','M00155':'Cytochrome c oxidase',
 'M00529':'Denitrification','M00006':'Oxidative pentose phosphate pathway',
 'M00417':'Cytochrome o ubiquinol oxidase','M01069':'Nitrogen cycle',
 'M00095':'Mevalonate pathway','M00157':'F-type ATPase','M00033':'Ectoine biosynthesis'}

fig=plt.figure(figsize=(7.3,7.5),facecolor='white')
norm=LogNorm(1,32)
for panel,data,key,names,rect,title in [
 ('a',ko,'term',ko_names,[.45,.535,.52,.345],'Enriched KEGG orthologues'),
 ('b',modules,'module',mod_names,[.45,.195,.52,.247],'Potential metabolic modules')]:
 ax=fig.add_axes(rect)
 order=list(names)
 n=len(order)
 for i in range(n):
  if i%2==0: ax.axhspan(i-.5,i+.5,color='#f3f6f9',zorder=0)
 for j in [1.5,3.5,5.5,7.5]: ax.axvline(j,color='#cbd5df',lw=.65)
 for _,r in data.iterrows():
  x=2*axes.index(r.axis)+(r.direction=='positive')
  y=order.index(r[key]); count=int(r.selected_term_genes)
  ax.scatter(x,y,s=22+count*1.28,c=[r.fold_enrichment],cmap=blue,norm=norm,
             edgecolors='#172431' if r.retained_gene_q005 else '#7c94a7',
             linewidths=1.25 if r.retained_gene_q005 else .35,zorder=3)
 ax.set_xlim(-.55,9.55);ax.set_ylim(n-.5,-.5)
 ax.set_yticks(range(n),[f'{k}  {names[k]}' for k in order],fontsize=7.3)
 ax.tick_params(axis='y',length=0,pad=6)
 ax.set_xticks(range(10),['−','+']*5,fontsize=8.4)
 ax.xaxis.tick_top();ax.tick_params(axis='x',length=0,pad=3)
 for j,label in enumerate(axis_labels):
  ax.text(2*j+.5,1.105,label,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=7.5,linespacing=1.15)
 for sp in ax.spines.values():sp.set_visible(False)
 fig.text(.018,rect[1]+rect[3]+.085,panel,fontweight='bold',fontsize=11)
 fig.text(.055,rect[1]+rect[3]+.085,title,fontweight='bold',fontsize=9)

cax=fig.add_axes([.105,.087,.30,.014])
cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=blue),cax=cax,orientation='horizontal',ticks=[1,2,4,8,16,32])
cb.set_ticklabels(['1','2','4','8','16','32']);cb.ax.tick_params(labelsize=7,length=2,pad=2);cb.outline.set_visible(False)
fig.text(.105,.117,'Fold enrichment',fontsize=8)
handles=[plt.scatter([],[],s=22+n*1.28,facecolor='#bfd7e7',edgecolor='#7c94a7',linewidth=.4,label=str(n)) for n in [5,25,100]]
fig.legend(handles=handles,title='Associated genes',loc='lower left',bbox_to_anchor=(.53,.066),ncol=3,
 frameon=False,fontsize=7.5,title_fontsize=8,handletextpad=.5,columnspacing=1.1)
fig.text(.105,.033,'− / + : association direction     Dark outline: retained with gene q ≤ 0.05',fontsize=7.5)
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'Figure5.{ext}',dpi=500,facecolor='white')
fig.savefig(OUT/'Figure5_preview.png',dpi=150);plt.close(fig)

print('Created Figure 5 in',OUT)
