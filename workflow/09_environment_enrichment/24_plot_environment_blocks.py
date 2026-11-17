"""Figure 5: preserve enrichment and expose descriptive component evidence."""
from pathlib import Path
import argparse,json,html,datetime
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,LogNorm
ap=argparse.ArgumentParser();ap.add_argument('--ora',type=Path,required=True);ap.add_argument('--blocks',type=Path,required=True);ap.add_argument('--out-root',type=Path,required=True);a=ap.parse_args();a.output=a.out_root/('24_plot_environment_blocks_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S'));a.output.mkdir(parents=True,exist_ok=False);(a.out_root/'LATEST_24.txt').write_text(str(a.output.resolve()))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
axes=['vegetation_landscape_PC1','water_inundation_PC1','moisture_stress_PC1','physicochemical_PC1','nutrients_redox_PC1'];axnames=['Vegetation /\nlandscape','Surface-water\nsignal','Moisture\n(NDMI)','Physico-\nchemistry','Nutrients /\nsulfide'];short=dict(zip(axes,['Vegetation','Water signal','NDMI','Physicochemistry','Nutrients–sulfide']))
kn={'K02456':'GspG secretion protein','K21572':'SusD starch-binding protein','K02027':'Multiple-sugar binding protein','K03088':'ECF sigma factor','K02026':'Multiple-sugar permease','K00854':'Xylulokinase','K01179':'Endoglucanase','K01715':'Enoyl-CoA hydratase','K02035':'Peptide / nickel binding protein','K00231':'Porphyrinogen oxidase','K03427':'Restriction-system methylase','K05889':'Polyvinyl alcohol dehydrogenase','K08738':'Cytochrome c','K00799':'Glutathione S-transferase'}
mn={'M00014':'Uronate pathway','M00001':'Glycolysis','M00155':'Cytochrome c oxidase','M00529':'Denitrification components','M00006':'Oxidative pentose phosphate','M00417':'Cytochrome o oxidase','M01069':'Nitrogen-cycle components','M00095':'Mevalonate pathway','M00157':'F-type ATPase','M00033':'Ectoine biosynthesis'}
t=pd.read_csv(a.ora/'primary_enriched_terms_named.tsv',sep='\t');t=t[t.ontology.eq('KO')];s=pd.read_csv(a.blocks/'block_contrast_summary.tsv',sep='\t');blue=LinearSegmentedColormap.from_list('blue',['#edf4fa','#9ecae1','#2171b5']);norm=LogNorm(1,32)
fig=plt.figure(figsize=(8.2,8.7));ax=fig.add_axes([.43,.59,.54,.30]);order=list(kn)
for i in range(14):
 if i%2==0:ax.axhspan(i-.5,i+.5,color='#f4f6f8')
for _,r in t.iterrows():ax.scatter(2*axes.index(r.axis)+int(r.direction=='positive'),order.index(r.term),s=18+1.1*r.selected_term_genes,c=[r.fold_enrichment],norm=norm,cmap=blue,edgecolors='#172431' if r.retained_gene_q005 else '#8d9aa5',linewidths=1.2 if r.retained_gene_q005 else .4,zorder=3)
ax.set_xlim(-.5,9.5);ax.set_ylim(13.5,-.5);ax.set_yticks(range(14),[k+'  '+kn[k] for k in order],fontsize=8.5);ax.tick_params(axis='y',length=0);ax.set_xticks(range(10),['−','+']*5);ax.xaxis.tick_top();ax.tick_params(axis='x',length=0)
for i,l in enumerate(axnames):ax.text(2*i+.5,1.08,l,transform=ax.get_xaxis_transform(),ha='center',fontsize=8.5,va='bottom')
for x in [1.5,3.5,5.5,7.5]:ax.axvline(x,color='#cbd5df',lw=.6)
for sp in ax.spines.values():sp.set_visible(False)
fig.text(.02,.975,'a  Enriched KEGG annotations',fontsize=11,weight='bold')
ca=fig.add_axes([.43,.545,.23,.012]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=blue),cax=ca,orientation='horizontal',ticks=[1,2,4,8,16,32]);cb.set_ticklabels(['1','2','4','8','16','32']);cb.ax.tick_params(labelsize=7);fig.text(.43,.57,'Fold enrichment',fontsize=8)
handles=[plt.scatter([],[],s=18+1.1*n,facecolor='#9ecae1',edgecolor='#8d9aa5',label=str(n)) for n in [5,25,100]];fig.legend(handles=handles,title='Associated genes',loc='lower left',bbox_to_anchor=(.70,.53),ncol=3,frameon=False,fontsize=7,title_fontsize=7)
fig.text(.02,.487,'b  Candidate modules: statistical and component evidence',fontsize=11,weight='bold')
bx=fig.add_axes([.43,.09,.54,.335]);bx.set_xlim(-.5,4.5);bx.set_ylim(10.5,-.5)
labels=[]
for i,r in s.iterrows():
 if i%2==0:bx.axhspan(i-.5,i+.5,color='#f4f6f8',zorder=0)
 labels.append(mn[r.module]+' · '+short[r.axis]+(' −' if r.direction=='negative' else ' +'))
 for x,v,n in [(2,r.all_assignments_associated_steps,r.steps),(3,r.single_KO_associated_steps,r.steps),(4,r.single_KO_bounds_met_samples,51)]:
  bx.add_patch(plt.Rectangle((x-.49,i-.48),.98,.96,facecolor=blue(.05+.75*v/n),edgecolor='white'))
  text=f'{int(v)}/{int(n)}'
  if x==4 and r.single_KO_bounds_unknown_samples:text+='\n'+str(int(r.single_KO_bounds_unknown_samples))+' unresolved'
  bx.text(x,i,text,ha='center',va='center',fontsize=8.5,color='white' if v/n>.85 else '#172431')
 bx.text(0,i,f'{r.q_ora:.4f}',ha='center',va='center',fontsize=8.6)
 bx.text(1,i,f'{r.q_single:.3f}' if pd.notna(r.q_single) else 'Ineligible',ha='center',va='center',fontsize=8.6,color='#555')
bx.set_yticks(range(11),labels,fontsize=8.8);bx.tick_params(length=0);bx.set_xticks(range(5),['ORA q\nprimary','ORA q\nsingle-KO','Associated\nblocks: all','Associated\nblocks: single','Samples meeting\nrule: single'],fontsize=8.4);bx.xaxis.tick_top()
for sp in bx.spines.values():sp.set_visible(False)
fig.text(.02,.046,'Panel b shading describes component coverage, not statistical significance.',fontsize=8)
fig.text(.02,.026,'Single = single-KO evidence. No module retains enrichment q ≤ 0.10 in the single-KO analysis.',fontsize=8)
for ext in ['png','pdf','svg']:fig.savefig(a.output/f'Figure5.{ext}',dpi=350)
plt.close(fig)
# Per-step descriptive atlas. No additional P values or artificial pathway abundances.
components=pd.read_csv(a.blocks/'component_associations.tsv',sep='\t');steps=pd.read_csv(a.blocks/'step_sample_evidence.tsv.gz',sep='\t');rules=pd.read_csv(a.blocks/'module_rules.tsv',sep='\t');assoc=pd.read_csv(a.blocks/'step_associations.tsv',sep='\t');parts=['<html><head><meta charset="utf-8"><title>Environmental component evidence</title><style>body{font:16px system-ui;max-width:1250px;margin:40px auto;color:#173047}table{border-collapse:collapse;font-size:13px;display:block;overflow:auto}td,th{padding:7px;border:1px solid #ccd5df}th{background:#e8f1f7}h2{margin-top:48px}code{white-space:pre-wrap}</style></head><body><h1>Environmental associations and functional components</h1><p>51 samples, 17 profiles. Existing gene estimates, enrichment and a descriptive component overlay. No new block hypothesis test. Lack of support is not equivalence; the small design does not prove that every unsupported pattern is a sampling limitation.</p><p>Rule representation is evaluated within each sample; it does not establish a pathway in one organism. Associated components need not form a coordinated response. Gene estimates include nonsignificant and unusable fits in the accompanying tables.</p>',s.to_html(index=False)]
for r in rules.itertuples():
 parts+=['<h2>'+html.escape(r.module+' — '+r.name)+'</h2><p><code>'+html.escape(r.definition)+'</code></p>']
 z=steps[(steps.module==r.module)&(steps.variant=='complete')];tab=z.groupby(['step','expression','evidence','state']).size().unstack(fill_value=0).reset_index();parts.append(tab.to_html(index=False))
 parts.append('<h3>Gene-association coverage</h3><p>Zero tested genes means not evaluated in the selected set, not absence. Median/range of gene coefficients are descriptive summaries, not a fitted KO or pathway effect.</p>');parts.append(components[components.module==r.module].to_html(index=False,float_format=lambda x:f'{x:.3g}'))
parts.append('</body></html>');(a.output/'Environmental_blocks_atlas.html').write_text('\n'.join(parts))
print(a.output)
