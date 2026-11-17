#!/usr/bin/env python3
"""Descriptive tables and figures. Does not change inference or select a new test.
Requires pandas, numpy, matplotlib; run after both analysis scripts.
"""
from pathlib import Path
from itertools import product
import argparse, io, json, gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

STAGES=['Degraded','Early restoration','Intermediate restoration','Advanced restoration','Conserved']
DEPTHS=[5,20,40]; LOCS=['Carmen','Cozumel','Tuxpan']
COLORS=['#0072B2','#D55E00','#009E73']
SHORT=['D','E','I','A','C']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
 'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})

def atomic(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp');tmp.write_bytes(data)
    assert tmp.read_bytes()==data;tmp.replace(path)

def table(d,path): atomic(path,d.to_csv(sep='\t',index=False,na_rep='NA',float_format='%.15g').encode())

def run(inputs,results):
    inputs=Path(inputs);results=Path(results);out=results/'figures';out.mkdir(exist_ok=True)
    tabs=results/'report_tables';tabs.mkdir(exist_ok=True)
    read=lambda f:pd.read_csv(results/f,sep='\t')
    m=read('611_main_results.tsv.gz');s=read('611_sensitivity_results.tsv.gz')
    l=read('611_LOLO_results.tsv.gz');sl=read('611_sensitivity_LOLO.tsv.gz')
    means=read('611_stage_means.tsv.gz')
    annotations=pd.read_csv(inputs/'feature_annotations.tsv',sep='\t')
    annotated=m.merge(annotations,on=['ontology','feature_id'],how='left')
    atomic(tabs/'611_all_results_annotated.tsv.gz',gzip.compress(annotated.to_csv(sep='\t',index=False,na_rep='NA',float_format='%.15g').encode(),mtime=0))
    wide=m[m.eligible].pivot(index=['ontology','feature_id'],columns='depth_cm',values='S').reindex(columns=DEPTHS)
    wide.columns=['S_5','S_20','S_40']
    changes=wide[(wide.gt(0).any(axis=1))&(wide.lt(0).any(axis=1))].reset_index()
    table(changes,tabs/'611_direction_changes_across_depths_descriptive.tsv')
    desc=m[m.eligible&m.monotone_observed&m.direction_stable].copy()
    desc=desc.merge(annotations,on=['ontology','feature_id'],how='left')
    table(desc,tabs/'611_descriptive_monotone_candidates.tsv')
    summary=read('611_summary.tsv')
    for name,mask in [('n_nominal',m.p_exact<=.05),('n_direction_stable',m.direction_stable),
                      ('n_descriptive_monotone',m.eligible&m.monotone_observed&m.direction_stable)]:
        counts=m.assign(v=mask).groupby(['ontology','depth_cm']).v.sum()
        summary[name]=[counts[(r.ontology,r.depth_cm)] for r in summary.itertuples()]
    table(summary,tabs/'611_summary_extended.tsv')
    ss=s.groupby(['ontology','variant']).agg(n_tests=('feature_id','size'),n_eligible=('eligible','sum'),
       min_p=('p_exact','min'),min_q=('q_BH','min'),n_BH=('q_BH',lambda x:int((x<=.05).sum())),
       n_BY=('q_BY',lambda x:int((x<=.05).sum()))).reset_index()
    table(ss,tabs/'611_sensitivity_summary.tsv')
    chosen=m[m.eligible & m.direction.ne(0)].sort_values(['p_exact','feature_id']).groupby(['ontology','depth_cm','direction']).head(2)
    chosen=chosen.merge(annotations,on=['ontology','feature_id'],how='left').sort_values(['ontology','depth_cm','direction','p_exact','feature_id'])
    table(chosen,tabs/'611_heatmap_examples_exploratory.tsv')
    top=m[m.eligible].sort_values(['p_exact','feature_id']).groupby(['ontology','depth_cm']).head(1).sort_values(['ontology','depth_cm'])
    table(top,tabs/'611_LOLO_examples_exploratory.tsv')
    all_members=[];sets=[];intersections=[]
    for criterion,frame in [('FDR_monotone_stable',m[m.stable_monotone_candidate]),('descriptive_monotone_stable',desc)]:
      for ont in ['KO','PFAM']:
       for direction in [1,-1]:
        z=frame[(frame.ontology==ont)&(frame.direction==direction)]
        ds={d:set(z.loc[z.depth_cm==d,'feature_id']) for d in DEPTHS}
        for dep in DEPTHS:sets.append([criterion,ont,direction,dep,len(ds[dep])])
        ids=set().union(*ds.values())
        for feature in sorted(ids):
            pattern=''.join('1' if feature in ds[d] else '0' for d in DEPTHS)
            all_members.append([criterion,ont,direction,feature,pattern])
        for pattern in ['100','010','001','110','101','011','111']:
            n=sum(''.join('1' if f in ds[d] else '0' for d in DEPTHS)==pattern for f in ids)
            intersections.append([criterion,ont,direction,pattern,n])
    setdf=pd.DataFrame(sets,columns=['criterion','ontology','direction','depth_cm','n_functions'])
    idf=pd.DataFrame(intersections,columns=['criterion','ontology','direction','pattern_5_20_40','n_functions'])
    table(setdf,tabs/'611_UpSet_sets.tsv');table(idf,tabs/'611_UpSet_intersections.tsv')
    table(pd.DataFrame(all_members,columns=['criterion','ontology','direction','feature_id','pattern_5_20_40']),tabs/'611_UpSet_members.tsv')
    meta=pd.read_csv(inputs/'metadata.tsv',sep='\t');x=pd.read_csv(inputs/'KO.tsv.gz',sep='\t',index_col=0).loc[meta.sample_id]
    p=x.div(x.sum(axis=1),axis=0);log=np.log2(p+5e-7);clr=log.sub(log.mean(axis=1),axis=0)
    obs=meta.copy();obs['relative_abundance']=p.K05825.to_numpy();obs['relative_percent']=100*obs.relative_abundance;obs['CLR']=clr.K05825.to_numpy()
    obs['stage_order']=obs.restoration4.map(dict(zip(STAGES,range(5))))
    table(obs,tabs/'611_K05825_observations.tsv')
    kstats=pd.concat([m[m.feature_id=='K05825'][['ontology','feature_id','depth_cm','variant','S','p_exact','q_BH','q_BY','monotone_observed']],s[s.feature_id=='K05825'][['ontology','feature_id','depth_cm','variant','S','p_exact','q_BH','q_BY','monotone_observed']]])
    table(kstats,tabs/'611_K05825_all_scales.tsv');table(sl[sl.feature_id=='K05825'],tabs/'611_K05825_sensitivity_LOLO.tsv')
    bookbuf=io.BytesIO();book=PdfPages(bookbuf)
    files=[]
    def save(fig,name):
        for ext in ['png','pdf','svg']:
            buf=io.BytesIO();fig.savefig(buf,format=ext,dpi=300,bbox_inches='tight');atomic(out/(name+'.'+ext),buf.getvalue())
        book.savefig(fig,bbox_inches='tight');plt.close(fig);files.append(name)
    for ont in ['KO','PFAM']:
        fig,axs=plt.subplots(3,1,figsize=(10.2,9.0));fig.subplots_adjust(left=.36,right=.91,top=.88,bottom=.08,hspace=.65)
        fig.suptitle(f'{ont} | Tendencias ilustrativas del análisis principal',x=.035,y=.985,ha='left',fontsize=16,fontweight='bold')
        fig.text(.035,.944,'Exploratorio: las 2 menores P por signo y profundidad. Ninguna supera FDR.',fontsize=11)
        fig.text(.035,.917,'Color: desviación de la media CLR de cada función; no es abundancia absoluta.',fontsize=10,color='#444444')
        for ax,dep in zip(axs,DEPTHS):
            z=chosen[(chosen.ontology==ont)&(chosen.depth_cm==dep)]
            mat=z[[f'mean_CLR_{j}' for j in range(5)]].to_numpy();mat-=mat.mean(axis=1,keepdims=True)
            lim=max(abs(mat).max(),.01)
            im=ax.imshow(mat,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto')
            ax.set_xticks(range(5),SHORT);ax.set_yticks(range(len(z)))
            labels=[]
            for r in z.itertuples():
                short=r.short_label if isinstance(r.short_label,str) else r.feature_id
                label=f'{r.feature_id}  {short}' if ont=='KO' else r.feature_id
                labels.append(label+(' *' if r.core_reference_9of9 else ''))
            ax.set_yticklabels(labels,fontsize=9)
            ax.set_title(f'{dep} cm | q mínima de estos ejemplos = {z.q_BH.min():.3f}',loc='left',fontsize=11,pad=8)
            for i,row in enumerate(mat):
                for j,v in enumerate(row):ax.text(j,i,f'{v:+.2f}',ha='center',va='center',fontsize=9,color='white' if abs(v)>.6*lim else '#111111')
            fig.colorbar(im,ax=ax,fraction=.025,pad=.025,label='Δ CLR (log₂)')
        fig.text(.035,.025,'D degradado · E temprano · I intermedio · A avanzado · C conservado · * core conservado 9/9\nPromedios por estadio con el mismo peso para cada localidad. Los perfiles son transversales.',fontsize=9)
        save(fig,f'611_01_heatmap_{ont}')
    fig,axs=plt.subplots(2,3,figsize=(11.4,7.4));fig.subplots_adjust(left=.08,right=.98,top=.84,bottom=.18,hspace=.6,wspace=.3)
    fig.suptitle('K05825 | La señal depende de la escala de abundancia',x=.04,y=.98,ha='left',fontsize=16,fontweight='bold')
    fig.text(.04,.929,'LYSN · 2-aminoadipate transaminase · potencial genético, no actividad medida',fontsize=11)
    for row,measure in enumerate(['relative_percent','CLR']):
      for col,dep in enumerate(DEPTHS):
        ax=axs[row,col];z=obs[obs.depth_cm==dep]
        for j,loc in enumerate(LOCS):
            zz=z[z.locality==loc];mn=zz.groupby('stage_order')[measure].mean()
            ax.plot(mn.index,mn.to_numpy(),color=COLORS[j],lw=1.5,label=loc)
            ax.scatter(zz.stage_order+(j-1)*.055,zz[measure],s=28,color=COLORS[j],edgecolor='white',lw=.5,zorder=3)
        variant='relative_abundance' if row==0 else 'main_CLR'
        st=kstats[(kstats.depth_cm==dep)&(kstats.variant==variant)].iloc[0]
        ax.set_title(f'{dep} cm · S = {st.S:.2f}\nP = {st.p_exact:.2g} · q = {st.q_BH:.3g}',fontsize=10,loc='left')
        ax.set_xticks(range(5),SHORT);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
        if col==0:ax.set_ylabel('Abundancia relativa (%)' if row==0 else 'CLR (log₂; δ = 0.5 ppm)')
    h,labs=axs[0,0].get_legend_handles_labels();fig.legend(h,labs,loc='lower center',bbox_to_anchor=(.5,.09),ncol=3,frameon=False)
    fig.text(.04,.025,'Arriba: sensibilidad con proporciones. Abajo: análisis principal. q de Benjamini–Hochberg.\nLíneas unen medias de estadios distintos; no representan seguimiento temporal de los mismos perfiles.',fontsize=9)
    save(fig,'611_02_K05825_scales')
    fig,axs=plt.subplots(2,3,figsize=(11.4,7.2),sharex=True);fig.subplots_adjust(left=.11,right=.98,top=.83,bottom=.16,wspace=.3,hspace=.5)
    fig.suptitle('LOLO | Influencia local y concordancia de dirección',x=.035,y=.98,ha='left',fontsize=16,fontweight='bold')
    fig.text(.035,.91,'Ejemplo con menor P por ontología y profundidad; todos exploratorios, q > 0.05.',fontsize=11)
    for ax,r in zip(axs.flat,top.itertuples()):
        ax.axvline(0,color='#999999',lw=.8);ax.axvline(r.S,color='#222222',lw=1.2,label='Completo')
        zz=l[(l.ontology==r.ontology)&(l.depth_cm==r.depth_cm)&(l.feature_id==r.feature_id)].set_index('omitted_locality')
        for j,loc in enumerate(LOCS):
            rr=zz.loc[loc];ax.plot([rr.S_training,rr.S_heldout],[j,j],color='#bbbbbb',zorder=1)
            ax.scatter(rr.S_training,j,marker='o',color=COLORS[j],s=45,zorder=3,label='Entrenamiento' if j==0 else None)
            ax.scatter(rr.S_heldout,j,marker='D',facecolor='white',edgecolor=COLORS[j],s=40,zorder=3,label='Localidad omitida' if j==0 else None)
        ax.set_yticks(range(3),LOCS,fontsize=9);ax.set_ylim(2.5,-.5);ax.set_xlim(-1.05,1.05)
        ax.set_title(f'{r.ontology} · {r.depth_cm} cm · {r.feature_id}\nq = {r.q_BH:.3f}',loc='left',fontsize=10)
        ax.set_xlabel('Concordancia ordinal S')
    h,labs=axs[0,0].get_legend_handles_labels();fig.legend(h,labs,loc='lower center',bbox_to_anchor=(.5,.065),ncol=3,frameon=False)
    fig.text(.035,.018,'Los tres pliegues se solapan. Concordancia de signo no equivale a selección FDR ni a validación externa.',fontsize=9)
    save(fig,'611_03_LOLO_examples')
    fig,axs=plt.subplots(1,2,figsize=(10.5,4.5));fig.subplots_adjust(top=.69,bottom=.17,wspace=.25)
    fig.suptitle('UpSet inferencial | No hay funciones seleccionadas',x=.04,y=.97,ha='left',fontsize=16,fontweight='bold')
    fig.text(.04,.86,'Criterio: q ≤ 0.05, medias monótonas y dirección estable en LOLO.\nLas intersecciones son vacías en ambas ontologías y ambos signos.',fontsize=11)
    for ax,ont in zip(axs,['KO','PFAM']):
        ax.barh(['5 cm','20 cm','40 cm'],[0,0,0],color='#0072B2');ax.set_xlim(0,1);ax.set_xticks([0]);ax.set_title(ont)
        for i in range(3):ax.text(.02,i,'0',va='center')
        ax.text(.5,1,'Sin intersecciones',ha='center',color='#555555');ax.set_xlabel('Funciones seleccionadas')
    save(fig,'611_04_UpSet_FDR_empty')
    fig=plt.figure(figsize=(11.4,8.6));outer=fig.add_gridspec(2,2,left=.07,right=.98,top=.85,bottom=.11,wspace=.25,hspace=.42)
    fig.suptitle('UpSet descriptivo | Coincidencias entre profundidades',x=.04,y=.98,ha='left',fontsize=16,fontweight='bold')
    fig.text(.04,.935,'Medias monótonas y dirección estable; sin filtro de P. Ninguna es un descubrimiento FDR.',fontsize=11)
    for spec,(ont,direction) in zip(outer,product(['KO','PFAM'],[1,-1])):
        gs=spec.subgridspec(2,2,height_ratios=[2.4,1.3],width_ratios=[1.15,3.5],hspace=.06,wspace=.08)
        bars=fig.add_subplot(gs[0,1]);mat=fig.add_subplot(gs[1,1]);sizes=fig.add_subplot(gs[1,0])
        ints=idf[(idf.criterion=='descriptive_monotone_stable')&(idf.ontology==ont)&(idf.direction==direction)]
        pats=ints.pattern_5_20_40.tolist();ns=ints.n_functions.to_numpy();bars.bar(range(7),ns,color='#0072B2' if direction==1 else '#D55E00')
        bars.set_title(f'{ont} · {"creciente" if direction==1 else "decreciente"}',loc='left',fontsize=11)
        bars.set_xticks([]);bars.set_ylabel('Intersección exclusiva',fontsize=9);bars.set_ylim(0,max(ns)*1.3+1)
        for j,n in enumerate(ns):bars.text(j,n+.3,str(n),ha='center',fontsize=9)
        for j,pat in enumerate(pats):
            active=[i for i,v in enumerate(pat) if v=='1'];mat.plot([j,j],[min(active),max(active)],color='#333333',lw=1.2)
            for i,v in enumerate(pat):mat.scatter(j,i,color='#333333' if v=='1' else '#dddddd',s=28)
        mat.set_xlim(-.5,6.5);mat.set_ylim(2.5,-.5);mat.set_xticks([]);mat.set_yticks([]);mat.set_frame_on(False)
        sz=setdf[(setdf.criterion=='descriptive_monotone_stable')&(setdf.ontology==ont)&(setdf.direction==direction)]
        ns2=[int(sz.loc[sz.depth_cm==d,'n_functions'].iloc[0]) for d in DEPTHS]
        sizes.barh(range(3),ns2,color='#777777');sizes.invert_xaxis();sizes.set_ylim(2.5,-.5)
        sizes.set_yticks(range(3),['5 cm','20 cm','40 cm'],fontsize=9);sizes.set_xlabel('Conjunto',fontsize=9)
        for i,n in enumerate(ns2):sizes.text(n+.7,i,str(n),ha='right',va='center',fontsize=8)
        sizes.set_xlim(max(ns2)*1.4,0)
    fig.text(.04,.035,'Barras superiores: funciones presentes exactamente en las profundidades marcadas con puntos oscuros.\nBarras laterales: tamaño de cada conjunto. Los signos se analizan por separado; no se agrupan cambios de dirección.',fontsize=9)
    save(fig,'611_05_UpSet_descriptive')
    book.close();atomic(out/'611_figures_all.pdf',bookbuf.getvalue())
    atomic(out/'figure_manifest.json',json.dumps({'figures':files,'inference_modified':False,'heatmap_rule':'two lowest exact P per sign, ontology and depth; exploratory only'},indent=2).encode())
    print(summary.to_string(index=False));print('FIGURES_COMPLETE',out)
    return out

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--inputs',required=True);a.add_argument('--results',required=True);o=a.parse_args();run(o.inputs,o.results)
