from pathlib import Path
import datetime,json,hashlib,shutil,base64,xml.etree.ElementTree as ET
import numpy as np,pandas as pd
from scipy.stats import binomtest
from stats_cns import depth_test,adjust,profile_matrix
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent;S=P;O=P/('19_nitrogen_depth_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S'));O.mkdir();(P/'LATEST.txt').write_text(str(O))
reg=pd.read_csv(S/'eligibility.tsv',sep='\t').set_index('KO');scope=pd.read_csv(S/'scope.tsv',sep='\t');keys=sorted(set(scope.loc[scope.cycle.eq('N'),'KO'])&set(reg.index[reg.eligible]));assert len(keys)==49
m=pd.read_csv(S/'samples.tsv',sep='\t');p=pd.read_csv(S/'profiles.tsv',sep='\t');assert len(p)==17
plan=dict(primary='paired depth 20/5 and 40/5; exact two-sided sign tests across 17 independent profiles',family='98 nitrogen depth tests, fixed eligible set of 49 KOs from the expanded panel',secondary='within-stage descriptive estimates only (3-5 profiles per stage); existing stage-by-depth tests retained as contextual evidence',interval='pointwise 95.0958% distribution-free common-median intervals; not simultaneous; heterogeneous profile distributions limit common-median interpretation',ties='absolute log2 difference <=1e-12 treated as tie, excluded from sign test; intervals use all 17',exploratory=True,global_context='also combine 98 new depth tests with the existing 1311 tests: Holm 1409; not a study-wide correction',variants=['primary','fractional','reference56','both_mates','supported_fractional'])
(O/'plan.json').write_text(json.dumps(plan,indent=2));tests=[];pairs=[];within=[];sensitivity=[]
for vn in plan['variants']:
 v=pd.read_csv(S/(vn+'_normalized.tsv'),sep='\t',index_col=0);v.loc[keys].to_csv(O/(vn+'_normalized.tsv'),sep='\t');rr=[]
 for k in keys:
  valid=(v.loc[k]>1e-7).all()
  arr=np.log2(profile_matrix(v.loc[[k]],m,p.profile_id)[0]) if valid else None
  for j,depth in [(1,20),(2,40)]:
   if not valid:rr.append(dict(variant=vn,KO=k,depth=depth,p=1,status='not_estimable_zero'));continue
   delta=arr[:,j]-arr[:,0];z=depth_test(delta);rr.append(dict(variant=vn,KO=k,depth=depth,status='estimated',ratio=2**z['median_log2'],ratio_low=2**z['ci_low'],ratio_high=2**z['ci_high'],**z))
   if vn=='primary':
    for i,pr in p.iterrows():pairs.append(dict(KO=k,depth=depth,profile_id=pr.profile_id,stage=pr.restoration4,locality=pr.locality,log2_change=delta[i],ratio=2**delta[i]))
    for stage,ix in p.groupby('restoration4').groups.items():
     d=delta[list(ix)];within.append(dict(KO=k,depth=depth,stage=stage,n_profiles=len(d),median_ratio=2**np.median(d),n_positive=int((d>0).sum()),n_negative=int((d<0).sum())))
    for typ,col in [('profile','profile_id'),('locality','locality')]:
     for val in p[col].unique():sensitivity.append(dict(KO=k,depth=depth,omission=typ+':'+str(val),median_log2=float(np.median(delta[p[col].ne(val)]))))
 for z,q in zip(rr,adjust([z['p'] for z in rr])):z['p_holm_98']=q
 old=pd.read_csv(S/'previous_tests.tsv',sep='\t');old=old[old.variant.eq(vn)];assert len(old)==1311
 q=adjust(old.p.tolist()+[z['p'] for z in rr])[-98:]
 for z,a in zip(rr,q):z['p_holm_1409']=a
 tests.extend(rr)
t=pd.DataFrame(tests);t.to_csv(O/'depth_tests.tsv',sep='\t',index=False);pd.DataFrame(pairs).to_csv(O/'paired_profiles.tsv',sep='\t',index=False);pd.DataFrame(within).to_csv(O/'within_stage.tsv',sep='\t',index=False);pd.DataFrame(sensitivity).to_csv(O/'leave_out.tsv',sep='\t',index=False)
primary=t[t.variant.eq('primary')]; summary=t.groupby(['variant','depth']).agg(tests=('p','size'),supported_98=('p_holm_98',lambda x:int((x<.05).sum())),supported_1409=('p_holm_1409',lambda x:int((x<.05).sum())));summary.to_csv(O/'summary.tsv',sep='\t')
rob=t.pivot(index=['KO','depth'],columns='variant',values='p_holm_98');rob['all_five_supported']=(rob[plan['variants']]<.05).all(axis=1);rob.to_csv(O/'robustness.tsv',sep='\t')
# All eligible genes, no selection by significance.
fig,ax=plt.subplots(figsize=(9,15));ys=np.arange(len(keys))
for d,offset,c in [(20,-.15,'#1984a3'),(40,.15,'#b54735')]:
 a=primary[primary.depth.eq(d)].set_index('KO').loc[keys];x=np.log2(a.ratio);lo=np.log2(a.ratio_low);hi=np.log2(a.ratio_high);ax.errorbar(x,ys+offset,xerr=np.array([x-lo,hi-x]),fmt='none',ecolor=c,alpha=.65,lw=1)
 ax.scatter(x,ys+offset,s=24,facecolors=[c if q<.05 else 'white' for q in a.p_holm_98],edgecolors=c,label=str(d)+'/5 cm',zorder=3)
ax.axvline(0,color='gray',lw=1);ax.set_yticks(ys,[k+(' †' if reg.loc[k,'quantification']!='single_KO' else '') for k in keys],fontsize=9);ax.invert_yaxis();ax.set_xlabel('log₂ de la razón pareada; 0 = sin cambio');ax.set_title('Nitrógeno: cambios con la profundidad\n17 perfiles independientes · 49 KOs · Holm sobre 98 pruebas');ax.legend();fig.text(.02,.01,'Intervalos puntuales 95.1%; puntos llenos: P Holm < 0.05. † Cuantificación provisional.',fontsize=9);fig.tight_layout(rect=[0,.025,1,1]);fig.savefig(O/'Nitrogeno_profundidad_completo.png',dpi=160);fig.savefig(O/'Nitrogeno_profundidad_completo.pdf');plt.close(fig)
selected={'K02588':'nifH · fijación de N','K15864':'nirS · nitrito/hidroxilamina','K00368':'nirK · reducción de nitrito','K04561':'norB · reducción de NO','K00376':'nosZ · reducción de N₂O','K03385':'nrfA · nitrito a amonio'}
a=pd.DataFrame(pairs);st=['Degraded','Early restoration','Intermediate restoration','Advanced restoration','Conserved'];fig,axes=plt.subplots(2,3,figsize=(12,7))
for ax,(k,label) in zip(axes.flat,selected.items()):
 for j,s in enumerate(st):
  v=a[(a.KO==k)&(a.depth==40)&(a.stage==s)].sort_values('profile_id');ax.scatter(j+np.linspace(-.12,.12,len(v)),v.ratio,s=25,color='#1984a3');ax.plot([j-.2,j+.2],[2**np.median(v.log2_change)]*2,color='black',lw=2)
 ax.set_yscale('log',base=2);ax.set_ylim(.025,2.5);ax.set_yticks([.03125,.0625,.125,.25,.5,1,2],['1/32','1/16','1/8','1/4','1/2','1','2']);ax.axhline(1,color='gray',lw=.8);ax.set_xticks(range(5),['D','E','I','A','C']);ax.set_title(label,fontsize=10);ax.set_ylabel('Razón 40/5 cm')
fig.suptitle('Profundidad dentro de cada estadio · cada punto es un perfil\nD degradado · E temprano · I intermedio · A avanzado · C conservado');fig.tight_layout();fig.savefig(O/'Nitrogeno_por_estadio.png',dpi=180);plt.close(fig)
selected_table=primary[primary.KO.isin(selected)].copy();selected_table['function']=selected_table.KO.map(selected);selected_table.to_csv(O/'selected_markers.tsv',sep='\t',index=False)
# Portable map of pooled paired changes; no stage-specific significance implied.
root=ET.parse(P/'ko00910.xml').getroot();entries=[]
for e in root.findall('entry'):
 if e.get('type')=='ortholog':
  g=e.find('graphics');entries.append(dict(kos=[k.replace('ko:','') for k in e.get('name').split()],**{x:float(g.get(x)) for x in ['x','y','width','height']}))
from PIL import Image
img=P/'ko00910.png';w,h=Image.open(img).size
payload={str(d):{z.KO:dict(ratio=z.ratio,p=z.p_holm_98,context_p=z.p_holm_1409,lo=z.ratio_low,hi=z.ratio_high) for z in primary[primary.depth.eq(d)].itertuples()} for d in [20,40]}
html='''<meta charset="utf-8"><title>Nitrógeno y profundidad</title><style>body{font:16px system-ui;margin:25px}#map{overflow:auto}img{max-width:100%}</style><h1>Nitrógeno: cambios pareados con profundidad</h1><p>17 perfiles · 51 muestras. Azul: menor; rojo: mayor representación que a 5 cm. Escala común ¼–4. Gris: sin estimación elegible. Cada franja es un KO, no una decisión sobre un complejo. Haz clic para consultar el cambio y las pruebas.</p><select id="d"><option>20</option><option selected>40</option></select> cm frente a 5 cm<div id="map"></div><p id="detail"></p><p>Los intervalos son puntuales, no simultáneos. Holm98 corresponde únicamente al piloto de profundidad; Holm1409 añade las pruebas previas de estadios e interacción. Es exploratorio. No identifica tasas ni demuestra recuperación. © KEGG/Kanehisa Laboratories.</p><script>const E=ENTRIES,D=DATA,IMG="IMAGE";const sel=document.getElementById('d');function color(v){let f=Math.min(1,Math.abs(Math.log2(v))/2),c=v<1?[33,102,172]:[178,24,43];return 'rgb('+c.map(x=>Math.round(255+(x-255)*f))+')'}function draw(){let s='<svg width="WIDTH" height="HEIGHT"><image width="WIDTH" height="HEIGHT" href="data:image/png;base64,'+IMG+'"/>';E.forEach(e=>e.kos.forEach((k,i)=>{let a=D[sel.value][k];s+='<rect data-k="'+k+'" x="'+(e.x-e.width/2+i*e.width/e.kos.length)+'" y="'+(e.y-e.height/2)+'" width="'+e.width/e.kos.length+'" height="'+e.height+'" fill="'+(a?color(a.ratio):'#aaa')+'" fill-opacity=".65" stroke="black" stroke-width=".3"/>'}));document.getElementById('map').innerHTML=s+'</svg>';document.querySelectorAll('[data-k]').forEach(x=>x.onclick=()=>{let k=x.dataset.k,a=D[sel.value][k];document.getElementById('detail').textContent=k+': '+(a?JSON.stringify(a):'Sin estimación elegible')})}sel.onchange=draw;draw();</script>'''
html=html.replace('ENTRIES',json.dumps(entries)).replace('DATA',json.dumps(payload)).replace('IMAGE',base64.b64encode(img.read_bytes()).decode()).replace('WIDTH',str(w)).replace('HEIGHT',str(h));(O/'Mapa_nitrogeno_profundidad.html').write_text(html)
for f in ['samples.tsv','profiles.tsv','eligibility.tsv'] :shutil.copy2(S/f,O/f)
shutil.copy2(__file__,O/Path(__file__).name);shutil.copy2(P/'stats_cns.py',O/'stats_cns.py')
# Independent sign probabilities and family-size checks.
for z in primary.itertuples():assert np.isclose(z.p,binomtest(z.n_positive,z.n_nonzero,.5).pvalue)
assert len(t)==490 and (t.groupby('variant').size()==98).all()
(O/'source_hashes.json').write_text(json.dumps({str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [S/(v+'_normalized.tsv') for v in plan['variants']]+[S/'scope.tsv',S/'previous_tests.tsv']},indent=2))
print(summary.to_string());print(selected_table[['KO','depth','ratio','ratio_low','ratio_high','n_positive','p_holm_98','p_holm_1409']].to_string(index=False));print('ALL_FIVE',rob[rob.all_five_supported].index.tolist());print('OUT',O)
