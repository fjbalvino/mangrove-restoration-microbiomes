from pathlib import Path
import urllib.request,time,hashlib,json,base64,xml.etree.ElementTree as ET
import pandas as pd,numpy as np
from PIL import Image
import matplotlib;matplotlib.use('Agg')
from matplotlib import colormaps,colors
P=Path(__file__).resolve().parent;O=Path((P/'LATEST.txt').read_text());C=P/'kegg';C.mkdir(exist_ok=True)
maps={'01200':'Carbon metabolism','00710':'Calvin cycle reference map','00720':'Other carbon fixation pathways','00680':'Methane metabolism','00910':'Nitrogen metabolism','00920':'Sulfur metabolism','00440':'Phosphonate and phosphinate metabolism'}
for mid in maps:
 for ext,endpoint in [('xml','kgml'),('png','image')]:
  dest=C/f'ko{mid}.{ext}'
  if not dest.exists():
   old=P.parent/'nitrogen_pilot'/dest.name
   if old.exists():dest.write_bytes(old.read_bytes())
   else:
    with urllib.request.urlopen(f'https://rest.kegg.jp/get/ko{mid}/{endpoint}',timeout=90) as r:dest.write_bytes(r.read())
    time.sleep(.4)
 print('Fetched',mid,flush=True)
if not (O/'effects.tsv').exists():print('Maps downloaded; rerun after statistics.');raise SystemExit
E=pd.read_csv(O/'effects.tsv',sep='\t');E=E[E.variant.eq('primary')];T=pd.read_csv(O/'tests.tsv',sep='\t');T=T[T.variant.eq('primary')];R=pd.read_csv(O/'eligibility.tsv',sep='\t').set_index('KO');data={}
for k,row in R.iterrows():
 data[k]=dict(available=bool(row.available),eligible=bool(row.eligible),provisional=row.quantification!='single_KO',effects={},p={})
 for z in E[E.KO.eq(k)].itertuples():data[k]['effects'][z.target+'|'+z.stage]=z.log2_ratio
 for z in T[T.KO.eq(k)].itertuples():data[k]['p'][z.target]=z.p_holm
payload=[]
for mid,title in maps.items():
 root=ET.parse(C/f'ko{mid}.xml').getroot();entries=[]
 for e in root.findall('entry'):
  if e.get('type')!='ortholog':continue
  for g in e.findall('graphics'):
   if g.get('type') not in ['rectangle','line']:continue
   entries.append(dict(type=g.get('type'),coords=g.get('coords',''),kos=[k.replace('ko:','') for k in e.get('name').split()],label=g.get('name'),**{k:float(g.get(k) or 0) for k in ['x','y','width','height']}))
 w,h=Image.open(C/f'ko{mid}.png').size
 payload.append(dict(id=mid,title=title,entries=entries,w=w,h=h,img=base64.b64encode((C/f'ko{mid}.png').read_bytes()).decode()))
# Opaque split cells make color scale exact; native labels redrawn above subdivisions.
html='''<!doctype html><meta charset="utf-8"><title>Restauración · Atlas funcional</title><style>body{font:16px system-ui;margin:24px;color:#233444}select{padding:9px;margin:8px}#map{overflow:auto;border:1px solid #aaa}#detail{padding:15px;background:#eef3f6}p{max-width:1100px}table{border-collapse:collapse}td,th{padding:5px;border:1px solid #ddd}</style><h1>¿Dónde cambia la representación genética entre estadios?</h1><p>51 muestras · 17 perfiles independientes · ajuste por localidad · referencia: conservado. Las magnitudes son estimaciones exploratorias. Un color intenso no significa una diferencia respaldada estadísticamente.</p><label>Ruta<select id="route"></select></label><label>Estadio<select id="stage"><option value="Degraded">Degradado</option><option value="Early restoration">Temprano</option><option value="Intermediate restoration">Intermedio</option><option value="Advanced restoration">Avanzado</option></select></label><label>Comparación<select id="target"><option value="restoration">Promedio de las tres profundidades</option><option value="depth20_by_stage">Cambio 20/5 frente al cambio del conservado</option><option value="depth40_by_stage">Cambio 40/5 frente al cambio del conservado</option></select></label><p>Azul: menor; blanco: igual; rojo: mayor. Escala común: ¼, ½, 1, 2, 4 veces (saturación fuera de ese rango). Gris: KO sin fila; beige: no elegible. Borde discontinuo: cuantificación fraccional provisional. Cada franja representa un KO; clic para valores exactos.</p><div id="map"></div><div id="detail">Selecciona un componente.</div><p>Las reglas de complejos y rutas se entregan por separado. Codetección comunitaria no demuestra una vía en un organismo. Nar/Nxr, Amo/Pmo y nirS requieren interpretación bioquímica adicional; dsr no resuelve dirección del azufre. El mapa de fósforo cubre fosfonatos; transporte y polifosfato están además en las tablas curadas. ADN indica potencial, no tasas ni recuperación demostrada.</p><p>Topología de referencia © Kanehisa Laboratories, KEGG. Las flechas no representan flujos medidos.</p><script>const M=MAPS,D=DATA;const r=document.getElementById('route'),s=document.getElementById('stage'),t=document.getElementById('target');M.forEach((m,i)=>r.add(new Option(m.title,i)));function col(v){let a=v<0?[33,102,172]:[178,24,43],f=Math.min(1,Math.abs(v)/2);return 'rgb('+a.map(x=>Math.round(255+(x-255)*f)).join(',')+')'}function draw(){let m=M[+r.value],svg='<svg xmlns="http://www.w3.org/2000/svg" width="'+m.w+'" height="'+m.h+'"><image href="data:image/png;base64,'+m.img+'" width="'+m.w+'" height="'+m.h+'"/>';m.entries.forEach(e=>{e.kos.forEach((k,i)=>{let d=D[k],v=d?.effects[t.value+'|'+s.value],c=!d?.available?'#aaa':v==null?'#e9dfc6':col(v);if(e.type==='line'){let cxy=e.coords.split(',').map(Number),off=(i-(e.kos.length-1)/2)*2;svg+='<polyline data-ko="'+k+'" points="'+cxy.reduce((a,v,j)=>{if(j%2===0)a.push((v+off)+','+(cxy[j+1]+off));return a},[]).join(' ')+'" fill="none" stroke="'+c+'" stroke-width="2" '+(d?.provisional?'stroke-dasharray="3,2"':'')+'><title>'+k+'</title></polyline>';return;}svg+='<rect data-ko="'+k+'" x="'+(e.x-e.width/2+i*e.width/e.kos.length)+'" y="'+(e.y-e.height/2)+'" width="'+e.width/e.kos.length+'" height="'+e.height+'" fill="'+c+'" stroke="#222" stroke-width="0.6" '+(d?.provisional?'stroke-dasharray="2,2"':'')+'><title>'+k+'</title></rect>'});if(e.type==='rectangle')svg+='<text pointer-events="none" x="'+e.x+'" y="'+(e.y+3)+'" text-anchor="middle" font-family="Arial" font-size="8" fill="black" stroke="white" stroke-width="1.3" paint-order="stroke">'+e.label.replaceAll('&','&amp;').replaceAll('<','&lt;')+'</text>'});svg+='</svg>';document.getElementById('map').innerHTML=svg;document.querySelectorAll('[data-ko]').forEach(x=>x.onclick=()=>{let k=x.dataset.ko,d=D[k],v=d?.effects[t.value+'|'+s.value];document.getElementById('detail').textContent=k+' | '+(v==null?'Sin estimación':(2**v).toPrecision(4)+' veces; log₂ = '+v.toFixed(3))+' | P Holm global entre estadios: '+(d?.p[t.value]??'NA')+' | '+(d?.provisional?'Cuantificación provisional':'Contribución single-KO recuperable')+'. La prueba global no es una prueba de este par específico.'})}r.onchange=s.onchange=t.onchange=draw;draw();</script>'''
(O/'Atlas_funcional_CF.html').write_text(html.replace('MAPS',json.dumps(payload)).replace('DATA',json.dumps(data)))
(O/'kegg_provenance.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in C.iterdir()},indent=2))
print('Atlas ready',flush=True)
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
m=next(x for x in payload if x['id']=='00910');fig,ax=plt.subplots(figsize=(m['w']/110,m['h']/110));ax.imshow(Image.open(C/'ko00910.png'));cm=colors.LinearSegmentedColormap.from_list('map',['#2166ac','white','#b2182b']);norm=colors.Normalize(-2,2)
for e in m['entries']:
 for i,k in enumerate(e['kos']):
  d=data.get(k,{});v=d.get('effects',{}).get('restoration|Degraded');color='#aaa' if not d.get('available') else '#e9dfc6' if v is None else cm(norm(v));ax.add_patch(Rectangle((e['x']-e['width']/2+i*e['width']/len(e['kos']),e['y']-e['height']/2),e['width']/len(e['kos']),e['height'],facecolor=color,edgecolor='black',lw=.3,linestyle='--' if d.get('provisional') else '-'))
 ax.text(e['x'],e['y'],e['label'],ha='center',va='center',fontsize=4,color='black')
ax.set_title('Nitrogen: degraded / conserved, adjusted for locality\nCommon colour scale ¼–4×; no supported primary stage difference',fontsize=12);ax.axis('off');fig.tight_layout();fig.savefig(O/'Nitrogen_stage_map.png',dpi=180);plt.close(fig)
rr=pd.read_csv(O/'rule_summary.tsv',sep='\t');rr=rr[(rr.variant=='complete')&(rr.evidence=='single_KO_bounds')][['cycle','rule_id','label','met','not_met','unknown','expression','limit']]
f=O/'Atlas_funcional_CF.html';ht=f.read_text();ht=ht.replace('<script>','<h2>Reglas por muestra: componentes requeridos y alternativas</h2><p>Frecuencias sobre 51 muestras, no 51 réplicas independientes. “met” indica evidencia comunitaria suficiente para la regla; “not_met”, regla no satisfecha por la evidencia; “unknown”, indeterminada. No son abundancias ni pruebas entre estadios.</p><div style="overflow:auto;max-height:500px">'+rr.to_html(index=False,escape=True)+'</div><script>',1);f.write_text(ht)
