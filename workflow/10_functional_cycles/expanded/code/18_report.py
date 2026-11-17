from pathlib import Path
import pandas as pd,json,hashlib,shutil,base64
P=Path(__file__).resolve().parent;O=Path((P/'LATEST.txt').read_text());v=json.load(open(O/'validation.json'));s=pd.read_csv(O/'test_summary.tsv',sep='\t');c=pd.read_csv(O/'cycle_coverage.tsv',sep='\t');t=pd.read_csv(O/'tests.tsv',sep='\t');e=pd.read_csv(O/'effects.tsv',sep='\t');a=e[(e.variant=='primary')&(e.target=='restoration')].pivot(index='KO',columns='stage',values='ratio')
def md(df):
 return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in df.itertuples(index=False,name=None))
report='''# Panel ampliado: restauración y potencial biogeoquímico

La comparación principal no identifica diferencias entre estadios respaldadas tras la corrección conjunta. Esto no demuestra igualdad: solo contamos con 17 perfiles independientes, entre tres y cinco por estadio, y los tamaños estimados varían entre componentes.

Se conservó la cuantificación de Ciencia de Frontera. Las reglas AND/OR de MAGENTA se usaron para interpretar componentes requeridos y alternativas, sin trasladar abundancias ni umbrales proteicos de ese proyecto.

## Qué recuperamos

Se buscó en las 8,416 filas de la matriz original. Siete mapas y las reglas de cinco ciclos dieron 748 KOs únicos: 497 disponibles, 437 elegibles para comparación y 251 sin fila. Otros 60 tenían fila pero no una contribución primaria positiva en todas las muestras. Entre los elegibles, 180 permiten recuperar una contribución single-KO y 257 conservan cuantificación fraccional provisional. Un KO puede pertenecer a varios ciclos: las filas de la tabla siguiente no se suman.

'''+md(c)+'''

## Dónde están las diferencias estimadas

La Figura 6 compara cada estadio con el conservado después de ajustar por localidad. Un valor 0.5 significa la mitad de representación estimada; 2 significa el doble. Los colores no representan pruebas por pares. Los ejemplos se eligieron por identidad del proceso y no por P.

'''
for k,name in [('K02588','nifH, componente de fijación de nitrógeno'),('K15864','nirS, transformaciones de nitrito/hidroxilamina'),('K11180','dsrA, transformaciones del sulfito'),('K02040','pstS, componente de captación de fosfato')]:
 if k in a.index:report+=f"- {name}: degradado/conservado = {a.loc[k,'Degraded']:.2f}; avanzado/conservado = {a.loc[k,'Advanced restoration']:.2f}.\n"
report+='''
Ninguno de estos cocientes demuestra una diferencia por pares o una trayectoria temporal de recuperación. La variabilidad entre perfiles se muestra en el análisis de profundidad y las tablas conservan todas las estimaciones y omisiones.

## Respaldo estadístico y sensibilidad

La familia principal contiene 1,311 pruebas: 437 comparaciones globales de estadios y 874 pruebas de si los cambios pareados 20/5 o 40/5 difieren entre estadios. Todas las P Holm primarias son 1. La mínima P sin corregir para estadio es 0.00818. El análisis secundario no prueba un efecto global de profundidad; responde si el patrón vertical cambia entre estadios.

'''+md(s)+'''

La única señal con P Holm < 0.05 en las alternativas es K01659 (prpC, 2-metilcitrato sintasa), para diferencias entre estadios en el cambio 40/5 cuando se usa toda la contribución fraccional: P Holm = 0.01311. En la cuantificación primaria, su P Holm es 1. Se conserva como resultado sensible a la asignación, sin presentarlo como recuperación ni hallazgo robusto. Identidad de referencia: https://www.kegg.jp/entry/K01659.

## Lo que dicen las reglas de rutas

Bajo el criterio de contribución single-KO con límite inferior positivo, nifHDK, pstSCAB y sat–apr–dsr satisfacen sus reglas en las 51 muestras. El conjunto Sox central las satisface en 44, con siete indeterminadas. Los conjuntos MCR y anammox permanecen indeterminados en las 51: falta evidencia obligatoria recuperable. No se interpreta esto como ausencia biológica.

La detección de componentes de una comunidad no establece una vía en un mismo organismo. La dirección del metabolismo del azufre, la especificidad Nar/Nxr y Amo/Pmo, y la amplitud funcional de nirS requieren evidencia adicional. Las reglas no se sumaron para fabricar una abundancia de ruta.

## Cómo abrir y reproducir

- `Atlas_funcional_CF.html`: siete mapas, selector de estadio y tipo de comparación; cada KO ofrece su estimación y P global. Incluye reglas y frecuencias por muestra. Mapa01200 colorea líneas; los demás, casillas subdivididas por KO. Escala común ¼–4, con saturación.
- `Figure6.png` y `Figure6.pdf`: figura integrada; `Depth_within_stages.png`: puntos individuales de cambios pareados.
- `Supplementary_Table16_tests.tsv`, `Supplementary_Table16_estimates.tsv`: pruebas globales y coeficientes descriptivos, respectivamente.
- `rule_sample_evidence.tsv.gz`: evidencia por muestra; `operational_rules.tsv`: expresiones AND/OR originales.
- `leave_out_stage_effects.tsv`: sensibilidad de coeficientes al omitir perfiles/localidades; no son intervalos de confianza.
- `validation.json`: verificaciones numéricas y límites. No se realizó nueva calibración de tamaño finito del bootstrap para el panel ampliado. La inferencia es exploratoria y no demuestra equivalencia.

El paquete conserva scripts, mapas, tablas, matrices normalizadas y los cinco insumos RPK necesarios. Para recalcular, entrar al directorio `reproducibilidad` y ejecutar: `OPENBLAS_NUM_THREADS=1 python3 expanded_analysis/15_compare_expanded.py`, después `17_present_expanded.py`, `16_build_atlas.py` y `18_report.py`. Se requieren Python, numpy, pandas, scipy, matplotlib y Pillow. Las descargas KEGG se reutilizan cuando están presentes. Las reglas por muestra son la capa auditada previamente y se incluyen con sus expresiones y procedencia. El cálculo de sus límites desde conteos puede reconstruirse con el paquete previo de evidencia ampliada.

Los resultados previos de nitrógeno se conservaron para comparación y sus P sin corregir coinciden exactamente. La familia de multiplicidad se recalculó para el universo ampliado; no se copiaron las P ajustadas del piloto. La figura antigua de 11 marcadores ya no corresponde a la pregunta principal actual.

Topología de mapas: KEGG, © Kanehisa Laboratories. Fuentes y hashes incluidos. Ninguna flecha representa un flujo medido. El atlas se verificó estructuralmente y con extracción/sintaxis de JavaScript; este entorno no permitió abrir Chromium, por lo que no se afirma validación interactiva en navegador.
'''
(O/'Informe.md').write_text(report)
html='<meta charset="utf-8"><title>Panel biogeoquímico CF</title><style>body{font:16px system-ui;max-width:1100px;margin:30px auto;line-height:1.5}img{max-width:100%}pre{white-space:pre-wrap}table{border-collapse:collapse}td,th{padding:5px;border:1px solid #ccc}</style><h1>Potencial biogeoquímico y restauración</h1><p>No se identificaron diferencias entre estadios respaldadas en el análisis principal. Las magnitudes observadas y la incertidumbre de asignación se conservan.</p>'
for name in ['Figure6.png','Depth_within_stages.png','Nitrogen_stage_map.png']:
 html+='<img src="data:image/png;base64,'+base64.b64encode((O/name).read_bytes()).decode()+'">'
import html as h
html+='<h2>Informe completo</h2><pre>'+h.escape(report)+'</pre>'
(O/'Informe.html').write_text(html)
t.to_csv(O/'Supplementary_Table16_tests.tsv',sep='\t',index=False);e.to_csv(O/'Supplementary_Table16_estimates.tsv',sep='\t',index=False)
shutil.copy2(P/'manuscript_edits.json',O/'manuscript_edits.json')
(O/'hashes.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.is_file() and p.name!='hashes.json'},indent=2))
