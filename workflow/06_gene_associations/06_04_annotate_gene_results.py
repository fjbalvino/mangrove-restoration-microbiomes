#!/usr/bin/env python3
# ============================================================
# 06_04_annotate_gene_results.py
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Annotates the full tested universe and selected hits; does not refit models or recompute FDR.
# Inputs (source expressions; complete list in docs/contracts/06_04_annotate_gene_results.json):
#   lines=p.read_text().splitlines()
#   manifest=json.loads((audit/'manifest.json').read_text())
# Outputs (source expressions; complete list in contract):
#   temp.write_text(str(run)+'\n')
#   (logs/'604c_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
# Algorithmic provenance:
# Join existing functional annotations to the tested gene universe and results; no model refit.
#   Bates et al. (2015), doi:10.18637/jss.v067.i01; Kuznetsova et al. (2017), doi:10.18637/jss.v082.i13; Benjamini & Hochberg (1995), doi:10.1111/j.2517-6161.1995.tb02031.x.
# Source SHA-256: e8827df8e53c0be02367ab3058bd0ec9a9b812b16835d3fa8b3bdf3bf1323bd1
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================

"""Annotate the audited 603c universe from the complete annotation table.

Python standard library only. No fitting, FDR recalculation or enrichment.
Counts/CLR of five previously specified candidates are extracted for review.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import shutil
import statistics
import tarfile
import time
from collections import Counter, defaultdict
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

NAME = '604c_anotar_universo_proteico_603c'
ROOT = '/home/fjbalvino/Tipping_points/resultados_finales'
ANNOTATION = '/data/Ciencia-Frontera/Results/04-assemblies/assemblies-annotations/final_tables/FINAL_ANNOTATED_ABUNDANCE_noMZ_bacteria.tsv.gz'
AXES = ('vegetation_landscape_PC1', 'water_inundation_PC1', 'moisture_stress_PC1',
        'physicochemical_PC1', 'nutrients_redox_PC1')
FOCAL = {
    'CZ5_2.contigs_935909': 'environment_depth_interaction_moisture',
    'CZ1_3.contigs_1220634': 'environment_depth_interaction_nutrients',
    'CC1_3.contigs_1003408': 'MHI_diagnostic_singular_not_primary',
    'CC2_3.contigs_980998': 'MHI_diagnostic_singular_not_primary',
    'CZ3_2.contigs_1844408': 'MHI_near_cutoff_not_discovery',
}
ANNOT_FIELDS = {'query','seed_ortholog','evalue','score','eggnog_ogs','max_annot_lvl',
    'cog_category','description','preferred_name','gos','ec','kegg_ko','kegg_pathway',
    'kegg_module','kegg_reaction','kegg_rclass','brite','kegg_tc','cazy','bigg_reaction','pfams'}
MISSING = {'','-','na','n/a','nan','none','null'}
csv.field_size_limit(10_000_000)


def need(condition, message):
    if not condition:
        raise ValueError(message)


def log(message):
    print(datetime.now().strftime('%H:%M:%S') + ' | ' + message, flush=True)


def opener(path, mode='rt'):
    if str(path).endswith('.gz'):
        return gzip.open(path, mode, newline='', encoding='utf-8-sig' if 'r' in mode else 'utf-8')
    return path.open(mode.replace('t',''), newline='', encoding='utf-8-sig' if 'r' in mode else 'utf-8')


def read_table(path, delimiter='\t'):
    with opener(path) as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        need(reader.fieldnames and len(set(reader.fieldnames))==len(reader.fieldnames), f'Cabecera invalida: {path}')
        rows = list(reader)
    need(all(None not in row and None not in row.values() for row in rows), f'Fila incompleta: {path}')
    return rows


def write_table(path, rows, fields):
    with opener(path, 'wt') as f:
        writer=csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def digest(path, algorithm='sha256'):
    h=hashlib.new(algorithm)
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def latest(root, name):
    p=root/('LATEST_'+name+'.txt')
    need(p.is_file(), f'Falta puntero: {p}')
    lines=p.read_text().splitlines()
    need(lines and lines[0].strip(), f'Puntero vacio: {p}')
    return Path(lines[0].strip()).resolve(strict=True)


def publish(root, run, status):
    suffix='' if status=='PASS' else '_review_required'
    dest=root/('LATEST_'+NAME+suffix+'.txt')
    temp=dest.with_name(dest.name+f'.{os.getpid()}.tmp')
    temp.write_text(str(run)+'\n')
    temp.replace(dest)
    return dest


class HashingReader(io.RawIOBase):
    def __init__(self, raw):
        super().__init__()
        self.raw=raw
        self.sha256=hashlib.sha256()
        self.bytes_read=0

    def readable(self):
        return True

    def readinto(self, buffer):
        n=self.raw.readinto(buffer)
        if n:
            self.sha256.update(memoryview(buffer)[:n])
            self.bytes_read+=n
        return n


@contextmanager
def hashed_text(path):
    with path.open('rb') as raw:
        tracker=HashingReader(raw)
        with io.BufferedReader(tracker, buffer_size=1024*1024) as buf:
            if str(path).endswith('.gz'):
                with gzip.GzipFile(fileobj=buf, mode='rb') as gz:
                    with io.TextIOWrapper(gz, encoding='utf-8-sig', newline='') as text:
                        yield text,tracker
            else:
                with io.TextIOWrapper(buf, encoding='utf-8-sig', newline='') as text:
                    yield text,tracker


def tokens(value, ontology):
    good, invalid=set(), set()
    for part in re.split(r'[,;]', value):
        term=part.strip()
        if term.lower() in MISSING:
            continue
        if ontology=='KO':
            match=re.fullmatch(r'(?:ko:)?(K\d{5})', term, flags=re.IGNORECASE)
            if match:
                good.add(match.group(1).upper())
            else:
                invalid.add(term)
        elif re.fullmatch(r'[A-Za-z0-9_.:+\-]+', term):
            # Preserve Pfam names/accessions exactly; do not invent accession mappings.
            good.add(term)
        else:
            invalid.add(term)
    return sorted(good), sorted(invalid)


def scan_annotations(path, targets):
    before=path.stat()
    hits={}
    duplicates=Counter()
    last=time.monotonic()
    n=0
    with hashed_text(path) as (handle, tracker):
        reader=csv.reader(handle, delimiter='\t')
        fields=None
        for row in reader:
            if not row or row[0].startswith('##'):
                continue
            fields=[v.strip().lstrip('#').strip().lower() for v in row]
            break
        need(fields and len(fields)==len(set(fields)), 'Cabecera duplicada/invalida en anotacion')
        need({'query','kegg_ko','pfams'}.issubset(fields), 'Se requieren query/#query, KEGG_ko y PFAMs')
        qi=fields.index('query')
        cols=[(i,name) for i,name in enumerate(fields) if name in ANNOT_FIELDS and name!='query']
        for row in reader:
            if not row:
                continue
            n+=1
            need(len(row)==len(fields), f'Anotacion: fila {n} tiene {len(row)} columnas, esperadas {len(fields)}')
            fid=row[qi]
            if fid in targets:
                if fid in hits:
                    duplicates[fid]+=1
                else:
                    hits[fid]={name:row[i] for i,name in cols}
            if time.monotonic()-last>=30:
                pct=100*tracker.bytes_read/before.st_size if before.st_size else 0
                log(f'Anotacion: {n:,} filas; {len(hits):,}/{len(targets):,} IDs; {pct:.1f}% de bytes leidos')
                last=time.monotonic()
        source_sha=tracker.sha256.hexdigest()
        source_bytes=tracker.bytes_read
    after=path.stat()
    need((before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns), 'La anotacion cambio durante la lectura')
    need(source_bytes==before.st_size, 'No se recorrio el archivo completo de anotacion')
    return hits,duplicates,[name for _,name in cols],dict(path=str(path),size=before.st_size,
        mtime_ns=before.st_mtime_ns,sha256=source_sha,rows_scanned=n)


def extract_matrix(path, target, samples):
    before=path.stat()
    found={}
    with opener(path) as h:
        reader=csv.reader(h, delimiter='\t')
        fields=next(reader)
        need(len(fields)==len(set(fields)) and fields[0]=='feature_id' and set(fields[1:])==set(samples), f'Columnas incompatibles: {path}')
        ii=[fields.index(s) for s in samples]
        for row in reader:
            need(len(row)==len(fields), f'Fila matriz incompleta: {path}')
            fid=row[0]
            if fid in target:
                need(fid not in found, f'ID duplicado: {fid}')
                found[fid]=[float(row[i]) for i in ii]
    need(set(found)==set(target), f'Faltan candidatos en {path}: {set(target)-set(found)}')
    after=path.stat()
    need((before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns), f'Matriz modificada: {path}')
    return found


def focal_data(source_run, meta, selection, tables):
    t=source_run/'tables'
    samples=[r['sample_id'] for r in meta]
    counts=extract_matrix(t/'601b_top_100000_counts_features_x_samples.tsv.gz', FOCAL, samples)
    clr=extract_matrix(t/'601b_top_100000_CLR_features_x_samples.tsv.gz', FOCAL, samples)
    cr=read_table(t/'601b_CLR_centers.tsv')
    need(len(cr)==51 and {r['sample_id'] for r in cr}==set(samples), 'Centros CLR incompatibles')
    centers={r['sample_id']:float(r['mean_log_count_plus1']) for r in cr}
    need(all(math.isfinite(v) for v in centers.values()), 'Centro CLR no finito')
    long, qc=[],[]
    meta_fields=['sample_id','profile_id','locality','restoration4','depth_cm','MHI_local',*AXES]
    for fid in FOCAL:
        cs,ys=counts[fid],clr[fid]
        need(all(math.isfinite(c) and c>=0 for c in cs) and all(math.isfinite(y) for y in ys), f'Valores invalidos: {fid}')
        need(max(abs(math.log1p(c)-centers[s]-y) for c,y,s in zip(cs,ys,samples))<1e-8, f'CLR no reproducible: {fid}')
        st=selection[fid]
        need(sum(c>0 for c in cs)==int(st['prevalence']), f'Prevalencia incompatible: {fid}')
        need(math.isclose(sum(cs),float(st['total_abundance']),rel_tol=1e-9,abs_tol=1e-8), f'Total incompatible: {fid}')
        need(math.isclose(statistics.stdev(ys),float(st['sd_clr']),rel_tol=1e-8,abs_tol=1e-9), f'SD incompatible: {fid}')
        for m,c,y in zip(meta,cs,ys):
            long.append(dict(feature_id=fid,candidate_role=FOCAL[fid],**{k:m[k] for k in meta_fields},count=c,CLR=y,detected=c>0))
        for loc in ('Carmen','Cozumel','Tuxpan'):
            for depth in (5,20,40):
                ii=[i for i,m in enumerate(meta) if m['locality']==loc and float(m['depth_cm'])==depth]
                qc.append(dict(feature_id=fid,candidate_role=FOCAL[fid],locality=loc,depth_cm=depth,
                    n_profiles=len(ii),detections=sum(cs[i]>0 for i in ii),mean_count=statistics.mean(cs[i] for i in ii),
                    mean_CLR=statistics.mean(ys[i] for i in ii)))
    write_table(tables/'604c_focal_profiles.tsv',long,list(long[0]))
    write_table(tables/'604c_focal_detection_by_locality_depth.tsv',qc,list(qc[0]))
    write_table(tables/'604c_focal_selection_stats.tsv',[dict(candidate_role=FOCAL[f],**selection[f]) for f in FOCAL],['candidate_role',*next(iter(selection.values())).keys()])
    return len(long)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out-root',default=ROOT)
    ap.add_argument('--run-603c',default='')
    ap.add_argument('--annotation',default=ANNOTATION)
    ap.add_argument('--audit-dir',default=str(Path(__file__).resolve().parent/'604c_inputs_verificados'))
    args=ap.parse_args()
    root=Path(args.out_root).resolve(strict=True)
    audit=Path(args.audit_dir).resolve(strict=True)
    annotation=Path(args.annotation).resolve(strict=True)
    manifest=json.loads((audit/'manifest.json').read_text())
    log('[1/5] Verificando corrida603c y entradas auditadas')
    for name,sha in manifest['files_sha256'].items():
        need(digest(audit/name)==sha, f'Archivo auditado modificado: {name}')
    run603=Path(args.run_603c).resolve(strict=True) if args.run_603c else latest(root,'603c_asociaciones_bloques_por_profundidad')
    summary=read_table(run603/'tables/603c_run_summary.tsv')
    expected={'status':'PASS','mode':'full','samples':'51','profiles':'17','axes':'5','features_per_axis':'100000',
        'pseudocount':'1','depth_averaged':'FALSE','axes_refitted':'FALSE','selection_uses_predictors':'FALSE',
        'HI_HFR_MHI_in_models':'FALSE','FDR_computed':'TRUE','FDR_scope':'joint_five_axes_within_each_test_family'}
    need(len(summary)==1 and all(summary[0].get(k)==v for k,v in expected.items()),'603c no coincide con la corrida full auditada')
    need(digest(run603/'tables/603c_all_results.tsv.gz')==manifest['source_results_sha256'],'603c difiere del resultado auditado; no mezclar corridas')
    need(digest(run603/'logs/603c_asociaciones_bloques_por_profundidad.R')==manifest['source_script_sha256'],'Script603c distinto del auditado')
    inputs=read_table(run603/'tables/603c_inputs.tsv')
    need(len(inputs)==4,'Procedencia603c inesperada')
    sources={Path(r['input']).name:Path(r['input']) for r in inputs}
    need(len(sources)==4,'Entradas repetidas')
    for row in inputs:
        need(digest(Path(row['input']),'md5')==row['md5'],f"Entrada603c modificada: {row['input']}")
    source_run=sources['601b_metadata_51samples.csv'].parent.parent
    need(all(p.parent.parent==source_run for p in sources.values()),'Entradas de distintas corridas601b')
    meta=read_table(sources['601b_metadata_51samples.csv'],',')
    need(len(meta)==51 and len({m['sample_id'] for m in meta})==51,'Metadata sin51IDs unicos')
    profiles=defaultdict(list)
    for row in meta:
        profiles[row['profile_id']].append(row)
        need(all(math.isfinite(float(row[k])) for k in ['MHI_local',*AXES]),'Predictor no finito')
    need(len(profiles)==17 and all(sorted(float(r['depth_cm']) for r in v)==[5,20,40]
         and len({r['locality'] for r in v})==len({r['restoration4'] for r in v})==1 for v in profiles.values()),'Perfiles incompletos/inconsistentes')
    u=read_table(audit/'603c_annotation_universe_100000.tsv.gz')
    universe={r['feature_id']:r for r in u}
    need(len(u)==len(universe)==100000,'Universo distinto de100000')
    selection_rows=read_table(sources['601b_selected_features_top_100000.tsv'])
    selection={r['feature_id']:r for r in selection_rows}
    need(len(selection_rows)==len(selection)==100000 and set(selection)==set(universe),'IDs no coinciden con seleccion601b')
    need(all(int(selection[f]['rank'])==int(universe[f]['selection_rank']) for f in universe),'Ranking601b discrepante')
    need(set(FOCAL).issubset(universe),'Candidato ausente del universo')
    primary=read_table(audit/'603c_primary_associations.tsv')
    unique=read_table(audit/'603c_primary_unique_proteins.tsv')
    need(len(primary)==1824 and len(unique)==1486 and len({r['feature_id'] for r in primary})==1486,'Lista principal inesperada')
    need({r['feature_id'] for r in unique}=={r['feature_id'] for r in primary},'Lista de proteinas unicas incompatible')
    need(all(r['feature_id'] in universe and r['status']=='ok' and 0<=float(r['q'])<=.10 for r in primary),'Asociacion principal incompatible')
    for row in primary:
        model='axis_common' if row['family']=='axis_common' else 'axis_by_depth'
        need(universe[row['feature_id']][row['axis']+'__'+model+'__status']=='ok','Hit fuera del fondo evaluable')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    run=root/(NAME+'_'+stamp)
    run.mkdir(exist_ok=False)
    tables,logs=run/'tables',run/'logs'
    tables.mkdir();logs.mkdir()
    log('Output: '+str(run))
    shutil.copy2(Path(__file__),logs/Path(__file__).name)
    shutil.copy2(audit/'manifest.json',logs/'604c_audit_manifest.json')
    for name in ['603c_run_summary.tsv','603c_inputs.tsv','603c_axis_audit.tsv','603c_metadata_used.tsv','603c_test_families.tsv']:
        shutil.copy2(run603/'tables'/name,tables/name)
    log('[2/5] Recuperando perfiles de los cinco candidatos; verificando CLR completo')
    focal_rows=focal_data(source_run,meta,selection,tables)
    log('[3/5] Lectura unica del catalogo completo de anotacion; se revisa hasta el final')
    hits,duplicates,rawfields,annotation_source=scan_annotations(annotation,universe)
    log(f'Lectura completada: {annotation_source["rows_scanned"]:,} filas; {len(hits):,}/100,000 IDs encontrados')
    annotated={}
    issues=[]
    mappings=[]
    ann_fields=['annotation_found','annotation_status','KO_ids','PFAM_terms','has_KO','has_PFAM',*['annotation_'+f for f in rawfields]]
    for fid in universe:
        raw=hits.get(fid,{})
        ko,badko=tokens(raw.get('kegg_ko',''),'KO')
        pf,badpf=tokens(raw.get('pfams',''),'PFAM')
        if fid not in hits:
            status='query_not_found'
            issues.append(dict(feature_id=fid,issue=status,detail=''))
        elif fid in duplicates:
            status='duplicate_query'
            issues.append(dict(feature_id=fid,issue=status,detail=str(duplicates[fid]+1)+' source rows'))
        elif badko or badpf:
            status='invalid_annotation_token'
        else:
            status='matched'
        for ont,bad in [('KO',badko),('PFAM',badpf)]:
            for value in bad:issues.append(dict(feature_id=fid,issue='invalid_'+ont+'_token',detail=value))
        # Ambiguous or malformed rows remain visible, but cannot populate term mappings.
        valid=status=='matched'
        annotated[fid]=dict(annotation_found=fid in hits,annotation_status=status,
            KO_ids=';'.join(ko),PFAM_terms=';'.join(pf),has_KO=valid and bool(ko),has_PFAM=valid and bool(pf),
            **{'annotation_'+f:raw.get(f,'') for f in rawfields})
        if valid:
            for ont,terms in [('KO',ko),('PFAM',pf)]:
                mappings.extend(dict(feature_id=fid,ontology=ont,term_id=term) for term in terms)
    log('[4/5] Exportando anotaciones, cobertura y fondos evaluables; sin enriquecimiento')
    write_table(tables/'604c_annotation_universe.tsv.gz',
        (dict(**row,**annotated[row['feature_id']]) for row in u),[*u[0],*ann_fields])
    write_table(tables/'604c_protein_to_KO_PFAM.tsv.gz',mappings,['feature_id','ontology','term_id'])
    write_table(tables/'604c_primary_associations_annotated.tsv',
        (dict(**row,**annotated[row['feature_id']]) for row in primary),[*primary[0],*ann_fields])
    write_table(tables/'604c_primary_proteins_annotated.tsv',
        (dict(**row,**annotated[row['feature_id']]) for row in unique),[*unique[0],*ann_fields])
    write_table(tables/'604c_focal_proteins_annotated.tsv',
        (dict(feature_id=f,candidate_role=FOCAL[f],**annotated[f]) for f in FOCAL),['feature_id','candidate_role',*ann_fields])
    write_table(tables/'604c_annotation_issues.tsv',issues,['feature_id','issue','detail'])
    background=[]
    eligibility={}
    for axis in AXES:
        for model in ('axis_common','axis_by_depth'):
            col=axis+'__'+model+'__status'
            need(all(col in row for row in u),f'Falta estado del modelo: {col}')
            eligible={f for f,row in universe.items() if row[col]=='ok'}
            for ont in ('KO','PFAM'):
                ann={f for f in eligible if annotated[f]['has_'+ont]}
                eligibility[(axis,model,ont)]=ann
                background.append(dict(axis=axis,model=model,ontology=ont,n_tested=100000,
                    n_primary_evaluable=len(eligible),n_evaluable_annotated=len(ann)))
    write_table(tables/'604c_background_coverage.tsv',background,list(background[0]))
    groups=defaultdict(set)
    for row in primary:
        direction='joint' if row['family']=='axis_by_depth' else ('positive' if float(row['estimate'])>0 else 'negative')
        depth=str(int(float(row['depth_cm']))) if row['family']=='axis_depth_slopes' else 'NA'
        groups[(row['axis'],row['family'],depth,direction)].add(row['feature_id'])
    coverage=[]
    for axis in AXES:
        for fam in ('axis_common','axis_by_depth','axis_depth_slopes'):
            for depth in (('5','20','40') if fam=='axis_depth_slopes' else ('NA',)):
                for direction in (('joint',) if fam=='axis_by_depth' else ('positive','negative')):
                    selected=groups[(axis,fam,depth,direction)]
                    model='axis_common' if fam=='axis_common' else 'axis_by_depth'
                    for ont in ('KO','PFAM'):
                        bg=eligibility[(axis,model,ont)]
                        coverage.append(dict(axis=axis,family=fam,depth_cm=depth,direction=direction,ontology=ont,
                            n_primary_proteins=len(selected),n_primary_annotated=len(selected & bg),n_background_annotated=len(bg)))
    write_table(tables/'604c_primary_annotation_coverage.tsv',coverage,list(coverage[0]))
    # Preserve source hashes and check small audited inputs once more before publishing.
    for row in inputs:
        if Path(row['input']).name!='601b_top_100000_CLR_features_x_samples.tsv.gz':
            need(digest(Path(row['input']),'md5')==row['md5'],'Entrada601b cambio durante ejecucion')
    status='PASS' if not issues else 'REVIEW_REQUIRED'
    result=dict(status=status,samples=51,profiles=17,universe_proteins=100000,primary_associations=1824,
        primary_unique_proteins=1486,annotation_queries_found=len(hits),queries_not_found=100000-len(hits),
        duplicate_queries=len(duplicates),invalid_token_count=sum(x['issue'].startswith('invalid_') for x in issues),
        proteins_with_KO=sum(v['has_KO'] for v in annotated.values()),proteins_with_PFAM=sum(v['has_PFAM'] for v in annotated.values()),
        focal_candidates=len(FOCAL),focal_sample_rows=focal_rows,CLR_verified=True,pseudocount=1,
        models_refitted=False,FDR_recalculated=False,enrichment_executed=False,
        interpretation='functional_gene_potential_not_activity')
    write_table(tables/'604c_run_summary.tsv',[result],list(result))
    provenance=dict(run603c=str(run603),run601b=str(source_run),audit_manifest=manifest,
        verified_601b_inputs=inputs,annotation_source=annotation_source,
        primary_selection='unchanged_603c_q_le_0.10_status_ok',diagnostic_singular_hits_merged=False,
        mapping='exact_query; KO prefix removed; Pfam identifiers preserved; unique feature-term edges',
        focal_roles=FOCAL,warning='Three MHI candidates are diagnostic cases, not three confirmed MHI associations.')
    (logs/'604c_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    archive=root/('figura4_anotacion_604c_resultados_'+stamp+'.tar.gz')
    log('[5/5] Preparando paquete de resultados')
    with tarfile.open(archive,'w:gz') as t:
        t.add(run,arcname=run.name)
    pointer=publish(root,run,status)
    log(json.dumps(result))
    print('Pointer: '+str(pointer),flush=True)
    print('Adjunta este archivo:\n'+str(archive),flush=True)
    if status!='PASS':
        print('Revisar 604c_annotation_issues.tsv; no se publica el puntero PASS ni se ejecuta enriquecimiento.',flush=True)


if __name__=='__main__':
    main()
