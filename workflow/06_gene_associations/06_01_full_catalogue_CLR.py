#!/usr/bin/env python3
# ============================================================
# 06_01_full_catalogue_CLR.py
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# CLR centring across 1,475,487 genes before selection of 100,000 variable genes.
# Inputs (source expressions; complete list in docs/contracts/06_01_full_catalogue_CLR.json):
#   Dynamic input discovery: see the contract and original source below.
# Outputs (source expressions; complete list in contract):
#   (out / "logs/601b_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
#   tmp.write_text(str(out) + "\n")
#   (out / "logs/601b_FAILED.txt").write_text(str(exc) + "\n")
# Algorithmic provenance:
# Full-catalogue log(count+1) centring before selection of 100,000 variable genes.
#   Bates et al. (2015), doi:10.18637/jss.v067.i01; Kuznetsova et al. (2017), doi:10.18637/jss.v082.i13; Benjamini & Hochberg (1995), doi:10.1111/j.2517-6161.1995.tb02031.x.
# Source SHA-256: 7ad2ebdaf09f5c28186d64ea4b9e2dc38daf48c0d75ed0aad4fe3b5004faa849
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================

"""Prepare protein CLR per sample; keep depths and select by variance only."""
import argparse
import csv
import gzip
import hashlib
import heapq
import json
import math
import os
import shutil
from array import array
from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path

NAME = "601b_preparar_CLR_proteico_por_muestra"
BASE = Path("/home/fjbalvino/Tipping_points")
CORE = ["sample_id", "profile_id", "locality", "restoration4", "depth_cm", "lat_block"]
AXES = ["vegetation_landscape_PC1", "water_inundation_PC1", "moisture_stress_PC1",
        "physicochemical_PC1", "nutrients_redox_PC1"]
HASH_MHI = "d9fb0bef5be343d76cdf0bd4ef5b4c7059386e9102e3a670dd67b34a1b3aaf1b"
HASH_CANON = "8232edd8008922d26f289ee8db9e273255bd3c1b360baf277fbe5430fac06993"


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path, sep="\t"):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=sep)
        fields = reader.fieldnames or []
        require(fields and len(fields) == len(set(fields)), f"Cabecera invalida: {path}")
        rows = list(reader)
    require(rows and all(set(r) == set(fields) and None not in r.values() for r in rows),
            f"Tabla vacia o filas incompletas: {path}")
    return rows


def keyed(rows, key):
    require(all(key in r and r[key].strip() for r in rows), f"Falta ID: {key}")
    require(len({r[key] for r in rows}) == len(rows), f"IDs duplicados: {key}")
    return {r[key]: r for r in rows}


def write_table(path, rows, fields=None, sep="\t"):
    fields = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=sep)
        writer.writeheader()
        writer.writerows(rows)


def metadata(canon_path, mhi_path, old_path):
    require(sha(canon_path) == HASH_CANON, "Canon202 distinto del auditado")
    require(sha(mhi_path) == HASH_MHI, "Metadata003b distinta de la validada")
    canon = read_table(canon_path, ",")
    current = keyed(read_table(mhi_path, ","), "sample_id")
    old = keyed(read_table(old_path, ","), "sample_id")
    ids = list(keyed(canon, "sample_id"))
    require(len(ids) == 51 and set(ids) == set(current) == set(old), "Universos incompatibles")
    profiles, result = {}, []
    for row in canon:
        require(all(row.get(k, "").strip() for k in CORE), "Canon con campos vacios")
        sid = row["sample_id"]
        for key in ["locality", "restoration4", "depth_cm", "lat_block"]:
            values = [r[key] for r in (row, current[sid], old[sid])]
            if key in ("depth_cm", "lat_block"):
                values = [Decimal(v) for v in values]
                require(all(v.is_finite() for v in values), f"Campo invalido: {sid}/{key}")
            require(len(set(values)) == 1, f"Metadata discrepante: {sid}/{key}")
        new = {k: row[k] for k in CORE}
        for key in ["MHI_local"] + AXES:
            value = current[sid].get(key, "")
            require(value and math.isfinite(float(value)), f"Predictor invalido: {sid}/{key}")
            new[key] = value
        require(-1 <= float(new["MHI_local"]) <= 1, f"MHI fuera de rango: {sid}")
        result.append(new)
        profiles.setdefault(row["profile_id"], []).append(row)
    require(len(profiles) == 17, "Se requieren 17 perfiles")
    audit = []
    for pid, rows in profiles.items():
        require(sorted(Decimal(r["depth_cm"]) for r in rows) == [5, 20, 40],
                f"Perfil incompleto: {pid}")
        require(all(len({r[k] for r in rows}) == 1 for k in
                    ["locality", "restoration4", "lat_block"]), f"Perfil inconsistente: {pid}")
        audit.append(dict(profile_id=pid, locality=rows[0]["locality"],
                          restoration4=rows[0]["restoration4"], n_samples=3,
                          depths="5|20|40", complete=True))
    return result, audit, ids


def matrix_rows(path, ids):
    with gzip.open(path, "rt", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        require(header and header[0] == "feature_id", "Falta columna feature_id")
        require(len(header) == 52 and len(set(header)) == 52 and set(header[1:]) == set(ids),
                "La matriz debe contener exactamente las 51 muestras canonicas")
        positions = [header.index(s) for s in ids]
        for i, row in enumerate(reader, 1):
            require(len(row) == len(header), f"Fila incompleta: {i + 1}")
            fid = row[0]
            require(fid and fid == fid.strip(), f"feature_id invalido: fila {i + 1}")
            try:
                values = [float(row[j]) for j in positions]
            except ValueError as exc:
                raise ValueError(f"Valor no numerico: {fid}") from exc
            require(all(math.isfinite(x) and x >= 0 for x in values),
                    f"Valor negativo/no finito: {fid}")
            yield i, fid, values


def first_pass(path, ids, expected_n):
    sums, detected, log_sums = [0.0] * 51, [0] * 51, [0.0] * 51
    seen, n = set(), 0
    for n, fid, values in matrix_rows(path, ids):
        require(fid not in seen, f"feature_id duplicado: {fid}")
        seen.add(fid)
        require(sum(v > 0 for v in values) >= 5 and math.fsum(values) >= 10 - 1e-9,
                f"Fila incompatible con filtro044: {fid}")
        for j, value in enumerate(values):
            sums[j] += value
            detected[j] += int(value > 0)
            log_sums[j] += math.log1p(value)
        if n % 100000 == 0:
            print(f"  Lectura1: {n:,} filas validadas", flush=True)
    require(n == expected_n, f"Filas observadas {n}; esperadas {expected_n}")
    return [x / n for x in log_sums], sums, detected


def select_features(path, ids, centers, top_n, expected_n):
    heap, positive, n = [], 0, 0
    for n, fid, values in matrix_rows(path, ids):
        clr = array("d", (math.log1p(v) - m for v, m in zip(values, centers)))
        mean = math.fsum(clr) / len(ids)
        variance = math.fsum((v - mean) ** 2 for v in clr) / (len(ids) - 1)
        require(math.isfinite(variance), f"Varianza invalida: {fid}")
        if variance > 0:
            positive += 1
            if len(heap) < top_n or (variance, -n) > heap[0][:2]:
                item = (variance, -n, fid, mean, sum(v > 0 for v in values),
                        math.fsum(values), clr, array("d", values))
                if len(heap) < top_n:
                    heapq.heappush(heap, item)
                else:
                    heapq.heapreplace(heap, item)
        if n % 100000 == 0:
            print(f"  Lectura2: {n:,} filas; seleccionadas {len(heap):,}", flush=True)
    require(n == expected_n and len(heap) == top_n, "Filas o seleccion final insuficientes")
    return sorted(heap, key=lambda x: (-x[0], -x[1])), positive


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run044", type=Path, default=BASE / "Results/044_build_protein_filtered_matrix_HI_20260529_000843")
    parser.add_argument("--run003b", type=Path, default=BASE / "resultados_finales/003b_integrar_MHI_ambiental_validado_20260915_192526")
    parser.add_argument("--canon", type=Path, default=BASE / "resultados_finales/202_construir_matriz_CLR_y_eje_bimodal_20260912_005021/tables/202_metadata_raw1031_alineada.csv")
    parser.add_argument("--out_root", type=Path, default=BASE / "resultados_finales")
    parser.add_argument("--top_n", type=int, default=100000)
    args = parser.parse_args()
    require(args.top_n > 0, "top_n debe ser positivo")
    t44, t3b = args.run044 / "tables", args.run003b / "tables"
    sources = {
        "canon202": args.canon,
        "metadata003b": t3b / "003_metadata_integrada_canon_51.csv",
        "summary003b": t3b / "003b_run_summary.tsv",
        "metadata044": t44 / "044_metadata_aligned_51samples.csv",
        "filter044": t44 / "044_filter_summary.tsv",
        "qc_filtered044": t44 / "044_sample_qc_filtered_matrix.tsv",
        "qc_raw044": t44 / "044_sample_qc_raw_common_samples.tsv",
    }
    matrix = t44 / "044_protein_filtered_counts_51samples.tsv.gz"
    for path in list(sources.values()) + [matrix]:
        require(path.is_file(), f"No existe: {path}")
    run3b = read_table(sources["summary003b"])
    require(len(run3b) == 1 and run3b[0]["status"] == "PASS", "003b no tiene PASS")
    require(int(run3b[0]["samples"]) == 51 and int(run3b[0]["profiles"]) == 17,
            "003b con universo inesperado")
    md, profiles, ids = metadata(sources["canon202"], sources["metadata003b"], sources["metadata044"])
    summary44 = {k: v["value"] for k, v in keyed(read_table(sources["filter044"]), "key").items()}
    require(float(summary44["min_prevalence"]) == 5 and float(summary44["min_total"]) == 10
            and float(summary44["min_sd"]) == 0, "Filtro044 distinto del auditado")
    expected_n = int(summary44["n_features_retained"])
    require(expected_n >= args.top_n, "Catalogo menor que top_n")
    qc_old = keyed(read_table(sources["qc_filtered044"]), "sample_id")
    qc_raw = keyed(read_table(sources["qc_raw044"]), "sample_id")
    require(set(qc_old) == set(qc_raw) == set(ids), "QC con muestras incompatibles")
    hashes = {key: sha(path) for key, path in sources.items()}
    now = datetime.now(timezone.utc)
    out = args.out_root.resolve() / f"{NAME}_{now:%Y%m%d_%H%M%S}"
    out.mkdir(parents=True, exist_ok=False)
    for name in ("tables", "inputs", "logs"):
        (out / name).mkdir()
    print(f"Output: {out}\n51 muestras / 17 perfiles; pc=1; top_n={args.top_n}", flush=True)
    try:
        print("[1/4] Verificando matriz y calculando centros CLR completos", flush=True)
        matrix_hash = sha(matrix)
        centers, totals, detected = first_pass(matrix, ids, expected_n)
        qc = []
        for j, sid in enumerate(ids):
            require(math.isclose(totals[j], float(qc_old[sid]["filtered_total_abundance"]),
                                 rel_tol=1e-9, abs_tol=1e-6), f"Suma distinta de044: {sid}")
            require(detected[j] == int(qc_old[sid]["filtered_detected_features"]),
                    f"Detecciones distintas de044: {sid}")
            raw_total = float(qc_raw[sid]["raw_total_abundance_common_samples"])
            require(math.isfinite(raw_total) and raw_total >= totals[j] > 0,
                    f"Totales invalidos: {sid}")
            qc.append(dict(sample_id=sid, filtered_total_abundance=totals[j],
                           filtered_detected_features=detected[j], raw_total_abundance=raw_total,
                           retained_abundance_fraction=totals[j] / raw_total))
        print("[2/4] Seleccion por varianza CLR; sin indices ni etiquetas", flush=True)
        selected, positive = select_features(matrix, ids, centers, args.top_n, expected_n)
        require(sha(matrix) == matrix_hash, "La matriz cambio durante la lectura")
        require(all(sha(path) == hashes[key] for key, path in sources.items()),
                "Una tabla de entrada cambio durante la lectura")
        print("[3/4] Exportando matrices y metadata por muestra", flush=True)
        tables = out / "tables"
        stem = f"601b_top_{args.top_n}"
        with gzip.open(tables / f"{stem}_CLR_features_x_samples.tsv.gz", "wt", newline="", encoding="utf-8") as ch, \
             gzip.open(tables / f"{stem}_counts_features_x_samples.tsv.gz", "wt", newline="", encoding="utf-8") as rh:
            cw, rw = csv.writer(ch, delimiter="\t"), csv.writer(rh, delimiter="\t")
            cw.writerow(["feature_id"] + ids)
            rw.writerow(["feature_id"] + ids)
            feature_rows, selected_detected = [], [0] * len(ids)
            for rank, (var, negrow, fid, mean, prev, total, clr, counts) in enumerate(selected, 1):
                cw.writerow([fid] + [format(v, ".17g") for v in clr])
                rw.writerow([fid] + [format(v, ".17g") for v in counts])
                for j, v in enumerate(counts):
                    selected_detected[j] += int(v > 0)
                feature_rows.append(dict(feature_id=fid, rank=rank, mean_clr=mean, var_clr=var,
                                         sd_clr=math.sqrt(var), prevalence=prev,
                                         total_abundance=total, source_row=-negrow))
        require(all(n > 0 for n in selected_detected), "Hay muestras vacias en el subconjunto seleccionado")
        for j, row in enumerate(qc):
            row["selected_detected_features"] = selected_detected[j]
        write_table(tables / f"601b_selected_features_top_{args.top_n}.tsv", feature_rows)
        write_table(tables / "601b_metadata_51samples.csv", md, sep=",")
        write_table(tables / "601b_profile_audit.tsv", profiles)
        write_table(tables / "601b_sample_QC.tsv", qc)
        write_table(tables / "601b_CLR_centers.tsv", [dict(sample_id=s, mean_log_count_plus1=m,
                    catalog_features=expected_n) for s, m in zip(ids, centers)])
        for key, path in sources.items():
            shutil.copy2(path, out / "inputs" / f"{key}{path.suffix}")
        shutil.copy2(Path(__file__), out / "logs" / Path(__file__).name)
        summary = dict(status="PASS", samples=51, profiles=17, catalog_features=expected_n,
                       positive_CLR_variance_features=positive, selected_features=len(selected),
                       pseudocount=1, depth_averaged=False, selection_uses_HI_MHI_HFR_stage=False,
                       CLR_center="full_filtered_catalog", selection="global_CLR_variance",
                       tie_break="original_matrix_row_order", models_executed=False)
        provenance = dict(started_utc=now.isoformat(), finished_utc=datetime.now(timezone.utc).isoformat(),
                          matrix=dict(path=str(matrix.resolve()), sha256=matrix_hash),
                          sources={k: dict(path=str(p.resolve()), sha256=hashes[k]) for k, p in sources.items()},
                          script_sha256=sha(Path(__file__)), summary=summary)
        (out / "logs/601b_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
        print("[4/4] Comprobaciones completas; publicando PASS", flush=True)
        write_table(tables / "601b_run_summary.tsv", [summary])
        pointer = args.out_root.resolve() / f"LATEST_{NAME}.txt"
        tmp = pointer.with_name(pointer.name + f".tmp_{os.getpid()}")
        tmp.write_text(str(out) + "\n")
        tmp.replace(pointer)
        print(json.dumps(summary), flush=True)
        print(f"Pointer: {pointer}", flush=True)
    except Exception as exc:
        (out / "logs/601b_FAILED.txt").write_text(str(exc) + "\n")
        raise


if __name__ == "__main__":
    main()
