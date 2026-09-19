#!/usr/bin/env python3
# ============================================================
# 02_02_filter_gene_counts.py
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Maintained successor of historical 044. The manuscript consumed archived 044 counts; equivalent output from this successor has not been demonstrated.
# Inputs (source expressions; complete list in docs/contracts/02_02_filter_gene_counts.json):
#   run_003 = Path(latest_003.read_text(encoding="utf-8").strip())
# Outputs (source expressions; complete list in contract):
#   def write_tsv(path, header, rows):
#   write_tsv(
# Algorithmic provenance:
# Stream gene counts; apply prevalence, total-abundance and variation filters.
#   Cantalapiedra et al. (2021), eggNOG-mapper v2, doi:10.1093/molbev/msab293; this is annotation provenance, not evidence that this script runs eggNOG.
# Source SHA-256: 9e5df3dd06b93c6101d62ce6f20df9a4efa01942d45c34b586258d56dba1d42b
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================

# -*- coding: utf-8 -*-

"""
502_construir_matriz_proteica_filtrada.py

Construye una matriz analítica proteína/ORF x muestra para análisis HI,
usando solo las muestras comunes entre la matriz proteica y la metadata HI.

Input principal:
  FINAL_ANNOTATED_ABUNDANCE_noMZ_bacteria_fran.tsv

Estructura esperada:
  columna 1 = #query / feature_id
  columnas restantes = muestras

Output:
  Results/502_construir_matriz_proteica_filtrada_YYYYMMDD_HHMMSS/
    tables/
      502_protein_filtered_counts_<N>samples.tsv.gz
      502_retained_feature_stats.tsv.gz
      502_metadata_aligned_<N>samples.csv
      502_sample_matching.tsv
      502_excluded_samples.tsv
      502_sample_qc_raw_common_samples.tsv
      502_sample_qc_filtered_matrix.tsv
      502_retained_prevalence_distribution.tsv
      502_retained_total_abundance_distribution.tsv
      502_retained_origin_sample_counts.tsv
      502_top_retained_features_by_total.tsv
      502_filter_summary.tsv
    logs/
      502_log.txt

Este script NO calcula CLR ni modelos. Eso queda para 601–603.
"""

import argparse
import csv
import gzip
import heapq
import math
import os
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path


def timestamp_now() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def log(msg: str, log_handle=None) -> None:
    text = str(msg)
    print(text, flush=True)
    if log_handle is not None:
        log_handle.write(text + "\n")
        log_handle.flush()


def clean_colname(x: str) -> str:
    x = re.sub(r"\s+", "_", x)
    x = re.sub(r"[^A-Za-z0-9_.]+", "_", x)
    x = re.sub(r"_+", "_", x)
    x = re.sub(r"^_|_$", "", x)
    return x.lower()


def detect_origin_sample(query_id: str) -> str:
    q = str(query_id)

    for pat in [r"\.contigs_.*$", r"\.contig_.*$", r"_contigs_.*$"]:
        if re.search(pat, q):
            return re.sub(pat, "", q)

    return re.sub(r"\..*$", "", q)


def parse_float_or_zero(x: str):
    try:
        val = float(x)
        if not math.isfinite(val):
            return 0.0, False
        return val, True
    except Exception:
        return 0.0, False


def format_number_for_matrix(x: str) -> str:
    """
    Mantiene el string original si es numérico legible.
    Si viene vacío o raro, lo convierte a 0.
    """
    if x is None:
        return "0"
    xs = str(x).strip()
    if xs == "":
        return "0"
    try:
        val = float(xs)
        if not math.isfinite(val):
            return "0"
        if val.is_integer():
            return str(int(val))
        return str(val)
    except Exception:
        return "0"


def total_bin(x: float) -> str:
    if x <= 0:
        return "0"
    if x <= 1:
        return "1"
    if x <= 5:
        return "2-5"
    if x <= 10:
        return "6-10"
    if x <= 20:
        return "11-20"
    if x <= 50:
        return "21-50"
    if x <= 100:
        return "51-100"
    if x <= 500:
        return "101-500"
    if x <= 1000:
        return "501-1000"
    if x <= 10000:
        return "1001-10000"
    return ">10000"


def read_metadata(meta_csv: str, id_col_requested: str):
    with open(meta_csv, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        original_fields = reader.fieldnames
        if original_fields is None:
            raise RuntimeError(f"Metadata vacía o inválida: {meta_csv}")

        clean_fields = [clean_colname(x) for x in original_fields]
        field_map = dict(zip(clean_fields, original_fields))

        rows = []
        for row in reader:
            clean_row = {}
            for clean_name, original_name in field_map.items():
                clean_row[clean_name] = row.get(original_name, "")
            rows.append(clean_row)

    candidate_cols = []
    for c in [id_col_requested, "sample_id", "id_metagenomics", "run_accession", "sample", "sample_name"]:
        cc = clean_colname(c)
        if cc in clean_fields and cc not in candidate_cols:
            candidate_cols.append(cc)

    if not candidate_cols:
        raise RuntimeError(
            "No encontré columnas ID candidatas en metadata. Columnas disponibles: "
            + ", ".join(original_fields)
        )

    return rows, candidate_cols, original_fields, clean_fields


def unique_nonempty(values):
    seen = set()
    out = []
    for v in values:
        if v is None:
            continue
        vv = str(v).strip()
        if vv == "":
            continue
        if vv not in seen:
            out.append(vv)
            seen.add(vv)
    return out


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def write_csv_rows(path, header, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in header})


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input_tsv",
        default="/home/fjbalvino/Tipping_points/FINAL_ANNOTATED_ABUNDANCE_noMZ_bacteria_fran.tsv",
    )
    parser.add_argument(
        "--meta_csv",
        default="",
    )
    parser.add_argument(
        "--out_root",
        default="/home/fjbalvino/Tipping_points/resultados_finales",
    )
    parser.add_argument("--id_col", default="sample_id")
    parser.add_argument("--min_prevalence", type=int, default=5)
    parser.add_argument("--min_total", type=float, default=10)
    parser.add_argument("--min_sd", type=float, default=0)
    parser.add_argument("--progress_every", type=int, default=100000)
    parser.add_argument("--top_n", type=int, default=10000)
    parser.add_argument("--expected_n", type=int, default=51)
    parser.add_argument("--expected_profiles", type=int, default=17)

    args = parser.parse_args()

    if args.expected_n != 51 or args.expected_profiles != 17:
        raise RuntimeError(
            "El canon esta congelado en 51 muestras y 17 perfiles; no se permite redefinirlo"
        )

    input_tsv = Path(args.input_tsv)
    out_root = Path(args.out_root)

    if args.meta_csv.strip():
        meta_csv = Path(args.meta_csv)
    else:
        latest_003 = out_root / "LATEST_003_integrar_ejes_ECI_HI_en_metadata.txt"
        if not latest_003.is_file():
            raise FileNotFoundError(f"No existe latest de metadata integrada: {latest_003}")
        run_003 = Path(latest_003.read_text(encoding="utf-8").strip())
        meta_csv = run_003 / "tables" / "003_metadata_integrada_canon_51.csv"

    if not input_tsv.exists():
        raise FileNotFoundError(f"No existe input_tsv: {input_tsv}")

    if not meta_csv.exists():
        raise FileNotFoundError(f"No existe meta_csv: {meta_csv}")

    if args.min_prevalence < 1:
        raise ValueError("--min_prevalence debe ser >= 1")

    if args.min_total < 0:
        raise ValueError("--min_total debe ser >= 0")

    if args.min_sd < 0:
        raise ValueError("--min_sd debe ser >= 0")

    out_dir = out_root / f"502_construir_matriz_proteica_filtrada_{timestamp_now()}"
    tables_dir = out_dir / "tables"
    logs_dir = out_dir / "logs"
    figures_dir = out_dir / "figures"

    tables_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    log_path = logs_dir / "502_log.txt"

    with open(log_path, "w") as log_handle:
        log("===== 502 construir matriz proteica filtrada =====", log_handle)
        log(f"Start: {datetime.now()}", log_handle)
        log(f"Input: {input_tsv}", log_handle)
        log(f"Metadata: {meta_csv}", log_handle)
        log(f"Out dir: {out_dir}", log_handle)
        log(f"min_prevalence: {args.min_prevalence}", log_handle)
        log(f"min_total: {args.min_total}", log_handle)
        log(f"min_sd: {args.min_sd}", log_handle)
        log(f"progress_every: {args.progress_every}", log_handle)
        log(f"top_n: {args.top_n}", log_handle)

        # ----------------------------------------------------------------------
        # Header matriz
        # ----------------------------------------------------------------------
        with open(input_tsv, "r", newline="") as fh:
            header_line = fh.readline().rstrip("\n")
        header = header_line.split("\t")

        if len(header) < 2:
            raise RuntimeError("La matriz tiene menos de 2 columnas.")

        query_col = header[0]
        matrix_samples = header[1:]
        matrix_sample_set = set(matrix_samples)

        sample_to_index = {s: i for i, s in enumerate(header)}

        log(f"Columns detected: {len(header)}", log_handle)
        log(f"Feature/query column: {query_col}", log_handle)
        log(f"Matrix sample columns: {len(matrix_samples)}", log_handle)
        log(f"First matrix samples: {', '.join(matrix_samples[:10])}", log_handle)

        # ----------------------------------------------------------------------
        # Metadata matching
        # ----------------------------------------------------------------------
        meta_rows, candidate_cols, original_fields, clean_fields = read_metadata(
            str(meta_csv),
            args.id_col,
        )

        matching_rows = []
        best_col = None
        best_common = -1

        for col in candidate_cols:
            meta_ids = unique_nonempty([r.get(col, "") for r in meta_rows])
            meta_set = set(meta_ids)
            common = [s for s in matrix_samples if s in meta_set]

            if len(common) > best_common:
                best_common = len(common)
                best_col = col

            matching_rows.append([
                col,
                len(meta_ids),
                len(matrix_samples),
                len(common),
                ";".join([s for s in matrix_samples if s not in meta_set]),
                ";".join([s for s in meta_ids if s not in matrix_sample_set]),
            ])

        write_tsv(
            tables_dir / "502_sample_matching.tsv",
            [
                "metadata_id_column",
                "n_metadata_ids",
                "n_matrix_samples",
                "n_common",
                "matrix_not_in_metadata",
                "metadata_not_in_matrix",
            ],
            matching_rows,
        )

        requested_clean_id = clean_colname(args.id_col)

        if requested_clean_id in candidate_cols:
            chosen_id_col = requested_clean_id
        else:
            chosen_id_col = best_col

        meta_ids_ordered = []
        meta_aligned_rows = []
        seen_meta = set()

        for row in meta_rows:
            sid = str(row.get(chosen_id_col, "")).strip()
            if sid == "":
                continue
            if sid in matrix_sample_set and sid not in seen_meta:
                rr = dict(row)
                rr["sample_id"] = sid
                rr["sample_id_for_matrix"] = sid
                meta_aligned_rows.append(rr)
                meta_ids_ordered.append(sid)
                seen_meta.add(sid)

        common_samples = meta_ids_ordered
        common_set = set(common_samples)
        common_indices = [sample_to_index[s] for s in common_samples]

        excluded_matrix_samples = [s for s in matrix_samples if s not in common_set]
        metadata_ids = unique_nonempty([r.get(chosen_id_col, "") for r in meta_rows])
        excluded_metadata_ids = [s for s in metadata_ids if s not in matrix_sample_set]

        canonical_ids = unique_nonempty([r.get(chosen_id_col, "") for r in meta_rows])
        if len(meta_rows) != args.expected_n or len(canonical_ids) != args.expected_n:
            raise RuntimeError(
                f"La metadata debe contener exactamente {args.expected_n} muestras canonicas unicas"
            )
        if len(common_samples) != args.expected_n:
            raise RuntimeError(
                f"La matriz debe contener las {args.expected_n} muestras canonicas; encontro {len(common_samples)}"
            )

        required_design = {"lat_block", "depth_cm", "locality", "restoration4"}
        missing_design = sorted(required_design.difference(clean_fields))
        if missing_design:
            raise RuntimeError("Faltan columnas de diseno: " + ", ".join(missing_design))
        profiles = {}
        for row in meta_aligned_rows:
            key = (str(row.get("locality", "")).strip(), str(row.get("lat_block", "")).strip())
            try:
                depth = int(float(str(row.get("depth_cm", "")).strip()))
            except ValueError as exc:
                raise RuntimeError("depth_cm invalido en metadata canonica") from exc
            profiles.setdefault(key, []).append(depth)
        if len(profiles) != args.expected_profiles or any(
            len(depths) != 3 or set(depths) != {5, 20, 40}
            for depths in profiles.values()
        ):
            raise RuntimeError(
                f"La metadata debe contener exactamente {args.expected_profiles} perfiles completos 5/20/40"
            )

        log(f"Chosen metadata ID column: {chosen_id_col}", log_handle)
        log(f"Common samples retained: {len(common_samples)}", log_handle)
        log(f"Excluded matrix samples: {len(excluded_matrix_samples)}", log_handle)
        log(f"Excluded matrix samples list: {', '.join(excluded_matrix_samples)}", log_handle)

        # Metadata alineada
        meta_header = ["sample_id", "sample_id_for_matrix"] + [
            c for c in clean_fields if c not in {"sample_id", "sample_id_for_matrix"}
        ]
        meta_header = list(dict.fromkeys(meta_header))

        write_csv_rows(
            tables_dir / f"502_metadata_aligned_{len(common_samples)}samples.csv",
            meta_header,
            meta_aligned_rows,
        )

        excluded_rows = []
        for s in excluded_matrix_samples:
            excluded_rows.append(["matrix_sample_excluded_no_metadata", s])
        for s in excluded_metadata_ids:
            excluded_rows.append(["metadata_sample_missing_from_matrix", s])

        write_tsv(
            tables_dir / "502_excluded_samples.tsv",
            ["category", "sample_id"],
            excluded_rows,
        )

        # ----------------------------------------------------------------------
        # Output writers
        # ----------------------------------------------------------------------
        filtered_matrix_path = tables_dir / f"502_protein_filtered_counts_{len(common_samples)}samples.tsv.gz"
        retained_stats_path = tables_dir / "502_retained_feature_stats.tsv.gz"

        n_common = len(common_samples)

        raw_col_sums = [0.0] * n_common
        raw_col_detected = [0] * n_common

        filt_col_sums = [0.0] * n_common
        filt_col_detected = [0] * n_common

        prev_counter = Counter()
        total_bin_counter = Counter()
        origin_counter = Counter()

        top_heap = []
        heap_counter = 0

        n_seen = 0
        n_retained = 0
        n_dropped = 0
        n_negative = 0
        n_nonfinite = 0

        log("Streaming input and writing filtered matrix...", log_handle)

        with open(input_tsv, "r", newline="") as fin, \
             gzip.open(filtered_matrix_path, "wt", newline="") as fout_matrix, \
             gzip.open(retained_stats_path, "wt", newline="") as fout_stats:

            reader = csv.reader(fin, delimiter="\t")
            next(reader)

            matrix_writer = csv.writer(fout_matrix, delimiter="\t", lineterminator="\n")
            stats_writer = csv.writer(fout_stats, delimiter="\t", lineterminator="\n")

            matrix_writer.writerow(["feature_id"] + common_samples)
            stats_writer.writerow([
                "feature_id",
                "origin_sample",
                "total_abundance",
                "prevalence",
                "mean_abundance",
                "sd_abundance",
                "max_abundance",
            ])

            for row in reader:
                n_seen += 1

                if len(row) < len(header):
                    row = row + ["0"] * (len(header) - len(row))

                feature_id = row[0]
                vals = []
                raw_vals = []

                for pos, idx in enumerate(common_indices):
                    raw = row[idx] if idx < len(row) else "0"
                    raw_clean = format_number_for_matrix(raw)
                    val, ok = parse_float_or_zero(raw_clean)

                    if not ok:
                        n_nonfinite += 1

                    if val < 0:
                        n_negative += 1

                    vals.append(val)
                    raw_vals.append(raw_clean)

                    raw_col_sums[pos] += val
                    if val > 0:
                        raw_col_detected[pos] += 1

                total = sum(vals)
                prevalence = sum(1 for v in vals if v > 0)
                mean = total / n_common

                if n_common > 1:
                    sumsq = sum(v * v for v in vals)
                    variance = (sumsq - (total * total / n_common)) / (n_common - 1)
                    if variance < 0 and variance > -1e-8:
                        variance = 0.0
                    sd = math.sqrt(variance) if variance >= 0 else 0.0
                else:
                    sd = 0.0

                max_val = max(vals) if vals else 0.0

                keep = (
                    prevalence >= args.min_prevalence
                    and total >= args.min_total
                    and math.isfinite(sd)
                    and sd >= args.min_sd
                )

                if keep:
                    n_retained += 1

                    matrix_writer.writerow([feature_id] + raw_vals)

                    origin = detect_origin_sample(feature_id)

                    stats_writer.writerow([
                        feature_id,
                        origin,
                        f"{total:.10g}",
                        prevalence,
                        f"{mean:.10g}",
                        f"{sd:.10g}",
                        f"{max_val:.10g}",
                    ])

                    for pos, v in enumerate(vals):
                        filt_col_sums[pos] += v
                        if v > 0:
                            filt_col_detected[pos] += 1

                    prev_counter[prevalence] += 1
                    total_bin_counter[total_bin(total)] += 1
                    origin_counter[origin] += 1

                    heap_counter += 1
                    heap_item = (total, prevalence, heap_counter, feature_id, origin, mean, sd, max_val)

                    if len(top_heap) < args.top_n:
                        heapq.heappush(top_heap, heap_item)
                    else:
                        if heap_item[:2] > top_heap[0][:2]:
                            heapq.heapreplace(top_heap, heap_item)
                else:
                    n_dropped += 1

                if n_seen % args.progress_every == 0:
                    log(
                        f"Processed: {n_seen:,} | retained: {n_retained:,} | dropped: {n_dropped:,}",
                        log_handle,
                    )

        log("Finished streaming.", log_handle)
        log(f"Features seen: {n_seen:,}", log_handle)
        log(f"Features retained: {n_retained:,}", log_handle)
        log(f"Features dropped: {n_dropped:,}", log_handle)

        # ----------------------------------------------------------------------
        # Tables
        # ----------------------------------------------------------------------
        sample_qc_raw_rows = []
        sample_qc_filtered_rows = []

        for i, s in enumerate(common_samples):
            sample_qc_raw_rows.append([
                s,
                f"{raw_col_sums[i]:.10g}",
                raw_col_detected[i],
            ])
            sample_qc_filtered_rows.append([
                s,
                f"{filt_col_sums[i]:.10g}",
                filt_col_detected[i],
                f"{(filt_col_detected[i] / n_retained) if n_retained else 0:.10g}",
            ])

        write_tsv(
            tables_dir / "502_sample_qc_raw_common_samples.tsv",
            [
                "sample_id",
                "raw_total_abundance_common_samples",
                "raw_detected_features_common_samples",
            ],
            sample_qc_raw_rows,
        )

        write_tsv(
            tables_dir / "502_sample_qc_filtered_matrix.tsv",
            [
                "sample_id",
                "filtered_total_abundance",
                "filtered_detected_features",
                "filtered_detected_fraction",
            ],
            sample_qc_filtered_rows,
        )

        prev_rows = []
        for k in sorted(prev_counter.keys()):
            n = prev_counter[k]
            prev_rows.append([k, n, f"{n / n_retained if n_retained else 0:.10g}"])

        write_tsv(
            tables_dir / "502_retained_prevalence_distribution.tsv",
            ["prevalence", "n_features", "fraction_features"],
            prev_rows,
        )

        bin_order = [
            "0",
            "1",
            "2-5",
            "6-10",
            "11-20",
            "21-50",
            "51-100",
            "101-500",
            "501-1000",
            "1001-10000",
            ">10000",
        ]

        bin_rows = []
        for b in bin_order:
            n = total_bin_counter.get(b, 0)
            if n > 0:
                bin_rows.append([b, n, f"{n / n_retained if n_retained else 0:.10g}"])

        write_tsv(
            tables_dir / "502_retained_total_abundance_distribution.tsv",
            ["total_abundance_bin", "n_features", "fraction_features"],
            bin_rows,
        )

        origin_rows = []
        for origin, n in origin_counter.most_common():
            origin_rows.append([origin, n])

        write_tsv(
            tables_dir / "502_retained_origin_sample_counts.tsv",
            ["origin_sample", "n_retained_features"],
            origin_rows,
        )

        top_items = sorted(top_heap, key=lambda x: (x[0], x[1]), reverse=True)

        top_rows = []
        for total, prevalence, _, feature_id, origin, mean, sd, max_val in top_items:
            top_rows.append([
                feature_id,
                origin,
                f"{total:.10g}",
                prevalence,
                f"{mean:.10g}",
                f"{sd:.10g}",
                f"{max_val:.10g}",
            ])

        write_tsv(
            tables_dir / "502_top_retained_features_by_total.tsv",
            [
                "feature_id",
                "origin_sample",
                "total_abundance",
                "prevalence",
                "mean_abundance",
                "sd_abundance",
                "max_abundance",
            ],
            top_rows,
        )

        filter_summary_rows = [
            ["script", "502_construir_matriz_proteica_filtrada.py"],
            ["timestamp", str(datetime.now())],
            ["input_tsv", str(input_tsv.resolve())],
            ["input_size_bytes", str(input_tsv.stat().st_size)],
            ["meta_csv", str(meta_csv.resolve())],
            ["out_dir", str(out_dir)],
            ["query_col", query_col],
            ["chosen_id_col", chosen_id_col],
            ["n_matrix_samples", len(matrix_samples)],
            ["n_metadata_rows", len(meta_rows)],
            ["n_common_samples", len(common_samples)],
            ["common_samples", ";".join(common_samples)],
            ["excluded_matrix_samples", ";".join(excluded_matrix_samples)],
            ["min_prevalence", args.min_prevalence],
            ["min_total", args.min_total],
            ["min_sd", args.min_sd],
            ["progress_every", args.progress_every],
            ["n_features_seen", n_seen],
            ["n_features_retained", n_retained],
            ["n_features_dropped", n_dropped],
            ["retained_fraction", f"{n_retained / n_seen if n_seen else 0:.10g}"],
            ["n_negative_values_common", n_negative],
            ["n_nonfinite_values_common", n_nonfinite],
            ["filtered_matrix_file", str(filtered_matrix_path)],
            ["retained_stats_file", str(retained_stats_path)],
        ]

        write_tsv(
            tables_dir / "502_filter_summary.tsv",
            ["key", "value"],
            filter_summary_rows,
        )

        # Latest
        latest_file = out_root / "LATEST_502_construir_matriz_proteica_filtrada.txt"
        with open(latest_file, "w") as fh:
            fh.write(str(out_dir) + "\n")

        log("===== DONE 502 =====", log_handle)
        log(f"Output: {out_dir}", log_handle)
        log(f"Latest: {latest_file}", log_handle)
        log("Main outputs:", log_handle)
        log(f"  {filtered_matrix_path}", log_handle)
        log(f"  {retained_stats_path}", log_handle)
        log(f"  {tables_dir / ('502_metadata_aligned_' + str(len(common_samples)) + 'samples.csv')}", log_handle)
        log(f"  {tables_dir / '502_filter_summary.tsv'}", log_handle)


if __name__ == "__main__":
    main()
