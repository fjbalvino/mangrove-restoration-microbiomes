#!/usr/bin/env python3
# ============================================================
# 02_03_aggregate_KO_PFAM.py
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Maintained successor of historical 082A. The manuscript consumed archived 082A matrices; a 504 rerun is not established by the capsule.
# Inputs (source expressions; complete list in docs/contracts/02_03_aggregate_KO_PFAM.json):
#   run_502 = Path(latest_502.read_text(encoding="utf-8").strip())
# Outputs (source expressions; complete list in contract):
#   temporary.write_text(text, encoding="utf-8")
#   def write_tsv(
#   write_tsv(
# Algorithmic provenance:
# Stream annotated gene abundances and aggregate separate KO and PFAM layers.
#   Cantalapiedra et al. (2021), eggNOG-mapper v2, doi:10.1093/molbev/msab293; this is annotation provenance, not evidence that this script runs eggNOG.
# Source SHA-256: 09773c8800b7da722b6840ca88f8f73fac8a0ce1490ba0db9582ed3de6e0a181
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================


"""
504_construir_matrices_completas_KO_PFAM.py

Reconstruct complete KEGG KO and PFAM abundance matrices by joining, in
lockstep by #query, the canonical bacteria-only/no-MZ protein abundance table
to its complete annotation table.

Scientific contract
-------------------
1. The analytical universe is fixed by the 51 samples in the 502 metadata.
2. KEGG_ko is the primary functional layer; PFAMs is an independent sensitivity
   layer. The two annotation types are never combined in one matrix.
3. Protein abundance and annotation rows must match exactly, in the same order,
   by #query. Any mismatch stops the run before a potentially invalid matrix is
   accepted.
4. Aggregation occurs before any protein-level prevalence/abundance filtering.
5. If a protein has multiple unique annotations within one layer, its abundance
   is divided equally among those annotations. This conserves annotated
   abundance and prevents multi-annotation inflation.
6. No HI, MHI_local, environmental axis, restoration state, locality, depth, or
   network result is used to select or rank functions.
7. This script creates complete matrices and audits their suitability for the
   downstream 801 contract. It does not select N=100/200/300
   functions, transform by CLR, infer networks, or test associations.

The implementation uses only the Python standard library and processes both
tables simultaneously without loading either protein-level table into memory.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import os
import re
import shutil
import statistics
import sys
import time
from array import array
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


SCRIPT_ID = "504"
SCRIPT_BASENAME = "construir_matrices_completas_KO_PFAM"
SCRIPT_NAME = f"{SCRIPT_ID}_{SCRIPT_BASENAME}"

DEFAULT_INPUT = (
    "/home/fjbalvino/Tipping_points/"
    "FINAL_ANNOTATED_ABUNDANCE_noMZ_bacteria_fran.tsv"
)
DEFAULT_ANNOTATIONS = (
    "/data/Ciencia-Frontera/Results/04-assemblies/"
    "assemblies-annotations/final_tables/"
    "FINAL_ANNOTATED_ABUNDANCE_noMZ_bacteria.tsv.gz"
)
DEFAULT_METADATA = ""
DEFAULT_OUT_ROOT = "/home/fjbalvino/Tipping_points/resultados_finales"

MISSING_ANNOTATIONS = {
    "",
    "-",
    ".",
    "na",
    "nan",
    "none",
    "null",
    "unannotated",
}
MISSING_NUMERIC = {"", "na", "nan", "none", "null", "."}
KO_PATTERN = re.compile(r"(?<![A-Za-z0-9])K[0-9]{5}(?![0-9])", re.IGNORECASE)
PFAM_SPLIT_PATTERN = re.compile(r"[,;]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build complete KEGG KO and PFAM abundance matrices by joining "
            "canonical protein abundances to complete annotations by #query."
        )
    )
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument("--annotations", default=DEFAULT_ANNOTATIONS)
    parser.add_argument("--protein_id_column", default="#query")
    parser.add_argument("--metadata", default=DEFAULT_METADATA)
    parser.add_argument("--historical_ko", default="")
    parser.add_argument("--out_root", default=DEFAULT_OUT_ROOT)
    parser.add_argument("--expected_rows", type=int, default=26_195_984)
    parser.add_argument("--expected_n", type=int, default=51)
    parser.add_argument("--expected_source_samples", type=int, default=60)
    parser.add_argument("--expected_lat_blocks", type=int, default=17)
    parser.add_argument("--expected_profile_size", type=int, default=3)
    parser.add_argument("--prevalence_threshold", type=float, default=0.50)
    parser.add_argument("--taxa_caps", default="100,200,300")
    parser.add_argument("--progress_every", type=int, default=100_000)
    parser.add_argument("--conservation_tolerance", type=float, default=1e-10)
    parser.add_argument("--gzip_level", type=int, default=6)
    parser.add_argument("--strict", choices=("0", "1"), default="1")
    return parser.parse_args()


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def timestamp_for_path() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def open_text(path: Path, mode: str = "rt"):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def format_number(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        if math.isinf(value):
            return "Inf" if value > 0 else "-Inf"
        return f"{value:.12g}"
    return str(value)


def write_tsv(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, object]],
    gzip_level: int = 6,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".gz":
        handle = gzip.open(
            path,
            "wt",
            encoding="utf-8",
            newline="",
            compresslevel=gzip_level,
        )
    else:
        handle = path.open("w", encoding="utf-8", newline="")
    with handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {field: format_number(row.get(field)) for field in fieldnames}
            )


def write_metadata_copy(
    path: Path,
    fieldnames: Sequence[str],
    rows: Sequence[Mapping[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_caps(raw: str) -> List[int]:
    values: List[int] = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        value = int(item)
        if value <= 0:
            raise ValueError("All taxa_caps values must be positive integers.")
        values.append(value)
    values = sorted(set(values))
    if not values:
        raise ValueError("taxa_caps cannot be empty.")
    return values


def resolve_metadata_id(fieldnames: Sequence[str]) -> str:
    for candidate in ("sample_id", ".sample_id"):
        if candidate in fieldnames:
            return candidate
    raise ValueError(
        "Metadata must contain the canonical sample identifier column "
        "'sample_id' or '.sample_id'."
    )


def read_metadata(
    path: Path,
) -> Tuple[List[Dict[str, str]], List[str], str, List[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("Metadata has no header.")
        fieldnames = [str(x).strip() for x in reader.fieldnames]
        if len(fieldnames) != len(set(fieldnames)):
            raise ValueError("Metadata contains duplicated column names.")
        id_column = resolve_metadata_id(fieldnames)
        required = {
            id_column,
            "restoration4",
            "locality",
            "lat_block",
            "depth_cm",
        }
        missing = sorted(required.difference(fieldnames))
        if missing:
            raise ValueError(
                "Metadata lacks required design columns: " + ", ".join(missing)
            )
        rows: List[Dict[str, str]] = []
        for raw_row in reader:
            row = {
                key: ("" if raw_row.get(key) is None else raw_row[key].strip())
                for key in fieldnames
            }
            rows.append(row)
    sample_ids = [row[id_column] for row in rows]
    if any(not sample_id for sample_id in sample_ids):
        raise ValueError("Metadata contains empty sample identifiers.")
    if len(sample_ids) != len(set(sample_ids)):
        duplicated = sorted(
            sample_id
            for sample_id, count in Counter(sample_ids).items()
            if count > 1
        )
        raise ValueError(
            "Metadata contains duplicated sample identifiers: "
            + ", ".join(duplicated)
        )
    return rows, fieldnames, id_column, sample_ids


def parse_ko_tokens(raw: str) -> Tuple[str, ...]:
    if raw.strip().lower() in MISSING_ANNOTATIONS:
        return ()
    tokens = {match.upper() for match in KO_PATTERN.findall(raw)}
    return tuple(sorted(tokens))


def parse_pfam_tokens(raw: str) -> Tuple[str, ...]:
    if raw.strip().lower() in MISSING_ANNOTATIONS:
        return ()
    tokens = set()
    for token in PFAM_SPLIT_PATTERN.split(raw):
        clean = " ".join(token.strip().split())
        if clean.lower() in MISSING_ANNOTATIONS:
            continue
        if clean.lower().startswith("pfam:"):
            clean = clean.split(":", 1)[1].strip()
        if clean:
            tokens.add(clean)
    return tuple(sorted(tokens, key=lambda item: (item.casefold(), item)))


class LayerAccumulator:
    def __init__(self, layer: str, n_samples: int) -> None:
        self.layer = layer
        self.n_samples = n_samples
        self.values: Dict[str, array] = {}
        self.rows_with_annotation = 0
        self.rows_with_annotation_and_abundance = 0
        self.rows_with_multiple_annotations = 0
        self.annotation_assignments = 0
        self.annotated_abundance_before_fractionation = 0.0

    def add(
        self,
        tokens: Sequence[str],
        nonzero_values: Sequence[Tuple[int, float]],
        row_total: float,
    ) -> None:
        if not tokens:
            return
        self.rows_with_annotation += 1
        self.annotation_assignments += len(tokens)
        if len(tokens) > 1:
            self.rows_with_multiple_annotations += 1
        if not nonzero_values:
            return
        self.rows_with_annotation_and_abundance += 1
        self.annotated_abundance_before_fractionation += row_total
        divisor = float(len(tokens))
        for token in tokens:
            target = self.values.get(token)
            if target is None:
                target = array("d", [0.0]) * self.n_samples
                self.values[token] = target
            for sample_index, abundance in nonzero_values:
                target[sample_index] += abundance / divisor


def read_source_header(path: Path, table_label: str) -> List[str]:
    with open_text(path, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t", quotechar='"')
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"{table_label} table is empty.") from exc
    header = [field.strip().lstrip("\ufeff") for field in header]
    if len(header) != len(set(header)):
        duplicates = sorted(
            field for field, count in Counter(header).items() if count > 1
        )
        raise ValueError(
            f"{table_label} table contains duplicated column names: "
            + ", ".join(duplicates)
        )
    return header


def stream_aggregate(
    abundance_path: Path,
    annotation_path: Path,
    abundance_header: Sequence[str],
    annotation_header: Sequence[str],
    sample_ids: Sequence[str],
    protein_id_column: str,
    progress_every: int,
) -> Tuple[LayerAccumulator, LayerAccumulator, Dict[str, object]]:
    abundance_id_index = abundance_header.index(protein_id_column)
    annotation_id_index = annotation_header.index(protein_id_column)
    ko_index = annotation_header.index("KEGG_ko")
    pfam_index = annotation_header.index("PFAMs")
    sample_indices = [
        abundance_header.index(sample_id) for sample_id in sample_ids
    ]
    n_abundance_columns = len(abundance_header)
    n_annotation_columns = len(annotation_header)

    ko = LayerAccumulator("KEGG_ko", len(sample_ids))
    pfam = LayerAccumulator("PFAM", len(sample_ids))

    counters: Dict[str, object] = {
        "rows_seen": 0,
        "abundance_rows_seen": 0,
        "annotation_rows_seen": 0,
        "protein_id_matches": 0,
        "empty_protein_ids": 0,
        "protein_id_mismatches": 0,
        "malformed_abundance_rows": 0,
        "malformed_annotation_rows": 0,
        "rows_without_KO_or_PFAM": 0,
        "rows_with_selected_sample_abundance": 0,
        "numeric_conversion_failures": 0,
        "nonfinite_values": 0,
        "negative_values": 0,
        "selected_abundance_all_rows": 0.0,
    }

    started = time.monotonic()
    with (
        open_text(abundance_path, "rt") as abundance_handle,
        open_text(annotation_path, "rt") as annotation_handle,
    ):
        abundance_reader = csv.reader(
            abundance_handle,
            delimiter="\t",
            quotechar='"',
        )
        annotation_reader = csv.reader(
            annotation_handle,
            delimiter="\t",
            quotechar='"',
        )
        next(abundance_reader)
        next(annotation_reader)

        paired_rows = itertools.zip_longest(
            abundance_reader,
            annotation_reader,
            fillvalue=None,
        )
        for line_number, pair in enumerate(paired_rows, start=2):
            abundance_row, annotation_row = pair

            if abundance_row is None:
                counters["annotation_rows_seen"] = (
                    int(counters["annotation_rows_seen"]) + 1
                )
                raise ValueError(
                    "Annotation table has more protein rows than the abundance "
                    f"table; first extra annotation row is line {line_number}."
                )
            if annotation_row is None:
                counters["abundance_rows_seen"] = (
                    int(counters["abundance_rows_seen"]) + 1
                )
                raise ValueError(
                    "Abundance table has more protein rows than the annotation "
                    f"table; first extra abundance row is line {line_number}."
                )

            counters["rows_seen"] = int(counters["rows_seen"]) + 1
            counters["abundance_rows_seen"] = (
                int(counters["abundance_rows_seen"]) + 1
            )
            counters["annotation_rows_seen"] = (
                int(counters["annotation_rows_seen"]) + 1
            )
            row_number = int(counters["rows_seen"])

            if len(abundance_row) != n_abundance_columns:
                counters["malformed_abundance_rows"] = (
                    int(counters["malformed_abundance_rows"]) + 1
                )
                raise ValueError(
                    "Malformed abundance row at line "
                    f"{line_number}: observed {len(abundance_row)} columns; "
                    f"expected {n_abundance_columns}."
                )
            if len(annotation_row) != n_annotation_columns:
                counters["malformed_annotation_rows"] = (
                    int(counters["malformed_annotation_rows"]) + 1
                )
                raise ValueError(
                    "Malformed annotation row at line "
                    f"{line_number}: observed {len(annotation_row)} columns; "
                    f"expected {n_annotation_columns}."
                )

            abundance_id = abundance_row[abundance_id_index].strip()
            annotation_id = annotation_row[annotation_id_index].strip()
            if not abundance_id or not annotation_id:
                counters["empty_protein_ids"] = (
                    int(counters["empty_protein_ids"]) + 1
                )
                raise ValueError(
                    "Empty protein identifier at line "
                    f"{line_number}: abundance={abundance_id!r}; "
                    f"annotation={annotation_id!r}."
                )
            if abundance_id != annotation_id:
                counters["protein_id_mismatches"] = (
                    int(counters["protein_id_mismatches"]) + 1
                )
                raise ValueError(
                    "Protein ID mismatch between abundance and annotation "
                    f"tables at line {line_number}: "
                    f"abundance={abundance_id!r}; "
                    f"annotation={annotation_id!r}. "
                    "The join is intentionally strict and requires identical "
                    "row order."
                )
            counters["protein_id_matches"] = (
                int(counters["protein_id_matches"]) + 1
            )

            ko_tokens = parse_ko_tokens(annotation_row[ko_index])
            pfam_tokens = parse_pfam_tokens(annotation_row[pfam_index])
            if not ko_tokens and not pfam_tokens:
                counters["rows_without_KO_or_PFAM"] = (
                    int(counters["rows_without_KO_or_PFAM"]) + 1
                )
                if progress_every > 0 and row_number % progress_every == 0:
                    elapsed = max(time.monotonic() - started, 1e-9)
                    print(
                        f"  processed {row_number:,} protein rows "
                        f"({row_number / elapsed:,.0f} rows/s)",
                        flush=True,
                    )
                continue

            nonzero_values: List[Tuple[int, float]] = []
            row_total = 0.0
            for output_index, source_index in enumerate(sample_indices):
                raw_value = abundance_row[source_index].strip()
                if raw_value.lower() in MISSING_NUMERIC:
                    if raw_value:
                        counters["numeric_conversion_failures"] = (
                            int(counters["numeric_conversion_failures"]) + 1
                        )
                    continue
                if raw_value in {"0", "0.0", "0.00", "0.000"}:
                    continue
                try:
                    value = float(raw_value)
                except ValueError:
                    counters["numeric_conversion_failures"] = (
                        int(counters["numeric_conversion_failures"]) + 1
                    )
                    continue
                if not math.isfinite(value):
                    counters["nonfinite_values"] = (
                        int(counters["nonfinite_values"]) + 1
                    )
                    continue
                if value < 0:
                    counters["negative_values"] = (
                        int(counters["negative_values"]) + 1
                    )
                    continue
                if value == 0:
                    continue
                nonzero_values.append((output_index, value))
                row_total += value

            if nonzero_values:
                counters["rows_with_selected_sample_abundance"] = (
                    int(counters["rows_with_selected_sample_abundance"]) + 1
                )
                counters["selected_abundance_all_rows"] = (
                    float(counters["selected_abundance_all_rows"]) + row_total
                )

            ko.add(ko_tokens, nonzero_values, row_total)
            pfam.add(pfam_tokens, nonzero_values, row_total)

            if progress_every > 0 and row_number % progress_every == 0:
                elapsed = max(time.monotonic() - started, 1e-9)
                print(
                    f"  processed {row_number:,} protein rows "
                    f"({row_number / elapsed:,.0f} rows/s)",
                    flush=True,
                )

    counters["elapsed_seconds_streaming"] = time.monotonic() - started
    return ko, pfam, counters


def sorted_features(layer: str, values: Mapping[str, array]) -> List[str]:
    if layer == "KEGG_ko":
        return sorted(values)
    return sorted(values, key=lambda item: (item.casefold(), item))


def compute_layer_qc(
    accumulator: LayerAccumulator,
    sample_ids: Sequence[str],
    min_prevalence_n: int,
) -> Tuple[List[str], List[Dict[str, object]], List[Dict[str, object]], Dict[str, object]]:
    features = sorted_features(accumulator.layer, accumulator.values)
    n_samples = len(sample_ids)
    sample_totals = [0.0] * n_samples
    for feature in features:
        values = accumulator.values[feature]
        for sample_index, value in enumerate(values):
            sample_totals[sample_index] += value

    feature_rows: List[Dict[str, object]] = []
    for feature in features:
        values = accumulator.values[feature]
        prevalence_n = sum(value > 0 for value in values)
        total_abundance = math.fsum(values)
        relative_values = [
            values[index] / sample_totals[index]
            for index in range(n_samples)
            if sample_totals[index] > 0
        ]
        feature_rows.append(
            {
                "annotation_layer": accumulator.layer,
                "feature_id": feature,
                "prevalence_n": prevalence_n,
                "prevalence_fraction": prevalence_n / n_samples,
                "total_abundance": total_abundance,
                "mean_relative_abundance": (
                    statistics.fmean(relative_values) if relative_values else 0.0
                ),
                "eligible_prevalence_ge_threshold": (
                    prevalence_n >= min_prevalence_n
                ),
            }
        )

    sample_rows = [
        {
            "annotation_layer": accumulator.layer,
            "sample_id": sample_id,
            "total_abundance": sample_totals[index],
            "n_nonzero_features": sum(
                accumulator.values[feature][index] > 0 for feature in features
            ),
        }
        for index, sample_id in enumerate(sample_ids)
    ]

    matrix_total = math.fsum(sample_totals)
    annotated_total = accumulator.annotated_abundance_before_fractionation
    absolute_difference = abs(matrix_total - annotated_total)
    relative_error = absolute_difference / max(abs(annotated_total), 1.0)
    eligible_n = sum(
        row["eligible_prevalence_ge_threshold"] for row in feature_rows
    )
    summary = {
        "annotation_layer": accumulator.layer,
        "n_samples": n_samples,
        "n_complete_features": len(features),
        "n_features_prevalence_ge_threshold": eligible_n,
        "min_prevalence_n": min_prevalence_n,
        "zero_fraction": (
            sum(
                value == 0
                for feature in features
                for value in accumulator.values[feature]
            )
            / (len(features) * n_samples)
            if features
            else math.nan
        ),
        "matrix_total_abundance": matrix_total,
        "annotated_abundance_before_fractionation": annotated_total,
        "conservation_absolute_difference": absolute_difference,
        "conservation_relative_error": relative_error,
        "rows_with_annotation": accumulator.rows_with_annotation,
        "rows_with_annotation_and_abundance": (
            accumulator.rows_with_annotation_and_abundance
        ),
        "rows_with_multiple_annotations": (
            accumulator.rows_with_multiple_annotations
        ),
        "annotation_assignments": accumulator.annotation_assignments,
        "median_sample_total": (
            statistics.median(sample_totals) if sample_totals else math.nan
        ),
        "minimum_sample_total": min(sample_totals) if sample_totals else math.nan,
        "maximum_sample_total": max(sample_totals) if sample_totals else math.nan,
    }
    return features, feature_rows, sample_rows, summary


def write_matrix_samples_x_features(
    path: Path,
    accumulator: LayerAccumulator,
    sample_ids: Sequence[str],
    features: Sequence[str],
    gzip_level: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(
        path,
        "wt",
        encoding="utf-8",
        newline="",
        compresslevel=gzip_level,
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", *features])
        for sample_index, sample_id in enumerate(sample_ids):
            writer.writerow(
                [
                    sample_id,
                    *(
                        format_number(accumulator.values[feature][sample_index])
                        for feature in features
                    ),
                ]
            )


def mean_rank(values: Sequence[float]) -> List[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + 1 + end) / 2.0
        for position in range(start, end):
            ranks[order[position]] = rank
        start = end
    return ranks


def pearson(values_x: Sequence[float], values_y: Sequence[float]) -> float:
    if len(values_x) != len(values_y) or len(values_x) < 2:
        return math.nan
    mean_x = statistics.fmean(values_x)
    mean_y = statistics.fmean(values_y)
    numerator = math.fsum(
        (x - mean_x) * (y - mean_y) for x, y in zip(values_x, values_y)
    )
    denominator_x = math.fsum((x - mean_x) ** 2 for x in values_x)
    denominator_y = math.fsum((y - mean_y) ** 2 for y in values_y)
    denominator = math.sqrt(denominator_x * denominator_y)
    return numerator / denominator if denominator > 0 else math.nan


def spearman(values_x: Sequence[float], values_y: Sequence[float]) -> float:
    if len(values_x) != len(values_y) or len(values_x) < 2:
        return math.nan
    return pearson(mean_rank(values_x), mean_rank(values_y))


def read_historical_ko(
    path: Path,
    sample_ids: Sequence[str],
) -> Dict[str, array]:
    with open_text(path, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t", quotechar='"')
        try:
            header = [item.strip().lstrip("\ufeff") for item in next(reader)]
        except StopIteration as exc:
            raise ValueError("Historical KO matrix is empty.") from exc

        if set(sample_ids).issubset(header[1:]):
            sample_positions = [header.index(sample_id) for sample_id in sample_ids]
            matrix: Dict[str, array] = {}
            for row in reader:
                if len(row) != len(header):
                    raise ValueError("Historical KO matrix has malformed rows.")
                feature = row[0].strip()
                if not feature:
                    continue
                values = array("d")
                for position in sample_positions:
                    value = float(row[position])
                    if not math.isfinite(value) or value < 0:
                        raise ValueError(
                            "Historical KO matrix contains invalid abundances."
                        )
                    values.append(value)
                matrix[feature] = values
            return matrix

        if header[0] in {"sample_id", ".sample_id"}:
            feature_names = header[1:]
            sample_rows: Dict[str, List[float]] = {}
            for row in reader:
                if len(row) != len(header):
                    raise ValueError("Historical KO matrix has malformed rows.")
                sample_id = row[0].strip()
                if sample_id in sample_ids:
                    values = [float(item) for item in row[1:]]
                    if any(
                        not math.isfinite(value) or value < 0 for value in values
                    ):
                        raise ValueError(
                            "Historical KO matrix contains invalid abundances."
                        )
                    sample_rows[sample_id] = values
            missing = sorted(set(sample_ids).difference(sample_rows))
            if missing:
                raise ValueError(
                    "Historical KO matrix lacks samples: " + ", ".join(missing)
                )
            matrix = {}
            for feature_index, feature in enumerate(feature_names):
                matrix[feature] = array(
                    "d",
                    [
                        sample_rows[sample_id][feature_index]
                        for sample_id in sample_ids
                    ],
                )
            return matrix

    raise ValueError(
        "Could not determine the orientation of the historical KO matrix."
    )


def compare_historical_ko(
    reconstructed: Mapping[str, array],
    historical: Mapping[str, array],
) -> Dict[str, object]:
    reconstructed_features = set(reconstructed)
    historical_features = set(historical)
    common = sorted(reconstructed_features.intersection(historical_features))
    union = reconstructed_features.union(historical_features)

    reconstructed_totals = {
        feature: math.fsum(values) for feature, values in reconstructed.items()
    }
    historical_totals = {
        feature: math.fsum(values) for feature, values in historical.items()
    }
    common_reconstructed = [
        math.log1p(reconstructed_totals[feature]) for feature in common
    ]
    common_historical = [
        math.log1p(historical_totals[feature]) for feature in common
    ]

    n_samples = len(next(iter(reconstructed.values()))) if reconstructed else 0
    reconstructed_sample_totals = [
        math.fsum(values[index] for values in reconstructed.values())
        for index in range(n_samples)
    ]
    historical_sample_totals = [
        math.fsum(values[index] for values in historical.values())
        for index in range(n_samples)
    ]

    return {
        "comparison_status": "COMPLETED_DESCRIPTIVE_ONLY",
        "n_reconstructed_KOs": len(reconstructed_features),
        "n_historical_KOs": len(historical_features),
        "n_common_KOs": len(common),
        "n_only_reconstructed_KOs": len(
            reconstructed_features.difference(historical_features)
        ),
        "n_only_historical_KOs": len(
            historical_features.difference(reconstructed_features)
        ),
        "feature_jaccard": len(common) / len(union) if union else math.nan,
        "common_feature_log1p_total_pearson": pearson(
            common_reconstructed, common_historical
        ),
        "common_feature_log1p_total_spearman": spearman(
            common_reconstructed, common_historical
        ),
        "sample_total_pearson": pearson(
            reconstructed_sample_totals, historical_sample_totals
        ),
        "sample_total_spearman": spearman(
            reconstructed_sample_totals, historical_sample_totals
        ),
        "reconstructed_total_abundance": math.fsum(
            reconstructed_sample_totals
        ),
        "historical_total_abundance": math.fsum(historical_sample_totals),
        "historical_is_acceptance_criterion": False,
        "interpretation": (
            "Descriptive control only: the historical matrix was built from a "
            "different protein universe and may also use a different "
            "multi-annotation allocation rule."
        ),
    }


def build_check(
    check: str,
    observed: object,
    expected: object,
    passed: bool,
    critical: bool = True,
) -> Dict[str, object]:
    return {
        "check": check,
        "observed": observed,
        "expected": expected,
        "pass": passed,
        "critical": critical,
    }


def main() -> int:
    args = parse_args()
    if (
        args.expected_n != 51
        or args.expected_lat_blocks != 17
        or args.expected_profile_size != 3
    ):
        raise RuntimeError(
            "The analytical canon is frozen at 51 samples, 17 profiles, and 3 depths."
        )
    if (
        args.expected_n != 51
        or args.expected_lat_blocks != 17
        or args.expected_profile_size != 3
    ):
        raise ValueError(
            "The analytical canon is frozen at 51 samples, 17 profiles, and 3 depths."
        )
    abundance_path = Path(args.input).expanduser().resolve()
    annotation_path = Path(args.annotations).expanduser().resolve()
    out_root = Path(args.out_root).expanduser().resolve()
    if args.metadata.strip():
        metadata_path = Path(args.metadata).expanduser().resolve()
    else:
        latest_502 = out_root / "LATEST_502_construir_matriz_proteica_filtrada.txt"
        if not latest_502.is_file():
            raise FileNotFoundError(
                "Missing --metadata and latest pointer: " + str(latest_502)
            )
        run_502 = Path(latest_502.read_text(encoding="utf-8").strip())
        candidates = sorted(
            (run_502 / "tables").glob("502_metadata_aligned_*samples.csv")
        )
        if len(candidates) != 1:
            raise RuntimeError(
                "Expected exactly one 502 aligned metadata file; found "
                + str(len(candidates))
            )
        metadata_path = candidates[0].resolve()
    historical_path = (
        Path(args.historical_ko).expanduser().resolve()
        if args.historical_ko.strip()
        else None
    )
    caps = parse_caps(args.taxa_caps)

    if not (0 < args.prevalence_threshold <= 1):
        raise ValueError("prevalence_threshold must be in (0, 1].")
    if args.progress_every < 0:
        raise ValueError("progress_every must be nonnegative.")
    if not abundance_path.is_file():
        raise FileNotFoundError(
            f"Protein abundance table not found: {abundance_path}"
        )
    if not annotation_path.is_file():
        raise FileNotFoundError(
            f"Protein annotation table not found: {annotation_path}"
        )
    if not metadata_path.is_file():
        raise FileNotFoundError(f"Metadata not found: {metadata_path}")

    stamp = timestamp_for_path()
    out_dir = out_root / f"{SCRIPT_NAME}_{stamp}"
    tables_dir = out_dir / "tables"
    provenance_dir = out_dir / "provenance"
    out_dir.mkdir(parents=True, exist_ok=False)
    tables_dir.mkdir(parents=True)
    provenance_dir.mkdir(parents=True)

    latest_path = out_root / f"LATEST_{SCRIPT_NAME}.txt"
    atomic_write_text(latest_path, str(out_dir) + "\n")

    start_time = utc_now_iso()
    wall_start = time.monotonic()
    print("=" * 72)
    print(SCRIPT_NAME)
    print(f"Started: {start_time}")
    print(f"Output: {out_dir}")
    print("=" * 72)

    print("[1/8] Loading and validating the frozen 51-sample metadata...", flush=True)
    metadata_rows, metadata_fields, metadata_id, sample_ids = read_metadata(
        metadata_path
    )
    profile_counts = Counter(row["lat_block"] for row in metadata_rows)
    n_lat_blocks = len(profile_counts)
    profile_sizes = sorted(profile_counts.values())

    print("[2/8] Auditing both protein-table schemas...", flush=True)
    abundance_header = read_source_header(abundance_path, "Abundance")
    annotation_header = read_source_header(annotation_path, "Annotation")
    if args.protein_id_column not in abundance_header:
        raise ValueError(
            "Abundance table lacks protein ID column: "
            f"{args.protein_id_column}"
        )
    if args.protein_id_column not in annotation_header:
        raise ValueError(
            "Annotation table lacks protein ID column: "
            f"{args.protein_id_column}"
        )

    required_annotation_columns = {"KEGG_ko", "PFAMs"}
    missing_annotation_columns = sorted(
        required_annotation_columns.difference(annotation_header)
    )
    if missing_annotation_columns:
        raise ValueError(
            "Annotation table lacks required columns: "
            + ", ".join(missing_annotation_columns)
        )

    sample_id_set = set(sample_ids)
    abundance_source_sample_columns = [
        column
        for column in abundance_header
        if column != args.protein_id_column
    ]
    annotation_pfam_index = annotation_header.index("PFAMs")
    annotation_source_sample_columns = annotation_header[
        annotation_pfam_index + 1 :
    ]

    missing_abundance_samples = sorted(
        sample_id_set.difference(abundance_source_sample_columns)
    )
    missing_annotation_samples = sorted(
        sample_id_set.difference(annotation_source_sample_columns)
    )
    metadata_outside_abundance_block = sorted(
        sample_id_set.difference(abundance_source_sample_columns)
    )
    metadata_outside_annotation_block = sorted(
        sample_id_set.difference(annotation_source_sample_columns)
    )
    excluded_source_samples = [
        sample_id
        for sample_id in abundance_source_sample_columns
        if sample_id not in sample_id_set
    ]
    annotation_excluded_source_samples = [
        sample_id
        for sample_id in annotation_source_sample_columns
        if sample_id not in sample_id_set
    ]
    if missing_abundance_samples or metadata_outside_abundance_block:
        raise ValueError(
            "The abundance table does not contain all metadata samples. "
            "Missing: "
            + ",".join(missing_abundance_samples)
            + "; outside abundance block: "
            + ",".join(metadata_outside_abundance_block)
        )
    if missing_annotation_samples or metadata_outside_annotation_block:
        raise ValueError(
            "The annotation table does not contain all metadata samples in "
            "its post-PFAM abundance block. Missing: "
            + ",".join(missing_annotation_samples)
            + "; outside annotation abundance block: "
            + ",".join(metadata_outside_annotation_block)
        )

    schema_rows: List[Dict[str, object]] = []
    for source_table, header in (
        ("abundance", abundance_header),
        ("annotation", annotation_header),
    ):
        for index, name in enumerate(header):
            if name == args.protein_id_column:
                role = "protein_id"
            elif source_table == "annotation" and name == "KEGG_ko":
                role = "KEGG_ko_annotation"
            elif source_table == "annotation" and name == "PFAMs":
                role = "PFAM_annotation"
            elif name in sample_id_set:
                role = (
                    "selected_sample_abundance_used"
                    if source_table == "abundance"
                    else "selected_sample_abundance_not_used"
                )
            elif (
                source_table == "abundance"
                or index > annotation_pfam_index
            ):
                role = (
                    "excluded_source_sample_abundance_used"
                    if source_table == "abundance"
                    else "excluded_source_sample_abundance_not_used"
                )
            else:
                role = "other_annotation"
            schema_rows.append(
                {
                    "source_table": source_table,
                    "column_index_1based": index + 1,
                    "column_name": name,
                    "column_role": role,
                }
            )
    write_tsv(
        tables_dir / "504_input_schema.tsv",
        [
            "source_table",
            "column_index_1based",
            "column_name",
            "column_role",
        ],
        schema_rows,
        args.gzip_level,
    )
    write_metadata_copy(
        tables_dir / "504_metadata_aligned_51samples.csv",
        metadata_fields,
        metadata_rows,
    )

    print(
        "[3/8] Joining protein rows by #query and aggregating KO/PFAM...",
        flush=True,
    )
    ko, pfam, stream_qc = stream_aggregate(
        abundance_path,
        annotation_path,
        abundance_header,
        annotation_header,
        sample_ids,
        args.protein_id_column,
        args.progress_every,
    )

    min_prevalence_n = math.ceil(
        args.prevalence_threshold * len(sample_ids) - 1e-12
    )
    print("[4/8] Computing matrix and feature-level QC...", flush=True)
    (
        ko_features,
        ko_feature_rows,
        ko_sample_rows,
        ko_summary,
    ) = compute_layer_qc(ko, sample_ids, min_prevalence_n)
    (
        pfam_features,
        pfam_feature_rows,
        pfam_sample_rows,
        pfam_summary,
    ) = compute_layer_qc(pfam, sample_ids, min_prevalence_n)
    layer_summaries = [ko_summary, pfam_summary]

    print("[5/8] Writing complete matrices and audit tables...", flush=True)
    write_matrix_samples_x_features(
        tables_dir
        / "504_KEGG_ko_abundance_samples_x_functions.tsv.gz",
        ko,
        sample_ids,
        ko_features,
        args.gzip_level,
    )
    write_matrix_samples_x_features(
        tables_dir / "504_PFAM_abundance_samples_x_functions.tsv.gz",
        pfam,
        sample_ids,
        pfam_features,
        args.gzip_level,
    )
    feature_fields = [
        "annotation_layer",
        "feature_id",
        "prevalence_n",
        "prevalence_fraction",
        "total_abundance",
        "mean_relative_abundance",
        "eligible_prevalence_ge_threshold",
    ]
    write_tsv(
        tables_dir / "504_KEGG_ko_feature_qc.tsv.gz",
        feature_fields,
        ko_feature_rows,
        args.gzip_level,
    )
    write_tsv(
        tables_dir / "504_PFAM_feature_qc.tsv.gz",
        feature_fields,
        pfam_feature_rows,
        args.gzip_level,
    )
    sample_fields = [
        "annotation_layer",
        "sample_id",
        "total_abundance",
        "n_nonzero_features",
    ]
    write_tsv(
        tables_dir / "504_sample_totals.tsv",
        sample_fields,
        [*ko_sample_rows, *pfam_sample_rows],
        args.gzip_level,
    )
    summary_fields = list(layer_summaries[0].keys())
    write_tsv(
        tables_dir / "504_functional_matrix_summary.tsv",
        summary_fields,
        layer_summaries,
        args.gzip_level,
    )

    annotation_rows = [
        {
            "annotation_layer": accumulator.layer,
            "rows_with_annotation": accumulator.rows_with_annotation,
            "rows_with_annotation_and_abundance": (
                accumulator.rows_with_annotation_and_abundance
            ),
            "rows_with_multiple_annotations": (
                accumulator.rows_with_multiple_annotations
            ),
            "annotation_assignments": accumulator.annotation_assignments,
            "multi_annotation_allocation": "equal_fraction_across_unique_tokens",
        }
        for accumulator in (ko, pfam)
    ]
    write_tsv(
        tables_dir / "504_annotation_parsing_summary.tsv",
        list(annotation_rows[0].keys()),
        annotation_rows,
        args.gzip_level,
    )

    print("[6/8] Comparing reconstructed KO with the optional historical control...", flush=True)
    historical_comparison: Dict[str, object]
    historical_warning = ""
    if historical_path is None:
        historical_comparison = {
            "comparison_status": "NOT_REQUESTED",
            "historical_path": "",
            "historical_is_acceptance_criterion": False,
        }
    elif not historical_path.is_file():
        historical_warning = f"Historical KO matrix not found: {historical_path}"
        historical_comparison = {
            "comparison_status": "NOT_FOUND",
            "historical_path": str(historical_path),
            "historical_is_acceptance_criterion": False,
            "warning": historical_warning,
        }
    else:
        try:
            historical_matrix = read_historical_ko(historical_path, sample_ids)
            historical_comparison = compare_historical_ko(
                ko.values,
                historical_matrix,
            )
            historical_comparison["historical_path"] = str(historical_path)
        except Exception as exc:
            historical_warning = (
                "Historical KO comparison could not be completed: "
                f"{type(exc).__name__}: {exc}"
            )
            historical_comparison = {
                "comparison_status": "ERROR",
                "historical_path": str(historical_path),
                "historical_is_acceptance_criterion": False,
                "warning": historical_warning,
            }
    write_tsv(
        tables_dir / "504_historical_KO_comparison.tsv",
        ["metric", "value"],
        [
            {"metric": key, "value": value}
            for key, value in historical_comparison.items()
        ],
        args.gzip_level,
    )

    print("[7/8] Evaluating acceptance criteria...", flush=True)
    max_cap = max(caps)
    checks = [
        build_check(
            "metadata_sample_count",
            len(sample_ids),
            args.expected_n,
            len(sample_ids) == args.expected_n,
        ),
        build_check(
            "metadata_unique_sample_ids",
            len(set(sample_ids)),
            len(sample_ids),
            len(set(sample_ids)) == len(sample_ids),
        ),
        build_check(
            "lat_block_count",
            n_lat_blocks,
            args.expected_lat_blocks,
            n_lat_blocks == args.expected_lat_blocks,
        ),
        build_check(
            "profile_sizes",
            ",".join(map(str, profile_sizes)),
            str(args.expected_profile_size),
            bool(profile_sizes)
            and all(size == args.expected_profile_size for size in profile_sizes),
        ),
        build_check(
            "abundance_source_sample_column_count",
            len(abundance_source_sample_columns),
            args.expected_source_samples,
            len(abundance_source_sample_columns)
            == args.expected_source_samples,
        ),
        build_check(
            "annotation_source_sample_column_count",
            len(annotation_source_sample_columns),
            args.expected_source_samples,
            len(annotation_source_sample_columns)
            == args.expected_source_samples,
        ),
        build_check(
            "all_metadata_samples_in_abundance_block",
            len(metadata_outside_abundance_block),
            0,
            not metadata_outside_abundance_block,
        ),
        build_check(
            "all_metadata_samples_in_annotation_block",
            len(metadata_outside_annotation_block),
            0,
            not metadata_outside_annotation_block,
        ),
        build_check(
            "paired_protein_rows_seen",
            stream_qc["rows_seen"],
            args.expected_rows,
            int(stream_qc["rows_seen"]) == args.expected_rows,
        ),
        build_check(
            "abundance_rows_seen",
            stream_qc["abundance_rows_seen"],
            args.expected_rows,
            int(stream_qc["abundance_rows_seen"]) == args.expected_rows,
        ),
        build_check(
            "annotation_rows_seen",
            stream_qc["annotation_rows_seen"],
            args.expected_rows,
            int(stream_qc["annotation_rows_seen"]) == args.expected_rows,
        ),
        build_check(
            "protein_IDs_matched_in_lockstep",
            stream_qc["protein_id_matches"],
            stream_qc["rows_seen"],
            int(stream_qc["protein_id_matches"])
            == int(stream_qc["rows_seen"]),
        ),
        build_check(
            "protein_ID_mismatches",
            stream_qc["protein_id_mismatches"],
            0,
            int(stream_qc["protein_id_mismatches"]) == 0,
        ),
        build_check(
            "empty_protein_IDs",
            stream_qc["empty_protein_ids"],
            0,
            int(stream_qc["empty_protein_ids"]) == 0,
        ),
        build_check(
            "malformed_abundance_rows",
            stream_qc["malformed_abundance_rows"],
            0,
            int(stream_qc["malformed_abundance_rows"]) == 0,
        ),
        build_check(
            "malformed_annotation_rows",
            stream_qc["malformed_annotation_rows"],
            0,
            int(stream_qc["malformed_annotation_rows"]) == 0,
        ),
        build_check(
            "numeric_conversion_failures",
            stream_qc["numeric_conversion_failures"],
            0,
            int(stream_qc["numeric_conversion_failures"]) == 0,
        ),
        build_check(
            "nonfinite_values",
            stream_qc["nonfinite_values"],
            0,
            int(stream_qc["nonfinite_values"]) == 0,
        ),
        build_check(
            "negative_values",
            stream_qc["negative_values"],
            0,
            int(stream_qc["negative_values"]) == 0,
        ),
        build_check(
            "KEGG_ko_abundance_conserved",
            ko_summary["conservation_relative_error"],
            f"<= {args.conservation_tolerance}",
            float(ko_summary["conservation_relative_error"])
            <= args.conservation_tolerance,
        ),
        build_check(
            "PFAM_abundance_conserved",
            pfam_summary["conservation_relative_error"],
            f"<= {args.conservation_tolerance}",
            float(pfam_summary["conservation_relative_error"])
            <= args.conservation_tolerance,
        ),
        build_check(
            "KEGG_ko_supports_max_cap_after_prevalence",
            ko_summary["n_features_prevalence_ge_threshold"],
            f">= {max_cap}",
            int(ko_summary["n_features_prevalence_ge_threshold"]) >= max_cap,
        ),
        build_check(
            "PFAM_supports_max_cap_after_prevalence",
            pfam_summary["n_features_prevalence_ge_threshold"],
            f">= {max_cap}",
            int(pfam_summary["n_features_prevalence_ge_threshold"]) >= max_cap,
        ),
        build_check(
            "outcome_independent_aggregation",
            "#query_join+sample_abundances+KEGG_ko+PFAMs_only",
            "no_HI_MHI_environment_or_network_results",
            True,
        ),
        build_check(
            "network_inference_performed",
            False,
            False,
            True,
        ),
    ]
    critical_failures = [
        check for check in checks if check["critical"] and not check["pass"]
    ]
    if critical_failures:
        overall_status = "FAIL"
    elif historical_warning:
        overall_status = "PASS_WITH_WARNINGS"
    else:
        overall_status = "PASS"

    write_tsv(
        tables_dir / "504_acceptance_checks.tsv",
        ["check", "observed", "expected", "pass", "critical"],
        checks,
        args.gzip_level,
    )

    finished_time = utc_now_iso()
    elapsed_total = time.monotonic() - wall_start
    run_summary = {
        "script": SCRIPT_NAME,
        "timestamp": stamp,
        "overall_status": overall_status,
        "input_protein_abundance_table": str(abundance_path),
        "input_protein_annotation_table": str(annotation_path),
        "protein_id_column": args.protein_id_column,
        "protein_join_mode": "strict_lockstep_exact_ID",
        "metadata": str(metadata_path),
        "historical_ko": str(historical_path) if historical_path else "",
        "n_protein_rows": stream_qc["rows_seen"],
        "n_abundance_source_sample_columns": len(
            abundance_source_sample_columns
        ),
        "n_annotation_source_sample_columns": len(
            annotation_source_sample_columns
        ),
        "n_samples": len(sample_ids),
        "n_lat_blocks": n_lat_blocks,
        "profile_sizes": ",".join(map(str, sorted(set(profile_sizes)))),
        "excluded_source_samples": ";".join(excluded_source_samples),
        "annotation_excluded_source_samples": ";".join(
            annotation_excluded_source_samples
        ),
        "primary_layer": "KEGG_ko",
        "sensitivity_layer": "PFAM",
        "multi_annotation_allocation": "equal_fraction_across_unique_tokens",
        "prevalence_threshold_for_future_selection": args.prevalence_threshold,
        "min_prevalence_n_for_future_selection": min_prevalence_n,
        "future_taxa_caps": ",".join(map(str, caps)),
        "n_complete_KEGG_ko": ko_summary["n_complete_features"],
        "n_eligible_KEGG_ko": ko_summary[
            "n_features_prevalence_ge_threshold"
        ],
        "n_complete_PFAM": pfam_summary["n_complete_features"],
        "n_eligible_PFAM": pfam_summary[
            "n_features_prevalence_ge_threshold"
        ],
        "indices_used_for_aggregation": False,
        "protein_level_filter_applied_before_aggregation": False,
        "functional_cap_applied": False,
        "transformation_applied": False,
        "network_inference_performed": False,
        "historical_comparison_status": historical_comparison.get(
            "comparison_status", ""
        ),
        "started_utc": start_time,
        "finished_utc": finished_time,
        "elapsed_seconds": elapsed_total,
    }
    write_tsv(
        tables_dir / "504_run_summary.tsv",
        list(run_summary.keys()),
        [run_summary],
        args.gzip_level,
    )

    print("[8/8] Writing report and provenance...", flush=True)
    source_script = Path(sys.argv[0]).expanduser().resolve()
    script_sha256 = ""
    if source_script.is_file():
        script_sha256 = sha256_file(source_script)
        shutil.copy2(
            source_script,
            provenance_dir / source_script.name,
        )

    provenance = {
        "script": SCRIPT_NAME,
        "script_path": str(source_script),
        "script_sha256": script_sha256,
        "command": [sys.executable, *sys.argv],
        "python_version": sys.version,
        "platform": sys.platform,
        "abundance_input": {
            "path": str(abundance_path),
            "size_bytes": abundance_path.stat().st_size,
            "mtime_utc": datetime.fromtimestamp(
                abundance_path.stat().st_mtime, tz=timezone.utc
            ).isoformat(timespec="seconds"),
            "expected_rows": args.expected_rows,
        },
        "annotation_input": {
            "path": str(annotation_path),
            "size_bytes": annotation_path.stat().st_size,
            "mtime_utc": datetime.fromtimestamp(
                annotation_path.stat().st_mtime, tz=timezone.utc
            ).isoformat(timespec="seconds"),
            "expected_rows": args.expected_rows,
        },
        "metadata": {
            "path": str(metadata_path),
            "sha256": sha256_file(metadata_path),
            "id_column": metadata_id,
        },
        "parameters": vars(args),
        "stream_qc": stream_qc,
        "run_summary": run_summary,
        "acceptance_checks": checks,
        "historical_comparison": historical_comparison,
    }
    atomic_write_text(
        provenance_dir / "504_provenance.json",
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n",
    )

    report_lines = [
        f"# {SCRIPT_NAME}",
        "",
        f"- Status: **{overall_status}**",
        f"- Protein rows processed: {int(stream_qc['rows_seen']):,}",
        (
            "- Protein-table join: exact lockstep match by "
            f"`{args.protein_id_column}`"
        ),
        f"- Samples: {len(sample_ids)} in {n_lat_blocks} lat_block profiles",
        (
            f"- KEGG KO: {int(ko_summary['n_complete_features']):,} complete; "
            f"{int(ko_summary['n_features_prevalence_ge_threshold']):,} "
            f"present in at least {min_prevalence_n}/{len(sample_ids)} samples"
        ),
        (
            f"- PFAM: {int(pfam_summary['n_complete_features']):,} complete; "
            f"{int(pfam_summary['n_features_prevalence_ge_threshold']):,} "
            f"present in at least {min_prevalence_n}/{len(sample_ids)} samples"
        ),
        (
            "- Multi-annotation allocation: equal fractions across unique "
            "tokens within each annotation layer"
        ),
        "- Protein-level filtering before aggregation: no",
        "- Functional truncation or transformation: no",
        "- HI, MHI_local, environmental axes, and network outcomes used: no",
        "- Network inference performed: no",
        "",
        "## Acceptance checks",
        "",
        "| Check | Observed | Expected | Pass |",
        "|---|---:|---:|:---:|",
    ]
    report_lines.extend(
        "| {check} | {observed} | {expected} | {passed} |".format(
            check=check["check"],
            observed=format_number(check["observed"]),
            expected=format_number(check["expected"]),
            passed="yes" if check["pass"] else "no",
        )
        for check in checks
    )
    if historical_warning:
        report_lines.extend(["", "## Warning", "", historical_warning])
    atomic_write_text(
        out_dir / "504_functional_matrix_build_report.md",
        "\n".join(report_lines) + "\n",
    )

    print("")
    print("=" * 72)
    print(f"Overall status: {overall_status}")
    print(f"Protein rows: {int(stream_qc['rows_seen']):,}")
    print(
        "KEGG KO complete/eligible: "
        f"{int(ko_summary['n_complete_features']):,}/"
        f"{int(ko_summary['n_features_prevalence_ge_threshold']):,}"
    )
    print(
        "PFAM complete/eligible: "
        f"{int(pfam_summary['n_complete_features']):,}/"
        f"{int(pfam_summary['n_features_prevalence_ge_threshold']):,}"
    )
    print(f"Output: {out_dir}")
    print(f"LATEST: {latest_path}")
    print(f"Finished: {finished_time}")
    print("=" * 72)

    if overall_status == "FAIL" and args.strict == "1":
        return 1
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\nInterrupted by user.", file=sys.stderr)
        sys.exit(130)
    except Exception as error:
        print(
            f"\nFATAL: {type(error).__name__}: {error}",
            file=sys.stderr,
        )
        sys.exit(1)
