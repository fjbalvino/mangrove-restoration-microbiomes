#!/usr/bin/env python3
# ============================================================
# 01_03_import_corrected_MHI.py
# Status: REVIEW CANDIDATE; see docs/REPRODUCIBILITY.md
# Imports the hash-verified 57-row historical environmental calibration; preserves the archived HI.
# Inputs (source expressions; complete list in docs/contracts/01_03_import_corrected_MHI.json):
#   return hashlib.sha256(path.read_bytes()).hexdigest()
#   def read_csv(path):
#   fields, base, base_ids = read_csv(sources["base_csv"])
#   _, mhi, mhi_ids = read_csv(sources["mhi_csv"])
#   _, canon, canon_ids = read_csv(sources["canon_csv"])
#   saved_fields, saved, _ = read_csv(target)
# Outputs (source expressions; complete list in contract):
#   (out / "logs/003b_provenance.json").write_text(
#   temporary.write_text(str(out) + "\n", encoding="utf-8")
# Algorithmic provenance:
# SHA-256-verified import of the fixed 57-observation environmental MHI; preserve HI.
#   Study-specific environmental calibration and data integration; see docs/METHODS.md.
# Source SHA-256: 3d7de01f8d2b0c9f5a1acc0fbe56eb957e28e691f519956394f33c0d7464f306
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Changes in this copy: descriptive filename and this documentation header only.
# Original body and internal output names are preserved exactly.
# ============================================================

"""Integrate audited environmental MHI without refitting or changing other axes."""
import argparse
import csv
import hashlib
import json
import math
import os
import shutil
from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path

NAME = "003b_integrar_MHI_ambiental_validado"
TP = Path("/home/fjbalvino/Tipping_points")
HASHES = {
    "base_csv": "a9fd144387f9cda33eae0fed23cdaefb110f2a75b8b0f58d4adbdfd5e0670c62",
    "mhi_csv": "ef48c2052d985b08893457dfb6f0249e1b0d209cda10aae3a1a749c4e23c413e",
    "canon_csv": "8232edd8008922d26f289ee8db9e273255bd3c1b360baf277fbe5430fac06993",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        rows = list(reader)
    require(len(fields) == len(set(fields)), f"Columnas duplicadas: {path}")
    require(rows and "sample_id" in fields, f"CSV sin muestras o sample_id: {path}")
    ids = [row["sample_id"] for row in rows]
    require(all(ids) and len(ids) == len(set(ids)), f"IDs vacios/duplicados: {path}")
    require(all(set(row) == set(fields) and None not in row.values() for row in rows),
            f"Filas incompletas: {path}")
    return fields, rows, {row["sample_id"]: row for row in rows}


def write_table(path, fields, rows, delimiter=","):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base_csv", type=Path, default=TP / "resultados_finales/003_integrar_ejes_ECI_HI_en_metadata_20260909_232743/tables/003_metadata_integrada_canon_51.csv")
    parser.add_argument("--mhi_csv", type=Path, default=Path("/home/fjbalvino/R/Ciencia-Frontera/resultados/46_merge_mhi_metadata_20260216_223800/metadata_with_MHI.csv"))
    parser.add_argument("--canon_csv", type=Path, default=TP / "resultados_finales/202_construir_matriz_CLR_y_eje_bimodal_20260912_005021/tables/202_metadata_raw1031_alineada.csv")
    parser.add_argument("--out_root", type=Path, default=TP / "resultados_finales")
    args = parser.parse_args()
    sources = {key: getattr(args, key).resolve() for key in HASHES}
    for key, path in sources.items():
        require(path.is_file(), f"No existe: {path}")
        require(sha(path) == HASHES[key], f"Archivo distinto del auditado: {path}")

    fields, base, base_ids = read_csv(sources["base_csv"])
    _, mhi, mhi_ids = read_csv(sources["mhi_csv"])
    _, canon, canon_ids = read_csv(sources["canon_csv"])
    require(len(base) == len(canon) == 51 and len(mhi) == 57, "Universo inesperado")
    require(set(base_ids) == set(canon_ids) <= set(mhi_ids), "Muestras incompatibles")
    require({"MHI_local", "MHI_global"} <= set(fields), "Faltan indices en 003")
    profiles = {}
    for row in canon:
        profiles.setdefault(row["profile_id"], []).append(row)
    require(len(profiles) == 17, "Se requieren 17 perfiles")
    for profile, rows in profiles.items():
        require(sorted(Decimal(r["depth_cm"]) for r in rows) == [5, 20, 40],
                f"Perfil incompleto: {profile}")
        require(all(len({r[key] for r in rows}) == 1
                    for key in ("locality", "restoration4", "lat_block")),
                f"Perfil inconsistente: {profile}")

    output_fields = [key for key in fields if key != "MHI_global"]
    unchanged = [key for key in fields if key not in ("MHI_local", "MHI_global")]
    corrected, audit = [], []
    for row in base:
        sid = row["sample_id"]
        reference, canonical = mhi_ids[sid], canon_ids[sid]
        for key in ("locality", "restoration4", "depth_cm", "lat_block"):
            vals = [obj[key] for obj in (row, reference, canonical)]
            if key in ("depth_cm", "lat_block"):
                vals = list(map(Decimal, vals))
            require(len(set(vals)) == 1, f"Metadata discrepante: {sid}, {key}")
        value = float(reference["MHI_local"])
        require(math.isfinite(value) and -1 <= value <= 1, f"MHI invalido: {sid}")
        new = {key: row[key] for key in output_fields}
        new["MHI_local"] = reference["MHI_local"]
        require(all(new[key] == row[key] for key in unchanged), "Cambio no permitido")
        corrected.append(new)
        audit.append({"sample_id": sid, "profile_id": canonical["profile_id"],
                      "MHI_local_47_anterior": row["MHI_local"],
                      "MHI_local_46_validado": new["MHI_local"],
                      "MHI_global_47_retirado": row["MHI_global"],
                      "diferencia": value - float(row["MHI_local"])})

    now = datetime.now(timezone.utc)
    out = args.out_root.resolve() / f"{NAME}_{now:%Y%m%d_%H%M%S}"
    out.mkdir(parents=True, exist_ok=False)
    for directory in ("tables", "inputs", "logs"):
        (out / directory).mkdir()
    target = out / "tables/003_metadata_integrada_canon_51.csv"
    write_table(target, output_fields, corrected)
    saved_fields, saved, _ = read_csv(target)
    require(saved_fields == output_fields and saved == corrected, "Error de exportacion")
    require([r["sample_id"] for r in saved] == [r["sample_id"] for r in base],
            "Cambio de orden")
    write_table(out / "tables/003b_cambios_MHI.tsv", list(audit[0]), audit, "\t")
    summary = {"status": "PASS", "samples": 51, "profiles": 17,
               "calibration_samples": 57, "calibration_profiles": 19,
               "unchanged_columns": len(unchanged), "MHI_global_removed": True,
               "MHI_local_changed": sum(abs(r["diferencia"]) > 1e-12 for r in audit),
               "MHI_refitted": False, "models_executed": False}
    for key, path in sources.items():
        shutil.copy2(path, out / "inputs" / f"{key}.csv")
    shutil.copy2(Path(__file__), out / "logs" / f"{NAME}.py")
    provenance = {"started_utc": now.isoformat(), "script_sha256": sha(Path(__file__)),
                  "sources": {key: {"path": str(path), "sha256": sha(path)}
                              for key, path in sources.items()},
                  "metadata_sha256": sha(target), "summary": summary}
    (out / "logs/003b_provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_table(out / "tables/003b_run_summary.tsv", list(summary), [summary], "\t")
    pointer = args.out_root.resolve() / f"LATEST_{NAME}.txt"
    temporary = pointer.with_name(pointer.name + f".tmp_{os.getpid()}")
    temporary.write_text(str(out) + "\n", encoding="utf-8")
    temporary.replace(pointer)
    print(json.dumps(summary, ensure_ascii=False), flush=True)
    print(f"Metadata: {target}\nPointer: {pointer}", flush=True)


if __name__ == "__main__":
    main()
