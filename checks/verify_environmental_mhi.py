#!/usr/bin/env python3
"""Reconstruct the fixed environmental MHI calibration without changing it.

Inputs: source_data/indices/{MHI_calibration_57,metadata_003b}.csv.
Output: validation/MHI_reconstruction_audit.json.
Algorithm: locality median/MAD scaling, Euclidean distance to the local
Conserved centroid (depths pooled), inverse locality min-max scaling.
AI-assisted implementation: OpenAI Codex (OpenAI, 2026).
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
VARIABLES = [
    "shdi", "ndvi_p50", "ndvi_sd", "ndwi_p50", "ndwi_sd", "mndwi_p50",
    "mndwi_sd", "water_prop_mndwi", "ndmi_p50", "ndmi_sd", "temperature_c",
    "salinity_ups", "ph", "redox_shallow_mv", "n_no2_umol_l", "n_no3_umol_l",
    "n_nh3_mg_l_campo", "n_nh4plus_umol_l", "p_po4_3_umol_l", "s_2_umol",
]


def main():
    calibration_path = ROOT / "source_data/indices/MHI_calibration_57.csv"
    current_path = ROOT / "source_data/indices/metadata_003b.csv"
    calibration = pd.read_csv(calibration_path).set_index("sample_id", verify_integrity=True)
    current = pd.read_csv(current_path).set_index("sample_id", verify_integrity=True)
    if len(calibration) != 57 or len(current) != 51:
        raise ValueError("Unexpected calibration/analysis cohort size")
    if not current.index.isin(calibration.index).all():
        raise ValueError("An analytical sample is missing from the calibration")
    reconstructed = pd.Series(index=calibration.index, dtype=float)
    raw = reconstructed.copy()
    records = []
    for locality, group in calibration.groupby("locality", sort=True):
        values = group[VARIABLES].to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f"Nonfinite calibration values: {locality}")
        medians = np.median(values, axis=0)
        mad = np.median(np.abs(values - medians), axis=0)
        if (mad <= 0).any():
            raise ValueError(f"Zero MAD in frozen calibration: {locality}")
        standardized = (values - medians) / mad
        conserved = group.restoration4.eq("Conserved").to_numpy()
        if not conserved.any():
            raise ValueError(f"No conserved calibration observations: {locality}")
        centroid = standardized[conserved].mean(axis=0)
        distances = np.linalg.norm(standardized - centroid, axis=1)
        span = np.ptp(distances)
        if span <= 0:
            raise ValueError(f"Degenerate distance range: {locality}")
        raw.loc[group.index] = -distances
        reconstructed.loc[group.index] = 1 - 2 * (distances - distances.min()) / span
        records.append({"locality": locality, "samples": len(group),
                        "conserved_observations": int(conserved.sum()),
                        "minimum_distance": float(distances.min()),
                        "maximum_distance": float(distances.max())})
    errors = {
        "raw_to_archived": float(np.max(np.abs(raw - calibration.MHI_local_raw))),
        "scaled_to_archived": float(np.max(np.abs(reconstructed - calibration.MHI_local))),
        "archived_to_current_51": float(np.max(np.abs(
            calibration.loc[current.index, "MHI_local"] - current.MHI_local))),
    }
    report = {
        "status": "PASS" if max(errors.values()) < 1e-9 else "FAIL",
        "calibration_samples": 57, "analytical_samples": 51,
        "variables": VARIABLES, "scaling": "within-locality median / unscaled MAD",
        "reference": "within-locality Conserved centroid; depths pooled",
        "calibration_refitted_to_51": False, "max_absolute_errors": errors,
        "localities": records,
        "input_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in [calibration_path, current_path]},
    }
    (ROOT / "validation").mkdir(parents=True,exist_ok=True)
    (ROOT / "validation/MHI_reconstruction_audit.json").write_text(
        json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ["status", "max_absolute_errors"]}, indent=2))
    if report["status"] != "PASS":
        raise ValueError("Reconstructed MHI does not match its fixed calibration")


if __name__ == "__main__":
    main()
