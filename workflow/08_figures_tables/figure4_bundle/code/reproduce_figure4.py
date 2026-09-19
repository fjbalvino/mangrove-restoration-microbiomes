#!/usr/bin/env python3
"""Reproduce Figure 4 A–D from frozen audited statistics, without fitting models.

All data values come from inputs/*.tsv. Graphical choices live in config/plot.json.
See README.md and docs/METHODS_AND_PROVENANCE.md for scope and upstream gaps.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import sys
import threading
import time
import traceback

NAME = "604d_reproducir_figura4_A_D"

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")

def build(args):
    source = args.source_dir.resolve()
    config_path = args.config.resolve()
    out = args.out_root.resolve() / (NAME + "_" + datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%fZ"))
    out.mkdir(parents=True, exist_ok=False)
    for folder in ("figures", "panels", "tables", "inputs", "code", "config"):
        (out / folder).mkdir()
    start = time.monotonic()
    lock = threading.Lock()
    stop = threading.Event()
    state = {"stage": "initialisation"}
    log_handle = (out / "run.log").open("w", buffering=1)

    def log(message):
        line = datetime.now(timezone.utc).strftime("%H:%M:%S UTC | ") + message
        with lock:
            print(line, flush=True)
            log_handle.write(line + "\n")

    def monitor():
        while not stop.wait(30):
            usage = resource.getrusage(resource.RUSAGE_SELF)
            log(f"MONITOR | pid={os.getpid()} | elapsed={time.monotonic()-start:.1f}s | "
                f"cpu={usage.ru_utime+usage.ru_stime:.1f}s | stage={state['stage']}")

    log(f"Output: {out}")
    log(f"MONITOR START | interval=30s | pid={os.getpid()}")
    worker = threading.Thread(target=monitor, daemon=True)
    worker.start()
    checks = []

    def check(name, condition):
        passed = bool(condition)
        checks.append({"check": name, "pass": passed})
        if not passed:
            raise ValueError(f"Validation failed: {name}")

    try:
        # Imports occur inside monitoring, so cold font-cache/dependency loading is logged.
        import numpy as np
        import pandas as pd
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch, Rectangle
        from matplotlib.colors import LinearSegmentedColormap

        state["stage"] = "input checks"
        cfg = json.loads(config_path.read_text())
        manifest_path = source / "SHA256SUMS.json"
        manifest = json.loads(manifest_path.read_text())
        for filename, expected in manifest.items():
            check("sha256:" + filename, sha256(source / filename) == expected)
            shutil.copy2(source / filename, out / "inputs" / filename)
        shutil.copy2(manifest_path, out / "inputs" / manifest_path.name)
        shutil.copy2(config_path, out / "config" / "plot.json")
        shutil.copy2(__file__, out / "code" / Path(__file__).name)
        counts = pd.read_csv(source / "Figure4_counts_q010_q005.tsv", sep="\t")
        stats = pd.read_csv(source / "Figure4_interactions_original_statistics.tsv", sep="\t")
        samples = pd.read_csv(source / "Figure4_two_candidates_sample_data.tsv", sep="\t")
        order, depths = cfg["axis_order"], cfg["depths_cm"]
        candidates = cfg["candidates"]
        check("counts_unique_keys", not counts.duplicated(["axis", "family", "depth_cm"]).any())
        check("six_axes", set(counts.axis) == set(order))
        check("counts_shape_30x10", counts.shape == (30, 10))
        check("nested_FDR_counts", (counts.n_q_0_05 <= counts.n_associations).all())
        directional = counts.family != "axis_by_depth"
        check("positive_plus_negative_equals_total", (counts.loc[directional, "n_positive"] + counts.loc[directional, "n_negative"] == counts.loc[directional, "n_associations"]).all())
        common = counts[counts.family == "axis_common"].set_index("axis").loc[order].reset_index()
        depth_counts = counts[counts.family == "axis_depth_slopes"].copy()
        check("common_total_461", common.n_associations.sum() == 461)
        check("depth_total_1361", depth_counts.n_associations.sum() == 1361)
        check("depth_q005_total_160", depth_counts.n_q_0_05.sum() == 160)
        check("interaction_total_two", counts.loc[counts.family == "axis_by_depth", "n_associations"].sum() == 2)
        check("MHI_zero_all_families", (counts.loc[counts.axis == "MHI", "n_associations"] == 0).all())
        check("sample_unique_keys", not samples.duplicated(["feature_id", "sample_id"]).any())
        check("102_sample_rows", len(samples) == 102)
        detected = samples.detected.astype(str).str.lower().map({"true": True, "false": False, "1": True, "0": False})
        check("detected_equals_count_positive", detected.notna().all() and (detected == (samples["count"] > 0)).all())
        detection = samples.assign(n_detected=(samples["count"] > 0).astype(int)).groupby(["feature_id", "depth_cm"]).agg(n_detected=("n_detected", "sum"), n_samples=("sample_id", "nunique")).reset_index()
        slopes = stats[stats.family == "axis_depth_slopes"].merge(detection, on=["feature_id", "depth_cm"], validate="one_to_one")
        check("six_slopes", len(slopes) == 6)
        check("17_samples_each_depth", (slopes.n_samples == 17).all())
        check("finite_effects_and_q", np.isfinite(slopes[["estimate_sd", "ci_low_sd", "ci_high_sd", "q"]]).all().all())
        check("CI_order", ((slopes.ci_low_sd <= slopes.estimate_sd) & (slopes.estimate_sd <= slopes.ci_high_sd)).all())
        for candidate in candidates:
            s = samples[samples.feature_id == candidate["feature_id"]]
            check(candidate["panel"] + "_17_complete_profiles", s.profile_id.nunique() == 17 and all(set(g.depth_cm) == set(depths) for _, g in s.groupby("profile_id")))
        common.to_csv(out / "tables" / "panel_A_counts.tsv", sep="\t", index=False)
        depth_counts.to_csv(out / "tables" / "panel_B_counts.tsv", sep="\t", index=False)
        slopes.to_csv(out / "tables" / "panels_C_D_effects_and_detection.tsv", sep="\t", index=False)
        stats[stats.family == "axis_by_depth"].to_csv(out / "tables" / "panels_C_D_interaction_tests.tsv", sep="\t", index=False)
        log(f"Input/data validation: PASS ({len(checks)} checks)")

        state["stage"] = "render A–D"
        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5,
            "axes.labelsize": 9, "axes.linewidth": .65, "svg.fonttype": "none",
            "pdf.fonttype": 42, "savefig.facecolor": "white", "axes.unicode_minus": True})
        dark, muted = cfg["dark"], cfg["muted"]
        fig = plt.figure(figsize=cfg["figsize_inches"])
        axes = {}
        width, height = np.array(cfg["figsize_inches"]) * 72
        for key, (left, top, right, bottom) in cfg["axes_rectangles_points"].items():
            ax = fig.add_axes([left / width, 1-bottom / height, (right-left) / width, (bottom-top) / height])
            for edge in ("left", "right", "top"):
                ax.spines[edge].set_visible(False)
            ax.tick_params(axis="y", length=0)
            axes[key] = ax
        ax = axes["A"]
        y = np.arange(len(order))
        ax.set(xlim=cfg["A_xlim"], ylim=(5.65, -.75))
        ax.set_axisbelow(True)
        ax.grid(axis="x", color="#e4e8eb", linewidth=.4)
        ax.axvline(0, color="#76838d", linewidth=.7)
        ax.barh(y, -common.n_negative, height=.55, color=cfg["negative_color"])
        ax.barh(y, common.n_positive, height=.55, color=cfg["positive_color"])
        for i, row in common.iterrows():
            if row.n_negative:
                ax.text(-row.n_negative-4, i, str(row.n_negative), va="center", ha="right", fontsize=8)
            if row.n_positive:
                ax.text(row.n_positive+4, i, str(row.n_positive), va="center", fontsize=8)
            if row.n_associations == 0:
                ax.plot(0, i, "o", color=dark, markersize=3)
            ax.text(217, i, str(row.n_associations), va="center", ha="right", weight="bold", fontsize=9)
        ax.set_yticks(y, common.label, fontsize=9)
        ax.tick_params(axis="y", pad=8)
        ticks = [-150, -100, -50, 0, 50, 100, 150]
        ax.set_xticks(ticks, [str(abs(x)) for x in ticks])
        ax.set_xlabel("Common-slope associations (q ≤ 0.10)")
        ax.text(.98, 1.01, "Total", transform=ax.transAxes, ha="right", color=muted, fontsize=8)
        ax.text(-.04, 1.086, "A", transform=ax.transAxes, fontsize=13, weight="bold", color=dark)
        ax.legend(handles=[Patch(color=cfg["negative_color"], label="Negative slope"), Patch(color=cfg["positive_color"], label="Positive slope")], loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, frameon=False, fontsize=8, handlelength=1.1, columnspacing=1.3, borderaxespad=0)
        ax.text(0, -.25, "MHI: 0 discoveries across all primary families.", transform=ax.transAxes, fontsize=8, color=muted)

        ax = axes["B"]
        ax.spines["bottom"].set_visible(False)
        palette = LinearSegmentedColormap.from_list("reference_teal", cfg["B_color_anchors"], N=4096)
        maximum = float(depth_counts.n_associations.max())
        for i, axis in enumerate(order):
            for j, depth in enumerate(depths):
                row = depth_counts[(depth_counts.axis == axis) & (depth_counts.depth_cm == depth)].iloc[0]
                n, small = int(row.n_associations), int(row.n_q_0_05)
                ax.add_patch(Rectangle((j-.5, i-.5), 1, 1, facecolor=palette(n/maximum), edgecolor="#b7c7cb", linewidth=.35))
                colour = "white" if n >= 160 else dark
                ax.text(j, i-.12, str(n), ha="center", va="center", fontsize=11, weight="bold", color=colour)
                ax.text(j, i+.23, f"({small})", ha="center", va="center", fontsize=8, color=colour)
        ax.set(xlim=(-.5, 2.5), ylim=(5.5, -.5))
        ax.set_yticks(y, common.label, fontsize=9)
        ax.tick_params(axis="y", pad=8)
        ax.set_xticks(range(3), [f"{d} cm" for d in depths])
        ax.tick_params(axis="x", length=0)
        ax.set_xlabel("Depth-specific associations")
        ax.text(0, 1.05, "Count at q ≤ 0.10  ·  (count at q ≤ 0.05)", transform=ax.transAxes, color=muted, fontsize=8)
        ax.text(-.04, 1.086, "B", transform=ax.transAxes, fontsize=13, weight="bold", color=dark)
        ax.text(0, -.25, "Counts do not test differences between slopes.", transform=ax.transAxes, fontsize=8, color=muted)

        for candidate in candidates:
            panel, fid, axis = candidate["panel"], candidate["feature_id"], candidate["axis"]
            ax = axes[panel]
            selected = slopes[(slopes.feature_id == fid) & (slopes.axis == axis)].set_index("depth_cm").loc[depths]
            interaction = stats[(stats.feature_id == fid) & (stats.axis == axis) & (stats.family == "axis_by_depth")].iloc[0]
            ax.set(xlim=cfg["CD_xlim"], ylim=(2.55, -.83))
            ax.axvline(0, color="#87929a", linestyle="--", linewidth=.7)
            for i, (depth, row) in enumerate(selected.iterrows()):
                colour = cfg["depth_colors"][str(int(depth))]
                filled = row.q <= cfg["primary_q"]
                ax.errorbar(row.estimate_sd, i, xerr=[[row.estimate_sd-row.ci_low_sd], [row.ci_high_sd-row.estimate_sd]], fmt="o", color=colour, markerfacecolor=colour if filled else "white", markeredgewidth=1.4, markersize=7, elinewidth=1.5, capsize=3, capthick=1.4)
                ax.text(2.2, i, f"{row.q:.4f}", ha="center", va="center")
                ax.text(2.96, i, f"{int(row.n_detected)}/{int(row.n_samples)}", ha="center", va="center")
            ax.text(2.2, -.57, "Slope q", ha="center", fontsize=8, color=muted)
            ax.text(2.96, -.57, "Detected", ha="center", fontsize=8, color=muted)
            ax.set_yticks(range(3), [f"{d} cm" for d in depths], fontsize=9)
            ax.set_xticks([-.5, 0, .5, 1, 1.5])
            ax.spines["bottom"].set_bounds(-.8, 1.7)
            ax.set_xlabel("Standardised slope (pointwise 95% CI)")
            ax.text(0, 1.2, candidate["title"], transform=ax.transAxes, weight="bold", fontsize=10, color=dark)
            ax.text(0, 1.08, f"{candidate['axis_label']}  |  interaction q = {interaction.q:.4f}", transform=ax.transAxes, color=muted)
            ax.text(-.10, 1.2073, panel, transform=ax.transAxes, fontsize=13, weight="bold", color=dark)
            ax.text(0, -.24, fid, transform=ax.transAxes, fontsize=8, color=muted)
        fig.legend(handles=[Line2D([], [], marker="o", linestyle="none", color=dark, label="Slope q ≤ 0.10", markersize=6), Line2D([], [], marker="o", linestyle="none", color=dark, markerfacecolor="white", label="Slope q > 0.10", markersize=6)], loc="lower center", bbox_to_anchor=(.58, .025), ncol=2, frameon=False, fontsize=8)
        for extension in ("png", "pdf", "svg"):
            fig.savefig(out / "figures" / f"Figure4_reproduced.{extension}", dpi=350)
        fig.savefig(out / "figures" / "Figure4_preview.png", dpi=150)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for key, ax in axes.items():
            box = ax.get_tightbbox(renderer).transformed(fig.dpi_scale_trans.inverted()).expanded(1.025, 1.04)
            for extension in ("pdf", "svg"):
                fig.savefig(out / "panels" / f"Figure4_panel_{key}.{extension}", bbox_inches=box)
        plt.close(fig)
        state["stage"] = "output checks"
        # Round-trip check that exported plotting data preserve all coefficients and q.
        reread = pd.read_csv(out / "tables" / "panels_C_D_effects_and_detection.tsv", sep="\t")
        check("export_numeric_roundtrip", np.allclose(reread[["estimate_sd", "ci_low_sd", "ci_high_sd", "q"]], slopes[["estimate_sd", "ci_low_sd", "ci_high_sd", "q"]], rtol=0, atol=1e-14))
        for extension in ("png", "pdf", "svg"):
            check("nonempty_figure_" + extension, (out / "figures" / f"Figure4_reproduced.{extension}").stat().st_size > 1000)
        write_json(out / "validation.json", {"status": "PASS", "scope": "frozen_statistics_to_graphics", "models_refitted": False, "FDR_recalculated": False, "HFR_computed": False, "checks": checks})
        write_json(out / "run_metadata.json", {"script": NAME, "started_utc": out.name.split(NAME + "_")[1], "command_argv": sys.argv, "python": sys.version, "platform": platform.platform(), "versions": {x: importlib.metadata.version(x) for x in ("numpy", "pandas", "matplotlib")}, "source_directory": str(source), "output_directory": str(out), "input_hashes": manifest, "plot_config_sha256": sha256(config_path), "script_sha256": sha256(__file__), "monitor_interval_seconds": 30, "elapsed_seconds": time.monotonic()-start, "meaning": "Graphical reconstruction; exact statistical values, approximately preserved layout. No upstream model rerun."})
        log(f"PASS | checks={len(checks)} | models_refitted=false | FDR_recalculated=false")
        return out
    except Exception:
        write_json(out / "validation.json", {"status": "FAIL", "checks": checks, "error": traceback.format_exc()})
        log(traceback.format_exc())
        raise
    finally:
        stop.set()
        worker.join(timeout=2)
        log(f"MONITOR END | elapsed={time.monotonic()-start:.1f}s")
        log_handle.close()
        # Includes the closed log, copied source code, frozen inputs, figures and metadata.
        write_json(out / "SHA256SUMS.json", {str(p.relative_to(out)): sha256(p) for p in sorted(out.rglob("*")) if p.is_file() and p != out / "SHA256SUMS.json"})

if __name__ == "__main__":
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=root / "inputs")
    parser.add_argument("--config", type=Path, default=root / "config" / "plot.json")
    parser.add_argument("--out-root", type=Path, default=root / "outputs")
    build(parser.parse_args())
