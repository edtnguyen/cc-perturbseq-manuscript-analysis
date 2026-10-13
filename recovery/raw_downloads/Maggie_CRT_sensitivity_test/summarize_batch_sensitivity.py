#!/usr/bin/env python3
"""Summarize NTC calibration across inferred-batch sensitivity variants."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kstest


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True, help="CRT_batch_sensitivity root directory")
    p.add_argument("--sel-thresh", default="0_4")
    p.add_argument("--output-dir", default=None)
    p.add_argument("--min-tests-per-program", type=int, default=30)
    return p.parse_args()


def qq_xy(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    p = np.asarray(p, float)
    p = p[np.isfinite(p)]
    p = np.clip(p, 1e-300, 1.0)
    p.sort()
    m = len(p)
    q = (np.arange(1, m + 1) - 0.5) / m
    return -np.log10(q), -np.log10(p)


def calibration_metrics(p: np.ndarray) -> dict[str, float]:
    p = np.asarray(p, float)
    p = p[np.isfinite(p)]
    if p.size == 0:
        return {"n_pvals": 0}
    ks = kstest(p, "uniform")
    out: dict[str, float] = {
        "n_pvals": int(p.size),
        "mean_p": float(np.mean(p)),
        "median_p": float(np.median(p)),
        "ks_stat": float(ks.statistic),
        "ks_pvalue": float(ks.pvalue),
    }
    for t in (0.05, 0.01, 0.001):
        observed = float(np.mean(p <= t))
        out[f"fraction_p_le_{t:g}"] = observed
        out[f"inflation_ratio_at_{t:g}"] = observed / t
    return out


def load_ntc_file(path: Path) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_csv(path, sep="\t")
    metadata = {"ensemble", "group_id", "condition", "batch"}
    programs = [
        c
        for c in df.columns
        if c not in metadata and pd.api.types.is_numeric_dtype(df[c])
    ]
    if not programs:
        raise ValueError(f"No numeric program p-value columns found in {path}")
    return df, programs


def main() -> None:
    args = parse_args()
    root = Path(args.root)
    out_dir = Path(args.output_dir) if args.output_dir else root / "comparison_summary"
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(root.glob(f"*/sel_thresh_{args.sel_thresh}/*_ntc_null_pvals.txt"))
    if not files:
        raise FileNotFoundError(
            f"No *_ntc_null_pvals.txt files found below {root} for sel_thresh_{args.sel_thresh}"
        )

    pooled_by_variant: dict[str, list[np.ndarray]] = {}
    overall_rows = []
    condition_rows = []
    program_rows = []

    for path in files:
        variant = path.relative_to(root).parts[0]
        # Filename form: 30_CRT_<condition>_ntc_null_pvals.txt
        prefix = path.name.removesuffix("_ntc_null_pvals.txt")
        condition = prefix.split("_CRT_", 1)[1] if "_CRT_" in prefix else prefix
        df, programs = load_ntc_file(path)
        values = df[programs].to_numpy(float).ravel()
        values = values[np.isfinite(values)]
        pooled_by_variant.setdefault(variant, []).append(values)

        row = {"variant": variant, "condition": condition, "file": str(path)}
        row.update(calibration_metrics(values))
        condition_rows.append(row)

        for program in programs:
            p = pd.to_numeric(df[program], errors="coerce").to_numpy(float)
            p = p[np.isfinite(p)]
            if len(p) < args.min_tests_per_program:
                continue
            prow = {
                "variant": variant,
                "condition": condition,
                "program": program,
            }
            prow.update(calibration_metrics(p))
            program_rows.append(prow)

    pooled_final: dict[str, np.ndarray] = {}
    for variant, chunks in pooled_by_variant.items():
        p = np.concatenate(chunks)
        pooled_final[variant] = p
        row = {"variant": variant}
        row.update(calibration_metrics(p))
        overall_rows.append(row)

    overall = pd.DataFrame(overall_rows).sort_values("ks_stat")
    conditions = pd.DataFrame(condition_rows).sort_values(["variant", "condition"])
    programs = pd.DataFrame(program_rows)
    if not programs.empty:
        programs["bad_at_0.01"] = programs["inflation_ratio_at_0.01"] >= 2.0
        programs["bad_ks"] = programs["ks_pvalue"] < 0.01

    overall.to_csv(out_dir / "ntc_calibration_overall.csv", index=False)
    conditions.to_csv(out_dir / "ntc_calibration_by_condition.csv", index=False)
    programs.to_csv(out_dir / "ntc_calibration_by_program.csv", index=False)

    print("\nOverall NTC calibration (sorted by KS statistic; lower is better):")
    cols = [
        "variant",
        "n_pvals",
        "ks_stat",
        "ks_pvalue",
        "fraction_p_le_0.05",
        "fraction_p_le_0.01",
        "fraction_p_le_0.001",
    ]
    print(overall[[c for c in cols if c in overall]].to_string(index=False))

    # Bulk QQ (raw p scale).
    plt.figure(figsize=(7, 6))
    for variant, p in pooled_final.items():
        p_sorted = np.sort(np.clip(p[np.isfinite(p)], 0, 1))
        q = (np.arange(1, len(p_sorted) + 1) - 0.5) / len(p_sorted)
        plt.plot(q, p_sorted, label=variant, linewidth=1.4)
    plt.plot([0, 1], [0, 1], color="black", linewidth=1)
    plt.xlabel("Expected p-value")
    plt.ylabel("Observed p-value")
    plt.title("NTC calibration sensitivity: bulk QQ")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(out_dir / "ntc_qq_bulk_all_variants.png", dpi=160)
    plt.close()

    # Tail QQ (-log10).
    plt.figure(figsize=(7, 6))
    max_lim = 1.0
    for variant, p in pooled_final.items():
        x, y = qq_xy(p)
        plt.plot(x, y, label=variant, linewidth=1.4)
        if x.size:
            max_lim = max(max_lim, float(np.nanmax(x)), float(np.nanmax(y)))
    plt.plot([0, max_lim], [0, max_lim], color="black", linewidth=1)
    plt.xlim(0, max_lim)
    plt.ylim(0, max_lim)
    plt.xlabel("Expected -log10(p)")
    plt.ylabel("Observed -log10(p)")
    plt.title("NTC calibration sensitivity: tail QQ")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(out_dir / "ntc_qq_tail_all_variants.png", dpi=160)
    plt.close()

    if not programs.empty:
        program_summary = (
            programs.groupby("variant")
            .agg(
                n_program_condition_tests=("program", "size"),
                fraction_bad_at_0_01=("bad_at_0.01", "mean"),
                fraction_bad_ks=("bad_ks", "mean"),
                median_program_ks=("ks_stat", "median"),
                max_program_inflation_0_01=("inflation_ratio_at_0.01", "max"),
            )
            .reset_index()
            .sort_values("median_program_ks")
        )
        program_summary.to_csv(out_dir / "program_calibration_summary.csv", index=False)
        print("\nProgram-specific calibration summary:")
        print(program_summary.to_string(index=False))

    print(f"\nWrote summaries and plots to {out_dir}")


if __name__ == "__main__":
    main()
