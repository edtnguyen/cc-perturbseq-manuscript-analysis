from __future__ import annotations

from pathlib import Path
import hashlib

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def sha256_file(path: str | Path) -> str:
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _save_figure(fig, outdir: Path, stem: str) -> dict[str, str]:
    outputs: dict[str, str] = {}
    for suffix in ("pdf", "svg", "png"):
        path = outdir / f"{stem}.{suffix}"
        fig.savefig(
            path,
            dpi=600 if suffix == "png" else None,
            bbox_inches="tight",
            pad_inches=0.03,
            facecolor="white",
        )
        outputs[suffix] = str(path)
    return outputs


# -----------------------------------------------------------------------------
# Figure 1c: validation screen guide-level EC percentages
# Source: validation_boxplots_verified bundle
# -----------------------------------------------------------------------------

PANEL_C_VALUE_COL = "estimated_EC_pct_control70"
PANEL_C_EXPECTED_SHA256 = (
    "e4d79dce678bde0a84dcaf6927573c24003f362aae3a02d11d4670a2516a5c1e"
)
PANEL_C_CONTROL_TYPES = {"NON-TARGETING", "SAFE"}
PANEL_C_GROUPS = [
    ("Control", ["NTC"]),
    ("BMP/TGFβ pathway", ["BMPR1A", "SMAD1", "SMAD5"]),
    ("Mesoderm regulators", ["MIXL1", "EOMES"]),
    ("EC + LM regulators", ["KDR", "KCTD10", "PRKCD", "ETV2", "PLCG1"]),
    ("Hippo pathway", ["TAOK1", "AMOTL2", "LATS2", "NF2"]),
]
PANEL_C_GROUP_LABELS = {
    "Control": "Control",
    "BMP/TGFβ pathway": "BMP/TGFβ\npathway",
    "Mesoderm regulators": "Mesoderm\nregulators",
    "EC + LM regulators": "EC + LM\nregulators",
    "Hippo pathway": "Hippo\npathway",
}


def panel_c_bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1 or not np.all(np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("BH adjustment requires finite p-values in [0, 1].")
    order = np.argsort(p, kind="stable")
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.empty(len(p), dtype=float)
    q[order] = np.clip(np.minimum.accumulate(ranked[::-1])[::-1], 0, 1)
    return q


def panel_c_stars(q: float) -> str:
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""


def read_panel_c_guides(path: str | Path):
    path = Path(path)
    table = pd.read_csv(path)
    required = {"guide_id", "target_gene", "guide_type", PANEL_C_VALUE_COL}
    if missing := required.difference(table.columns):
        raise ValueError(f"Missing validation-table columns: {sorted(missing)}")
    if table["guide_id"].isna().any() or table["guide_id"].duplicated().any():
        raise ValueError("Guide IDs must be nonempty and unique.")

    table = table.copy()
    table["guide_id"] = table["guide_id"].astype(str).str.strip()
    table["target_gene"] = table["target_gene"].astype(str).str.strip().str.upper()
    table["guide_type"] = table["guide_type"].astype(str).str.strip().str.upper()
    table[PANEL_C_VALUE_COL] = pd.to_numeric(
        table[PANEL_C_VALUE_COL], errors="coerce"
    )
    if table[PANEL_C_VALUE_COL].isna().any():
        raise ValueError("The verified input contains a nonfinite EC estimate.")
    if not table[PANEL_C_VALUE_COL].between(0, 100).all():
        raise ValueError("EC percentages must be in [0, 100].")

    allowed = PANEL_C_CONTROL_TYPES | {"TARGETING"}
    bad_types = sorted(set(table["guide_type"]).difference(allowed))
    if bad_types:
        raise ValueError(f"Unrecognized guide types: {bad_types}")

    controls = table.loc[
        table["guide_type"].isin(PANEL_C_CONTROL_TYPES)
    ].sort_values("guide_id")
    targeting = {
        gene: group.sort_values("guide_id").copy()
        for gene, group in table.loc[
            table["guide_type"].eq("TARGETING")
        ].groupby("target_gene", sort=True)
    }
    for _, genes in PANEL_C_GROUPS[1:]:
        for gene in genes:
            if gene not in targeting:
                raise ValueError(f"Requested Figure 1c gene is missing: {gene}")
    return table, controls.copy(), targeting


def panel_c_uniform_subsets(
    n_control: int,
    n_target: int,
    permutations: int,
    seed: int,
) -> np.ndarray:
    total = n_control + n_target
    rng = np.random.default_rng(np.random.SeedSequence([seed, n_control, n_target]))
    dtype = np.uint16 if total < 65536 else np.uint32
    output = np.empty((permutations, n_target), dtype=dtype)
    for lo in range(0, permutations, 5000):
        hi = min(lo + 5000, permutations)
        keys = rng.random((hi - lo, total))
        chosen = np.argpartition(keys, n_target - 1, axis=1)[:, :n_target]
        chosen.sort(axis=1)
        output[lo:hi] = chosen
    return output


def compute_panel_c_statistics(
    controls: pd.DataFrame,
    targeting: dict[str, pd.DataFrame],
    outdir: str | Path,
    *,
    permutations: int = 100_000,
    seed: int = 20260908,
) -> pd.DataFrame:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    control_values = controls[PANEL_C_VALUE_COL].to_numpy(float)
    control_mean = float(control_values.mean())
    subset_cache = {
        n: panel_c_uniform_subsets(len(control_values), n, permutations, seed)
        for n in sorted({len(group) for group in targeting.values()})
    }
    np.savez_compressed(
        outdir / "permutation_indices_100k.npz",
        **{f"n_target_{n}": indices for n, indices in subset_cache.items()},
    )

    rows = []
    for gene, group in targeting.items():
        target_values = group[PANEL_C_VALUE_COL].to_numpy(float)
        n_target = len(target_values)
        target_mean = float(target_values.mean())
        observed = target_mean - control_mean
        pooled = np.concatenate([target_values, control_values])
        centered = pooled - pooled.mean(dtype=np.float64)
        centered_total = centered.sum(dtype=np.float64)
        tolerance = 100 * np.finfo(float).eps * max(1.0, abs(observed))
        extreme = 0
        indices = subset_cache[n_target]
        for lo in range(0, permutations, 10_000):
            chunk = indices[lo : lo + 10_000]
            selected_sum = centered[chunk].sum(axis=1, dtype=np.float64)
            permuted = selected_sum / n_target - (
                centered_total - selected_sum
            ) / len(control_values)
            extreme += int(
                np.count_nonzero(np.abs(permuted) >= abs(observed) - tolerance)
            )
        rows.append(
            {
                "gene": gene,
                "n_targeting_guides": n_target,
                "n_control_guides": len(control_values),
                "mean_target_EC_pct": target_mean,
                "median_target_EC_pct": float(np.median(target_values)),
                "mean_control_EC_pct": control_mean,
                "observed_difference_target_minus_control_pp": observed,
                "effect_log2_control_over_target": float(
                    np.log2(control_mean / target_mean)
                ),
                "n_extreme_permutations": extreme,
                "n_permutations": permutations,
                "permutation_pvalue_two_sided_b_over_B": extreme / permutations,
                "permutation_pvalue_plus_one_sensitivity": (extreme + 1)
                / (permutations + 1),
                "seed": seed,
            }
        )

    stats = pd.DataFrame(rows)
    stats["FDR_BH"] = panel_c_bh_adjust(
        stats["permutation_pvalue_two_sided_b_over_B"].to_numpy()
    )
    stats["FDR_BH_plus_one_sensitivity"] = panel_c_bh_adjust(
        stats["permutation_pvalue_plus_one_sensitivity"].to_numpy()
    )
    stats["significance"] = stats["FDR_BH"].map(panel_c_stars)
    stats["BH_number_genes"] = len(stats)
    return stats


def plot_panel_c(
    controls: pd.DataFrame,
    targeting: dict[str, pd.DataFrame],
    stats: pd.DataFrame,
):
    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.family": "DejaVu Sans",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.transparent": False,
        }
    )
    stats_index = stats.set_index("gene")
    datasets: list[np.ndarray] = []
    labels: list[str] = []
    positions: list[float] = []
    spans: list[tuple[str, float, float]] = []
    names: list[str] = []
    jitter_rows: list[dict[str, object]] = []
    x = 0.0
    for group_name, genes in PANEL_C_GROUPS:
        start = x
        for gene in genes:
            group = controls if gene == "NTC" else targeting[gene]
            values = group[PANEL_C_VALUE_COL].to_numpy(float)
            positions.append(x)
            datasets.append(values)
            names.append(gene)
            labels.append(f"{gene}\nn = {len(values)}")
            x += 1.0
        spans.append((group_name, start, x - 1.0))
        x += 0.90

    x_min, x_max = -0.6, positions[-1] + 0.6
    width_per_position = 0.58
    ax_width = (x_max - x_min) * width_per_position
    ax_height, left, right, bottom, top = 3.85, 1.15, 0.15, 1.80, 0.62
    fig_w, fig_h = left + ax_width + right, bottom + ax_height + top
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes(
        [left / fig_w, bottom / fig_h, ax_width / fig_w, ax_height / fig_h]
    )
    ax.axhline(
        np.mean(datasets[0]),
        color="#9EA3A8",
        linestyle="--",
        linewidth=1.6,
        zorder=1,
    )
    boxplot = ax.boxplot(
        datasets,
        positions=positions,
        widths=0.42,
        whis=1.5,
        patch_artist=True,
        showfliers=False,
        manage_ticks=False,
        medianprops={"color": "black", "linewidth": 2.0, "zorder": 5},
        whiskerprops={"color": "black", "linewidth": 1.8, "zorder": 5},
        capprops={"color": "black", "linewidth": 1.8, "zorder": 5},
        boxprops={"edgecolor": "black", "linewidth": 1.8, "zorder": 3},
    )
    for box, gene in zip(boxplot["boxes"], names):
        box.set_facecolor("white" if gene == "NTC" else "#E3B51E")

    for gene, position, values in zip(names, positions, datasets):
        digest = int.from_bytes(hashlib.sha256(gene.encode()).digest()[:4], "little")
        rng = np.random.default_rng(np.random.SeedSequence([7, digest]))
        jitter = np.clip(
            rng.normal(0, 0.09 if gene == "NTC" else 0.055, len(values)),
            -0.16,
            0.16,
        )
        ax.scatter(
            np.full(len(values), position) + jitter,
            values,
            s=16 if gene == "NTC" else 18,
            c="#7A7A7A",
            alpha=0.48,
            linewidths=0,
            edgecolors="none",
            zorder=6,
            rasterized=False,
        )
        group = controls if gene == "NTC" else targeting[gene]
        for (_, guide), offset in zip(group.iterrows(), jitter):
            jitter_rows.append(
                {
                    "plot_gene": gene,
                    "guide_id": guide["guide_id"],
                    PANEL_C_VALUE_COL: guide[PANEL_C_VALUE_COL],
                    "x_center": position,
                    "x_jitter": float(offset),
                }
            )
        if gene != "NTC":
            annotation = stats_index.loc[gene, "significance"]
            if annotation:
                ax.text(
                    position,
                    95,
                    annotation,
                    ha="center",
                    va="bottom",
                    fontsize=12,
                    fontweight="bold",
                    zorder=8,
                )

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0, 105)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.set_xticks(positions)
    ax.set_xticklabels(
        labels,
        rotation=58,
        ha="right",
        rotation_mode="anchor",
        fontsize=12,
    )
    ax.set_ylabel("Estimated D3 EC cells (%)", fontsize=18, labelpad=8)
    ax.tick_params(axis="y", labelsize=12, width=2.0, length=6)
    ax.tick_params(axis="x", width=2.0, length=6, pad=6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(2.2)

    transform = ax.get_xaxis_transform()
    for group_name, start, end in spans:
        ax.plot(
            [start - 0.32, end + 0.32],
            [-0.290, -0.290],
            transform=transform,
            clip_on=False,
            color="black",
            linewidth=1,
        )
        ax.text(
            (start + end) / 2,
            -0.323,
            PANEL_C_GROUP_LABELS[group_name],
            transform=transform,
            ha="center",
            va="top",
            fontsize=10.5,
            linespacing=1.15,
            fontweight="bold",
            clip_on=False,
        )

    fig.text(
        left / fig_w,
        (bottom + ax_height + 0.24) / fig_h,
        "100,000 two-sided guide-label permutations; BH correction across all tested genes",
        ha="left",
        va="bottom",
        fontsize=9.5,
    )
    fig.text(
        (left + ax_width) / fig_w,
        (bottom + ax_height + 0.24) / fig_h,
        "* FDR < 0.05   ** < 0.01   *** < 0.001",
        ha="right",
        va="bottom",
        fontsize=9.5,
    )
    return fig, pd.DataFrame(jitter_rows)


def run_figure_01c(
    validation_path: str | Path,
    outdir: str | Path,
    *,
    verify_hash: bool = True,
    permutations: int = 100_000,
    seed: int = 20260908,
):
    validation_path = Path(validation_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    actual_hash = sha256_file(validation_path)
    if verify_hash and actual_hash != PANEL_C_EXPECTED_SHA256:
        raise ValueError(
            "Figure 1c input does not match the verified bundle. "
            f"Expected {PANEL_C_EXPECTED_SHA256}, found {actual_hash}."
        )

    table, controls, targeting = read_panel_c_guides(validation_path)
    if len(controls) != 120 or len(targeting) != 297:
        raise ValueError(
            "Expected 120 control guides and 297 target genes; "
            f"found {len(controls)} and {len(targeting)}."
        )
    stats = compute_panel_c_statistics(
        controls,
        targeting,
        outdir,
        permutations=permutations,
        seed=seed,
    )
    fig, jitter = plot_panel_c(controls, targeting, stats)
    outputs = _save_figure(
        fig, outdir, "Figure_01c_validation_boxplots_verified"
    )

    selected_genes = [gene for _, genes in PANEL_C_GROUPS[1:] for gene in genes]
    selected_stats = stats.set_index("gene").loc[selected_genes].reset_index()
    stats.to_csv(outdir / "permutation_statistics_all_297_genes.csv", index=False)
    selected_stats.to_csv(
        outdir / "permutation_statistics_plotted_genes.csv", index=False
    )
    table.to_csv(outdir / "guide_values_used.csv", index=False)
    jitter.to_csv(outdir / "plot_jitter_coordinates.csv", index=False)
    return fig, stats, selected_stats, outputs, actual_hash
