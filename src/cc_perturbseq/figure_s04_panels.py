"""Ported S4 Q-Q and biological-replicate effect producers.
The original cells 13 and 21 are scientific provenance, not runtime inputs.
"""
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from IPython.display import display
from plotnine import (
    ggplot, aes, geom_abline, geom_point,
    scale_color_manual, scale_x_continuous, scale_y_continuous,
    coord_equal, labs, guides, guide_legend,
    theme_bw, theme,
    element_text, element_line, element_rect, element_blank,
)

# Keep text editable in Illustrator when possible.
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"

QQ_REAL_COLOR = "#17365D"
QQ_NTC_COLOR = "#B22222"
QQ_DIAGONAL_COLOR = "#AEB7C2"
QQ_GRID_COLOR = "#E8EDF3"
QQ_BORDER_COLOR = "#182635"


def _qq_columns_for_one_series(d: pd.DataFrame) -> pd.DataFrame:
    """Use only if input has raw p-values but not precomputed QQ x/y columns."""
    d = d.copy()
    d["pvalue"] = pd.to_numeric(d["pvalue"], errors="coerce")
    d = d[np.isfinite(d["pvalue"])].copy()
    d["pvalue"] = d["pvalue"].clip(1e-300, 1.0)

    d = d.sort_values("pvalue", kind="mergesort").reset_index(drop=True)
    n = len(d)
    d["rank"] = np.arange(1, n + 1)
    d["expected_pvalue"] = (d["rank"] - 0.5) / n
    d["expected_neglog10_p"] = -np.log10(d["expected_pvalue"])
    d["observed_neglog10_p"] = -np.log10(d["pvalue"])
    return d


def plot_ipsc_ec_qq(
    points,
    *,
    day=None,
    pvalue_type="empirical",
    title=None,
    subtitle="RNA QC, guide count and guide UMI burden",
    output_stem=None,
    limit=None,
    width=4.6,
    height=4.9,
    dpi=300,
    formats=("pdf", "svg", "eps", "png"),
):
    """Plot exported iPSC-EC QQ points with an above-panel legend."""
    df = pd.read_csv(points) if isinstance(points, (str, Path)) else points.copy()

    if day is not None and "day" in df.columns:
        df = df[df["day"].astype(str).eq(str(day))].copy()

    if pvalue_type is not None and "pvalue_type" in df.columns:
        df = df[df["pvalue_type"].astype(str).eq(str(pvalue_type))].copy()

    series_map = {
        "observed_target": "Real targets",
        "ntc_pseudotarget": "NTC pseudo-targets",
    }

    if "series_type" not in df.columns:
        raise ValueError("Expected a 'series_type' column in the QQ points table.")

    df = df[df["series_type"].isin(series_map)].copy()
    df["qq_series"] = df["series_type"].map(series_map)

    if not {"expected_neglog10_p", "observed_neglog10_p"}.issubset(df.columns):
        if "pvalue" not in df.columns:
            raise ValueError("Need either QQ coordinate columns or a raw 'pvalue' column.")
        df = pd.concat(
            [_qq_columns_for_one_series(g) for _, g in df.groupby("qq_series", sort=False)],
            ignore_index=True,
            sort=False,
        )

    df = df[
        np.isfinite(df["expected_neglog10_p"])
        & np.isfinite(df["observed_neglog10_p"])
    ].copy()

    if df.empty:
        raise ValueError("No points to plot after filtering.")

    ntc_df = df[df["qq_series"].eq("NTC pseudo-targets")].copy()
    real_df = df[df["qq_series"].eq("Real targets")].copy()

    if limit is None:
        limit = float(
            np.nanmax(
                [
                    df["expected_neglog10_p"].max(),
                    df["observed_neglog10_p"].max(),
                    4.0,
                ]
            )
        )
        limit = np.ceil(limit * 10) / 10

    breaks = np.arange(0, np.floor(limit) + 1, 1)

    if title is None:
        title = "Real target and NTC empirical p-values"
        if day is not None:
            title = f"{day}: {title}"

    p = (
        ggplot()
        + geom_abline(
            intercept=0,
            slope=1,
            linetype="dashed",
            color=QQ_DIAGONAL_COLOR,
            size=0.95,
        )
        + geom_point(
            data=ntc_df,
            mapping=aes(
                "expected_neglog10_p",
                "observed_neglog10_p",
                color='"NTC pseudo-targets"',
            ),
            size=1.55,
            alpha=0.78,
            stroke=0,
        )
        + geom_point(
            data=real_df,
            mapping=aes(
                "expected_neglog10_p",
                "observed_neglog10_p",
                color='"Real targets"',
            ),
            size=2.05,
            alpha=0.96,
            stroke=0,
        )
        + scale_color_manual(
            values={
                "Real targets": QQ_REAL_COLOR,
                "NTC pseudo-targets": QQ_NTC_COLOR,
            },
            breaks=["Real targets", "NTC pseudo-targets"],
            name="",
        )
        + scale_x_continuous(
            limits=(0, limit),
            breaks=breaks,
            expand=(0.01, 0.01),
        )
        + scale_y_continuous(
            limits=(0, limit),
            breaks=breaks,
            expand=(0.01, 0.01),
        )
        + coord_equal()
        + labs(
            x="Expected -log10(p)",
            y="Observed -log10(p)",
            title=title,
            subtitle=subtitle,
        )
        + guides(
            color=guide_legend(
                nrow=1,
                override_aes={"size": 3.8, "alpha": 1},
            )
        )
        + theme_bw(base_size=10)
        + theme(
            figure_size=(width, height),
            plot_title=element_text(size=12, weight="bold", ha="center"),
            plot_subtitle=element_text(size=8.5, weight="bold", ha="center"),
            axis_title=element_text(size=10),
            axis_text=element_text(size=9),
            panel_grid_major=element_line(color=QQ_GRID_COLOR, size=0.55),
            panel_grid_minor=element_blank(),
            panel_border=element_rect(color=QQ_BORDER_COLOR, size=0.8),
            legend_position="top",
            legend_direction="horizontal",
            legend_background=element_blank(),
            legend_key=element_blank(),
            legend_title=element_blank(),
            legend_text=element_text(size=9),
        )
    )

    saved = {}
    if output_stem is not None:
        output_stem = Path(output_stem)
        output_stem.parent.mkdir(parents=True, exist_ok=True)

        for fmt in formats:
            out = output_stem.with_suffix(f".{fmt}")
            kwargs = dict(width=width, height=height, units="in", verbose=False)
            if fmt.lower() == "png":
                kwargs["dpi"] = dpi
            p.save(out, **kwargs)
            saved[fmt] = out

    return p, saved

from pathlib import Path
import itertools
import matplotlib as mpl
import numpy as np
import pandas as pd
from IPython.display import display
from plotnine import (
    ggplot, aes,
    geom_abline, geom_hline, geom_vline, geom_point,
    facet_wrap,
    scale_x_continuous, scale_y_continuous,
    coord_equal,
    labs,
    theme_bw, theme,
    element_text, element_line, element_rect, element_blank,
)

NTC_LABELS = ["SAFE", "non-targeting", "NTC"]
NTC_KEY_SET = {str(x).strip().lower() for x in NTC_LABELS}
def read_beta_wide_as_long(beta_path, *, day, biological_sample):
    beta_wide = pd.read_csv(beta_path, sep="\t", index_col=0)
    beta_wide.index = beta_wide.index.astype(str)
    beta_wide.columns = beta_wide.columns.astype(str)
    beta_wide.index.name = "target_name"

    beta_long = (
        beta_wide
        .reset_index()
        .melt(
            id_vars="target_name",
            var_name="program_name",
            value_name="beta",
        )
    )

    beta_long["day"] = day
    beta_long["biological_sample"] = biological_sample
    beta_long["beta"] = pd.to_numeric(beta_long["beta"], errors="coerce")
    return beta_long

def format_correlation(value):
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:.2f}"

def save_plotnine_formats(
    p,
    output_stem,
    *,
    width,
    height,
    dpi=300,
    formats=("pdf", "svg", "eps", "png"),
):
    output_stem = Path(output_stem)
    output_stem.parent.mkdir(parents=True, exist_ok=True)

    saved = {}
    for fmt in formats:
        out = output_stem.with_suffix(f".{fmt}")
        kwargs = dict(width=width, height=height, units="in", verbose=False)
        if fmt.lower() == "png":
            kwargs["dpi"] = dpi
        try:
            p.save(out, **kwargs)
            saved[fmt] = out
        except Exception as exc:
            print(f"Could not save {fmt.upper()} for {output_stem.name}: {exc}")
    return saved

def build_effect_pair_table(manifest):
    pair_frames = []

    for day, day_manifest in manifest.groupby("day", sort=False):
        records = day_manifest.to_dict("records")

        if len(records) < 2:
            print(f"{day}: fewer than 2 biological samples with beta files; skipping.")
            continue

        for rec_a, rec_b in itertools.combinations(records, 2):
            sample_a = str(rec_a["biological_sample"])
            sample_b = str(rec_b["biological_sample"])

            beta_a = read_beta_wide_as_long(
                rec_a["beta_path"],
                day=day,
                biological_sample=sample_a,
            ).rename(columns={"beta": "beta_x"})

            beta_b = read_beta_wide_as_long(
                rec_b["beta_path"],
                day=day,
                biological_sample=sample_b,
            ).rename(columns={"beta": "beta_y"})

            merged = beta_a[["target_name", "program_name", "beta_x"]].merge(
                beta_b[["target_name", "program_name", "beta_y"]],
                on=["target_name", "program_name"],
                how="inner",
            )

            # Remove NTC target rows; these are not real target effects.
            merged = merged[
                ~merged["target_name"].astype(str).str.strip().str.lower().isin(NTC_KEY_SET)
            ].copy()

            merged = merged[
                np.isfinite(merged["beta_x"])
                & np.isfinite(merged["beta_y"])
            ].copy()

            if merged.empty:
                print(f"{day}: {sample_a} vs {sample_b} produced no matched finite beta pairs.")
                continue

            merged["day"] = day
            merged["replicate_x"] = sample_a
            merged["replicate_y"] = sample_b
            merged["replicate_pair"] = f"{sample_a} vs {sample_b}"
            pair_frames.append(merged)

    if not pair_frames:
        raise ValueError("No replicate pairs available for plotting.")

    return pd.concat(pair_frames, ignore_index=True, sort=False)

def summarize_effect_pairs(effects):
    summary_rows = []
    for (day, pair), group in effects.groupby(["day", "replicate_pair"], sort=False):
        pearson_r = group["beta_x"].corr(group["beta_y"], method="pearson")
        spearman_r = group["beta_x"].corr(group["beta_y"], method="spearman")

        summary_rows.append({
            "day": day,
            "replicate_pair": pair,
            "n_target_program_pairs": int(len(group)),
            "pearson_r": pearson_r,
            "spearman_r": spearman_r,
        })

    summary = pd.DataFrame(summary_rows)
    summary["facet_label"] = summary.apply(
        lambda r: (
            f"{r['day']}: {r['replicate_pair']}\n"
            f"Pearson r={format_correlation(r['pearson_r'])}; "
            f"Spearman rho={format_correlation(r['spearman_r'])}; "
            f"n={int(r['n_target_program_pairs']):,}"
        ),
        axis=1,
    )
    return summary

def plot_effect_correlation(plot_df, *, title, subtitle, ncol=2):
    axis_limit = float(
        np.nanmax(np.abs(plot_df[["beta_x", "beta_y"]].to_numpy(dtype=float)))
    )
    axis_limit = axis_limit * 1.05 if axis_limit > 0 else 1.0

    n_facets = int(plot_df["facet_label"].nunique())
    height = max(4.5, 3.15 * np.ceil(n_facets / ncol))
    width = 8.5 if ncol == 2 else 5.2

    p = (
        ggplot(plot_df, aes("beta_x", "beta_y"))
        + geom_hline(yintercept=0, color="#D5DCE3", size=0.35)
        + geom_vline(xintercept=0, color="#D5DCE3", size=0.35)
        + geom_abline(
            intercept=0,
            slope=1,
            linetype="dashed",
            color="#AEB7C2",
            size=0.7,
        )
        + geom_point(
            size=0.75,
            alpha=0.45,
            color="#17365D",
            stroke=0,
        )
        + facet_wrap("~facet_label", ncol=ncol)
        + scale_x_continuous(limits=(-axis_limit, axis_limit))
        + scale_y_continuous(limits=(-axis_limit, axis_limit))
        + coord_equal()
        + labs(
            title=title,
            subtitle=subtitle,
            x="DE program beta effect in biological replicate 1",
            y="DE program beta effect in biological replicate 2",
        )
        + theme_bw(base_size=10)
        + theme(
            figure_size=(width, height),
            plot_title=element_text(size=13, weight="bold", ha="center"),
            plot_subtitle=element_text(size=9, ha="center"),
            strip_text=element_text(size=8.5),
            axis_title=element_text(size=10),
            axis_text=element_text(size=8),
            panel_grid_major=element_line(color="#E8EDF3", size=0.45),
            panel_grid_minor=element_blank(),
            panel_border=element_rect(color="#182635", size=0.65),
        )
    )
    return p, width, float(height)
