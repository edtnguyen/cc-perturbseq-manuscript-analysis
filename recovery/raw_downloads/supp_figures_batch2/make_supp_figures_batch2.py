
import os
import zipfile
from pathlib import Path

os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib"

import numpy as np
import pandas as pd
from plotnine import (
    ggplot, aes, geom_tile, geom_text, geom_hline, geom_vline,
    geom_point, geom_line,
    scale_fill_gradient, scale_fill_gradient2, scale_color_manual, scale_y_discrete,
    theme_minimal, theme, labs, element_text, element_blank, element_rect
)

BASE = Path("/mnt/data")
OUTDIR = BASE / "supp_figures_batch2"
OUTDIR.mkdir(exist_ok=True, parents=True)

program_long = pd.read_csv(BASE / "program_activity_long_v4_updated_colors_labels.csv")
d2 = pd.read_csv(BASE / "target_D2_state_counts_and_fractions.csv")
d3 = pd.read_csv(BASE / "target_D3_state_counts_and_fractions.csv")
qc = pd.read_csv(BASE / "target_qc_counts.csv")

def save_plot(p, stem, width, height, dpi=320):
    png = OUTDIR / f"{stem}.png"
    svg = OUTDIR / f"{stem}.svg"
    pdf = OUTDIR / f"{stem}.pdf"
    p.save(filename=str(png), width=width, height=height, dpi=dpi, verbose=False)
    p.save(filename=str(svg), width=width, height=height, verbose=False)
    p.save(filename=str(pdf), width=width, height=height, verbose=False)
    return png, svg, pdf

def pub_theme(base_size=10, fig_w=6, fig_h=5, rotate_x=False):
    return theme_minimal(base_size=base_size) + theme(
        figure_size=(fig_w, fig_h),
        panel_grid_major=element_blank(),
        panel_grid_minor=element_blank(),
        panel_background=element_rect(fill="white", color="white"),
        plot_background=element_rect(fill="white", color="white"),
        axis_title=element_blank(),
        axis_text_x=element_text(
            size=base_size,
            rotation=45 if rotate_x else 0,
            ha="right" if rotate_x else "center",
            va="top" if rotate_x else "center",
            weight="bold"
        ),
        axis_text_y=element_text(size=base_size, ha="right"),
        axis_ticks=element_blank(),
        legend_title=element_text(size=base_size, weight="bold"),
        legend_text=element_text(size=base_size - 1),
        plot_title=element_text(size=base_size + 1, weight="bold"),
        plot_subtitle=element_text(size=base_size - 1),
        strip_text=element_text(size=base_size, weight="bold"),
    )

HIPPO = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]
TGFB = ["BMPR1A", "ACVR2A", "ACVR2B", "SMAD1", "SMAD4", "SMAD5"]
TARGET_ORDER = HIPPO + TGFB

# Panel 1: temporal landscape
program_df = program_long.copy()
program_df["Day"] = pd.Categorical(program_df["Day"], categories=["D0", "D1", "D2", "D3"], ordered=True)
program_df["program_id"] = program_df["program_name"].astype(str)
program_df["program_order"] = pd.to_numeric(program_df["ypos"], errors="coerce")
program_df["program_label_clean"] = program_df["program_label"].astype(str)

row_df = (
    program_df[["program_id", "program_label_clean", "program_order"]]
    .drop_duplicates()
    .sort_values(["program_order", "program_id"])
)
row_order = row_df["program_id"].tolist()
label_map = dict(zip(row_df["program_id"], row_df["program_label_clean"]))

plot_df = program_df[["program_id", "program_label_clean", "program_order", "Day", "mean_score"]].copy()
plot_df["program_id"] = pd.Categorical(plot_df["program_id"].astype(str), categories=row_order[::-1], ordered=True)

if plot_df["mean_score"].min() < 0:
    land_scale = scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C", midpoint=0, name="Mean score")
else:
    land_scale = scale_fill_gradient(low="white", high="#1F4E79", name="Mean score")

p_land = (
    ggplot(plot_df, aes("Day", "program_id", fill="mean_score"))
    + geom_tile(color="white", size=0.30)
    + land_scale
    + scale_y_discrete(labels=label_map)
    + labs(title="Temporal activity landscape of all cNMF programs")
    + pub_theme(base_size=9, fig_w=5.3, fig_h=max(9.5, 0.18 * len(row_order) + 1.2), rotate_x=False)
)
save_plot(
    p_land,
    "supp_panel_temporal_program_landscape_all50",
    width=5.3,
    height=max(9.5, 0.18 * len(row_order) + 1.2),
)
plot_df.to_csv(OUTDIR / "source_temporal_program_landscape_all50.csv", index=False)

# Panels 2-4: composition heatmaps
def make_composition_table(df, phase):
    frac_col = "fraction_among_all_D2" if phase == "D2" else "fraction_among_all_D3"
    state_order = ["Pre-EC / P15", "Lat-Meso / P17", "Other D2"] if phase == "D2" else ["artEC", "Fibro-Mesen", "Other D3"]

    x = df.copy()
    x = x[(x["control_mode"] == "strict_nt") & (x["target"].isin(TARGET_ORDER)) & (x["state"].isin(state_order))].copy()

    target_df = x[x["group"] == "target"][["target", "state", frac_col]].rename(columns={frac_col: "target_fraction"})
    ctrl_df = x[x["group"] == "control"][["target", "state", frac_col]].rename(columns={frac_col: "control_fraction"})
    merged = target_df.merge(ctrl_df, on=["target", "state"], how="left")
    merged["phase"] = phase
    merged["target_pct"] = 100 * merged["target_fraction"]
    merged["control_pct"] = 100 * merged["control_fraction"]
    merged["delta_pp"] = merged["target_pct"] - merged["control_pct"]
    merged["label_text"] = merged["target_pct"].map(lambda v: f"{v:.1f}")
    merged["category"] = np.where(merged["target"].isin(HIPPO), "Hippo", "TGFβ/BMP")

    present_targets = [t for t in TARGET_ORDER if t in merged["target"].unique()]
    merged["target"] = pd.Categorical(merged["target"], categories=present_targets[::-1], ordered=True)
    merged["state"] = pd.Categorical(merged["state"], categories=state_order, ordered=True)
    return merged

d2_comp = make_composition_table(d2, "D2")
d3_comp = make_composition_table(d3, "D3")

def hline_positions(targets_present):
    current = [t for t in TARGET_ORDER if t in targets_present]
    current_rev = current[::-1]
    pos = []
    last_cat = None
    for i, t in enumerate(current_rev):
        cat = "Hippo" if t in HIPPO else "TGFβ/BMP"
        if last_cat is not None and cat != last_cat:
            pos.append(i + 0.5)
        last_cat = cat
    return pos

def comp_heatmap(df, title, stem):
    if df.empty:
        return
    mx = max(2.0, float(np.nanmax(np.abs(df["delta_pp"]))))
    width = 4.8
    height = max(4.5, 0.38 * df["target"].nunique() + 1.2)
    p = (
        ggplot(df, aes("state", "target", fill="delta_pp"))
        + geom_tile(color="white", size=0.45)
        + geom_text(aes(label="label_text"), size=8, color="black")
        + scale_fill_gradient2(
            low="#2C7BB6", mid="white", high="#D7191C",
            midpoint=0, limits=(-mx, mx), name="Δ vs control\n(pp)"
        )
        + labs(title=title)
        + pub_theme(base_size=10, fig_w=width, fig_h=height, rotate_x=True)
    )
    for y in hline_positions(set(df["target"].astype(str))):
        p = p + geom_hline(yintercept=y, color="#666666", size=0.55)
    save_plot(p, stem, width=width, height=height)

comp_heatmap(d2_comp, "Selected target composition shifts at D2", "supp_panel_selected_targets_D2_state_composition")
comp_heatmap(d3_comp, "Selected target composition shifts at D3", "supp_panel_selected_targets_D3_state_composition")
d2_comp.to_csv(OUTDIR / "source_selected_targets_D2_state_composition.csv", index=False)
d3_comp.to_csv(OUTDIR / "source_selected_targets_D3_state_composition.csv", index=False)

combo = pd.concat([d2_comp, d3_comp], ignore_index=True)
combo["state_display"] = combo["phase"].astype(str) + " | " + combo["state"].astype(str)
state_order_combo = ["D2 | Pre-EC / P15", "D2 | Lat-Meso / P17", "D2 | Other D2",
                     "D3 | artEC", "D3 | Fibro-Mesen", "D3 | Other D3"]
combo["state_display"] = pd.Categorical(combo["state_display"], categories=state_order_combo, ordered=True)
combo_targets = [t for t in TARGET_ORDER if t in combo["target"].astype(str).unique()]
combo["target"] = pd.Categorical(combo["target"], categories=combo_targets[::-1], ordered=True)

mx = max(2.0, float(np.nanmax(np.abs(combo["delta_pp"]))))
p_combo = (
    ggplot(combo, aes("state_display", "target", fill="delta_pp"))
    + geom_tile(color="white", size=0.45)
    + geom_text(aes(label="label_text"), size=7.5, color="black")
    + scale_fill_gradient2(
        low="#2C7BB6", mid="white", high="#D7191C",
        midpoint=0, limits=(-mx, mx), name="Δ vs control\n(pp)"
    )
    + geom_vline(xintercept=3.5, color="#444444", size=0.75)
    + labs(title="Selected target composition shifts across D2 and D3")
    + pub_theme(
        base_size=10,
        fig_w=7.2,
        fig_h=max(4.8, 0.38 * combo["target"].nunique() + 1.2),
        rotate_x=True
    )
)
for y in hline_positions(set(combo["target"].astype(str))):
    p_combo = p_combo + geom_hline(yintercept=y, color="#666666", size=0.55)
save_plot(
    p_combo,
    "supp_panel_selected_targets_D2_D3_state_composition_combined",
    width=7.2,
    height=max(4.8, 0.38 * combo["target"].nunique() + 1.2),
)
combo.to_csv(OUTDIR / "source_selected_targets_D2_D3_state_composition_combined.csv", index=False)

# Panel 5: counts
qc_plot = qc[qc["target"].isin(TARGET_ORDER)].copy()
qc_long = pd.DataFrame({
    "target": list(qc_plot["target"]) * 2,
    "phase": ["D2"] * len(qc_plot) + ["D3"] * len(qc_plot),
    "n_cells": list(qc_plot["n_D2_carrier_cells"]) + list(qc_plot["n_D3_carrier_cells"])
})
qc_targets_present = [t for t in TARGET_ORDER if t in qc_long["target"].unique()]
qc_long["target"] = pd.Categorical(qc_long["target"], categories=qc_targets_present[::-1], ordered=True)

p_qc = (
    ggplot(qc_long, aes("n_cells", "target", color="phase"))
    + geom_line(aes(group="target"), color="#A0A0A0", size=0.5)
    + geom_point(size=2.6)
    + scale_color_manual(values={"D2": "#3F7FBF", "D3": "#C85A54"})
    + labs(title="Carrier-cell counts for selected targets", color="Day")
    + pub_theme(base_size=10, fig_w=5.4, fig_h=max(4.5, 0.38 * qc_long["target"].nunique() + 1.1), rotate_x=False)
)
for y in hline_positions(set(qc_long["target"].astype(str))):
    p_qc = p_qc + geom_hline(yintercept=y, color="#666666", size=0.55)
save_plot(
    p_qc,
    "supp_panel_selected_targets_carrier_cell_counts",
    width=5.4,
    height=max(4.5, 0.38 * qc_long["target"].nunique() + 1.1),
)
qc_long.to_csv(OUTDIR / "source_selected_targets_carrier_cell_counts.csv", index=False)
