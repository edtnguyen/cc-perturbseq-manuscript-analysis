# %% [markdown]
# # Fig. 3b and Fig. 4a Hippo/program panels
#
# This notebook/script builds a set of Illustrator-editable SVG panels for the revised manuscript outline:
#
# **Fig. 3b**: regulator-effect correlations linking D2 branch programs to D3 EC/alternative endpoint programs.  
# **Fig. 4a**: Hippo-specific perturb-seq KD and small-molecule effects on selected D2/D3 programs, plus collapsed marker category effects.
#
# The plotting code uses source tables only. It does not require loading large `.h5ad` files.

# %%
from __future__ import annotations

from pathlib import Path
import math
import re
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

try:
    from openpyxl import load_workbook
except Exception:
    load_workbook = None

mpl.rcParams.update({
    "font.family": "DejaVu Sans",  # change to Arial on Sherlock if available
    "font.size": 7,
    "figure.dpi": 160,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",  # keep text editable in Illustrator
})

# %% [markdown]
# ## Inputs and selected program sets

# %%
# Output directory
OUT_DIR = Path("/scratch/users/edtng/20260325_Hippo_timecourse/fig3b_fig4a_hippo_panels")
LOCAL_OUT_DIR = Path("/mnt/data/fig3b_fig4a_hippo_panels")
if Path("/scratch/users/edtng").exists():
    FIG_DIR = OUT_DIR
else:
    FIG_DIR = LOCAL_OUT_DIR
FIG_DIR.mkdir(parents=True, exist_ok=True)

# cNMF workbook with Program Loadings and Perturbation Association 1 sheets.
CNMF_WORKBOOK_CANDIDATES = [
    Path("/scratch/users/edtng/cNMF_50_0_2_032526_6co.xlsx"),
    Path("/mnt/data/cNMF_50_0_2_032526_6co.xlsx"),
]

# D3 EC endpoint percentage by target.
D3_EC_PERCENT_CANDIDATES = [
    Path("/scratch/users/edtng/20260325_Hippo_timecourse/d3_EC_percent_by_target_biorep_wide.csv"),
    Path("/mnt/data/d3_EC_percent_by_target_biorep_wide.csv"),
]

# Small-molecule program effect tables. Prefer CLR-matched effects if present.
SMALL_MOL_CLR_EFFECT_CANDIDATES = [
    Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_matched_DE/program_score_DE_clr/small_molecule_projected_cNMF_CLR_matched_effects.tsv"),
    Path("/mnt/data/small_molecule_projected_cNMF_CLR_matched_effects.tsv"),
]

# Fallback normalized-usage effect table if CLR table is not available.
SMALL_MOL_NORM_EFFECT_CANDIDATES = [
    Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_matched_DE/program_score_DE/projected_cNMF_program_score_matched_contrasts.tsv"),
    Path("/mnt/data/tmp_fig_inputs/hippo_publication_plots_package/source_tables/projected_cNMF_program_score_matched_contrasts.tsv"),
]

# Marker collapsed category tables from the marker-gene analysis.
MARKER_CATEGORY_SUMMARY_CANDIDATES = [
    Path("/scratch/users/edtng/20260325_Hippo_timecourse/marker_gene_control_centered_plots/F_latest_collapsed_marker_category_delta_summary.tsv"),
    Path("/mnt/data/F_latest_collapsed_marker_category_delta_summary.tsv"),
]
MARKER_CATEGORY_REPLICATES_CANDIDATES = [
    Path("/scratch/users/edtng/20260325_Hippo_timecourse/marker_gene_control_centered_plots/F_latest_collapsed_marker_category_delta_replicates.tsv"),
    Path("/mnt/data/F_latest_collapsed_marker_category_delta_replicates.tsv"),
]

# Program groups for the panels. Adjust these as annotations evolve.
D2_BRANCH_PROGRAMS = ["P13", "P30", "P14", "P15", "P17", "P19", "P23"]
D2_CORE_PROGRAMS_FOR_HEATMAP = ["P14", "P15", "P17", "P19", "P23"]
D3_EC_PROGRAMS = ["P2", "P3", "P8", "P11", "P42"]
D3_ALT_PROGRAMS = ["P16", "P24", "P27", "P33", "P34", "P35"]
SELECTED_PROGRAMS = D2_CORE_PROGRAMS_FOR_HEATMAP + D3_EC_PROGRAMS + D3_ALT_PROGRAMS
KEY_PROGRAMS = ["P15", "P17", "P34", "P2", "P23", "P24"]

HIPPO_KD_TARGETS = ["LATS2", "NF2", "SAV1", "TAOK1", "WWC3", "AMOTL2"]
PERTURBSEQ_EFFECT_COL = "approx_log2FC"  # use beta_hat instead if that is your final programDE effect column
CORRELATION_METHOD = "spearman"  # "spearman" or "pearson"

# Small molecule contrasts shown in Fig. 4a.
SMALL_MOL_CONTRAST_ORDER = [
    "D2_LM+TDI",
    "D2_LM+IAG",
    "D3_AE+TDI",
    "D3_AE+IAG",
    "D2_LM+TDI->D3_AE+TDI",
    "D2_LM+IAG->D3_AE+IAG",
]

print("Output directory:", FIG_DIR)

# %% [markdown]
# ## Helper functions

# %%
def find_existing(candidates: Iterable[Path], required: bool = True, label: str = "file") -> Path | None:
    for p in candidates:
        if p.exists():
            return p
    if required:
        raise FileNotFoundError(f"Could not find {label}. Tried: {[str(p) for p in candidates]}")
    return None


def read_xlsx_sheet_fast(path: Path, sheet_name: str) -> pd.DataFrame:
    """Read one xlsx sheet using openpyxl read-only mode; drops the leading unnamed index column."""
    if load_workbook is None:
        return pd.read_excel(path, sheet_name=sheet_name)
    wb = load_workbook(path, read_only=True, data_only=True)
    try:
        ws = wb[sheet_name]
        rows = ws.iter_rows(values_only=True)
        header = list(next(rows))
        data = list(rows)
    finally:
        wb.close()
    df = pd.DataFrame(data, columns=header)
    # Drop completely empty rows and columns, and the leading unnamed index column if present.
    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    unnamed = [c for c in df.columns if c is None or str(c).startswith("Unnamed")]
    if unnamed:
        df = df.drop(columns=unnamed)
    return df


def read_perturbation_assoc_subset(path: Path, sheet_name: str, samples, programs, effect_col: str) -> pd.DataFrame:
    """Stream-filter Perturbation Association sheet to selected samples/programs."""
    if load_workbook is None:
        df = pd.read_excel(path, sheet_name=sheet_name)
        df = df.copy()
        unnamed = [c for c in df.columns if c is None or str(c).startswith("Unnamed")]
        if unnamed:
            df = df.drop(columns=unnamed)
        df["program"] = "P" + df["program_name"].map(program_to_int).astype(str)
        df[effect_col] = pd.to_numeric(df[effect_col], errors="coerce")
        return df[df["Sample"].astype(str).isin(samples) & df["program"].isin(programs)]

    sample_set = set(map(str, samples))
    program_set = set(map(str, programs))
    wb = load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet_name]
    rows = ws.iter_rows(values_only=True)
    header = list(next(rows))
    col = {str(name): idx for idx, name in enumerate(header) if name is not None}
    required = ["target_name", "program_name", "Sample", effect_col]
    missing = [x for x in required if x not in col]
    if missing:
        raise ValueError(f"Missing required columns {missing} from {sheet_name}; header={header}")

    out = []
    for row in rows:
        if row is None or all(v is None for v in row):
            continue
        sample = str(row[col["Sample"]])
        if sample not in sample_set:
            continue
        try:
            pid = program_to_int(row[col["program_name"]])
        except Exception:
            continue
        program = f"P{pid}"
        if program not in program_set:
            continue
        out.append({
            "target_name": str(row[col["target_name"]]),
            "program_name": str(row[col["program_name"]]),
            "program": program,
            "program_id": pid,
            "Sample": sample,
            effect_col: row[col[effect_col]],
        })
    df = pd.DataFrame(out)
    if len(df):
        df[effect_col] = pd.to_numeric(df[effect_col], errors="coerce")
    return df


def program_to_int(program: str | int | float) -> int:
    s = str(program).replace("P", "").strip()
    return int(float(s))


def program_sort_key(program: str) -> int:
    return program_to_int(program)


def label_program(program: str, labels: dict[str, str] | None = None, max_chars: int = 38) -> str:
    if labels is None or program not in labels or not labels[program]:
        return program
    label = str(labels[program]).strip()
    out = f"{program}: {label}"
    if len(out) > max_chars:
        out = out[: max_chars - 1].rstrip() + "…"
    return out


def add_program_groups_to_ylabels(programs, labels, group_map):
    out = []
    last_group = None
    for p in programs:
        group = group_map.get(p, "")
        prefix = ""
        if group != last_group:
            prefix = f"{group}\n" if group else ""
        out.append(prefix + label_program(p, labels, max_chars=34))
        last_group = group
    return out


def corr_pair(a: pd.Series, b: pd.Series, method: str = "spearman"):
    m = pd.concat([a, b], axis=1, join="inner").dropna()
    if m.shape[0] < 3:
        return np.nan, np.nan, m.shape[0]
    if method == "spearman":
        r, p = spearmanr(m.iloc[:, 0], m.iloc[:, 1])
    elif method == "pearson":
        r, p = pearsonr(m.iloc[:, 0], m.iloc[:, 1])
    else:
        raise ValueError(method)
    return float(r), float(p), int(m.shape[0])


def effect_vector(prog_df: pd.DataFrame, program: str, day: str, effect_col: str = PERTURBSEQ_EFFECT_COL) -> pd.Series:
    pid = program_to_int(program)
    sub = prog_df.loc[(prog_df["program"] == f"P{pid}") & (prog_df["Sample"].astype(str) == day), ["target_name", effect_col]]
    return sub.dropna().groupby("target_name")[effect_col].mean().rename(f"{program}_{day}")


def plot_matrix_heatmap(df: pd.DataFrame, title: str, out_prefix: Path, cmap="RdBu_r", center_zero=True,
                        annotate=True, x_rotation=45, y_fontsize=7, vlim=None, figsize=None):
    mat = df.to_numpy(dtype=float)
    n_rows, n_cols = mat.shape
    if vlim is not None:
        vmin, vmax = (-vlim, vlim) if center_zero else (np.nanmin(mat), vlim)
    elif center_zero:
        vmax = np.nanmax(np.abs(mat))
        vmin = -vmax
    else:
        vmin = np.nanmin(mat)
        vmax = np.nanmax(mat)
        if np.isclose(vmin, vmax):
            vmax = vmin + 1e-9

    if figsize is None:
        figsize = (max(5, 0.45 * n_cols + 2.5), max(4, 0.28 * n_rows + 1.5))

    fig, ax = plt.subplots(figsize=figsize)
    mesh = ax.pcolormesh(
        np.arange(n_cols + 1), np.arange(n_rows + 1), mat,
        cmap=cmap, vmin=vmin, vmax=vmax, shading="flat",
        linewidth=0, antialiased=False,
    )
    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_rows, 0)
    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_xticklabels(df.columns, rotation=x_rotation, ha="right" if x_rotation else "center")
    ax.set_yticklabels(df.index, fontsize=y_fontsize)
    ax.set_title(title)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if annotate:
        for i in range(n_rows):
            for j in range(n_cols):
                val = mat[i, j]
                if np.isfinite(val):
                    ax.text(j + 0.5, i + 0.5, f"{val:.2f}", ha="center", va="center", fontsize=5)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.ax.tick_params(labelsize=6)
    fig.tight_layout()
    for ext in ["svg", "pdf", "png"]:
        fig.savefig(out_prefix.with_suffix(f".{ext}"), bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close(fig)


def save_table(df: pd.DataFrame, name: str):
    path = FIG_DIR / name
    df.to_csv(path, sep="\t", index=True)
    return path

# %% [markdown]
# ## Load cNMF programDE workbook and endpoint data

# %%
cnmf_workbook = find_existing(CNMF_WORKBOOK_CANDIDATES, label="cNMF workbook")
d3_ec_path = find_existing(D3_EC_PERCENT_CANDIDATES, required=False, label="D3 EC endpoint table")
print("cNMF workbook:", cnmf_workbook)
print("D3 EC endpoint table:", d3_ec_path)

# Only D2/D3 and selected programs are needed for the panels.
assoc = read_perturbation_assoc_subset(
    cnmf_workbook,
    "Perturbation Association 1",
    samples=["D2", "D3"],
    programs=SELECTED_PROGRAMS + D2_BRANCH_PROGRAMS,
    effect_col=PERTURBSEQ_EFFECT_COL,
)
assoc = assoc.copy()
assoc["target_name"] = assoc["target_name"].astype(str)
assoc["Sample"] = assoc["Sample"].astype(str)
print("assoc subset shape:", assoc.shape)
print("samples:", sorted(assoc["Sample"].unique()))
print("programs:", sorted(assoc["program"].unique(), key=program_sort_key)[:5], "...", sorted(assoc["program"].unique(), key=program_sort_key)[-5:])

# Program labels from Program Loadings are not ideal; use manual descriptions here.
PROGRAM_LABELS = {
    "P13": "early/branch program",
    "P30": "early/branch program",
    "P14": "lateral mesoderm",
    "P15": "vascular progenitor",
    "P17": "epithelial/alt LM",
    "P19": "posterior mesoderm",
    "P23": "vascular LM / EC induction",
    "P2": "arterial endothelium",
    "P3": "activated endothelium",
    "P8": "EC adhesion/Notch",
    "P11": "angiogenic EC",
    "P42": "vascular signaling",
    "P16": "cytoskeletal/alt",
    "P24": "GATA4/RSPO2 alt",
    "P27": "mesenchymal ECM",
    "P33": "anterior-like mesoderm",
    "P34": "YAP/mesenchymal",
    "P35": "EMT/ECM",
}
PROGRAM_GROUP = {p: "D2 branch" for p in D2_CORE_PROGRAMS_FOR_HEATMAP}
PROGRAM_GROUP.update({p: "D3 EC" for p in D3_EC_PROGRAMS})
PROGRAM_GROUP.update({p: "D3 alternative" for p in D3_ALT_PROGRAMS})

# %% [markdown]
# ## Fig. 3b: D2 branch programs vs D3 EC/alternative programs and endpoint

# %%
# Build a D3 EC endpoint vector, if the endpoint table is available.
ec_endpoint = None
if d3_ec_path is not None:
    ec = pd.read_csv(d3_ec_path)
    # Weighted EC percentage across available D3 bioreps.
    ec["n_carrier_total"] = ec[[c for c in ec.columns if c.endswith("n_target_carrier_cells")]].sum(axis=1)
    ec["n_ec_total"] = ec[[c for c in ec.columns if c.endswith("n_EC_target_carrier_cells")]].sum(axis=1)
    ec["d3_ec_percent_weighted"] = 100 * ec["n_ec_total"] / ec["n_carrier_total"].replace(0, np.nan)
    ec_endpoint = ec.set_index("target")["d3_ec_percent_weighted"].rename("D3 EC %")
    ec[["target", "d3_ec_percent_weighted", "n_carrier_total", "n_ec_total"]].to_csv(
        FIG_DIR / "D3_EC_percent_endpoint_weighted_by_target.tsv", sep="\t", index=False
    )
    print("EC endpoint targets:", ec_endpoint.shape[0])

rows = []
fig3b_mat = pd.DataFrame(index=D2_BRANCH_PROGRAMS)
fig3b_p_mat = pd.DataFrame(index=D2_BRANCH_PROGRAMS)

fig3b_cols = D3_EC_PROGRAMS + D3_ALT_PROGRAMS + (["D3 EC %"] if ec_endpoint is not None else [])

for p_d2 in D2_BRANCH_PROGRAMS:
    a = effect_vector(assoc, p_d2, "D2")
    vals = []
    pvals = []
    for p_d3 in D3_EC_PROGRAMS + D3_ALT_PROGRAMS:
        b = effect_vector(assoc, p_d3, "D3")
        r, p, n = corr_pair(a, b, CORRELATION_METHOD)
        vals.append(r)
        pvals.append(p)
        rows.append({
            "source_day": "D2",
            "source_program": p_d2,
            "target_day": "D3",
            "target_program": p_d3,
            "target_type": "D3 EC program" if p_d3 in D3_EC_PROGRAMS else "D3 alternative program",
            "method": CORRELATION_METHOD,
            "correlation": r,
            "pvalue": p,
            "n_targets": n,
        })
    if ec_endpoint is not None:
        r, p, n = corr_pair(a, ec_endpoint, CORRELATION_METHOD)
        vals.append(r)
        pvals.append(p)
        rows.append({
            "source_day": "D2",
            "source_program": p_d2,
            "target_day": "D3",
            "target_program": "D3 EC %",
            "target_type": "D3 endpoint",
            "method": CORRELATION_METHOD,
            "correlation": r,
            "pvalue": p,
            "n_targets": n,
        })
    fig3b_mat.loc[p_d2, fig3b_cols] = vals
    fig3b_p_mat.loc[p_d2, fig3b_cols] = pvals

fig3b_edges = pd.DataFrame(rows)
fig3b_edges.to_csv(FIG_DIR / "Fig3b_D2_to_D3_program_endpoint_correlations.tsv", sep="\t", index=False)
fig3b_mat.to_csv(FIG_DIR / "Fig3b_D2_to_D3_program_endpoint_correlation_matrix.tsv", sep="\t")

# Heatmap labels.
fig3b_mat_labeled = fig3b_mat.copy()
fig3b_mat_labeled.index = [label_program(p, PROGRAM_LABELS, 32) for p in fig3b_mat.index]
fig3b_mat_labeled.columns = [label_program(c, PROGRAM_LABELS, 30) if c.startswith("P") else c for c in fig3b_mat.columns]

plot_matrix_heatmap(
    fig3b_mat_labeled,
    title=f"Fig. 3b candidate: D2 program effects correlated with D3 programs/endpoints\n{CORRELATION_METHOD.capitalize()} correlation across perturbation targets",
    out_prefix=FIG_DIR / "Fig3b_D2_branch_vs_D3_program_endpoint_correlation_heatmap",
    cmap="RdBu_r",
    center_zero=True,
    annotate=True,
    x_rotation=45,
    y_fontsize=7,
    figsize=(9.5, 4.2),
)

# %% [markdown]
# ## Optional Fig. 3b network/schematic from the same edge table

# %%
def plot_transition_network(edge_df: pd.DataFrame, out_prefix: Path, threshold=0.45):
    # Draw all biologically selected edges among P15/P17/P2/P34/D3 EC%, plus thresholded edges.
    forced = (
        edge_df["source_program"].isin(["P15", "P17"]) &
        edge_df["target_program"].isin(["P2", "P34", "D3 EC %"])
    ) | (
        edge_df["target_program"].isin(["P2", "P34", "D3 EC %"]) &
        (edge_df["correlation"].abs() >= 0.35)
    )
    plot_edges = edge_df.loc[(edge_df["correlation"].abs() >= threshold) | forced].copy()

    d2_nodes = D2_BRANCH_PROGRAMS
    d3_nodes = D3_EC_PROGRAMS + D3_ALT_PROGRAMS
    if "D3 EC %" in edge_df["target_program"].values:
        endpoint_nodes = ["D3 EC %"]
    else:
        endpoint_nodes = []

    pos = {}
    for i, node in enumerate(d2_nodes):
        pos[node] = (0, i)
    for i, node in enumerate(d3_nodes):
        pos[node] = (1.7, i * (len(d2_nodes)-1) / max(1, len(d3_nodes)-1))
    for i, node in enumerate(endpoint_nodes):
        pos[node] = (3.4, (len(d2_nodes)-1)/2)

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    # Edges
    for _, row in plot_edges.iterrows():
        src = row["source_program"]
        tgt = row["target_program"]
        if src not in pos or tgt not in pos:
            continue
        x1, y1 = pos[src]
        x2, y2 = pos[tgt]
        r = row["correlation"]
        color = "#c0392b" if r >= 0 else "#2166ac"
        lw = 0.5 + 3.0 * min(1, abs(r))
        ax.plot([x1, x2], [y1, y2], color=color, linewidth=lw, alpha=0.7)
        if src in ["P15", "P17"] and tgt in ["P2", "P34", "D3 EC %"]:
            ax.text((x1+x2)/2, (y1+y2)/2, f"r={r:.2f}", fontsize=5, ha="center", va="center")

    # Nodes
    for node, (x, y) in pos.items():
        if node in D2_BRANCH_PROGRAMS:
            face = "#f2f2f2"
        elif node in D3_EC_PROGRAMS:
            face = "#e6f2ff"
        elif node in D3_ALT_PROGRAMS:
            face = "#fff0e6"
        else:
            face = "#e8f5e9"
        circ = Circle((x, y), 0.12, facecolor=face, edgecolor="black", linewidth=0.8)
        ax.add_patch(circ)
        ax.text(x, y+0.19, label_program(node, PROGRAM_LABELS, 18), fontsize=6, ha="center", va="bottom")

    ax.text(0, -0.8, "D2 branch programs", ha="center", fontsize=8, weight="bold")
    ax.text(1.7, -0.8, "D3 programs", ha="center", fontsize=8, weight="bold")
    if endpoint_nodes:
        ax.text(3.4, -0.8, "Endpoint", ha="center", fontsize=8, weight="bold")
    ax.set_xlim(-0.5, 3.9)
    ax.set_ylim(len(d2_nodes)-0.2, -1.1)
    ax.axis("off")
    fig.tight_layout()
    for ext in ["svg", "pdf", "png"]:
        fig.savefig(out_prefix.with_suffix(f".{ext}"), dpi=300 if ext == "png" else None)
    plt.close(fig)

plot_transition_network(fig3b_edges, FIG_DIR / "Fig3b_selected_program_transition_network", threshold=0.50)

# %% [markdown]
# ## Fig. 4a-1: Hippo KD effects on selected programs in perturb-seq

# %%
hippo = assoc.loc[assoc["target_name"].isin(HIPPO_KD_TARGETS)].copy()
hippo = hippo.loc[hippo["Sample"].isin(["D2", "D3"])]

# Individual Hippo KD target heatmap: rows target/day, columns selected programs.
hippo_rows = []
for day in ["D2", "D3"]:
    for target in HIPPO_KD_TARGETS:
        sub = hippo.loc[(hippo["Sample"] == day) & (hippo["target_name"] == target)]
        vals = []
        for p in SELECTED_PROGRAMS:
            v = sub.loc[sub["program"] == p, PERTURBSEQ_EFFECT_COL]
            vals.append(float(v.mean()) if len(v) else np.nan)
        hippo_rows.append(pd.Series(vals, index=SELECTED_PROGRAMS, name=f"{day} {target}"))
hippo_individual_mat = pd.DataFrame(hippo_rows)
hippo_individual_mat.to_csv(FIG_DIR / "Fig4a_HippoKD_individual_targets_selected_program_effects.tsv", sep="\t")

# Mean Hippo KD effect across targets.
hippo_mean_rows = []
for day in ["D2", "D3"]:
    sub = hippo.loc[hippo["Sample"] == day]
    means = sub.groupby("program")[PERTURBSEQ_EFFECT_COL].mean()
    hippo_mean_rows.append(means.reindex(SELECTED_PROGRAMS).rename(f"Hippo KD mean {day}"))
hippo_mean_mat = pd.DataFrame(hippo_mean_rows)
hippo_mean_mat.to_csv(FIG_DIR / "Fig4a_HippoKD_mean_selected_program_effects.tsv", sep="\t")

selected_display_cols = [label_program(p, PROGRAM_LABELS, 26) for p in SELECTED_PROGRAMS]
hippo_indiv_plot = hippo_individual_mat.copy()
hippo_indiv_plot.columns = selected_display_cols

plot_matrix_heatmap(
    hippo_indiv_plot,
    title=f"Fig. 4a candidate: perturb-seq Hippo KD effects on selected programs\nEffect column: {PERTURBSEQ_EFFECT_COL}",
    out_prefix=FIG_DIR / "Fig4a_HippoKD_individual_targets_selected_program_effects_heatmap",
    cmap="RdBu_r",
    center_zero=True,
    annotate=False,
    x_rotation=50,
    y_fontsize=6,
    figsize=(11, 5.5),
)

hippo_mean_plot = hippo_mean_mat.copy()
hippo_mean_plot.columns = selected_display_cols
plot_matrix_heatmap(
    hippo_mean_plot,
    title=f"Fig. 4a candidate: mean Hippo KD effects on selected programs\nAveraged across {len(HIPPO_KD_TARGETS)} Hippo targets",
    out_prefix=FIG_DIR / "Fig4a_HippoKD_mean_selected_program_effects_heatmap",
    cmap="RdBu_r",
    center_zero=True,
    annotate=True,
    x_rotation=50,
    y_fontsize=7,
    figsize=(11, 2.6),
)

# %% [markdown]
# ## Fig. 4a-2: Small-molecule effects on the same selected programs

# %%
sm_effect_path = find_existing(SMALL_MOL_CLR_EFFECT_CANDIDATES, required=False, label="small-molecule CLR effects")
sm_effect_source = "CLR"
if sm_effect_path is None:
    sm_effect_path = find_existing(SMALL_MOL_NORM_EFFECT_CANDIDATES, required=False, label="small-molecule normalized usage effects")
    sm_effect_source = "normalized_usage_difference"
print("Small-molecule effect source:", sm_effect_source if sm_effect_path is not None else "not found")
print("Small-molecule effect table:", sm_effect_path)

if sm_effect_path is not None:
    sm = pd.read_csv(sm_effect_path, sep="\t")
    # Harmonize columns.
    if sm_effect_source == "CLR":
        effect_col = "mean_diff_clr"
        if "treated_condition" not in sm.columns and "contrast" in sm.columns:
            sm["treated_condition"] = sm["contrast"]
    else:
        candidates = ["diff_treated_minus_control", "mean_diff", "mean_diff_usage", "mean_paired_diff"]
        effect_col = next((c for c in candidates if c in sm.columns), None)
        if effect_col is None:
            raise ValueError(f"Could not infer small-molecule effect column. Columns: {list(sm.columns)}")

    if effect_col not in sm.columns:
        raise ValueError(f"Expected effect column {effect_col!r} in {sm_effect_path}. Columns: {list(sm.columns)}")

    if "treated_condition" not in sm.columns:
        if "contrast" in sm.columns:
            sm["treated_condition"] = sm["contrast"]
        elif "contrast_name" in sm.columns:
            sm["treated_condition"] = sm["contrast_name"]
        else:
            raise ValueError("Could not identify treated condition column in small-molecule effect table")

    sm["treated_condition"] = sm["treated_condition"].astype(str).map(lambda x: {
        "TDI_at_D2_LM": "D2_LM+TDI",
        "IAG_at_D2_LM": "D2_LM+IAG",
        "TDI_at_D3_AE": "D3_AE+TDI",
        "IAG_at_D3_AE": "D3_AE+IAG",
        "TDI_D2_to_D3_vs_D3_AE": "D2_LM+TDI->D3_AE+TDI",
        "IAG_D2_to_D3_vs_D3_AE": "D2_LM+IAG->D3_AE+IAG",
    }.get(x, x))

    sm_rows = []
    for treated in SMALL_MOL_CONTRAST_ORDER:
        sub = sm.loc[sm["treated_condition"] == treated]
        vals = []
        for p in SELECTED_PROGRAMS:
            v = sub.loc[sub["program"].astype(str) == p, effect_col]
            vals.append(float(v.mean()) if len(v) else np.nan)
        sm_rows.append(pd.Series(vals, index=SELECTED_PROGRAMS, name=treated))
    sm_mat = pd.DataFrame(sm_rows)
    sm_mat.to_csv(FIG_DIR / f"Fig4a_small_molecule_{sm_effect_source}_selected_program_effects.tsv", sep="\t")

    sm_plot = sm_mat.copy()
    sm_plot.columns = selected_display_cols
    plot_matrix_heatmap(
        sm_plot,
        title=f"Fig. 4a candidate: TDI/IAG effects on selected projected programs\nScale: {sm_effect_source}; values are treated - matched control",
        out_prefix=FIG_DIR / f"Fig4a_small_molecule_{sm_effect_source}_selected_program_effects_heatmap",
        cmap="RdBu_r",
        center_zero=True,
        annotate=True,
        x_rotation=50,
        y_fontsize=7,
        figsize=(11, 4.2),
    )
else:
    print("Small-molecule program-effect table not found; skipped Fig4a TDI/IAG selected-program heatmap.")

# %% [markdown]
# ## Fig. 4a-3: Collapsed marker category treatment effect plot
#
# This uses the already computed control-centered marker module category table if present.

# %%
marker_summary_path = find_existing(MARKER_CATEGORY_SUMMARY_CANDIDATES, required=False, label="marker category summary")
marker_rep_path = find_existing(MARKER_CATEGORY_REPLICATES_CANDIDATES, required=False, label="marker category replicate table")
print("marker summary:", marker_summary_path)
print("marker replicate:", marker_rep_path)

if marker_summary_path is not None and marker_rep_path is not None:
    marker_summary = pd.read_csv(marker_summary_path, sep="\t")
    marker_rep = pd.read_csv(marker_rep_path, sep="\t")

    # Normalize contrast names for order.
    contrast_order_short = ["D2 TDI", "D2 IAG", "D3 TDI", "D3 IAG", "D2→D3 TDI", "D2→D3 IAG"]
    category_order = ["EC/progenitor markers", "Alternative/YAP markers"]
    category_short = {"EC/progenitor markers": "EC / progenitor", "Alternative/YAP markers": "Alternative / YAP"}

    marker_rep = marker_rep.loc[marker_rep["contrast"].isin(contrast_order_short)].copy()
    marker_summary = marker_summary.loc[marker_summary["contrast"].isin(contrast_order_short)].copy()

    x_positions = {}
    xs = []
    labels = []
    spacing = 1.1
    dx = 0.18
    for i, contrast in enumerate(contrast_order_short):
        for cat_j, cat in enumerate(category_order):
            x = i * spacing + (-dx if cat_j == 0 else dx)
            x_positions[(contrast, cat)] = x
        xs.append(i * spacing)
        labels.append(contrast)

    fig, ax = plt.subplots(figsize=(7.5, 3.2))
    # replicate points
    for _, row in marker_rep.iterrows():
        key = (row["contrast"], row["category"])
        if key not in x_positions:
            continue
        x = x_positions[key]
        jitter = (int(row["replicate"]) - 2) * 0.025 if "replicate" in row and pd.notna(row["replicate"]) else 0
        marker = "o" if row["category"] == "EC/progenitor markers" else "s"
        ax.scatter(x + jitter, row["mean_delta_log2cpm"], s=18, marker=marker, facecolors="white", edgecolors="black", linewidth=0.6, zorder=3)

    # means and connecting lines between categories within each contrast.
    for contrast in contrast_order_short:
        sub = marker_summary.loc[marker_summary["contrast"] == contrast]
        vals = {}
        for _, row in sub.iterrows():
            vals[row["category"]] = row["mean"]
        if all(cat in vals for cat in category_order):
            x1 = x_positions[(contrast, category_order[0])]
            x2 = x_positions[(contrast, category_order[1])]
            y1 = vals[category_order[0]]
            y2 = vals[category_order[1]]
            ax.plot([x1, x2], [y1, y2], color="0.35", linewidth=1.0, zorder=2)
            ax.scatter([x1, x2], [y1, y2], s=50, marker="D", color="black", zorder=4)

    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Matched-control Δ marker module log2(CPM+1)")
    ax.set_title("Fig. 4a candidate: collapsed marker categories\nEC/progenitor vs Alternative/YAP, matched-control centered")
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="white", markeredgecolor="black", label="EC / progenitor"),
        Line2D([0], [0], marker="s", linestyle="none", markerfacecolor="white", markeredgecolor="black", label="Alternative / YAP"),
        Line2D([0], [0], marker="D", linestyle="none", color="black", label="mean"),
    ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=6, loc="best")
    fig.tight_layout()
    for ext in ["svg", "pdf", "png"]:
        fig.savefig(FIG_DIR / f"Fig4a_marker_category_collapsed_control_centered.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
else:
    print("Skipping marker category plot because marker summary/replicate tables are not available.")

# %% [markdown]
# ## Optional combined preview: Fig. 4a triptych
#
# This combines mean Hippo KD, small-molecule effect, and marker category summary. It is intended as a draft layout only.

# %%
# Save a README/source summary.
summary = {
    "CNMF workbook": str(cnmf_workbook),
    "D3 EC endpoint table": str(d3_ec_path) if d3_ec_path is not None else "not found",
    "Small molecule effect source": sm_effect_source,
    "Small molecule effect table": str(sm_effect_path),
    "Marker summary": str(marker_summary_path) if marker_summary_path is not None else "not found",
    "Marker replicate": str(marker_rep_path) if marker_rep_path is not None else "not found",
    "Perturbseq effect column": PERTURBSEQ_EFFECT_COL,
    "Hippo targets": ";".join(HIPPO_KD_TARGETS),
    "Selected programs": ";".join(SELECTED_PROGRAMS),
}
with open(FIG_DIR / "figure_panel_generation_summary.txt", "w") as f:
    for k, v in summary.items():
        f.write(f"{k}: {v}\n")

print("Saved figure panels to:", FIG_DIR)
for p in sorted(FIG_DIR.glob("*")):
    print(p.name)
