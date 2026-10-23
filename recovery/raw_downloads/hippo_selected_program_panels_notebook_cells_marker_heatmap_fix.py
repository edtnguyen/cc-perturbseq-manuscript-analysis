# %% [markdown]
# # Hippo selected-program panels: D2/D3 programs, perturb-seq correlations, and marker genes
# 
# This notebook builds the additional Hippo validation panels discussed:
# 
# 1. **Selected-program effect heatmap** for D2 and D3 programs using matched small-molecule contrasts.
# 2. **Paired replicate plots** for selected day-relevant programs.
# 3. **Perturb-seq D2 × D3 program correlation matrix** using the `Perturbation Association 1` tab.
# 4. **Marker gene candidates** selected from top cNMF loading genes and concordant PyDESeq2 DE.
# 5. Optional **marker gene pseudobulk expression heatmap** if pseudobulk count files are available.
# 
# All SVG outputs keep text editable in Illustrator.

# %%
from pathlib import Path
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr, fisher_exact
from statsmodels.stats.multitest import multipletests

try:
    from IPython.display import display
except Exception:
    def display(x):
        try:
            print(x.to_string())
        except Exception:
            print(x)

warnings.filterwarnings("ignore")

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",  # Illustrator-editable text
})

# ---------- paths ----------
SEARCH_DIRS = [
    Path.cwd(),
    Path("/mnt/data"),
    Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_matched_DE/program_score_DE"),
    Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_CNMf_matched_DE/program_score_DE"),
    Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_matched_DE/pydeseq2_gene_DE"),
]

def find_file(filename, extra_dirs=None, required=True):
    dirs = list(SEARCH_DIRS)
    if extra_dirs:
        dirs = list(extra_dirs) + dirs
    for d in dirs:
        p = d / filename
        if p.exists():
            return p
    if required:
        raise FileNotFoundError(f"Could not find {filename} in: {[str(d) for d in dirs]}")
    return None

# Uploaded/current files or Sherlock paths.
PROGRAM_EFFECTS_PATH = find_file("projected_cNMF_program_score_matched_contrasts.tsv")
PAIRED_DIFFS_PATH = find_file("projected_cNMF_program_score_paired_replicate_differences.tsv")
CNMF_XLSX_PATH = find_file("cNMF_50_0_2_032526_6co.xlsx")

# Optional pseudobulk expression files for marker heatmap.
PSEUDOBULK_COUNTS_PATH = Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_matched_DE/pseudobulk_gene_counts.filtered.tsv.gz")
PSEUDOBULK_META_PATH = Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_matched_DE/pseudobulk_sample_metadata.filtered.tsv")

# Output folder. Use Sherlock folder when available, otherwise local folder.
DEFAULT_OUT = Path("/scratch/users/edtng/20260325_Hippo_timecourse/hippo_selected_program_panels")
OUT_DIR = DEFAULT_OUT if DEFAULT_OUT.parent.exists() else Path("hippo_selected_program_panels")
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("PROGRAM_EFFECTS_PATH:", PROGRAM_EFFECTS_PATH)
print("PAIRED_DIFFS_PATH:", PAIRED_DIFFS_PATH)
print("CNMF_XLSX_PATH:", CNMF_XLSX_PATH)
print("OUT_DIR:", OUT_DIR.resolve())

# %%
# ---------- program sets ----------
# Anchor programs remain the central claim.
ANCHOR_PROGRAMS = ["P15", "P17", "P34"]

# Day-relevant supporting programs from peak-day/activity classification.
D2_PROGRAMS = ["P14", "P15", "P17", "P19", "P23"]
D3_EC_PROGRAMS = ["P2", "P3", "P8", "P11", "P42"]
D3_ALT_PROGRAMS = ["P16", "P24", "P27", "P33", "P34", "P35"]

SELECTED_PROGRAMS = D2_PROGRAMS + D3_EC_PROGRAMS + D3_ALT_PROGRAMS
SELECTED_PROGRAMS = list(dict.fromkeys(SELECTED_PROGRAMS))

PROGRAM_SET_MAP = {}
for p in D2_PROGRAMS: PROGRAM_SET_MAP[p] = "D2 active"
for p in D3_EC_PROGRAMS: PROGRAM_SET_MAP[p] = "D3 EC / vascular"
for p in D3_ALT_PROGRAMS: PROGRAM_SET_MAP[p] = "D3 alternative / mesenchymal"

# Fallback labels from your existing peak-day program map.
FALLBACK_PROGRAM_LABELS = {
    "P14": "lateral mesoderm",
    "P15": "endothelial / pre-EC",
    "P17": "epithelial-like",
    "P19": "posterior mesoderm",
    "P23": "EC induction",
    "P2": "EC migratory",
    "P3": "EC activated",
    "P8": "EC adhesion/Notch",
    "P11": "EC angiogenic",
    "P42": "vascular signaling",
    "P16": "cytoskeletal",
    "P24": "HOX patterning",
    "P27": "mesenchymal ECM",
    "P33": "mesenchymal-guidance",
    "P34": "contractile/YAP",
    "P35": "EMT/ECM",
}

# Contrasts for D2/D3-focused Hippo panel, using original names.
FOCAL_CONTRASTS = [
    "TDI_at_D2_LM",
    "IAG_at_D2_LM",
    "TDI_at_D3_AE",
    "IAG_at_D3_AE",
    "TDI_D2_to_D3_vs_D3_AE",
    "IAG_D2_to_D3_vs_D3_AE",
]

FOCAL_CONTRAST_LABELS = {
    "TDI_at_D2_LM": "D2_LM+TDI",
    "IAG_at_D2_LM": "D2_LM+IAG",
    "TDI_at_D3_AE": "D3_AE+TDI",
    "IAG_at_D3_AE": "D3_AE+IAG",
    "TDI_D2_to_D3_vs_D3_AE": "D2_LM+TDI→D3_AE+TDI",
    "IAG_D2_to_D3_vs_D3_AE": "D2_LM+IAG→D3_AE+IAG",
}

# For marker gene selection: expected directions for TDI claim.
# Direction refers to treated vs matched control.
MARKER_RULES = [
    {"program": "P15", "contrast_name": "TDI_at_D2_LM", "direction": "down", "label": "P15 down in D2 TDI"},
    {"program": "P17", "contrast_name": "TDI_at_D2_LM", "direction": "up", "label": "P17 up in D2 TDI"},
    {"program": "P34", "contrast_name": "TDI_at_D3_AE", "direction": "up", "label": "P34 up in D3 TDI"},
    {"program": "P34", "contrast_name": "TDI_D2_to_D3_vs_D3_AE", "direction": "up", "label": "P34 up in D2→D3 TDI"},
]

# Optional broader rules for D3 EC/alternative sets.
# Uncomment/add rules if you want marker genes for additional programs.
# for p in D3_EC_PROGRAMS:
#     MARKER_RULES.append({"program": p, "contrast_name": "TDI_at_D3_AE", "direction": "down", "label": f"{p} down in D3 TDI"})
# for p in D3_ALT_PROGRAMS:
#     MARKER_RULES.append({"program": p, "contrast_name": "TDI_at_D3_AE", "direction": "up", "label": f"{p} up in D3 TDI"})

print("Selected programs:", SELECTED_PROGRAMS)
print("Focal contrasts:", FOCAL_CONTRASTS)

# %%
# ---------- helpers ----------
def program_num(program):
    m = re.fullmatch(r"P(\d+)", str(program))
    if m is None:
        raise ValueError(f"Program must be P<number>, got {program!r}")
    return int(m.group(1))


def compact_label(program):
    label = FALLBACK_PROGRAM_LABELS.get(program, program)
    return f"{program} — {label}"


def save_figure(fig, out_prefix: Path):
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for ext in ["svg", "pdf", "png"]:
        kwargs = {"bbox_inches": "tight"}
        if ext == "png":
            kwargs["dpi"] = 300
        fig.savefig(out_prefix.with_suffix(f".{ext}"), **kwargs)
    print("Saved:", out_prefix.with_suffix(".svg"))


def bh_adjust(pvals):
    pvals = pd.Series(pvals, dtype=float)
    ok = pvals.notna()
    out = pd.Series(np.nan, index=pvals.index, dtype=float)
    if ok.sum() > 0:
        out.loc[ok] = multipletests(pvals.loc[ok], method="fdr_bh")[1]
    return out


def zscore_rows(df):
    z = df.sub(df.mean(axis=1), axis=0)
    z = z.div(df.std(axis=1).replace(0, np.nan), axis=0)
    return z.fillna(0)


def plot_matrix_heatmap(df, title, out_prefix, cmap="RdBu_r", center_zero=True, annotate=False, vlim=None, x_rotation=45, y_fontsize=7):
    mat = df.to_numpy(dtype=float)
    n_rows, n_cols = mat.shape
    if vlim is not None:
        vmin, vmax = -vlim, vlim if center_zero else (np.nanmin(mat), vlim)
    elif center_zero:
        vmax = np.nanmax(np.abs(mat))
        vmin = -vmax
    else:
        vmin, vmax = np.nanmin(mat), np.nanmax(mat)
        if np.isclose(vmin, vmax):
            vmax = vmin + 1e-9

    fig_w = max(7, 0.75 * n_cols + 4)
    fig_h = max(5, 0.33 * n_rows + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    mesh = ax.pcolormesh(
        np.arange(n_cols + 1), np.arange(n_rows + 1), mat,
        cmap=cmap, vmin=vmin, vmax=vmax, shading="flat", linewidth=0.3,
        edgecolors="white", antialiased=False,
    )
    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_rows, 0)
    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_xticklabels(df.columns.tolist(), rotation=x_rotation, ha="right")
    ax.set_yticklabels(df.index.tolist(), fontsize=y_fontsize)
    ax.set_title(title)
    ax.tick_params(length=0)
    for spine in ax.spines.values(): spine.set_visible(False)

    if annotate:
        for i in range(n_rows):
            for j in range(n_cols):
                val = mat[i, j]
                if np.isfinite(val):
                    ax.text(j + 0.5, i + 0.5, f"{val:.2f}", ha="center", va="center", fontsize=6)

    cbar = fig.colorbar(mesh, ax=ax, pad=0.02, fraction=0.035)
    fig.tight_layout()
    save_figure(fig, out_prefix)
    plt.show()
    plt.close(fig)

# %% [markdown]
# ## 1. Selected-program effect heatmap
# 
# This panel focuses on day-relevant programs rather than all 50. It uses the sample-level matched effect table, so values are treatment minus matched control after sample-level aggregation.

# %%
effects = pd.read_csv(PROGRAM_EFFECTS_PATH, sep="\t")
paired = pd.read_csv(PAIRED_DIFFS_PATH, sep="\t")

print("effects columns:", effects.columns.tolist())
print("paired columns:", paired.columns.tolist())
print(effects.head())

# Filter focal programs/contrasts.
sel = effects[
    effects["program"].isin(SELECTED_PROGRAMS) &
    effects["contrast_name"].isin(FOCAL_CONTRASTS)
].copy()

# Pivot to selected programs x focal contrasts.
effect_col = "diff_treated_minus_control"
effect_mat = sel.pivot_table(index="program", columns="contrast_name", values=effect_col, aggfunc="first")
effect_mat = effect_mat.reindex(index=SELECTED_PROGRAMS, columns=FOCAL_CONTRASTS)
effect_mat.columns = [FOCAL_CONTRAST_LABELS[c] for c in effect_mat.columns]
effect_mat.index = [compact_label(p) for p in effect_mat.index]

# Save source table.
effect_mat.to_csv(OUT_DIR / "selected_D2_D3_program_effect_matrix.tsv", sep="\t")

# Plot.
plot_matrix_heatmap(
    effect_mat,
    title="Selected D2/D3 projected cNMF program effects\nsmall molecule treated − matched control",
    out_prefix=OUT_DIR / "selected_D2_D3_program_effect_heatmap",
    cmap="RdBu_r",
    center_zero=True,
    annotate=True,
    x_rotation=45,
    y_fontsize=7,
)

display(effect_mat)

# %% [markdown]
# ## 2. Paired replicate plots for selected programs
# 
# These panels show biological-replicate pairing. They are useful for checking whether effects are consistent or driven by one replicate.

# %%
def plot_paired_lines_for_programs(programs, contrasts, out_prefix, title):
    n_rows = len(programs)
    n_cols = len(contrasts)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(max(10, 1.9*n_cols), max(5, 1.7*n_rows)), squeeze=False)

    for i, program in enumerate(programs):
        for j, contrast in enumerate(contrasts):
            ax = axes[i, j]
            sub = paired[(paired["program"] == program) & (paired["contrast_name"] == contrast)].copy()
            if sub.empty:
                ax.axis("off")
                continue

            for _, row in sub.iterrows():
                ax.plot([0, 1], [row["control_value"], row["treated_value"]], marker="o", linewidth=0.9, alpha=0.8)
                ax.text(1.03, row["treated_value"], str(row["replicate_de"]).replace("rep", "r"), fontsize=5, va="center")

            mean_control = sub["control_value"].mean()
            mean_treated = sub["treated_value"].mean()
            ax.plot([0, 1], [mean_control, mean_treated], marker="s", linewidth=2.0, color="black")
            ax.set_xlim(-0.25, 1.35)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["ctrl", "treated"], rotation=45, ha="right")
            if i == 0:
                ax.set_title(FOCAL_CONTRAST_LABELS.get(contrast, contrast), fontsize=7)
            if j == 0:
                ax.set_ylabel(compact_label(program), fontsize=7)
            ax.grid(axis="y", linewidth=0.3, alpha=0.3)
            for spine in ["top", "right"]: ax.spines[spine].set_visible(False)

    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    save_figure(fig, out_prefix)
    plt.show()
    plt.close(fig)

plot_paired_lines_for_programs(
    programs=SELECTED_PROGRAMS,
    contrasts=FOCAL_CONTRASTS,
    out_prefix=OUT_DIR / "selected_D2_D3_program_paired_replicate_lines",
    title="Selected D2/D3 projected cNMF program usage by matched replicate",
)

# %% [markdown]
# ## 3. Expanded perturb-seq D2 × D3 program correlation matrix
# 
# This tests the developmental relationship between D2-active programs and D3 endpoint programs using perturb-seq KD effects. Each correlation is computed across perturbation targets in the `Perturbation Association 1` tab.
# 
# Rows: D2 programs. Columns: D3 EC and D3 alternative programs.

# %%
PERT_SHEET = "Perturbation Association 1"
PERT_EFFECT_COL = "approx_log2FC"  # use beta_hat if you want natural-log beta scale
PERT_TARGET_COL = "target_name"

pert = pd.read_excel(CNMF_XLSX_PATH, sheet_name=PERT_SHEET)
required = [PERT_TARGET_COL, "program_name", "Sample", PERT_EFFECT_COL]
missing = [c for c in required if c not in pert.columns]
if missing:
    raise ValueError(f"Missing required columns in {PERT_SHEET}: {missing}")

pert = pert.copy()
pert["program"] = "P" + pert["program_name"].astype(int).astype(str)
pert[PERT_EFFECT_COL] = pd.to_numeric(pert[PERT_EFFECT_COL], errors="coerce")

D2_CORR_PROGRAMS = D2_PROGRAMS
D3_CORR_PROGRAMS = D3_EC_PROGRAMS + D3_ALT_PROGRAMS

corr_records = []
for p2 in D2_CORR_PROGRAMS:
    d2 = pert[(pert["Sample"].astype(str) == "D2") & (pert["program"] == p2)][[PERT_TARGET_COL, PERT_EFFECT_COL]].dropna()
    d2 = d2.rename(columns={PERT_EFFECT_COL: "D2_effect"}).drop_duplicates(PERT_TARGET_COL)
    for p3 in D3_CORR_PROGRAMS:
        d3 = pert[(pert["Sample"].astype(str) == "D3") & (pert["program"] == p3)][[PERT_TARGET_COL, PERT_EFFECT_COL]].dropna()
        d3 = d3.rename(columns={PERT_EFFECT_COL: "D3_effect"}).drop_duplicates(PERT_TARGET_COL)
        m = d2.merge(d3, on=PERT_TARGET_COL, how="inner")
        if m.shape[0] >= 5:
            rho, pval = spearmanr(m["D2_effect"], m["D3_effect"])
            pr, pp = pearsonr(m["D2_effect"], m["D3_effect"])
        else:
            rho, pval, pr, pp = np.nan, np.nan, np.nan, np.nan
        corr_records.append({
            "D2_program": p2,
            "D3_program": p3,
            "D3_program_class": PROGRAM_SET_MAP.get(p3, ""),
            "n_targets": m.shape[0],
            "spearman_r": rho,
            "spearman_pvalue": pval,
            "pearson_r": pr,
            "pearson_pvalue": pp,
        })

corr_df = pd.DataFrame(corr_records)
corr_df["spearman_padj_bh"] = bh_adjust(corr_df["spearman_pvalue"])
corr_df.to_csv(OUT_DIR / "perturbseq_D2_vs_D3_program_correlation_table.tsv", sep="\t", index=False)

corr_mat = corr_df.pivot(index="D2_program", columns="D3_program", values="spearman_r").reindex(index=D2_CORR_PROGRAMS, columns=D3_CORR_PROGRAMS)
corr_mat_labeled = corr_mat.copy()
corr_mat_labeled.index = [compact_label(p) for p in corr_mat_labeled.index]
corr_mat_labeled.columns = [compact_label(p) for p in corr_mat_labeled.columns]

plot_matrix_heatmap(
    corr_mat_labeled,
    title="Perturb-seq program effect correlations\nD2 active programs vs D3 endpoint programs across perturbations",
    out_prefix=OUT_DIR / "perturbseq_D2_vs_D3_program_correlation_heatmap",
    cmap="RdBu_r",
    center_zero=True,
    annotate=True,
    x_rotation=45,
    y_fontsize=7,
)

display(corr_df.sort_values("spearman_r", ascending=False).head(15))
display(corr_df.sort_values("spearman_r").head(15))

# %% [markdown]
# ## 4. Marker gene candidates from top cNMF loadings and concordant DE
# 
# This selects genes that are both top-loading in a program and DE in the expected direction in a small-molecule contrast. It creates a reproducible gene list for marker-expression panels.

# %%
LOADINGS_SHEET = "Program Loadings"
TOP_N_LOADING = 300
MAX_MARKERS_PER_RULE = 8
DE_PADJ_CUTOFF = 0.05
DE_ABS_LFC_MIN = 0.25

loadings = pd.read_excel(CNMF_XLSX_PATH, sheet_name=LOADINGS_SHEET)
required = ["Program", "Rank", "Gene"]
missing = [c for c in required if c not in loadings.columns]
if missing:
    raise ValueError(f"Missing columns in Program Loadings: {missing}")

loadings = loadings.copy()
loadings["program"] = "P" + loadings["Program"].astype(int).astype(str)
loadings["Rank"] = pd.to_numeric(loadings["Rank"], errors="coerce")
loadings["Gene"] = loadings["Gene"].astype(str)
loadings_top = loadings[loadings["Rank"] <= TOP_N_LOADING].copy()

# Locate PyDESeq2 result table for each contrast.
def find_deseq_file(contrast_name):
    filename = f"PyDESeq2_{contrast_name}.tsv"
    return find_file(filename, required=True)

candidate_records = []
for rule in MARKER_RULES:
    program = rule["program"]
    contrast_name = rule["contrast_name"]
    direction = rule["direction"]
    label = rule["label"]

    de_path = find_deseq_file(contrast_name)
    de = pd.read_csv(de_path, sep="\t")
    if "gene_symbol" not in de.columns or "log2FoldChange" not in de.columns or "padj" not in de.columns:
        raise ValueError(f"DE table lacks required columns: {de_path}")

    de = de.copy()
    de["gene_symbol"] = de["gene_symbol"].astype(str)
    de["log2FoldChange"] = pd.to_numeric(de["log2FoldChange"], errors="coerce")
    de["padj"] = pd.to_numeric(de["padj"], errors="coerce")

    top = loadings_top[loadings_top["program"] == program][["program", "Rank", "Gene"]].copy()
    merged = top.merge(de, left_on="Gene", right_on="gene_symbol", how="inner")

    if direction == "up":
        merged = merged[(merged["padj"] < DE_PADJ_CUTOFF) & (merged["log2FoldChange"] >= DE_ABS_LFC_MIN)]
        merged = merged.sort_values(["Rank", "padj", "log2FoldChange"], ascending=[True, True, False])
    elif direction == "down":
        merged = merged[(merged["padj"] < DE_PADJ_CUTOFF) & (merged["log2FoldChange"] <= -DE_ABS_LFC_MIN)]
        merged = merged.sort_values(["Rank", "padj", "log2FoldChange"], ascending=[True, True, True])
    else:
        raise ValueError("direction must be up or down")

    merged = merged.head(MAX_MARKERS_PER_RULE).copy()
    merged["marker_rule"] = label
    merged["program"] = program
    merged["contrast_name"] = contrast_name
    merged["expected_direction"] = direction
    candidate_records.append(merged)

marker_candidates = pd.concat(candidate_records, ignore_index=True) if candidate_records else pd.DataFrame()
keep_cols = [
    "marker_rule", "program", "Gene", "Rank", "contrast_name", "expected_direction",
    "gene_id", "gene_symbol", "log2FoldChange", "padj", "baseMean", "stat"
]
keep_cols = [c for c in keep_cols if c in marker_candidates.columns]
marker_candidates = marker_candidates[keep_cols]
marker_candidates.to_csv(OUT_DIR / "selected_program_marker_gene_candidates.tsv", sep="\t", index=False)

display(marker_candidates)
print("Saved:", OUT_DIR / "selected_program_marker_gene_candidates.tsv")

# %% [markdown]
# ## 5. Optional marker gene pseudobulk expression heatmap
# 
# This requires pseudobulk counts and sample metadata files. If those files are not present in the current environment, this cell will skip plotting after saving the marker candidate table.

# %%
# Robust marker gene pseudobulk expression heatmap.
# This version handles the common case where:
#   - marker_candidates["Gene"] contains gene symbols
#   - pseudobulk counts rows are Ensembl gene IDs
#   - gene_annotation.tsv maps gene_id -> gene_symbol
GENE_ANNOTATION_PATH = PSEUDOBULK_COUNTS_PATH.parent / "gene_annotation.tsv"

if not PSEUDOBULK_COUNTS_PATH.exists() or not PSEUDOBULK_META_PATH.exists():
    print("Pseudobulk counts/meta not found. Skipping marker expression heatmap.")
    print("Expected counts:", PSEUDOBULK_COUNTS_PATH)
    print("Expected metadata:", PSEUDOBULK_META_PATH)
else:
    counts = pd.read_csv(PSEUDOBULK_COUNTS_PATH, sep="\t", index_col=0)
    sample_meta = pd.read_csv(PSEUDOBULK_META_PATH, sep="\t")

    print("counts shape, genes x samples:", counts.shape)
    print("sample_meta shape:", sample_meta.shape)
    print("counts first rows:", counts.index[:5].tolist())
    print("sample_meta columns:", sample_meta.columns.tolist())

    # Find sample column robustly.
    sample_col = None
    for candidate in ["sample_name_de", "sample_name", "sample"]:
        if candidate in sample_meta.columns:
            sample_col = candidate
            break
    if sample_col is None:
        sample_col = sample_meta.columns[0]
        print(f"WARNING: could not find sample column; using first column {sample_col!r}")

    sample_meta = sample_meta.set_index(sample_col, drop=False)
    common_samples = [s for s in counts.columns if s in sample_meta.index]
    if len(common_samples) == 0:
        raise ValueError(
            "No overlapping samples between pseudobulk count columns and metadata index.\n"
            f"counts columns first 10: {counts.columns[:10].tolist()}\n"
            f"metadata index first 10: {sample_meta.index[:10].tolist()}"
        )

    counts = counts[common_samples]
    sample_meta = sample_meta.loc[common_samples]

    # If no DE-concordant marker candidates were found, fall back to top loading genes.
    marker_table = marker_candidates.copy()
    if marker_table.empty:
        print("WARNING: marker_candidates is empty. Falling back to top loading genes from MARKER_RULES.")
        fallback_records = []
        for rule in MARKER_RULES:
            program = rule["program"]
            label = rule["label"]
            top = loadings_top.loc[loadings_top["program"] == program].sort_values("Rank").head(MAX_MARKERS_PER_RULE)
            for _, row in top.iterrows():
                fallback_records.append({
                    "marker_rule": label,
                    "program": program,
                    "Gene": row["Gene"],
                    "Rank": row["Rank"],
                    "contrast_name": rule["contrast_name"],
                    "expected_direction": rule["direction"],
                    "fallback_reason": "top_loading_no_DE_filter",
                })
        marker_table = pd.DataFrame(fallback_records)

    if marker_table.empty:
        raise ValueError("No marker genes available after fallback. Check Program Loadings and MARKER_RULES.")

    marker_genes = marker_table["Gene"].astype(str).drop_duplicates().tolist()
    print(f"Candidate marker genes: {len(marker_genes)}")
    print("first candidate marker genes:", marker_genes[:20])

    # Build counts matrix indexed by marker gene SYMBOL.
    # Case 1: counts rows are already symbols.
    direct_present = [g for g in marker_genes if g in counts.index]
    print(f"Direct symbol matches in counts index: {len(direct_present)} / {len(marker_genes)}")

    marker_counts_by_symbol = {}
    gene_symbol_to_count_ids = {}

    if len(direct_present) > 0:
        for gene in direct_present:
            marker_counts_by_symbol[gene] = counts.loc[[gene]].sum(axis=0)
            gene_symbol_to_count_ids[gene] = [gene]

    # Case 2: counts rows are gene IDs; map marker symbols through gene_annotation.tsv.
    if GENE_ANNOTATION_PATH.exists():
        gene_annot = pd.read_csv(GENE_ANNOTATION_PATH, sep="\t")
        if {"gene_id", "gene_symbol"}.issubset(gene_annot.columns):
            gene_annot["gene_id"] = gene_annot["gene_id"].astype(str)
            gene_annot["gene_symbol"] = gene_annot["gene_symbol"].astype(str)
            gene_annot = gene_annot[gene_annot["gene_id"].isin(counts.index)]

            for gene in marker_genes:
                if gene in marker_counts_by_symbol:
                    continue
                ids = gene_annot.loc[gene_annot["gene_symbol"] == gene, "gene_id"].drop_duplicates().tolist()
                ids = [x for x in ids if x in counts.index]
                if ids:
                    # Sum duplicated gene IDs for the same symbol.
                    marker_counts_by_symbol[gene] = counts.loc[ids].sum(axis=0)
                    gene_symbol_to_count_ids[gene] = ids
        else:
            print("WARNING: gene_annotation.tsv lacks gene_id/gene_symbol columns:", gene_annot.columns.tolist())
    else:
        print("WARNING: gene_annotation.tsv not found:", GENE_ANNOTATION_PATH)

    marker_genes_present = [g for g in marker_genes if g in marker_counts_by_symbol]
    missing_marker_genes = [g for g in marker_genes if g not in marker_counts_by_symbol]
    print(f"Marker genes mappable to pseudobulk counts: {len(marker_genes_present)} / {len(marker_genes)}")
    if missing_marker_genes:
        print("First missing marker genes:", missing_marker_genes[:30])

    if len(marker_genes_present) == 0:
        raise ValueError(
            "No marker genes could be mapped to pseudobulk counts. This usually means counts rows are gene IDs "
            "and gene_annotation.tsv is missing or has mismatched symbols."
        )

    marker_counts = pd.DataFrame({g: marker_counts_by_symbol[g] for g in marker_genes_present}).T
    marker_counts = marker_counts[common_samples]

    # CPM normalization and log transform.
    lib_size = counts.sum(axis=0)
    cpm = marker_counts.div(lib_size, axis=1) * 1e6
    expr = np.log2(cpm + 1)

    # Sample order: D2/D3 focused, using the sample names as they appear.
    preferred_conditions = [
        "D2_LM", "D2_LM+TDI", "D2_LM+IAG",
        "D3_AE", "D3_AE+TDI", "D3_AE+IAG",
        "D2_LM+TDI->D3_AE+TDI", "D2_LM+IAG->D3_AE+IAG",
    ]
    ordered_samples = []
    for cond in preferred_conditions:
        ordered_samples += [s for s in expr.columns if str(s).startswith(cond + "_") or str(s) == cond]
    ordered_samples += [s for s in expr.columns if s not in ordered_samples]
    expr = expr[ordered_samples]

    # Order genes by marker rule/program, but keep only genes that mapped to expression.
    gene_order = []
    gene_labels = []
    seen = set()
    for _, row in marker_table.iterrows():
        gene = str(row["Gene"])
        if gene in expr.index and gene not in seen:
            gene_order.append(gene)
            seen.add(gene)
            rule_label = row.get("marker_rule", row.get("program", "marker"))
            gene_labels.append(f"{row['program']} {gene}")

    if len(gene_order) == 0:
        raise ValueError(
            "Marker genes mapped to counts, but none survived gene_order construction. "
            "Check marker_table['Gene'] and expr.index."
        )

    expr = expr.loc[gene_order]
    expr.index = gene_labels

    if expr.shape[0] == 0 or expr.shape[1] == 0:
        raise ValueError(
            f"Marker expression matrix is empty: shape={expr.shape}. "
            "Cannot plot a heatmap."
        )

    expr_z = zscore_rows(expr)
    expr.to_csv(OUT_DIR / "selected_program_marker_gene_log2CPM_by_sample.tsv", sep="\t")
    expr_z.to_csv(OUT_DIR / "selected_program_marker_gene_log2CPM_zscore_by_gene.tsv", sep="\t")

    # Save mapping table for transparency.
    mapping_rows = []
    for gene in marker_genes_present:
        mapping_rows.append({
            "Gene": gene,
            "count_row_ids": ";".join(gene_symbol_to_count_ids.get(gene, [])),
            "n_count_rows_summed": len(gene_symbol_to_count_ids.get(gene, [])),
        })
    pd.DataFrame(mapping_rows).to_csv(OUT_DIR / "selected_marker_gene_symbol_to_count_rows.tsv", sep="\t", index=False)

    print("Final expression heatmap matrix shape:", expr_z.shape)
    display(expr_z.iloc[:5, :5])

    plot_matrix_heatmap(
        expr_z,
        title="Selected marker genes\npseudobulk log2(CPM+1), z-scored by gene",
        out_prefix=OUT_DIR / "selected_program_marker_gene_expression_heatmap",
        cmap="RdBu_r",
        center_zero=True,
        annotate=False,
        x_rotation=60,
        y_fontsize=7,
    )

display(marker_candidates)

# %% [markdown]
# ## Interpretation guide
# 
# Use these outputs as follows:
# 
# - `selected_D2_D3_program_effect_heatmap.svg`: compact program effect panel for D2/D3 programs.
# - `selected_D2_D3_program_paired_replicate_lines.svg`: replicate-aware program score evidence.
# - `perturbseq_D2_vs_D3_program_correlation_heatmap.svg`: expanded D2-to-D3 program relationship map from perturb-seq effects.
# - `selected_program_marker_gene_candidates.tsv`: reproducible marker-gene set based on cNMF loadings and concordant DE.
# - `selected_program_marker_gene_expression_heatmap.svg`: marker expression heatmap if pseudobulk counts are available.
# 
# For the main figure, I would use the selected-program effect heatmap + a marker gene heatmap. For the supplement, use the paired replicate lines and full D2×D3 correlation matrix.