# %% [markdown]
# # Marker gene heatmaps from Kyle/Hunter program annotation notes
#
# This notebook uses the representative genes from the first table of
# `Program annotation notes (Kyle+Hunter) (1).docx`, not the top-loading gene list.
# It makes pseudobulk marker-gene expression heatmaps across the Hippo time-course samples.

# %%
from pathlib import Path
import re
import textwrap
import warnings

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

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

# %% [markdown]
# ## Paths and settings

# %%
BASE = Path("/scratch/users/edtng/20260325_Hippo_timecourse")
DE_DIR = BASE / "perturbseq_cNMF_matched_DE"
OUT_DIR = BASE / "hippo_marker_gene_panels_from_annotation_notes"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PSEUDOBULK_COUNTS_PATH = DE_DIR / "pseudobulk_gene_counts.filtered.tsv.gz"
PSEUDOBULK_META_PATH = DE_DIR / "pseudobulk_sample_metadata.filtered.tsv"
GENE_ANNOTATION_PATH = DE_DIR / "gene_annotation.tsv"

# If you copy the annotation notes docx onto Sherlock, set this path.
DOCX_PATH = Path("/scratch/users/edtng/20260325_Hippo_timecourse/Program annotation notes (Kyle+Hunter) (1).docx")

SAMPLE_COL = "sample_name"
CONDITION_COL = "condition"

# Set whether to plot biological replicate samples or condition means.
MAKE_SAMPLE_REPLICATE_HEATMAP = True
MAKE_CONDITION_MEAN_HEATMAP = True

# Representative programs to show in the Hippo panel.
# Edit if you want a broader/all-program marker panel.
PROGRAMS_TO_PLOT = [
    "P15", "P17", "P23",   # D2 vascular / lateral mesoderm programs
    "P2", "P33", "P34", "P24",  # D3 EC / alternative programs with representative genes
]

# Exact condition order for plotting.
CONDITION_TIMECOURSE_ORDER = [
    "D0_hPSC",
    "D1_MPS", "D1_MPS+TDI", "D1_MPS+IAG",
    "D1.33_LM", "D1.33_LM+TDI", "D1.33_LM+IAG",
    "D1.67_LM", "D1.67_LM+TDI", "D1.67_LM+IAG",
    "D2_LM", "D2_LM+TDI", "D2_LM+IAG",
    "D2_LM+TDI->D3_AE+TDI", "D2_LM+IAG->D3_AE+IAG",
    "D3_AE", "D3_AE+TDI", "D3_AE+IAG",
]

# Optional focused conditions for a smaller figure panel.
FOCUSED_CONDITIONS = [
    "D2_LM", "D2_LM+TDI", "D2_LM+IAG",
    "D3_AE", "D3_AE+TDI", "D3_AE+IAG",
    "D2_LM+TDI->D3_AE+TDI", "D2_LM+IAG->D3_AE+IAG",
]

# %% [markdown]
# ## Representative genes from first table of annotation notes
#
# These were extracted from the first table in the annotation note document. Programs marked
# "to eliminate" or rows without concrete gene symbols are omitted by default from the heatmap.

# %%
PROGRAM_REPRESENTATIVE_GENES = {
    "P7":  ["L1TD1", "DPPA4", "TERF1", "SOX2"],
    "P22": ["POU5F1", "TDGF1", "DNMT3B"],
    "P29": ["ESRG", "DPPA4", "NANOG"],
    "P28": ["MIXL1", "TBXT", "SP5"],
    "P13": ["NKD1", "CDX2", "MSX1"],
    "P30": ["DKK1", "NKD1", "TBXT"],
    "P48": ["NKD1", "TBXT", "WLS"],
    "P15": ["TAL1", "FLI1", "SOX7"],
    "P23": ["KDR", "ETV2"],
    "P17": ["KRT8", "KRT18"],
    "P33": ["CER1", "LHX1", "OTX2", "CYP26A1"],
    "P34": ["ANKRD1", "KRT8", "KRT18"],
    "P24": ["GATA4", "RSPO2"],
    "P2":  ["MECOM", "NRP1", "EFNB2"],
}

PROGRAM_ANNOTATIONS = {
    "P7":  "Pluripotency",
    "P22": "Pluripotency",
    "P29": "Pluripotency",
    "P28": "Primitive streak",
    "P13": "Posterior primitive streak",
    "P30": "Primitive streak, BMP-repressed 1",
    "P48": "Primitive streak, BMP-repressed 2",
    "P15": "Vascular progenitor",
    "P23": "Vascular lateral mesoderm",
    "P17": "Epithelial (?) lateral mesoderm",
    "P33": "Anterior-like mesoderm",
    "P34": "Mesenchyme, YAP/TAZ-induced",
    "P24": "GATA4/RSPO2 program",
    "P2":  "Artery endothelium",
}

records = []
for program, genes in PROGRAM_REPRESENTATIVE_GENES.items():
    for rank, gene in enumerate(genes, start=1):
        records.append({
            "program": program,
            "program_annotation": PROGRAM_ANNOTATIONS.get(program, ""),
            "gene": gene,
            "representative_gene_rank": rank,
        })
representative_gene_df = pd.DataFrame(records)
representative_gene_df.to_csv(OUT_DIR / "representative_genes_from_annotation_notes.tsv", sep="\t", index=False)
display(representative_gene_df)

# %% [markdown]
# ## Load pseudobulk counts and metadata

# %%
counts = pd.read_csv(PSEUDOBULK_COUNTS_PATH, sep="\t", index_col=0)
meta = pd.read_csv(PSEUDOBULK_META_PATH, sep="\t")

# The expected count table is genes x samples.
print("counts shape genes x samples:", counts.shape)
print("counts first row IDs:", counts.index[:5].tolist())
print("meta shape:", meta.shape)
print("meta columns:", list(meta.columns))
display(meta.head())

# Infer sample_name column if needed.
if SAMPLE_COL not in meta.columns:
    for candidate in ["sample", "sample_id", "pseudobulk_sample", "sample_name"]:
        if candidate in meta.columns:
            meta[SAMPLE_COL] = meta[candidate].astype(str)
            break
if SAMPLE_COL not in meta.columns:
    # fallback: use index if metadata rows match columns
    if meta.shape[0] == counts.shape[1]:
        meta[SAMPLE_COL] = counts.columns.astype(str)
    else:
        raise ValueError(f"Could not find or infer {SAMPLE_COL} in metadata.")

# Align metadata to count columns.
meta[SAMPLE_COL] = meta[SAMPLE_COL].astype(str)
meta = meta.drop_duplicates(SAMPLE_COL).set_index(SAMPLE_COL).reindex(counts.columns.astype(str))
if meta.isna().all(axis=1).any():
    missing = meta.index[meta.isna().all(axis=1)].tolist()
    raise ValueError(f"Metadata missing for count columns: {missing[:10]}")

# Infer condition if needed.
def infer_condition_from_sample_name(sample_name: str):
    sample_name = str(sample_name)
    m = re.match(r"(.+)_rep\d+$", sample_name)
    if m:
        return m.group(1)
    return sample_name

if CONDITION_COL not in meta.columns:
    meta[CONDITION_COL] = [infer_condition_from_sample_name(x) for x in meta.index]
meta[CONDITION_COL] = meta[CONDITION_COL].astype(str)

display(meta[[CONDITION_COL]].head())

# %% [markdown]
# ## Map representative gene symbols to count rows

# %%
def infer_gene_annotation_columns(annotation_df):
    cols = list(annotation_df.columns)
    gene_id_candidates = ["gene_id", "gene_ids", "ensembl", "ensembl_id", "var_names", "id"]
    symbol_candidates = ["gene", "gene_symbol", "symbol", "gene_name", "hgnc_symbol", "external_gene_name"]

    gene_id_col = next((c for c in gene_id_candidates if c in cols), None)
    symbol_col = next((c for c in symbol_candidates if c in cols), None)

    # Common case: first column was saved as unnamed gene id index
    if gene_id_col is None:
        unnamed = [c for c in cols if str(c).startswith("Unnamed")]
        if unnamed:
            gene_id_col = unnamed[0]

    if gene_id_col is None or symbol_col is None:
        raise ValueError(
            f"Could not infer gene annotation columns. Columns: {cols}. "
            "Expected something like gene_id and gene/gene_symbol."
        )
    return gene_id_col, symbol_col

marker_symbols = representative_gene_df.loc[
    representative_gene_df["program"].isin(PROGRAMS_TO_PLOT),
    "gene",
].drop_duplicates().astype(str).tolist()

print("Marker symbols requested:", marker_symbols)

# Direct symbol matches in count index.
direct_matches = {g: [g] for g in marker_symbols if g in counts.index}

symbol_to_count_rows = dict(direct_matches)

if len(symbol_to_count_rows) < len(marker_symbols) and GENE_ANNOTATION_PATH.exists():
    gene_anno = pd.read_csv(GENE_ANNOTATION_PATH, sep="\t")
    gene_id_col, symbol_col = infer_gene_annotation_columns(gene_anno)
    gene_anno[gene_id_col] = gene_anno[gene_id_col].astype(str)
    gene_anno[symbol_col] = gene_anno[symbol_col].astype(str)

    available_rows = set(counts.index.astype(str))
    for symbol in marker_symbols:
        if symbol in symbol_to_count_rows:
            continue
        rows = gene_anno.loc[gene_anno[symbol_col] == symbol, gene_id_col].astype(str).tolist()
        rows = [r for r in rows if r in available_rows]
        if rows:
            symbol_to_count_rows[symbol] = rows

mapping_records = []
for symbol in marker_symbols:
    rows = symbol_to_count_rows.get(symbol, [])
    mapping_records.append({
        "gene_symbol": symbol,
        "n_count_rows": len(rows),
        "count_rows": ",".join(rows),
        "found": len(rows) > 0,
    })
marker_map_df = pd.DataFrame(mapping_records)
marker_map_df.to_csv(OUT_DIR / "representative_gene_symbol_to_count_rows.tsv", sep="\t", index=False)
display(marker_map_df)

missing_symbols = marker_map_df.loc[~marker_map_df["found"], "gene_symbol"].tolist()
if missing_symbols:
    print("WARNING: missing marker symbols:", missing_symbols)

if marker_map_df["found"].sum() == 0:
    raise ValueError("No representative marker genes found in pseudobulk counts after mapping.")

# %% [markdown]
# ## Build marker expression matrix: log2(CPM+1)

# %%
# Sum counts by marker gene symbol if multiple count rows map to one symbol.
marker_count_rows = []
for symbol in marker_symbols:
    rows = symbol_to_count_rows.get(symbol, [])
    if not rows:
        continue
    summed = counts.loc[rows].sum(axis=0)
    summed.name = symbol
    marker_count_rows.append(summed)

marker_counts = pd.DataFrame(marker_count_rows)
print("marker_counts shape genes x samples:", marker_counts.shape)
display(marker_counts.iloc[:5, :5])

# CPM using total pseudobulk library size across all genes.
lib_size = counts.sum(axis=0)
marker_cpm = marker_counts.div(lib_size, axis=1) * 1e6
marker_log2cpm = np.log2(marker_cpm + 1)

marker_log2cpm.to_csv(OUT_DIR / "representative_gene_log2CPM_by_sample.tsv", sep="\t")

# %% [markdown]
# ## Create gene display labels grouped by program

# %%
# Build row order: programs in PROGRAMS_TO_PLOT, genes in table order.
row_gene_order = []
row_display_labels = []
row_programs = []
for program in PROGRAMS_TO_PLOT:
    genes = representative_gene_df.loc[representative_gene_df["program"] == program, "gene"].astype(str).tolist()
    for gene in genes:
        if gene in marker_log2cpm.index:
            row_gene_order.append(gene)
            row_programs.append(program)
            row_display_labels.append(f"{program} {gene}")

expr = marker_log2cpm.loc[row_gene_order].copy()
expr.index = row_display_labels

row_anno = pd.DataFrame({
    "display_label": row_display_labels,
    "program": row_programs,
    "gene": row_gene_order,
    "program_annotation": [PROGRAM_ANNOTATIONS.get(p, "") for p in row_programs],
})
row_anno.to_csv(OUT_DIR / "representative_gene_heatmap_rows.tsv", sep="\t", index=False)
display(row_anno)

# %% [markdown]
# ## Order samples by time course

# %%
def sample_order_from_timecourse(meta, focus_conditions=None):
    # metadata index is sample_name
    rows = meta.reset_index().rename(columns={meta.index.name or "index": SAMPLE_COL})
    if SAMPLE_COL not in rows.columns:
        rows[SAMPLE_COL] = meta.index.astype(str)

    if focus_conditions is not None:
        rows = rows[rows[CONDITION_COL].isin(focus_conditions)].copy()

    ordered_samples = []
    for cond in CONDITION_TIMECOURSE_ORDER:
        sub = rows.loc[rows[CONDITION_COL] == cond, SAMPLE_COL].astype(str).tolist()
        # sort rep1, rep2, rep3 naturally when possible
        sub = sorted(sub, key=lambda x: int(re.search(r"_rep(\d+)$", x).group(1)) if re.search(r"_rep(\d+)$", x) else 999)
        ordered_samples.extend(sub)

    remaining = [s for s in rows[SAMPLE_COL].astype(str).tolist() if s not in ordered_samples]
    ordered_samples += remaining
    return [s for s in ordered_samples if s in expr.columns]

all_sample_order = sample_order_from_timecourse(meta)
focused_sample_order = sample_order_from_timecourse(meta, FOCUSED_CONDITIONS)

pd.DataFrame({"sample_name": all_sample_order}).to_csv(OUT_DIR / "all_sample_order_used.tsv", sep="\t", index=False)
pd.DataFrame({"sample_name": focused_sample_order}).to_csv(OUT_DIR / "focused_sample_order_used.tsv", sep="\t", index=False)

print("all samples:", len(all_sample_order))
print(all_sample_order)
print("focused samples:", len(focused_sample_order))
print(focused_sample_order)

# %% [markdown]
# ## Condition means and z-scores

# %%
# Sample-level replicate heatmap matrix.
expr_sample = expr.loc[:, all_sample_order]
expr_sample_focused = expr.loc[:, focused_sample_order]

# Condition mean matrix.
condition_to_samples = {
    cond: [s for s in all_sample_order if meta.loc[s, CONDITION_COL] == cond]
    for cond in CONDITION_TIMECOURSE_ORDER
}
condition_to_samples = {k: v for k, v in condition_to_samples.items() if v}
expr_condition_mean = pd.DataFrame({
    cond: expr.loc[:, samples].mean(axis=1)
    for cond, samples in condition_to_samples.items()
})
expr_condition_mean_focused = expr_condition_mean.loc[:, [c for c in FOCUSED_CONDITIONS if c in expr_condition_mean.columns]]

def row_zscore(df):
    z = df.sub(df.mean(axis=1), axis=0)
    z = z.div(df.std(axis=1).replace(0, np.nan), axis=0)
    return z.fillna(0)

expr_sample_z = row_zscore(expr_sample)
expr_sample_focused_z = row_zscore(expr_sample_focused)
expr_condition_mean_z = row_zscore(expr_condition_mean)
expr_condition_mean_focused_z = row_zscore(expr_condition_mean_focused)

# Save matrices.
expr_sample.to_csv(OUT_DIR / "representative_gene_log2CPM_by_sample.replicates.tsv", sep="\t")
expr_sample_z.to_csv(OUT_DIR / "representative_gene_log2CPM_by_sample.replicates.zscore_by_gene.tsv", sep="\t")
expr_condition_mean.to_csv(OUT_DIR / "representative_gene_log2CPM_by_condition.mean.tsv", sep="\t")
expr_condition_mean_z.to_csv(OUT_DIR / "representative_gene_log2CPM_by_condition.mean.zscore_by_gene.tsv", sep="\t")

# %% [markdown]
# ## Plot Illustrator-editable red-blue heatmaps

# %%
def plot_heatmap(df, title, out_prefix, center_zero=True, cmap="RdBu_r", x_rotation=55):
    if df.shape[0] == 0 or df.shape[1] == 0:
        raise ValueError(f"Cannot plot empty matrix with shape {df.shape}")

    mat = df.to_numpy(dtype=float)
    n_rows, n_cols = mat.shape

    if center_zero:
        vmax = np.nanmax(np.abs(mat))
        vmin = -vmax
    else:
        vmin = np.nanmin(mat)
        vmax = np.nanmax(mat)

    fig_w = max(8, 0.42 * n_cols + 4)
    fig_h = max(4, 0.35 * n_rows + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    mesh = ax.pcolormesh(
        np.arange(n_cols + 1),
        np.arange(n_rows + 1),
        mat,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading="flat",
        linewidth=0,
        antialiased=False,
    )

    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_rows, 0)
    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_xticklabels(df.columns.tolist(), rotation=x_rotation, ha="right", fontsize=7)
    ax.set_yticklabels(df.index.tolist(), fontsize=7)
    ax.set_title(title)
    ax.tick_params(length=0)

    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label("row z-score" if center_zero else "log2(CPM + 1)")

    fig.tight_layout()
    for ext in ["svg", "pdf", "png"]:
        fig.savefig(OUT_DIR / f"{out_prefix}.{ext}", bbox_inches="tight")
    plt.show()
    plt.close(fig)

# Focused condition mean is probably best for a compact main figure panel.
plot_heatmap(
    expr_condition_mean_focused_z,
    "Representative genes from program annotation notes\ncondition mean log2(CPM+1), z-scored by gene",
    "representative_marker_genes_focused_condition_mean_zscore_redblue_editable",
    center_zero=True,
)

# Replicate-level focused heatmap for seeing biological replicate variation.
plot_heatmap(
    expr_sample_focused_z,
    "Representative genes from program annotation notes\nreplicate log2(CPM+1), z-scored by gene",
    "representative_marker_genes_focused_replicates_zscore_redblue_editable",
    center_zero=True,
)

# Full time-course condition mean.
plot_heatmap(
    expr_condition_mean_z,
    "Representative genes from program annotation notes\nfull time course, condition mean log2(CPM+1), z-scored by gene",
    "representative_marker_genes_full_timecourse_condition_mean_zscore_redblue_editable",
    center_zero=True,
)

print("Saved outputs to:", OUT_DIR)
for p in sorted(OUT_DIR.glob("representative_marker_genes*.svg")):
    print(p)
