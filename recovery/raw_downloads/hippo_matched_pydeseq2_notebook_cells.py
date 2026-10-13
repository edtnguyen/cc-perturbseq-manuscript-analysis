# %% [markdown]
# # Hippo time-course matched pseudobulk DE + projected cNMF program contrasts
#
# Copy these cells into a Jupyter notebook and run in order.
#
# This workflow:
# 1. Loads the zero-based projected AnnData.
# 2. Applies the condition-specific mito/ribo filters.
# 3. Creates sample-level pseudobulk raw counts.
# 4. Runs PyDESeq2 for matched treated-vs-control gene DE.
# 5. Runs matched sample-level tests on projected cNMF program scores.
# 6. Optionally scores cNMF program enrichment in DE gene log2FC.

# %%
# CELL 1 — Imports, paths, sample metadata, QC cutoffs, and contrasts

from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
import matplotlib as mpl
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

mpl.rcParams.update({
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})

# Input/output paths
H5AD_PATH = Path(
    "/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_projection/"
    "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5ad"
)

BASE_OUT = Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_matched_DE")
BASE_OUT.mkdir(parents=True, exist_ok=True)

SPECTRA_TPM_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50."
    "gene_spectra_tpm.k_50.dt_0_2.txt"
)

# Use raw counts for DE. Change to layer:mature or layer:nascent if you want those layers.
COUNTS_SOURCE = "X"  # "X" or "layer:<layer_name>"
PROGRAM_SCORE_KEY = "X_perturbseq_cNMF_usage"  # normalized projected usage
SAMPLE_COL = "sample_name"
CHUNK_SIZE = 4000
GENE_SYMBOL_COL = "gene"

# PyDESeq2 settings
N_CPUS = 8
MIN_SAMPLES_PER_GROUP = 2
MIN_COUNT = 10
MIN_SAMPLES_EXPRESSED = 2
SIZE_FACTORS_FIT_TYPE = "ratio"  # ratio / poscounts / iterative
USE_PAIRED_REPLICATE_DESIGN = True
QUIET_PYDESEQ2 = False

# Programs to highlight in the program-score contrast output. Zero-based labels.
KEY_PROGRAMS = ["P15", "P17", "P34"]

sample_mapping_dict = {
    '1-1': 'D0_hPSC_rep1', '1-2': 'D0_hPSC_rep2', '1-3': 'D0_hPSC_rep3',
    '1-4': 'D1_MPS_rep1', '1-5': 'D1_MPS_rep2', '1-6': 'D1_MPS_rep3',
    '1-7': 'D1_MPS+TDI_rep1', '1-8': 'D1_MPS+TDI_rep2', '1-9': 'D1_MPS+TDI_rep3',
    '1-10': 'D1_MPS+IAG_rep1', '1-11': 'D1_MPS+IAG_rep2', '1-12': 'D1_MPS+IAG_rep3',
    '2-1': 'D1.33_LM_rep1', '2-2': 'D1.33_LM_rep2', '2-3': 'D1.33_LM_rep3',
    '2-4': 'D1.33_LM+TDI_rep1', '2-5': 'D1.33_LM+TDI_rep2', '2-6': 'D1.33_LM+TDI_rep3',
    '2-7': 'D1.33_LM+IAG_rep1', '2-8': 'D1.33_LM+IAG_rep2', '2-9': 'D1.33_LM+IAG_rep3',
    '3-1': 'D1.67_LM_rep1', '3-2': 'D1.67_LM_rep2', '3-3': 'D1.67_LM_rep3',
    '3-4': 'D1.67_LM+TDI_rep1', '3-5': 'D1.67_LM+TDI_rep2', '3-6': 'D1.67_LM+TDI_rep3',
    '3-7': 'D1.67_LM+IAG_rep1', '3-8': 'D1.67_LM+IAG_rep2', '3-9': 'D1.67_LM+IAG_rep3',
    '4-1': 'D2_LM_rep1', '4-2': 'D2_LM_rep2', '4-3': 'D2_LM_rep3',
    '4-4': 'D2_LM+TDI_rep1', '4-5': 'D2_LM+TDI_rep2', '4-6': 'D2_LM+TDI_rep3',
    '4-7': 'D2_LM+IAG_rep1', '4-8': 'D2_LM+IAG_rep2', '4-9': 'D2_LM+IAG_rep3',
    '5-1': 'D3_AE_rep1', '5-2': 'D3_AE_rep2', '5-3': 'D3_AE_rep3',
    '5-4': 'D3_AE+TDI_rep1', '5-5': 'D3_AE+TDI_rep2', '5-6': 'D3_AE+TDI_rep3',
    '5-7': 'D3_AE+IAG_rep1', '5-8': 'D3_AE+IAG_rep2', '5-9': 'D3_AE+IAG_rep3',
    '5-10': 'D2_LM+TDI->D3_AE+TDI_rep1', '5-11': 'D2_LM+TDI->D3_AE+TDI_rep2',
    '5-12': 'D2_LM+TDI->D3_AE+TDI_rep3', '5-13': 'D2_LM+IAG->D3_AE+IAG_rep1',
    '5-14': 'D2_LM+IAG->D3_AE+IAG_rep2', '5-15': 'D2_LM+IAG->D3_AE+IAG_rep3',
}

cond_pct_ribo_cutoffs = {
    'D0_hPSC': 2.4,
    'D1_MPS': 2.3,
    'D1_MPS+TDI': 2.5,
    'D1_MPS+IAG': 2.2,
    'D1.33_LM': 1.9,
    'D1.33_LM+TDI': 2.1,
    'D1.33_LM+IAG': 2.2,
    'D1.67_LM': 2.0,
    'D1.67_LM+TDI': 2.2,
    'D1.67_LM+IAG': 1.9,
    'D2_LM': 2.0,
    'D2_LM+TDI': 2.4,
    'D2_LM+IAG': 1.9,
    'D3_AE': 1.0,
    'D3_AE+TDI': 2.0,
    'D3_AE+IAG': 1.0,
    'D2_LM+TDI->D3_AE+TDI': 2.5,
    'D2_LM+IAG->D3_AE+IAG': 1.0,
}

cond_pct_mt_cutoffs = {
    'D0_hPSC': 5.5,
    'D1_MPS': 3.0,
    'D1_MPS+TDI': 3.0,
    'D1_MPS+IAG': 4,
    'D1.33_LM': 3.5,
    'D1.33_LM+TDI': 3.75,
    'D1.33_LM+IAG': 5,
    'D1.67_LM': 4.0,
    'D1.67_LM+TDI': 5.5,
    'D1.67_LM+IAG': 3.5,
    'D2_LM': 6.5,
    'D2_LM+TDI': 8,
    'D2_LM+IAG': 3.5,
    'D3_AE': 3.25,
    'D3_AE+TDI': 7.0,
    'D3_AE+IAG': 3.75,
    'D2_LM+TDI->D3_AE+TDI': 7,
    'D2_LM+IAG->D3_AE+IAG': 3.75,
}

# Same-timepoint treated vs control contrasts.
contrasts = pd.DataFrame([
    ('D1_MPS+TDI', 'D1_MPS', 'TDI_at_D1_MPS'),
    ('D1_MPS+IAG', 'D1_MPS', 'IAG_at_D1_MPS'),
    ('D1.33_LM+TDI', 'D1.33_LM', 'TDI_at_D1.33_LM'),
    ('D1.33_LM+IAG', 'D1.33_LM', 'IAG_at_D1.33_LM'),
    ('D1.67_LM+TDI', 'D1.67_LM', 'TDI_at_D1.67_LM'),
    ('D1.67_LM+IAG', 'D1.67_LM', 'IAG_at_D1.67_LM'),
    ('D2_LM+TDI', 'D2_LM', 'TDI_at_D2_LM'),
    ('D2_LM+IAG', 'D2_LM', 'IAG_at_D2_LM'),
    ('D3_AE+TDI', 'D3_AE', 'TDI_at_D3_AE'),
    ('D3_AE+IAG', 'D3_AE', 'IAG_at_D3_AE'),
    ('D2_LM+TDI->D3_AE+TDI', 'D3_AE', 'TDI_D2_to_D3_vs_D3_AE'),
    ('D2_LM+IAG->D3_AE+IAG', 'D3_AE', 'IAG_D2_to_D3_vs_D3_AE'),
], columns=['treated_condition', 'control_condition', 'contrast_name'])

print("Base output:", BASE_OUT)
print("Contrasts:")
display(contrasts)

# %%
# CELL 2 — Helper functions

def infer_condition(sample_name: str) -> str:
    mapped = sample_mapping_dict.get(str(sample_name), str(sample_name))
    return re.sub(r"_rep\d+$", "", mapped)


def infer_replicate(sample_name: str) -> str:
    mapped = sample_mapping_dict.get(str(sample_name), str(sample_name))
    match = re.search(r"_rep(\d+)$", mapped)
    return f"rep{match.group(1)}" if match else "repNA"


def infer_treatment(condition: str) -> str:
    condition = str(condition)
    if "+TDI" in condition:
        return "TDI"
    if "+IAG" in condition:
        return "IAG"
    return "Control"


def get_counts_chunk(adata, start: int, end: int, source: str):
    sub = adata[start:end, :]
    if source == "X":
        return sub.X
    if source.startswith("layer:"):
        layer = source.split(":", 1)[1]
        if layer not in sub.layers:
            raise KeyError(f"Layer {layer!r} not found. Available layers: {list(sub.layers.keys())}")
        return sub.layers[layer]
    raise ValueError("COUNTS_SOURCE must be 'X' or 'layer:<layer_name>'")


def as_1d_sum(x):
    if sparse.issparse(x):
        return np.asarray(x.sum(axis=0)).ravel()
    return np.asarray(x).sum(axis=0).ravel()


def bh_adjust(pvals):
    pvals = np.asarray(pvals, dtype=float)
    out = np.full_like(pvals, np.nan, dtype=float)
    valid = np.isfinite(pvals)
    p = pvals[valid]
    if p.size == 0:
        return out
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * p.size / np.arange(1, p.size + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    tmp = np.empty_like(adj)
    tmp[order] = adj
    out[valid] = tmp
    return out


def safe_name(x):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x)).strip("_")

# %%
# CELL 3 — Load AnnData, QC-filter cells, export pseudobulk gene counts and sample-level cNMF scores

adata = sc.read_h5ad(H5AD_PATH, backed='r')
obs = adata.obs.copy()

if SAMPLE_COL not in obs.columns:
    raise KeyError(f"{SAMPLE_COL!r} not found in adata.obs. Available: {list(obs.columns)}")
for col in ["pct_counts_mt", "pct_counts_ribo"]:
    if col not in obs.columns:
        raise KeyError(f"Required QC column {col!r} not found in adata.obs")

raw_sample = obs[SAMPLE_COL].astype(str)
obs['sample_name_de'] = raw_sample.map(lambda x: sample_mapping_dict.get(x, x))
obs['condition_de'] = obs['sample_name_de'].map(infer_condition)
obs['replicate_de'] = obs['sample_name_de'].map(infer_replicate)
obs['treatment_de'] = obs['condition_de'].map(infer_treatment)

conditions = sorted(obs['condition_de'].unique())
missing_cutoffs = sorted(set(conditions) - set(cond_pct_mt_cutoffs))
if missing_cutoffs:
    raise ValueError(f"Conditions missing QC cutoffs: {missing_cutoffs}")

mt = pd.to_numeric(obs['pct_counts_mt'], errors='coerce').to_numpy()
ribo = pd.to_numeric(obs['pct_counts_ribo'], errors='coerce').to_numpy()
cond_array = obs['condition_de'].astype(str).to_numpy()

keep = np.ones(adata.n_obs, dtype=bool)
filter_records = []
for cond in conditions:
    cond_mask = cond_array == cond
    fail_mt = cond_mask & (mt >= cond_pct_mt_cutoffs[cond])
    fail_ribo = cond_mask & (ribo <= cond_pct_ribo_cutoffs[cond])
    fail_any = fail_mt | fail_ribo
    keep &= ~fail_any
    n = int(cond_mask.sum())
    filter_records.append({
        'condition_de': cond,
        'n_cells_before': n,
        'n_fail_mt': int(fail_mt.sum()),
        'n_fail_low_ribo': int(fail_ribo.sum()),
        'n_removed_any': int(fail_any.sum()),
        'n_cells_after': int(n - fail_any.sum()),
        'fraction_removed': float(fail_any.sum() / n) if n else np.nan,
        'pct_counts_mt_cutoff': cond_pct_mt_cutoffs[cond],
        'pct_counts_ribo_cutoff': cond_pct_ribo_cutoffs[cond],
    })

filter_summary = pd.DataFrame(filter_records)
filter_summary.to_csv(BASE_OUT / 'qc_filter_summary_by_condition.tsv', sep='\t', index=False)
print(f"Cells before filter: {adata.n_obs:,}")
print(f"Cells after filter:  {keep.sum():,}")
print(f"Removed:             {(~keep).sum():,} ({(~keep).mean():.2%})")
display(filter_summary)

obs_f = obs.loc[keep].copy()

# Sample order follows the provided mapping dictionary.
mapped_order = list(sample_mapping_dict.values())
sample_order = [x for x in mapped_order if x in set(obs_f['sample_name_de'])]
sample_order += [x for x in obs_f['sample_name_de'].drop_duplicates() if x not in sample_order]
n_samples = len(sample_order)
sample_to_idx = {s: i for i, s in enumerate(sample_order)}
sample_codes_all = obs['sample_name_de'].map(sample_to_idx).fillna(-1).astype(int).to_numpy()

sample_meta_records = []
for sample in sample_order:
    mask_sample = obs_f['sample_name_de'] == sample
    cond = infer_condition(sample)
    sample_meta_records.append({
        'sample_name_de': sample,
        'condition_de': cond,
        'replicate_de': infer_replicate(sample),
        'treatment_de': infer_treatment(cond),
        'n_cells_after_filter': int(mask_sample.sum()),
        'mean_pct_counts_mt': float(pd.to_numeric(obs_f.loc[mask_sample, 'pct_counts_mt'], errors='coerce').mean()),
        'mean_pct_counts_ribo': float(pd.to_numeric(obs_f.loc[mask_sample, 'pct_counts_ribo'], errors='coerce').mean()),
    })

sample_meta = pd.DataFrame(sample_meta_records).set_index('sample_name_de')
sample_meta.to_csv(BASE_OUT / 'pseudobulk_sample_metadata.filtered.tsv', sep='\t')
contrasts.to_csv(BASE_OUT / 'matched_contrasts.tsv', sep='\t', index=False)

display(sample_meta)

# Gene IDs and annotations.
gene_id = adata.var_names.astype(str).to_numpy()
gene_annot = pd.DataFrame({'gene_id': gene_id})
if GENE_SYMBOL_COL in adata.var.columns:
    gene_annot['gene_symbol'] = adata.var[GENE_SYMBOL_COL].astype(str).to_numpy()
gene_annot.to_csv(BASE_OUT / 'gene_annotation.tsv', sep='\t', index=False)

# Pseudobulk raw counts: genes x samples output.
n_genes = adata.n_vars
counts = np.zeros((n_samples, n_genes), dtype=np.float64)

print(f"Aggregating pseudobulk counts from {COUNTS_SOURCE!r} in chunks of {CHUNK_SIZE} cells...")
for start in range(0, adata.n_obs, CHUNK_SIZE):
    end = min(start + CHUNK_SIZE, adata.n_obs)
    chunk_keep = keep[start:end]
    if not np.any(chunk_keep):
        continue
    X = get_counts_chunk(adata, start, end, COUNTS_SOURCE)
    if sparse.issparse(X):
        X = X.tocsr()[chunk_keep, :]
    else:
        X = np.asarray(X)[chunk_keep, :]
    codes = sample_codes_all[start:end][chunk_keep]
    for code in np.unique(codes):
        if code < 0 or code >= n_samples:
            continue
        counts[code, :] += as_1d_sum(X[codes == code, :])
    if start % (CHUNK_SIZE * 10) == 0:
        print(f"  processed {end:,}/{adata.n_obs:,} cells", flush=True)

max_abs_round_error = float(np.max(np.abs(counts - np.round(counts))))
print(f"Max abs count rounding error: {max_abs_round_error:.3g}")
if max_abs_round_error > 1e-4:
    print("WARNING: counts are not integer-like. Rounding for PyDESeq2; verify COUNTS_SOURCE.")

counts_int = np.rint(counts).astype(np.int64)
counts_df = pd.DataFrame(counts_int.T, index=gene_id, columns=sample_order)
counts_path = BASE_OUT / 'pseudobulk_gene_counts.filtered.tsv.gz'
counts_df.to_csv(counts_path, sep='\t', compression='gzip')
print("Wrote pseudobulk counts:", counts_path)

# Sample-level projected cNMF program scores.
print(f"Aggregating projected scores from obsm[{PROGRAM_SCORE_KEY!r}]...")
if PROGRAM_SCORE_KEY not in adata.obsm:
    raise KeyError(f"{PROGRAM_SCORE_KEY!r} not found in adata.obsm")
usage = np.asarray(adata.obsm[PROGRAM_SCORE_KEY])
projection_uns = adata.uns.get('perturbseq_cNMF_projection', {})
program_names = list(projection_uns.get('program_names', [])) if isinstance(projection_uns, dict) else []
if len(program_names) != usage.shape[1]:
    program_names = [f'P{i}' for i in range(usage.shape[1])]

score_mean = np.zeros((n_samples, usage.shape[1]), dtype=np.float64)
score_median = np.zeros((n_samples, usage.shape[1]), dtype=np.float64)
for sample, idx in sample_to_idx.items():
    cell_mask = keep & (obs['sample_name_de'].to_numpy() == sample)
    if cell_mask.sum() == 0:
        score_mean[idx, :] = np.nan
        score_median[idx, :] = np.nan
    else:
        score_mean[idx, :] = usage[cell_mask, :].mean(axis=0)
        score_median[idx, :] = np.median(usage[cell_mask, :], axis=0)

score_mean_df = pd.DataFrame(score_mean, index=sample_order, columns=program_names)
score_median_df = pd.DataFrame(score_median, index=sample_order, columns=program_names)
score_mean_df.to_csv(BASE_OUT / 'pseudobulk_projected_cNMF_mean_usage.filtered.tsv', sep='\t')
score_median_df.to_csv(BASE_OUT / 'pseudobulk_projected_cNMF_median_usage.filtered.tsv', sep='\t')
print("Wrote pseudobulk inputs to:", BASE_OUT)

adata.file.close()

# %%
# CELL 4 — Run matched pseudobulk gene DE with PyDESeq2
#
# Install first if needed:
#     python -m pip install "pydeseq2>=0.5.4"

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

DE_OUT = BASE_OUT / 'pydeseq2_gene_DE'
DE_OUT.mkdir(parents=True, exist_ok=True)

counts_gene_by_sample = pd.read_csv(BASE_OUT / 'pseudobulk_gene_counts.filtered.tsv.gz', sep='\t', index_col=0)
meta_all = pd.read_csv(BASE_OUT / 'pseudobulk_sample_metadata.filtered.tsv', sep='\t', index_col=0)
contrasts_df = pd.read_csv(BASE_OUT / 'matched_contrasts.tsv', sep='\t')
gene_annot = pd.read_csv(BASE_OUT / 'gene_annotation.tsv', sep='\t')

# Align samples and transpose to samples x genes, as expected by PyDESeq2.
common_samples = [s for s in counts_gene_by_sample.columns if s in meta_all.index]
counts_sample_by_gene = counts_gene_by_sample.loc[:, common_samples].T.round().astype(np.int64)
meta_all = meta_all.loc[common_samples, :].copy()

print("counts_sample_by_gene:", counts_sample_by_gene.shape)
print("meta_all:", meta_all.shape)

def run_one_pydeseq2_contrast(treated_condition, control_condition, contrast_name):
    print(f"\n=== {contrast_name}: {treated_condition} vs {control_condition} ===")

    samples = meta_all.index[meta_all['condition_de'].isin([control_condition, treated_condition])].tolist()
    meta_sub = meta_all.loc[samples, :].copy()

    use_paired = False
    if USE_PAIRED_REPLICATE_DESIGN:
        reps_control = set(meta_sub.loc[meta_sub['condition_de'] == control_condition, 'replicate_de'].astype(str))
        reps_treated = set(meta_sub.loc[meta_sub['condition_de'] == treated_condition, 'replicate_de'].astype(str))
        common_reps = sorted(reps_control & reps_treated)
        if len(common_reps) >= MIN_SAMPLES_PER_GROUP:
            meta_sub = meta_sub.loc[meta_sub['replicate_de'].astype(str).isin(common_reps), :].copy()
            use_paired = True

    n_control = int((meta_sub['condition_de'] == control_condition).sum())
    n_treated = int((meta_sub['condition_de'] == treated_condition).sum())
    print(f"samples control={n_control}, treated={n_treated}, paired_design={use_paired}")

    if n_control < MIN_SAMPLES_PER_GROUP or n_treated < MIN_SAMPLES_PER_GROUP:
        print("Skipping: not enough samples")
        return None, {
            'contrast_name': contrast_name,
            'treated_condition': treated_condition,
            'control_condition': control_condition,
            'n_control_samples': n_control,
            'n_treated_samples': n_treated,
            'paired_design': use_paired,
            'skipped': True,
            'skip_reason': 'not_enough_samples',
        }

    sample_ids = meta_sub.index.tolist()
    counts_sub = counts_sample_by_gene.loc[sample_ids, :].copy()

    keep_genes = ((counts_sub >= MIN_COUNT).sum(axis=0) >= MIN_SAMPLES_EXPRESSED) & (counts_sub.sum(axis=0) >= MIN_COUNT)
    print(f"genes before filter={counts_sub.shape[1]}, after filter={int(keep_genes.sum())}")
    counts_sub = counts_sub.loc[:, keep_genes]

    coldata = meta_sub.copy()
    coldata['condition_group'] = np.where(coldata['condition_de'] == treated_condition, 'treated', 'control')
    coldata['condition_group'] = pd.Categorical(coldata['condition_group'], categories=['control', 'treated'])
    coldata['replicate_de'] = pd.Categorical(coldata['replicate_de'].astype(str))

    design = '~ replicate_de + condition_group' if use_paired else '~ condition_group'
    print("design:", design)

    inference = DefaultInference(n_cpus=N_CPUS)
    dds = DeseqDataSet(
        counts=counts_sub,
        metadata=coldata,
        design=design,
        refit_cooks=False,
        size_factors_fit_type=SIZE_FACTORS_FIT_TYPE,
        inference=inference,
        quiet=QUIET_PYDESEQ2,
        low_memory=True,
    )
    dds.deseq2()

    stat_res = DeseqStats(
        dds,
        contrast=['condition_group', 'treated', 'control'],
        alpha=0.05,
        cooks_filter=True,
        independent_filter=True,
        inference=inference,
        quiet=QUIET_PYDESEQ2,
    )
    stat_res.summary()

    res_df = stat_res.results_df.copy().reset_index().rename(columns={'index': 'gene_id'})
    if 'gene_id' not in res_df.columns:
        res_df = res_df.rename(columns={res_df.columns[0]: 'gene_id'})

    res_df['contrast_name'] = contrast_name
    res_df['treated_condition'] = treated_condition
    res_df['control_condition'] = control_condition
    res_df['n_control_samples'] = n_control
    res_df['n_treated_samples'] = n_treated
    res_df['paired_design'] = use_paired
    res_df['design'] = design
    res_df['method'] = 'PyDESeq2'

    if 'gene_id' in gene_annot.columns:
        res_df = res_df.merge(gene_annot, on='gene_id', how='left', sort=False)

    preferred_cols = [
        'contrast_name', 'treated_condition', 'control_condition', 'gene_id', 'gene_symbol',
        'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj',
        'n_control_samples', 'n_treated_samples', 'paired_design', 'design', 'method',
    ]
    extra_cols = [c for c in res_df.columns if c not in preferred_cols]
    res_df = res_df[[c for c in preferred_cols if c in res_df.columns] + extra_cols]

    sig = int(((res_df['padj'] < 0.05) & res_df['padj'].notna()).sum()) if 'padj' in res_df else 0
    sig_up = int(((res_df['padj'] < 0.05) & (res_df['log2FoldChange'] > 0) & res_df['padj'].notna()).sum())
    sig_down = int(((res_df['padj'] < 0.05) & (res_df['log2FoldChange'] < 0) & res_df['padj'].notna()).sum())

    summary = {
        'contrast_name': contrast_name,
        'treated_condition': treated_condition,
        'control_condition': control_condition,
        'n_control_samples': n_control,
        'n_treated_samples': n_treated,
        'paired_design': use_paired,
        'design': design,
        'n_genes_tested': int(res_df.shape[0]),
        'n_sig_padj_0_05': sig,
        'n_sig_up': sig_up,
        'n_sig_down': sig_down,
        'skipped': False,
        'skip_reason': '',
    }
    return res_df, summary

all_de_summaries = []
for _, row in contrasts_df.iterrows():
    res_df, summary = run_one_pydeseq2_contrast(
        treated_condition=str(row['treated_condition']),
        control_condition=str(row['control_condition']),
        contrast_name=str(row['contrast_name']),
    )
    all_de_summaries.append(summary)
    if res_df is not None:
        out_path = DE_OUT / f"PyDESeq2_{row['contrast_name']}.tsv"
        res_df.to_csv(out_path, sep='\t', index=False)
        print("wrote", out_path)

de_summary_df = pd.DataFrame(all_de_summaries)
de_summary_df.to_csv(DE_OUT / 'PyDESeq2_contrast_summary.tsv', sep='\t', index=False)
display(de_summary_df)

# %%
# CELL 5 — Matched treatment-vs-control tests for projected cNMF program scores

PROGRAM_OUT = BASE_OUT / 'program_score_DE'
PROGRAM_OUT.mkdir(parents=True, exist_ok=True)

scores = pd.read_csv(BASE_OUT / 'pseudobulk_projected_cNMF_mean_usage.filtered.tsv', sep='\t', index_col=0)
meta = pd.read_csv(BASE_OUT / 'pseudobulk_sample_metadata.filtered.tsv', sep='\t', index_col=0)
contrasts_df = pd.read_csv(BASE_OUT / 'matched_contrasts.tsv', sep='\t')

scores = scores.loc[meta.index]
program_names = list(scores.columns)

records = []
paired_records = []
EPS = 1e-8

for _, contrast in contrasts_df.iterrows():
    cname = contrast['contrast_name']
    trt = contrast['treated_condition']
    ctrl = contrast['control_condition']

    ctrl_samples = meta.index[meta['condition_de'] == ctrl].tolist()
    trt_samples = meta.index[meta['condition_de'] == trt].tolist()
    reps_ctrl = set(meta.loc[ctrl_samples, 'replicate_de'])
    reps_trt = set(meta.loc[trt_samples, 'replicate_de'])
    common_reps = sorted(reps_ctrl & reps_trt)
    use_paired = len(common_reps) >= 2

    for program in program_names:
        ctrl_vals = scores.loc[ctrl_samples, program].astype(float)
        trt_vals = scores.loc[trt_samples, program].astype(float)

        mean_ctrl = float(ctrl_vals.mean())
        mean_trt = float(trt_vals.mean())
        diff = mean_trt - mean_ctrl
        log2_ratio = float(np.log2((mean_trt + EPS) / (mean_ctrl + EPS)))

        if use_paired:
            paired_diffs = []
            for rep in common_reps:
                c_s = meta.index[(meta['condition_de'] == ctrl) & (meta['replicate_de'] == rep)].tolist()
                t_s = meta.index[(meta['condition_de'] == trt) & (meta['replicate_de'] == rep)].tolist()
                if len(c_s) == 1 and len(t_s) == 1:
                    d = float(scores.loc[t_s[0], program] - scores.loc[c_s[0], program])
                    paired_diffs.append(d)
                    paired_records.append({
                        'contrast_name': cname,
                        'program': program,
                        'replicate_de': rep,
                        'control_sample': c_s[0],
                        'treated_sample': t_s[0],
                        'control_value': float(scores.loc[c_s[0], program]),
                        'treated_value': float(scores.loc[t_s[0], program]),
                        'paired_diff': d,
                    })
            paired_diffs = np.asarray(paired_diffs)
            if paired_diffs.size >= 2 and np.std(paired_diffs) > 0:
                stat, pval = stats.ttest_1samp(paired_diffs, popmean=0.0)
            else:
                stat, pval = np.nan, np.nan
            test = 'paired_ttest_on_replicate_differences'
            n_pairs = int(paired_diffs.size)
        else:
            if len(ctrl_vals) >= 2 and len(trt_vals) >= 2:
                stat, pval = stats.ttest_ind(trt_vals, ctrl_vals, equal_var=False, nan_policy='omit')
            else:
                stat, pval = np.nan, np.nan
            test = 'welch_ttest'
            n_pairs = 0

        records.append({
            'contrast_name': cname,
            'treated_condition': trt,
            'control_condition': ctrl,
            'program': program,
            'mean_control': mean_ctrl,
            'mean_treated': mean_trt,
            'diff_treated_minus_control': diff,
            'log2_ratio_treated_vs_control': log2_ratio,
            'statistic': stat,
            'pvalue': pval,
            'test': test,
            'n_control_samples': len(ctrl_samples),
            'n_treated_samples': len(trt_samples),
            'n_matched_pairs': n_pairs,
        })

program_res = pd.DataFrame(records)
program_res['padj_bh_within_contrast'] = np.nan
for cname, idx in program_res.groupby('contrast_name').groups.items():
    program_res.loc[idx, 'padj_bh_within_contrast'] = bh_adjust(program_res.loc[idx, 'pvalue'].to_numpy())

program_res = program_res.sort_values(['contrast_name', 'padj_bh_within_contrast', 'pvalue', 'program'])
program_res.to_csv(PROGRAM_OUT / 'projected_cNMF_program_score_matched_contrasts.tsv', sep='\t', index=False)

if paired_records:
    pd.DataFrame(paired_records).to_csv(PROGRAM_OUT / 'projected_cNMF_program_score_paired_replicate_differences.tsv', sep='\t', index=False)

key_program_res = program_res[program_res['program'].isin(KEY_PROGRAMS)].copy()
key_program_res.to_csv(PROGRAM_OUT / 'key_programs_P15_P17_P34_matched_contrasts.tsv', sep='\t', index=False)
display(key_program_res)

# Plot key programs.
for program in KEY_PROGRAMS:
    if program not in program_names:
        print(f"Skipping {program}: not in program score table")
        continue
    prog_res = program_res[program_res['program'] == program].set_index('contrast_name')
    fig, ax = plt.subplots(figsize=(max(8, 0.5 * len(prog_res)), 3.5))
    x = np.arange(len(prog_res))
    ax.axhline(0, linewidth=0.8)
    ax.bar(x, prog_res['diff_treated_minus_control'].to_numpy())
    ax.set_xticks(x)
    ax.set_xticklabels(prog_res.index.tolist(), rotation=90)
    ax.set_ylabel(f'{program} mean usage difference\n(treated - control)')
    ax.set_title(f'{program} matched treatment effects')
    fig.tight_layout()
    fig.savefig(PROGRAM_OUT / f'{program}_matched_treatment_effects.svg')
    fig.savefig(PROGRAM_OUT / f'{program}_matched_treatment_effects.png')
    plt.show()
    plt.close(fig)

print('Wrote program-score contrast outputs to:', PROGRAM_OUT)

# %%
# CELL 6 — Optional: cNMF program enrichment from PyDESeq2 gene log2FC
#
# This checks whether genes with high weights in a perturb-seq cNMF program tend
# to increase or decrease in each DE contrast.

ENRICH_OUT = BASE_OUT / 'program_enrichment_from_PyDESeq2_gene_DE'
ENRICH_OUT.mkdir(parents=True, exist_ok=True)
TOP_N = 200

spectra = pd.read_csv(SPECTRA_TPM_PATH, sep='\t', index_col=0)
spectra.index = [f'P{int(x) - 1}' for x in spectra.index.astype(str)]  # source rows 1..50 -> P0..P49
spectra.columns = spectra.columns.astype(str)
spectra = spectra.apply(pd.to_numeric, errors='coerce').fillna(0)
spectra[spectra < 0] = 0

result_files = sorted(DE_OUT.glob('PyDESeq2_*.tsv'))
if not result_files:
    raise FileNotFoundError(f"No PyDESeq2_*.tsv files found in {DE_OUT}")

all_records = []
for path in result_files:
    if path.name == 'PyDESeq2_contrast_summary.tsv':
        continue
    de = pd.read_csv(path, sep='\t')
    contrast = path.stem.replace('PyDESeq2_', '')
    print('Reading:', path.name)

    if 'gene_symbol' in de.columns:
        gene_col = 'gene_symbol'
    elif 'gene_id' in de.columns:
        gene_col = 'gene_id'
    else:
        print(f'Skipping {path.name}: no gene_symbol or gene_id column')
        continue

    if 'log2FoldChange' not in de.columns:
        print(f'Skipping {path.name}: no log2FoldChange column')
        continue

    de = de.dropna(subset=['log2FoldChange']).copy()
    de[gene_col] = de[gene_col].astype(str)
    de = de.drop_duplicates(gene_col, keep='first').set_index(gene_col)

    common = de.index.intersection(spectra.columns)
    print(f"  shared genes with spectra: {len(common):,}")
    if len(common) == 0:
        continue

    lfc = de.loc[common, 'log2FoldChange'].astype(float)
    stat = de.loc[common, 'stat'].astype(float) if 'stat' in de.columns else lfc

    for program in spectra.index:
        weights = spectra.loc[program, common].astype(float)
        weight_sum = weights.sum()
        if weight_sum <= 0:
            continue

        weighted_lfc = float(np.sum(weights * lfc) / weight_sum)
        weighted_stat = float(np.sum(weights * stat) / weight_sum)

        top_genes = weights.sort_values(ascending=False).head(TOP_N).index
        top_lfc_mean = float(lfc.loc[top_genes].mean())
        top_stat_mean = float(stat.loc[top_genes].mean())
        rho, pval = stats.spearmanr(weights.to_numpy(), lfc.to_numpy(), nan_policy='omit')

        all_records.append({
            'method': 'PyDESeq2',
            'contrast_name': contrast,
            'program': program,
            'n_shared_genes': len(common),
            'weighted_mean_log2FC_by_spectra': weighted_lfc,
            'weighted_mean_stat_by_spectra': weighted_stat,
            f'top_{TOP_N}_genes_mean_log2FC': top_lfc_mean,
            f'top_{TOP_N}_genes_mean_stat': top_stat_mean,
            'spearman_weight_vs_log2FC': rho,
            'spearman_pvalue': pval,
            f'top_{TOP_N}_genes': ','.join(top_genes.tolist()),
        })

enrich_res = pd.DataFrame(all_records)
if enrich_res.empty:
    raise RuntimeError('No enrichment records were generated')

enrich_res['spearman_padj_bh_within_contrast'] = np.nan
for contrast, idx in enrich_res.groupby('contrast_name').groups.items():
    enrich_res.loc[idx, 'spearman_padj_bh_within_contrast'] = bh_adjust(enrich_res.loc[idx, 'spearman_pvalue'].to_numpy())

enrich_res = enrich_res.sort_values(['contrast_name', 'weighted_mean_log2FC_by_spectra'], ascending=[True, False])
enrich_res.to_csv(ENRICH_OUT / 'program_enrichment_from_PyDESeq2_log2FC.tsv', sep='\t', index=False)

key_enrich = enrich_res[enrich_res['program'].isin(KEY_PROGRAMS)].copy()
key_enrich.to_csv(ENRICH_OUT / 'key_programs_P15_P17_P34_enrichment_from_PyDESeq2.tsv', sep='\t', index=False)

display(key_enrich)
print('Wrote enrichment outputs to:', ENRICH_OUT)
