# %% [markdown]
# # Relabel projected perturb-seq cNMF programs from 1-based to 0-based
#
# Use this after the projection has already run and the notebook kernel still has
# `adata`, `usage_raw`, and `usage_norm` in memory.
#
# This does not recompute the projection. It only changes labels:
#
# - old `P1` -> new `P0`
# - old `P2` -> new `P1`
# - ...
# - old `P50` -> new `P49`

# %%
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

try:
    import muon as mu
except Exception as e:
    mu = None
    print("WARNING: muon could not be imported; .h5mu saving will be skipped.")
    print(e)

OUT_DIR = Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_projection")
ADATA_PATH = Path("/scratch/users/edtng/20260325_Hippo_timecourse/Hippo_timecourse_combined.h5ad")

SPECTRA_TPM_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50."
    "gene_spectra_tpm.k_50.dt_0_2.txt"
)

PROGRAM_META_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/data/program_activity_long_v4.csv"
)

OUT_H5AD_ZERO = OUT_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5ad"
OUT_H5MU_ZERO = OUT_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5mu"
OUT_SCORES_TSV_ZERO = OUT_DIR / "Hippo_timecourse_combined.perturbseq_cNMF_50_dt_0_2.program_scores.zero_based.tsv.gz"

# %%
required_vars = ["adata", "usage_raw", "usage_norm"]
missing = [v for v in required_vars if v not in globals()]
if missing:
    raise RuntimeError(
        "This relabel cell expects the projection notebook to still have these variables in memory: "
        + ", ".join(required_vars)
        + ". Missing: "
        + ", ".join(missing)
        + ". If the kernel was restarted, use the zero-based full projection notebook/script instead."
    )

assert usage_norm.shape[1] == 50, usage_norm.shape
assert usage_raw.shape[1] == 50, usage_raw.shape

program_ids_zero = list(range(50))
program_names_zero = [f"P{i}" for i in program_ids_zero]
spectra_source_programs_one_based = [str(i) for i in range(1, 51)]

print("Zero-based program names:", program_names_zero[:10], "...", program_names_zero[-5:])

# %%
# Rename old obs columns P1..P50 to P0..P49.
old_obs_cols = [f"perturbseq_cNMF_P{i}_usage" for i in range(1, 51)]
new_obs_cols = [f"perturbseq_cNMF_P{i}_usage" for i in range(0, 50)]

rename_map = {
    old: new
    for old, new in zip(old_obs_cols, new_obs_cols)
    if old in adata.obs.columns
}

print("Obs score columns to rename:", len(rename_map))
print(list(rename_map.items())[:5], "...", list(rename_map.items())[-5:])

adata.obs.rename(columns=rename_map, inplace=True)

# Guarantee that the zero-based columns match the normalized usage matrix.
for j, col in enumerate(new_obs_cols):
    adata.obs[col] = usage_norm[:, j]

# Drop old columns if any remain.
remaining_old = [c for c in old_obs_cols if c in adata.obs.columns and c not in new_obs_cols]
if remaining_old:
    adata.obs.drop(columns=remaining_old, inplace=True)

print("New obs score columns:", sum(c.startswith("perturbseq_cNMF_P") and c.endswith("_usage") for c in adata.obs.columns))
print([c for c in adata.obs.columns if c.startswith("perturbseq_cNMF_P")][:10])

# %%
# Build program label map from 0-based program metadata.
program_label_map = {p: p for p in program_names_zero}

if PROGRAM_META_PATH.exists():
    program_meta = pd.read_csv(PROGRAM_META_PATH).drop_duplicates("program_name")
    if "program_label" in program_meta.columns:
        for _, row in program_meta.iterrows():
            if pd.isna(row["program_name"]):
                continue
            pid = int(row["program_name"])
            pname = f"P{pid}"
            program_label_map[pname] = str(row["program_label"])
    print("Loaded metadata labels from:", PROGRAM_META_PATH)
else:
    print("Program metadata not found:", PROGRAM_META_PATH)

# Update AnnData metadata.
if "perturbseq_cNMF_projection" not in adata.uns:
    adata.uns["perturbseq_cNMF_projection"] = {}

adata.uns["perturbseq_cNMF_projection"]["program_names"] = program_names_zero
adata.uns["perturbseq_cNMF_projection"]["program_ids_zero_based"] = program_ids_zero
adata.uns["perturbseq_cNMF_projection"]["spectra_source_programs_one_based"] = spectra_source_programs_one_based
adata.uns["perturbseq_cNMF_projection"]["program_labels"] = [
    program_label_map.get(p, p) for p in program_names_zero
]
adata.uns["perturbseq_cNMF_projection"]["program_index_note"] = (
    "Programs are labeled with 0-based IDs P0..P49. "
    "The cNMF spectra_tpm file rows were originally labeled 1..50; source row 1 maps to P0."
)

# %%
# Save zero-based score table.
program_score_df_zero = pd.DataFrame(
    usage_norm,
    index=adata.obs_names.astype(str),
    columns=program_names_zero,
)

program_score_df_zero.to_csv(OUT_SCORES_TSV_ZERO, sep="\t", compression="gzip")
print("Saved zero-based score table:", OUT_SCORES_TSV_ZERO)
display(program_score_df_zero.iloc[:5, :10])

# %%
# Save zero-based projected AnnData.
# This writes a new file, leaving the original 1-based output untouched.
adata.write_h5ad(OUT_H5AD_ZERO, compression="gzip")
print("Saved zero-based AnnData:", OUT_H5AD_ZERO)

# %%
# Save zero-based MuData.
if mu is None:
    print("muon is not installed/importable; skipping .h5mu save.")
else:
    program_var = pd.DataFrame(index=pd.Index(program_names_zero, name="program"))
    program_var["program_name"] = program_names_zero
    program_var["program_id_zero_based"] = program_ids_zero
    program_var["spectra_source_program_one_based"] = spectra_source_programs_one_based
    program_var["program_label"] = [program_label_map.get(p, p) for p in program_names_zero]
    program_var["reference"] = "perturbseq_cNMF_50_dt_0_2"

    program_obs = pd.DataFrame(index=adata.obs_names.copy())
    for col in ["day", "time", "timepoint", "treatment", "condition", "sample_name", "sample_id", "biological_sample"]:
        if col in adata.obs.columns:
            program_obs[col] = adata.obs[col].values

    prog_adata = sc.AnnData(
        X=usage_norm.astype(np.float32),
        obs=program_obs,
        var=program_var,
    )
    prog_adata.layers["usage_raw"] = usage_raw.astype(np.float32)
    prog_adata.layers["usage_sum_to_1"] = usage_norm.astype(np.float32)
    prog_adata.uns["projection"] = adata.uns["perturbseq_cNMF_projection"]

    for key in ["X_umap", "X_pca"]:
        if key in adata.obsm:
            prog_adata.obsm[key] = adata.obsm[key].copy()

    mdata = mu.MuData({
        "rna": adata,
        "perturbseq_cNMF": prog_adata,
    })
    mdata.uns["description"] = (
        "Hippo inhibitor time-course scRNA-seq projected onto perturb-seq cNMF programs. "
        "Program labels are zero-based P0..P49; spectra source rows 1..50 map to P0..P49."
    )

    mdata.write_h5mu(OUT_H5MU_ZERO)
    print("Saved zero-based MuData:", OUT_H5MU_ZERO)

# %%
# Reload sanity checks.
adata_check = sc.read_h5ad(OUT_H5AD_ZERO, backed="r")
print(adata_check)
score_cols = [c for c in adata_check.obs.columns if c.startswith("perturbseq_cNMF_P") and c.endswith("_usage")]
print("score columns:", len(score_cols))
print(score_cols[:10], "...", score_cols[-5:])
print("program names in uns:", list(adata_check.uns["perturbseq_cNMF_projection"]["program_names"])[:10])
adata_check.file.close()

if OUT_H5MU_ZERO.exists() and mu is not None:
    mdata_check = mu.read_h5mu(OUT_H5MU_ZERO, backed="r")
    print(mdata_check)
    print("program modality var_names:", mdata_check.mod["perturbseq_cNMF"].var_names[:10].tolist())
    mdata_check.file.close()
