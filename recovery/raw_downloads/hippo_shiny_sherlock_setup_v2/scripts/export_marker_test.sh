#!/usr/bin/env bash
set -euo pipefail

BASE=${BASE:-/scratch/users/edtng/Hippo_Shiny_adata}
REPO_ROOT=${REPO_ROOT:-/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE}
MDATA=${MDATA:-${BASE}/raw/cNMF_50_0_2.h5mu}
if [[ ! -f "${MDATA}" && -f "${BASE}/cNMF_50_0_2.h5mu" ]]; then
  MDATA="${BASE}/cNMF_50_0_2.h5mu"
fi
ADATA_INT=${ADATA_INT:-/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad}
OUT=${OUT:-${BASE}/data/umap_gene_program_data.h5}
GENE_LIST=${GENE_LIST:-${BASE}/markers_for_shiny.example.txt}
PYTHON=${PYTHON:-python}

mkdir -p "${BASE}/data" "${BASE}/logs"

echo "Python:    ${PYTHON}"
echo "mdata:     ${MDATA}"
echo "adata_int: ${ADATA_INT}"
echo "repo:      ${REPO_ROOT}"
echo "out:       ${OUT}"
echo "genes:     ${GENE_LIST}"

"${PYTHON}" "${BASE}/scripts/export_umap_shiny_h5.py" \
  --repo-root "${REPO_ROOT}" \
  --mdata "${MDATA}" \
  --adata-int "${ADATA_INT}" \
  --allow-extra-adata-cells \
  --gene-mode list \
  --gene-list "${GENE_LIST}" \
  --out "${OUT}"

ls -lh "${OUT}"
