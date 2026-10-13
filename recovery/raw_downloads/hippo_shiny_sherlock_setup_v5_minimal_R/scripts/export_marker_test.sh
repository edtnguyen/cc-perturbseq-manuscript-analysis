#!/usr/bin/env bash
set -euo pipefail

BASE=${BASE:-/scratch/users/edtng/Hippo_Shiny_adata}
PYTHON=${PYTHON:-python}
REPO_ROOT=${REPO_ROOT:-/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE}
MDATA=${MDATA:-${BASE}/raw/cNMF_50_0_2.h5mu}
ADATA_INT=${ADATA_INT:-/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad}
OUT=${OUT:-${BASE}/data/umap_gene_program_data.h5}
GENE_LIST=${GENE_LIST:-${BASE}/scripts/markers_for_shiny.txt}
OBS_COLUMNS=${OBS_COLUMNS:-day,cell_type}
REQUIRE_OBS_COLUMNS=${REQUIRE_OBS_COLUMNS:-day}
LOG=${LOG:-${BASE}/logs/export_marker_test.$(date +%Y%m%d_%H%M%S).log}

mkdir -p "${BASE}/data" "${BASE}/logs" "${BASE}/scripts"

if [[ ! -f "${GENE_LIST}" ]]; then
  cat > "${GENE_LIST}" <<'GENES'
KDR
PECAM1
VWF
CDH5
SOX17
FLT1
ERG
TEK
GATA2
APLNR
GENES
fi

{
  echo "Python:    ${PYTHON}"
  echo "mdata:     ${MDATA}"
  echo "adata_int: ${ADATA_INT}"
  echo "repo:      ${REPO_ROOT}"
  echo "out:       ${OUT}"
  echo "genes:     ${GENE_LIST}"
  echo "log:       ${LOG}"
  echo
} | tee "${LOG}"

PYTHONUNBUFFERED=1 "${PYTHON}" -u "${BASE}/scripts/export_umap_shiny_h5.py" \
  --repo-root "${REPO_ROOT}" \
  --mdata "${MDATA}" \
  --adata-int "${ADATA_INT}" \
  --allow-extra-adata-cells \
  --obs-columns "${OBS_COLUMNS}" \
  --required-obs-columns "${REQUIRE_OBS_COLUMNS}" \
  --gene-mode list \
  --gene-list "${GENE_LIST}" \
  --out "${OUT}" 2>&1 | tee -a "${LOG}"

ls -lh "${OUT}" | tee -a "${LOG}"
