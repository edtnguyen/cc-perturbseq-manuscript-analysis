#!/bin/bash
#SBATCH --job-name=CRT_true_strat
#SBATCH --output=/oak/stanford/groups/engreitz/Users/tri/tmp/Maggie_CRT_sensitivity_test/logs/true_v2_%A_%a.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/tri/tmp/Maggie_CRT_sensitivity_test/logs/true_v2_%A_%a.err
#SBATCH --partition=bigmem,engreitz,owners
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=512G
#SBATCH --array=0-2

set -euo pipefail

BASE_DIR="/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/9.cnmf/20260505_k30_300kvicchd_d9_torch"
SENS_ROOT="/oak/stanford/groups/engreitz/Users/tri/tmp/Maggie_CRT_sensitivity_test"
GUIDE_H5AD="$SENS_ROOT/rna_with_guides_batch_sensitivity.h5ad"
CRT_SCRIPT="$SENS_ROOT/CRT_true_stratified_calibration_floored_v2.py"

# IMPORTANT: this is the patched checkout containing stratified_perm.
PROGRAMDE_DIR="/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE"
CONDA_ENV="/oak/stanford/groups/engreitz/Users/msychen/conda/envs/NMF_Benchmarking"
PYTHON="$CONDA_ENV/bin/python"

B="${B:-1024}"
N_ENSEMBLES="${N_ENSEMBLES:-20}"
N_MAX_GROUPS="${N_MAX_GROUPS:-9}"
PROPENSITY_BINS="${PROPENSITY_BINS:-20}"
SEL_THRESH="${SEL_THRESH:-0.4}"
USAGE_EPS="${USAGE_EPS:-1e-8}"

VARIANTS=(propensity_only inferred2 inferred4)
IDX="${SLURM_ARRAY_TASK_ID:-0}"
VARIANT="${VARIANTS[$IDX]}"
SAVE_DIR="$SENS_ROOT/true_stratified_${VARIANT}/sel_thresh_${SEL_THRESH//./_}"

mkdir -p "$SENS_ROOT/logs" "$SAVE_DIR"

for required in "$GUIDE_H5AD" "$CRT_SCRIPT" "$PYTHON"; do
    [[ -f "$required" ]] || { echo "ERROR: missing $required" >&2; exit 1; }
done
[[ -d "$PROGRAMDE_DIR/src/sceptre" ]] || {
    echo "ERROR: patched programDE checkout missing: $PROGRAMDE_DIR" >&2
    exit 1
}

export PATH="$CONDA_ENV/bin:$PATH"
export PYTHONPATH="$PROGRAMDE_DIR:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-20}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-20}"
export PYTHONHASHSEED=0

# Fail early if the batch job imports the wrong checkout.
"$PYTHON" - <<'PY'
import inspect
from src.sceptre.pipeline_helpers import _sample_crt_indices
from src.sceptre.ntc_groups import crt_pvals_for_ntc_groups_ensemble
print("_sample_crt_indices:", inspect.signature(_sample_crt_indices))
print("ntc ensemble:", inspect.signature(crt_pvals_for_ntc_groups_ensemble))
assert "resampling_method" in inspect.signature(_sample_crt_indices).parameters
assert "resampling_method" in inspect.signature(crt_pvals_for_ntc_groups_ensemble).parameters
PY

CMD=(
    "$PYTHON" "$CRT_SCRIPT"
    --out-dir "$BASE_DIR"
    --run-name Inference_guideumis
    --guide-h5ad "$GUIDE_H5AD"
    --save-dir "$SAVE_DIR"
    --component 30
    --sel-thresh "$SEL_THRESH"
    --outer-key sample
    --variant "$VARIANT"
    --covariates
        sublib
        pct_counts_mt
        pct_counts_ribo
        S_score
        G2M_score
        leiden_res0_25
    --log-covariates
        total_counts
        n_genes_by_counts
        guides_per_cell
        guide_umis_per_cell
    --ntc-labels non-targeting
    --guides-per-unit 15
    --B "$B"
    --propensity-bins "$PROPENSITY_BINS"
    --ntc-ensembles "$N_ENSEMBLES"
    --ntc-max-groups "$N_MAX_GROUPS"
    --usage-eps "$USAGE_EPS"
)

printf 'Running:'
printf ' %q' "${CMD[@]}"
printf '\n'
"${CMD[@]}"
