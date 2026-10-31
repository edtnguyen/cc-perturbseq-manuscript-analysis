#!/bin/bash

#SBATCH --job-name=CRT_batch_sens
#SBATCH --output=/oak/stanford/groups/engreitz/Users/tri/tmp/Maggie_CRT_sensitivity_test/logs/%A_%a.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/tri/tmp/Maggie_CRT_sensitivity_test/logs/%A_%a.err
#SBATCH --partition=bigmem,engreitz,owners
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=512G
#SBATCH --array=0-7
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
BASE_DIR="/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/9.cnmf/20260505_k30_300kvicchd_d9_torch"
SENS_ROOT="/oak/stanford/groups/engreitz/Users/tri/tmp/Maggie_CRT_sensitivity_test"
GUIDE_H5AD="$SENS_ROOT/rna_with_guides_batch_sensitivity.h5ad"
CRT_SCRIPT="$SENS_ROOT/CRT_batch_sensitivity.py"
PROGRAMDE_DIR="/oak/stanford/groups/engreitz/Users/msychen/cNMF_benchmarking/programDE"
CONDA_ENV="/oak/stanford/groups/engreitz/Users/msychen/conda/envs/NMF_Benchmarking"
PYTHON="$CONDA_ENV/bin/python"

# -----------------------------------------------------------------------------
# Settings. These can be overridden at submission time with --export.
# -----------------------------------------------------------------------------
CALIBRATION_ONLY="${CALIBRATION_ONLY:-1}"
SEL_THRESH="${SEL_THRESH:-0.4}"
N_PERM="${N_PERM:-1024}"
N_ENSEMBLES="${N_ENSEMBLES:-20}"
N_MAX_GROUPS="${N_MAX_GROUPS:-9}"
N_JOBS="${SLURM_CPUS_PER_TASK:-20}"

VARIANTS=(
    baseline
    inferred2_covariate
    inferred4_covariate
    technical_pcs
    inferred2_plus_pcs
    inferred4_plus_pcs
    inferred2_stratified
    inferred4_stratified
)

IDX="${SLURM_ARRAY_TASK_ID:-0}"
if (( IDX < 0 || IDX >= ${#VARIANTS[@]} )); then
    echo "ERROR: invalid SLURM_ARRAY_TASK_ID=$IDX" >&2
    exit 2
fi
VARIANT="${VARIANTS[$IDX]}"

# Baseline values. Append directly to ALL_COVS in each case branch so the
# baseline never expands an empty array under `set -u`.
CAT_KEY="sample"
ALL_COVS=(
    sublib
    pct_counts_mt
    pct_counts_ribo
    S_score
    G2M_score
    leiden_res0_25
)
LOG_COVS=(
    total_counts
    n_genes_by_counts
    guides_per_cell
    guide_umis_per_cell
)

case "$VARIANT" in
    baseline)
        ;;
    inferred2_covariate)
        ALL_COVS+=(inferred_batch_2)
        ;;
    inferred4_covariate)
        ALL_COVS+=(inferred_batch_4)
        ;;
    technical_pcs)
        ALL_COVS+=(technical_pc1 technical_pc2 technical_pc3 technical_pc4)
        ;;
    inferred2_plus_pcs)
        ALL_COVS+=(
            inferred_batch_2
            technical_pc1 technical_pc2 technical_pc3 technical_pc4
        )
        ;;
    inferred4_plus_pcs)
        ALL_COVS+=(
            inferred_batch_4
            technical_pc1 technical_pc2 technical_pc3 technical_pc4
        )
        ;;
    inferred2_stratified)
        CAT_KEY="sample_x_inferred_batch_2"
        ;;
    inferred4_stratified)
        CAT_KEY="sample_x_inferred_batch_4"
        ;;
    *)
        echo "ERROR: unknown variant: $VARIANT" >&2
        exit 2
        ;;
esac

SAVE_DIR="$SENS_ROOT/$VARIANT/sel_thresh_${SEL_THRESH//./_}"
mkdir -p "$SAVE_DIR" "$SENS_ROOT/logs"

# -----------------------------------------------------------------------------
# Validate paths before starting the expensive analysis.
# -----------------------------------------------------------------------------
for REQUIRED_FILE in "$GUIDE_H5AD" "$CRT_SCRIPT" "$PYTHON"; do
    if [[ ! -f "$REQUIRED_FILE" ]]; then
        echo "ERROR: required file is missing: $REQUIRED_FILE" >&2
        exit 1
    fi
done

if [[ ! -x "$PYTHON" ]]; then
    echo "ERROR: Python is not executable: $PYTHON" >&2
    exit 1
fi

if [[ ! -d "$PROGRAMDE_DIR" ]]; then
    echo "ERROR: programDE directory is missing: $PROGRAMDE_DIR" >&2
    exit 1
fi

CNMF_INPUT="$BASE_DIR/Inference_guideumis/adata/cNMF_30_${SEL_THRESH//./_}.h5mu"
if [[ ! -f "$CNMF_INPUT" ]]; then
    echo "ERROR: expected cNMF input is missing: $CNMF_INPUT" >&2
    exit 1
fi

export PATH="$CONDA_ENV/bin:$PATH"
export PYTHONPATH="$PROGRAMDE_DIR:${PYTHONPATH:-}"
export OMP_NUM_THREADS="$N_JOBS"
export MKL_NUM_THREADS="$N_JOBS"

# Build one command array. Optional arguments are appended conditionally, so no
# empty array is expanded under `set -u`.
CMD=(
    "$PYTHON" "$CRT_SCRIPT"
    --out_dir "$BASE_DIR"
    --run_name "Inference_guideumis"
    --mdata_guide_path "$GUIDE_H5AD"
    --categorical_key "$CAT_KEY"
    --covariates "${ALL_COVS[@]}"
    --log_covariates "${LOG_COVS[@]}"
    --components 30
    --sel_thresh "$SEL_THRESH"
    --guide_annotation_key non-targeting
    --number_permutations "$N_PERM"
    --number_guide 15
    --FDR_method StoreyQ
    --ntc-ensembles "$N_ENSEMBLES"
    --ntc-max-groups "$N_MAX_GROUPS"
    --n_jobs "$N_JOBS"
    --save_dir "$SAVE_DIR"
)

if [[ "$CALIBRATION_ONLY" == "1" ]]; then
    CMD+=(--calibration-only)
elif [[ "$CALIBRATION_ONLY" != "0" ]]; then
    echo "ERROR: CALIBRATION_ONLY must be 0 or 1; got $CALIBRATION_ONLY" >&2
    exit 2
fi

START_TIME=$(date +%s)
echo "Started: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Array task: $IDX"
echo "Variant: $VARIANT"
echo "Categorical key: $CAT_KEY"
echo "Covariates: ${ALL_COVS[*]}"
echo "Log covariates: ${LOG_COVS[*]}"
echo "Calibration only: $CALIBRATION_ONLY"
echo "Output: $SAVE_DIR"
printf 'Command:'
printf ' %q' "${CMD[@]}"
printf '\n'

"${CMD[@]}"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo "Completed: $(date)"
echo "Elapsed: $((ELAPSED / 3600))h $(((ELAPSED % 3600) / 60))m $((ELAPSED % 60))s"
