#!/bin/bash

#SBATCH --job-name=CRT_batch_sens
#SBATCH --output=/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/9.cnmf/20260505_k30_300kvicchd_d9_torch/CRT_batch_sensitivity/logs/%A_%a.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/9.cnmf/20260505_k30_300kvicchd_d9_torch/CRT_batch_sensitivity/logs/%A_%a.err
#SBATCH --partition=bigmem,engreitz,owners
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=512G
#SBATCH --array=0-7
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=msychen@stanford.edu

set -euo pipefail

BASE_DIR="/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/9.cnmf/20260505_k30_300kvicchd_d9_torch"
SENS_ROOT="$BASE_DIR/CRT_batch_sensitivity"
GUIDE_H5AD="$SENS_ROOT/rna_with_guides_batch_sensitivity.h5ad"

# Copy CRT_batch_sensitivity.py to this location or edit the variable.
CRT_SCRIPT="/oak/stanford/groups/engreitz/Users/msychen/cNMF_benchmarking/cNMF_benchmarking_pipeline/Perturbation_association_calibration/Slurm_version/CRT/CRT_batch_sensitivity.py"
PROGRAMDE_DIR="/oak/stanford/groups/engreitz/Users/msychen/cNMF_benchmarking/programDE"

# First-pass screen: one cNMF solution and only NTC calibration tests.
# Set CALIBRATION_ONLY=0 before submission to run full discovery for selected variants.
CALIBRATION_ONLY="${CALIBRATION_ONLY:-1}"
SEL_THRESH="${SEL_THRESH:-0.4}"
N_PERM="${N_PERM:-1024}"
N_ENSEMBLES="${N_ENSEMBLES:-20}"
N_MAX_GROUPS="${N_MAX_GROUPS:-9}"

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
VARIANT="${VARIANTS[$IDX]}"
CAT_KEY="sample"
EXTRA_COVS=()

case "$VARIANT" in
    baseline)
        ;;
    inferred2_covariate)
        EXTRA_COVS+=(inferred_batch_2)
        ;;
    inferred4_covariate)
        EXTRA_COVS+=(inferred_batch_4)
        ;;
    technical_pcs)
        EXTRA_COVS+=(technical_pc1 technical_pc2 technical_pc3 technical_pc4)
        ;;
    inferred2_plus_pcs)
        EXTRA_COVS+=(inferred_batch_2 technical_pc1 technical_pc2 technical_pc3 technical_pc4)
        ;;
    inferred4_plus_pcs)
        EXTRA_COVS+=(inferred_batch_4 technical_pc1 technical_pc2 technical_pc3 technical_pc4)
        ;;
    inferred2_stratified)
        CAT_KEY="sample_x_inferred_batch_2"
        ;;
    inferred4_stratified)
        CAT_KEY="sample_x_inferred_batch_4"
        ;;
    *)
        echo "Unknown variant: $VARIANT" >&2
        exit 2
        ;;
esac

SAVE_DIR="$SENS_ROOT/$VARIANT/sel_thresh_${SEL_THRESH//./_}"
mkdir -p "$SAVE_DIR" "$SENS_ROOT/logs"

if [[ ! -f "$GUIDE_H5AD" ]]; then
    echo "ERROR: prepared guide file missing: $GUIDE_H5AD" >&2
    echo "Run prepare_batch_sensitivity_covariates.py first." >&2
    exit 1
fi
if [[ ! -f "$CRT_SCRIPT" ]]; then
    echo "ERROR: modified CRT wrapper missing: $CRT_SCRIPT" >&2
    exit 1
fi

BASE_COVS=(
    sublib
    pct_counts_mt
    pct_counts_ribo
    S_score
    G2M_score
    leiden_res0_25
)
ALL_COVS=("${BASE_COVS[@]}" "${EXTRA_COVS[@]}")
LOG_COVS=(total_counts n_genes_by_counts guides_per_cell guide_umis_per_cell)

CALIBRATION_FLAG=()
if [[ "$CALIBRATION_ONLY" == "1" ]]; then
    CALIBRATION_FLAG+=(--calibration-only)
fi

START_TIME=$(date +%s)
echo "Started: $(date)"
echo "Variant: $VARIANT"
echo "Categorical key: $CAT_KEY"
echo "Additional covariates: ${EXTRA_COVS[*]:-(none)}"
echo "Output: $SAVE_DIR"

eval "$(conda shell.bash hook)"
conda activate /oak/stanford/groups/engreitz/Users/msychen/conda/envs/NMF_Benchmarking
export PYTHONPATH="$PROGRAMDE_DIR:$PYTHONPATH"

python3 "$CRT_SCRIPT" \
    --out_dir "$BASE_DIR" \
    --run_name "Inference_guideumis" \
    --mdata_guide_path "$GUIDE_H5AD" \
    --categorical_key "$CAT_KEY" \
    --covariates "${ALL_COVS[@]}" \
    --log_covariates "${LOG_COVS[@]}" \
    --components 30 \
    --sel_thresh "$SEL_THRESH" \
    --guide_annotation_key non-targeting \
    --number_permutations "$N_PERM" \
    --number_guide 15 \
    --FDR_method StoreyQ \
    --ntc-ensembles "$N_ENSEMBLES" \
    --ntc-max-groups "$N_MAX_GROUPS" \
    --n_jobs "$SLURM_CPUS_PER_TASK" \
    --save_dir "$SAVE_DIR" \
    "${CALIBRATION_FLAG[@]}"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo "Completed: $(date)"
echo "Elapsed: $((ELAPSED/3600))h $(((ELAPSED%3600)/60))m $((ELAPSED%60))s"
