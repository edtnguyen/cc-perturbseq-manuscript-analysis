#!/bin/bash
set -euo pipefail

PACKAGE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=/dev/null
source "$PACKAGE_DIR/config.env"
MAPPING_CSV="$PACKAGE_DIR/tf_file_set_mapping.csv"

fail() { echo "ERROR: $*" >&2; exit 2; }
[[ -d "$DATA_DIR" ]] || fail "DATA_DIR does not exist: $DATA_DIR"
[[ -f "$MAPPING_CSV" ]] || fail "Mapping CSV missing: $MAPPING_CSV"
[[ -x "$IGVF_PYTHON" ]] || fail "IGVF_PYTHON is not executable: $IGVF_PYTHON"
[[ -f "$IU_REGISTER" ]] || fail "iu_register.py missing: $IU_REGISTER"
[[ "$SEQUENCING_PLATFORM" != *REPLACE_ME* ]] || fail "Set SEQUENCING_PLATFORM in config.env"
[[ -n "${IGVF_API_KEY:-}" ]] || fail "IGVF_API_KEY is not exported in this shell"
[[ -n "${IGVF_SECRET_KEY:-}" ]] || fail "IGVF_SECRET_KEY is not exported in this shell"

mkdir -p "$WORK_DIR/logs"
rm -f "$WORK_DIR/DRYRUN_OK"

# Lightweight metadata/file-name preflight on the login node; no FASTQ contents are read here.
"$IGVF_PYTHON" "$PACKAGE_DIR/prepare_tf_upload.py" \
  --data-dir "$DATA_DIR" \
  --mapping "$MAPPING_CSV" \
  --work-dir "$WORK_DIR" \
  --modality gex

N=$("$IGVF_PYTHON" - "$WORK_DIR/inventory.csv" <<'PY'
import csv, sys
with open(sys.argv[1], newline='') as f:
    print(sum(1 for _ in csv.DictReader(f)))
PY
)
[[ "$N" -gt 0 ]] || fail "Inventory is empty"

export PACKAGE_DIR WORK_DIR IGVF_PYTHON IU_REGISTER SEQUENCING_PLATFORM SEQUENCING_KIT SEQUENCING_RUN
export IGVF_API_KEY IGVF_SECRET_KEY

common=(--parsable --partition="$SLURM_PARTITION" --cpus-per-task=1 --mem="$MEMORY" --export=ALL)
if [[ -n "$SLURM_ACCOUNT" ]]; then
  common+=(--account="$SLURM_ACCOUNT")
fi

validation_job=$(sbatch "${common[@]}" \
  --job-name=tf_fastq_check \
  --time="$VALIDATE_TIME" \
  --array="1-${N}%${MAX_VALIDATE}" \
  --output="$WORK_DIR/logs/validate_%A_%a.out" \
  "$PACKAGE_DIR/validate_fastq_array.sbatch")
validation_job=${validation_job%%;*}

finalize_job=$(sbatch "${common[@]}" \
  --job-name=tf_make_tsv \
  --time="$FINALIZE_TIME" \
  --dependency="afterok:${validation_job}" \
  --output="$WORK_DIR/logs/finalize_%j.out" \
  "$PACKAGE_DIR/finalize_sequence_file.sbatch")
finalize_job=${finalize_job%%;*}

dryrun_job=$(sbatch "${common[@]}" \
  --job-name=tf_igvf_dry \
  --time="$DRYRUN_TIME" \
  --dependency="afterok:${finalize_job}" \
  --output="$WORK_DIR/logs/dryrun_%j.out" \
  "$PACKAGE_DIR/dryrun_sequence_file.sbatch")
dryrun_job=${dryrun_job%%;*}

cat > "$WORK_DIR/pipeline_job_ids.env" <<EOF
VALIDATION_JOB=$validation_job
FINALIZE_JOB=$finalize_job
DRYRUN_JOB=$dryrun_job
EOF

cat <<EOF
Submitted validation and dry-run pipeline.
FASTQ tasks:       $N
Validation array:  $validation_job
Finalize job:      $finalize_job
Dry-run job:       $dryrun_job

Monitor:
  squeue -j $validation_job,$finalize_job,$dryrun_job

After completion, inspect:
  cat "$WORK_DIR/preflight_report.txt"
  cat "$WORK_DIR/finalization_report.txt"
  cat "$WORK_DIR/DRYRUN_OK"
  less "$WORK_DIR/logs/sequence_file_dryrun_${dryrun_job}.log"

Do not run stage 2 unless DRYRUN_OK exists and the dry-run log is clean.
EOF
