#!/bin/bash
set -euo pipefail

PACKAGE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=/dev/null
source "$PACKAGE_DIR/config.env"

fail() { echo "ERROR: $*" >&2; exit 2; }
[[ "${1:-}" == "--confirm" ]] || fail "Run as: bash 02_submit_upload_array.sh --confirm"
[[ -s "$WORK_DIR/DRYRUN_OK" ]] || fail "DRYRUN_OK is missing; finish and review stage 1 first"
[[ -s "$WORK_DIR/sequence_file.tsv" ]] || fail "sequence_file.tsv is missing"
[[ -s "$WORK_DIR/upload_manifest.csv" ]] || fail "upload_manifest.csv is missing"
[[ -x "$IGVF_PYTHON" ]] || fail "IGVF_PYTHON is not executable: $IGVF_PYTHON"
[[ -f "$IU_REGISTER" ]] || fail "iu_register.py missing: $IU_REGISTER"
[[ -n "${IGVF_API_KEY:-}" ]] || fail "IGVF_API_KEY is not exported in this shell"
[[ -n "${IGVF_SECRET_KEY:-}" ]] || fail "IGVF_SECRET_KEY is not exported in this shell"

N=$("$IGVF_PYTHON" - "$WORK_DIR/upload_manifest.csv" <<'PY'
import csv, sys
with open(sys.argv[1], newline='') as f:
    print(sum(1 for _ in csv.DictReader(f)))
PY
)
[[ "$N" -gt 0 ]] || fail "Upload manifest is empty"

export PACKAGE_DIR WORK_DIR IGVF_PYTHON IU_REGISTER IGVF_API_KEY IGVF_SECRET_KEY
mkdir -p "$WORK_DIR/logs" "$WORK_DIR/upload_logs" "$WORK_DIR/upload_done"

common=(--parsable --partition="$SLURM_PARTITION" --cpus-per-task=1 --mem="$MEMORY" --export=ALL)
if [[ -n "$SLURM_ACCOUNT" ]]; then
  common+=(--account="$SLURM_ACCOUNT")
fi

upload_job=$(sbatch "${common[@]}" \
  --job-name=tf_igvf_upload \
  --time="$UPLOAD_TIME" \
  --array="1-${N}%${MAX_UPLOAD}" \
  --output="$WORK_DIR/logs/upload_%A_%a.out" \
  "$PACKAGE_DIR/upload_sequence_file_array.sbatch")
upload_job=${upload_job%%;*}

collect_job=$(sbatch "${common[@]}" \
  --job-name=tf_igvf_collect \
  --time="00:20:00" \
  --dependency="afterany:${upload_job}" \
  --output="$WORK_DIR/logs/collect_%j.out" \
  "$PACKAGE_DIR/collect_upload_results.sbatch")
collect_job=${collect_job%%;*}

cat >> "$WORK_DIR/pipeline_job_ids.env" <<EOF
UPLOAD_JOB=$upload_job
COLLECT_JOB=$collect_job
EOF

cat <<EOF
Submitted IGVF upload array.
Upload tasks:  $N
Upload job:    $upload_job (maximum $MAX_UPLOAD simultaneous tasks)
Collector job: $collect_job

Monitor:
  squeue -j $upload_job,$collect_job
  sacct -j $upload_job --format=JobID,State,ExitCode,Elapsed

After completion:
  cat "$WORK_DIR/upload_summary.txt"
  column -s, -t < "$WORK_DIR/sequence_file_upload_results.csv" | less -S
EOF
