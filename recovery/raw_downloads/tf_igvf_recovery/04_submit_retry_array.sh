#!/bin/bash
set -euo pipefail
if [[ "${1:-}" != "--confirm" ]]; then
  echo "Refusing real upload without --confirm" >&2
  exit 2
fi
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/config.env"
: "${IGVF_API_KEY:?Export IGVF_API_KEY first}"
: "${IGVF_SECRET_KEY:?Export IGVF_SECRET_KEY first}"
N=$(find "$RECOVERY_DIR/retry_rows" -maxdepth 1 -name '*.tsv' -type f | wc -l)
[[ "$N" -eq 289 ]] || { echo "Expected 289 row shards, found $N" >&2; exit 2; }
mkdir -p "$RECOVERY_DIR/slurm" "$RECOVERY_DIR/upload_logs" "$RECOVERY_DIR/upload_done" "$RECOVERY_DIR/portal_after_attempt"
job=$(sbatch --parsable --partition="$SLURM_PARTITION" --array="1-${N}%${MAX_UPLOAD}" \
  --time="$UPLOAD_TIME" --mem="$MEMORY" --export=ALL,RECOVERY_DIR="$RECOVERY_DIR",IGVF_PYTHON="$IGVF_PYTHON",IU_REGISTER="$IU_REGISTER" \
  "$RECOVERY_DIR/retry_upload_array.sbatch")
echo "$job" | tee "$RECOVERY_DIR/retry_job_id.txt"
echo "Submitted retry array job $job with concurrency $MAX_UPLOAD"
