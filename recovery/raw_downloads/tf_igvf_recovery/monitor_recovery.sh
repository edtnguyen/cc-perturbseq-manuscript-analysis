#!/bin/bash
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/config.env"
job=$(cat "$RECOVERY_DIR/retry_job_id.txt" 2>/dev/null || true)
TOTAL=289
OK=$(find "$RECOVERY_DIR/upload_done" -maxdepth 1 -type f -name '*.json' 2>/dev/null | wc -l)
PARTIAL=$(grep -l '"status": "FOUND"' "$RECOVERY_DIR"/portal_after_attempt/*.json 2>/dev/null | wc -l)
NOTFOUND=$(grep -l '"status": "NOT_FOUND"' "$RECOVERY_DIR"/portal_after_attempt/*.json 2>/dev/null | wc -l)
echo "Total=$TOTAL Success=$OK Portal_found_after_failure=$PARTIAL Portal_not_found_after_failure=$NOTFOUND Remaining=$((TOTAL-OK-PARTIAL-NOTFOUND))"
if [[ -n "$job" ]]; then squeue -j "$job" -o '%.22i %.10P %.10T %.10M %.40R' || true; fi
