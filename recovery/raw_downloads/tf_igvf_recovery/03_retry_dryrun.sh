#!/bin/bash
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/config.env"
: "${IGVF_API_KEY:?Export IGVF_API_KEY first}"
: "${IGVF_SECRET_KEY:?Export IGVF_SECRET_KEY first}"
mkdir -p "$RECOVERY_DIR/logs"
LOG="$RECOVERY_DIR/logs/retry_not_found_dryrun_$(date +%Y%m%d_%H%M%S).log"
"$IGVF_PYTHON" "$IU_REGISTER" -m prod -p sequence_file -d \
  -i "$RECOVERY_DIR/retry_289_not_found_sequence_files.tsv" 2>&1 | tee "$LOG"
echo "Dry-run log: $LOG"
