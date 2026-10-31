#!/bin/bash
set -euo pipefail
if [[ "${1:-}" != "--confirm" ]]; then
  echo "Refusing real PATCH without --confirm" >&2
  exit 2
fi
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/config.env"
: "${IGVF_API_KEY:?Export IGVF_API_KEY first}"
: "${IGVF_SECRET_KEY:?Export IGVF_SECRET_KEY first}"
mkdir -p "$RECOVERY_DIR/logs"
LOG="$RECOVERY_DIR/logs/patch_wrong_grna_submit_$(date +%Y%m%d_%H%M%S).log"
"$IGVF_PYTHON" "$IU_REGISTER" -m prod -p sequence_file --patch \
  --tries 5 --delay 10 --backoff 2 \
  -i "$RECOVERY_DIR/patch_55_wrong_grna_file_sets.tsv" 2>&1 | tee "$LOG"
echo "PATCH log: $LOG"
