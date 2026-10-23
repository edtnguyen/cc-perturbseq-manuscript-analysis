#!/usr/bin/env bash
set -euo pipefail

BASE=${BASE:-/scratch/users/edtng/Hippo_Shiny_adata}
COLLAB=${COLLAB:-edtng}

# Directories need x permission for traversal. Files need r permission.
setfacl -m "u:${COLLAB}:x" "$(dirname "${BASE}")" || true
find "${BASE}" -type d -exec setfacl -m "u:${COLLAB}:rx" {} \;
find "${BASE}" -type f -exec setfacl -m "u:${COLLAB}:r" {} \;
find "${BASE}" -type d -exec setfacl -d -m "u:${COLLAB}:rx" {} \;

echo "ACLs for ${BASE}:"
getfacl "${BASE}" | sed -n '1,60p'
