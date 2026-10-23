#!/usr/bin/env bash
set -euo pipefail

TARGET=${1:-/scratch/users/edtng/Hippo_Shiny_adata}
SRC_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

mkdir -p "${TARGET}"/{app,scripts,data,raw,logs}
cp -f "${SRC_DIR}/app/app.R" "${TARGET}/app/app.R"
cp -f "${SRC_DIR}/run_app.R" "${TARGET}/run_app.R"
cp -f "${SRC_DIR}/scripts/"*.R "${TARGET}/scripts/"
cp -f "${SRC_DIR}/scripts/"*.sh "${TARGET}/scripts/"
cp -f "${SRC_DIR}/scripts/export_umap_shiny_h5.py" "${TARGET}/scripts/"
cp -f "${SRC_DIR}/scripts/check_python_packages.py" "${TARGET}/scripts/"
cp -f "${SRC_DIR}/markers_for_shiny.example.txt" "${TARGET}/markers_for_shiny.example.txt"
chmod +x "${TARGET}/run_app.R" "${TARGET}/scripts/"*.R "${TARGET}/scripts/"*.sh "${TARGET}/scripts/"*.py

echo "Installed app files under ${TARGET}"
find "${TARGET}" -maxdepth 2 -type f | sort
