#!/usr/bin/env bash
set -euo pipefail

EXPECTED_SHA256="d2fe4834b1efd3c00e01062b56b907a6178bfbc48a40beb5ea7b347538bfc824"
DEFAULT_SOURCE="/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/data/cNMF_50_0_2_032526_6co.xlsx"
SOURCE="$DEFAULT_SOURCE"
PUSH=false

for arg in "$@"; do
  case "$arg" in
    --push) PUSH=true ;;
    --help|-h)
      cat <<'EOF'
Usage: bash scripts/install_figure2_cnmf_workbook.sh [SOURCE_XLSX] [--push]

Copies the verified Figure 2 cNMF workbook into data/figure_02/.
With --push, the script also commits and pushes the workbook on the current branch.
EOF
      exit 0
      ;;
    *) SOURCE="$arg" ;;
  esac
done

REPO_ROOT="$(git rev-parse --show-toplevel)"
DEST_DIR="$REPO_ROOT/data/figure_02"
DEST="$DEST_DIR/cNMF_50_0_2_032526_6co.xlsx"

if [[ ! -f "$SOURCE" ]]; then
  echo "ERROR: source workbook not found: $SOURCE" >&2
  echo "Pass an alternate source path as the first argument." >&2
  exit 1
fi

ACTUAL_SHA256="$(sha256sum "$SOURCE" | awk '{print $1}')"
if [[ "$ACTUAL_SHA256" != "$EXPECTED_SHA256" ]]; then
  echo "ERROR: source workbook SHA-256 does not match the supplied workbook." >&2
  echo "Expected: $EXPECTED_SHA256" >&2
  echo "Observed: $ACTUAL_SHA256" >&2
  exit 1
fi

mkdir -p "$DEST_DIR"
cp -f "$SOURCE" "$DEST"
test "$(sha256sum "$DEST" | awk '{print $1}')" = "$EXPECTED_SHA256"
echo "Installed verified workbook: $DEST"

if [[ "$PUSH" == true ]]; then
  branch="$(git branch --show-current)"
  if [[ -z "$branch" ]]; then
    echo "ERROR: detached HEAD; check out the intended branch first." >&2
    exit 1
  fi
  git add "$DEST"
  if git diff --cached --quiet; then
    echo "Workbook is already committed; nothing to commit."
  else
    git commit -m "Add Figure 2 cNMF source workbook"
  fi
  git push origin "$branch"
else
  echo "Next: git add '$DEST' && git commit -m 'Add Figure 2 cNMF source workbook' && git push"
fi
