#!/usr/bin/env bash
set -euo pipefail

REPO="${1:-/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction}"
KEY="/home/users/edtng/.ssh/id_ed25519_github_edtnguyen"

NOTEBOOKS=(
  "notebooks/main/Figure_01.ipynb"
  "notebooks/main/Figure_02.ipynb"
  "notebooks/main/Figure_03.ipynb"
  "notebooks/main/Figure_04.ipynb"
  "notebooks/supplementary/Figure_S02.ipynb"
  "notebooks/supplementary/Figure_S04.ipynb"
  "notebooks/supplementary/Figure_S06.ipynb"
)

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

[ -d "$REPO/.git" ] || fail "Repository not found: $REPO"
cd "$REPO"

[ -n "${SSH_AUTH_SOCK:-}" ] || fail "SSH_AUTH_SOCK is not set. Start/source the SSH agent first."
[ -S "$SSH_AUTH_SOCK" ] || fail "SSH agent socket is unavailable: $SSH_AUTH_SOCK"
ssh-add -l >/dev/null 2>&1 || fail "No SSH identity is loaded. Run ssh-add $KEY first."

# Use only the loaded agent. BatchMode prevents another passphrase prompt.
export GIT_SSH_COMMAND="ssh -o BatchMode=yes -o IdentitiesOnly=no -o IdentityAgent=$SSH_AUTH_SOCK"

# Verify GitHub authentication without blocking for input. GitHub returns exit 1
# after a successful identity greeting because it does not provide shell access.
AUTH_TEXT="$(ssh -o BatchMode=yes -o IdentitiesOnly=no -o IdentityAgent="$SSH_AUTH_SOCK" -T git@github.com 2>&1 || true)"
echo "$AUTH_TEXT"
echo "$AUTH_TEXT" | grep -q "successfully authenticated" || fail "The loaded SSH agent did not authenticate to GitHub."

echo
echo "Checking saved notebook outputs..."
python - "${NOTEBOOKS[@]}" <<'PY'
import json
import sys
from pathlib import Path

visual = {"image/png", "image/jpeg", "image/svg+xml", "application/pdf"}
problems = []

for name in sys.argv[1:]:
    path = Path(name)
    if not path.is_file():
        problems.append(f"{name}: missing")
        continue
    size = path.stat().st_size
    if size >= 95 * 1024 * 1024:
        problems.append(f"{name}: {size / 1024**2:.1f} MiB, too close to GitHub's 100 MiB limit")
    try:
        nb = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        problems.append(f"{name}: invalid JSON ({exc})")
        continue
    figures = errors = executed = 0
    for cell in nb.get("cells", []):
        if cell.get("cell_type") == "code" and cell.get("execution_count") is not None:
            executed += 1
        for output in cell.get("outputs", []):
            if output.get("output_type") == "error":
                errors += 1
            data = output.get("data", {})
            if any(mime in data for mime in visual):
                figures += 1
    print(f"{name}: {size/1024**2:.2f} MiB; executed={executed}; figures={figures}; errors={errors}")
    if figures == 0:
        problems.append(f"{name}: no embedded figure output")
    if errors:
        problems.append(f"{name}: contains {errors} saved error output(s)")

if problems:
    print("\nNot ready to publish:")
    for item in problems:
        print(" -", item)
    raise SystemExit(1)
PY

STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="$HOME/executed_figure_notebooks_${STAMP}.tar.gz"
WORKTREE_ROOT="${SCRATCH:-/tmp}"
WORKTREE="$WORKTREE_ROOT/publish_executed_notebooks_${STAMP}"

tar -czf "$BACKUP" "${NOTEBOOKS[@]}"
echo "Backup: $BACKUP"

echo
echo "Fetching current GitHub main through the loaded agent..."
git fetch origin +refs/heads/main:refs/remotes/origin/main
MAIN_BEFORE="$(git rev-parse origin/main)"
echo "origin/main before publication: $MAIN_BEFORE"

[ ! -e "$WORKTREE" ] || fail "Temporary worktree already exists: $WORKTREE"
git worktree add --detach "$WORKTREE" "$MAIN_BEFORE"

cleanup() {
  cd "$REPO"
  git worktree remove --force "$WORKTREE" >/dev/null 2>&1 || true
}
trap cleanup EXIT

for relpath in "${NOTEBOOKS[@]}"; do
  mkdir -p "$WORKTREE/$(dirname "$relpath")"
  cp -p "$REPO/$relpath" "$WORKTREE/$relpath"
  blob="$(git -C "$WORKTREE" hash-object -w --no-filters "$WORKTREE/$relpath")"
  git -C "$WORKTREE" update-index --add --cacheinfo 100644 "$blob" "$relpath"
done

echo
echo "Staged changes:"
git -C "$WORKTREE" diff --cached --name-status

python - "$WORKTREE" "${NOTEBOOKS[@]}" <<'PY'
import json
import subprocess
import sys
from pathlib import Path

repo = Path(sys.argv[1])
allowed = set(sys.argv[2:])
actual = set(subprocess.check_output(
    ["git", "-C", str(repo), "diff", "--cached", "--name-only"],
    text=True,
).splitlines())

unexpected = actual - allowed
if unexpected:
    raise SystemExit("Unexpected staged paths: " + ", ".join(sorted(unexpected)))
if not actual:
    print("No notebook changes relative to GitHub main.")
    raise SystemExit(0)

visual = {"image/png", "image/jpeg", "image/svg+xml", "application/pdf"}
for relpath in sorted(actual):
    raw = subprocess.check_output(["git", "-C", str(repo), "show", f":{relpath}"])
    nb = json.loads(raw.decode("utf-8"))
    figures = errors = 0
    for cell in nb.get("cells", []):
        for output in cell.get("outputs", []):
            if output.get("output_type") == "error":
                errors += 1
            if any(m in output.get("data", {}) for m in visual):
                figures += 1
    print(f"staged {relpath}: figures={figures}; errors={errors}")
    if figures == 0 or errors:
        raise SystemExit(f"Invalid staged notebook: {relpath}")
PY

if git -C "$WORKTREE" diff --cached --quiet; then
  echo "Nothing to commit. GitHub main already has these exact notebook bytes."
  FINAL_SHA="$MAIN_BEFORE"
else
  git -C "$WORKTREE" diff --cached --check
  git -C "$WORKTREE" \
    -c user.name="Tri Nguyen" \
    -c user.email="72826864+edtnguyen@users.noreply.github.com" \
    commit \
    -m "Publish executed figure notebooks with embedded outputs" \
    -m "Update Figures 1-4, S2, S4, and S6 with saved Sherlock output cells so generated panels render directly on GitHub."

  git fetch origin +refs/heads/main:refs/remotes/origin/main
  MAIN_NOW="$(git rev-parse origin/main)"
  [ "$MAIN_NOW" = "$MAIN_BEFORE" ] || fail "GitHub main changed during preparation. Rerun this script."

  git -C "$WORKTREE" push origin HEAD:main
  git fetch origin +refs/heads/main:refs/remotes/origin/main
  FINAL_SHA="$(git rev-parse origin/main)"
  LOCAL_SHA="$(git -C "$WORKTREE" rev-parse HEAD)"
  [ "$FINAL_SHA" = "$LOCAL_SHA" ] || fail "Remote main does not match the published commit."
fi

echo
echo "SUCCESS"
echo "Current GitHub main commit: $FINAL_SHA"
echo "Backup: $BACKUP"
