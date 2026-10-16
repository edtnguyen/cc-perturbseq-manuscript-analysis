#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Copy every Jupyter notebook from one or more Sherlock source trees into the
repository and stage the complete archive with Git LFS.

No notebook is excluded because of file size. Size determines neither selection
nor copying. Git LFS is used for the full-fidelity Sherlock archive so notebooks
larger than ordinary GitHub blob limits can be pushed.

Usage:
  scripts/sync_all_sherlock_notebooks.sh [options] --source ROOT [--source ROOT ...]

Options:
  --repo PATH          Repository clone. Default: current Git repository.
  --dest PATH          Repository-relative archive directory.
                       Default: recovery/sherlock_notebooks
  --source ROOT        Source directory to scan recursively. Repeatable.
  --source-list FILE   Text file containing one source directory per line.
                       Blank lines and lines beginning with # are ignored.
  --commit             Commit staged changes after copying and verification.
  --push               Commit, upload all LFS objects, and push the current branch.
  --dry-run            Inventory matching notebooks without copying or staging.
  -h, --help           Show this help text.

Examples:
  scripts/sync_all_sherlock_notebooks.sh \
    --source /oak/stanford/groups/engreitz/Users/tri/CodeBase \
    --source /oak/stanford/groups/engreitz/Users/ymo/ipsc_ec \
    --commit

  scripts/sync_all_sherlock_notebooks.sh \
    --source-list config/sherlock_notebook_roots.local.txt \
    --push
EOF
}

repo_root=""
dest_rel="recovery/sherlock_notebooks"
do_commit=0
do_push=0
dry_run=0
declare -a source_roots=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo)
      repo_root=${2:?"--repo requires a path"}
      shift 2
      ;;
    --dest)
      dest_rel=${2:?"--dest requires a repository-relative path"}
      shift 2
      ;;
    --source)
      source_roots+=("${2:?"--source requires a directory"}")
      shift 2
      ;;
    --source-list)
      list_file=${2:?"--source-list requires a file"}
      [[ -f "$list_file" ]] || { echo "Source-list file not found: $list_file" >&2; exit 2; }
      while IFS= read -r line || [[ -n "$line" ]]; do
        line=${line%$'\r'}
        [[ -z "${line//[[:space:]]/}" ]] && continue
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        source_roots+=("$line")
      done < "$list_file"
      shift 2
      ;;
    --commit)
      do_commit=1
      shift
      ;;
    --push)
      do_commit=1
      do_push=1
      shift
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

[[ ${#source_roots[@]} -gt 0 ]] || {
  echo "At least one --source or --source-list entry is required." >&2
  usage >&2
  exit 2
}

if [[ -z "$repo_root" ]]; then
  repo_root=$(git rev-parse --show-toplevel 2>/dev/null || true)
fi
[[ -n "$repo_root" ]] || { echo "Run inside a Git clone or pass --repo." >&2; exit 2; }
repo_root=$(realpath -e -- "$repo_root")
cd "$repo_root"
git rev-parse --is-inside-work-tree >/dev/null

if [[ "$dest_rel" = /* || "$dest_rel" == ".." || "$dest_rel" == ../* || "$dest_rel" == */../* ]]; then
  echo "--dest must be a safe repository-relative path: $dest_rel" >&2
  exit 2
fi

if ! git lfs version >/dev/null 2>&1; then
  cat >&2 <<'EOF'
Git LFS is required but is not available in this Sherlock environment.
Install it in a conda/mamba environment, for example:

  mamba install -c conda-forge git-lfs
  git lfs install

Then rerun this script. No notebook should be omitted as a workaround.
EOF
  exit 3
fi

git lfs install --local >/dev/null

dest_abs="$repo_root/$dest_rel"
manifest_rel="$dest_rel/MANIFEST.tsv"
manifest_abs="$repo_root/$manifest_rel"
mkdir -p "$dest_abs"

manifest_tmp=$(mktemp)
trap 'rm -f "$manifest_tmp"' EXIT
printf 'source_absolute_path\trepository_path\tsize_bytes\tsha256\tstorage\n' > "$manifest_tmp"

declare -A seen_paths=()
declare -a very_large_files=()
notebook_count=0
total_bytes=0
copied_count=0
unchanged_count=0

for source_root in "${source_roots[@]}"; do
  [[ -d "$source_root" ]] || { echo "Source directory not found: $source_root" >&2; exit 2; }
  source_root=$(realpath -e -- "$source_root")
  echo "Scanning: $source_root" >&2

  while IFS= read -r -d '' source_file; do
    source_file=$(realpath -e -- "$source_file")
    if [[ -n "${seen_paths[$source_file]:-}" ]]; then
      continue
    fi
    seen_paths[$source_file]=1

    # Preserve the absolute Sherlock path below the archive root. This avoids
    # basename collisions and makes provenance recoverable without a side map.
    repository_path="$dest_rel/${source_file#/}"
    destination="$repo_root/$repository_path"
    size_bytes=$(stat -c '%s' -- "$source_file")
    sha256=$(sha256sum -- "$source_file" | awk '{print $1}')

    notebook_count=$((notebook_count + 1))
    total_bytes=$((total_bytes + size_bytes))
    if (( size_bytes > 2147483648 )); then
      very_large_files+=("$source_file")
    fi

    printf '%s\t%s\t%s\t%s\tgit-lfs\n' \
      "$source_file" "$repository_path" "$size_bytes" "$sha256" >> "$manifest_tmp"

    if (( dry_run )); then
      continue
    fi

    mkdir -p -- "$(dirname -- "$destination")"
    if [[ -f "$destination" ]] && cmp -s -- "$source_file" "$destination"; then
      unchanged_count=$((unchanged_count + 1))
    else
      cp -p -- "$source_file" "$destination"
      copied_count=$((copied_count + 1))
    fi
  done < <(
    find "$source_root" \
      \( -type d \( -name .git -o -name .ipynb_checkpoints \) -prune \) -o \
      \( -type f -name '*.ipynb' -print0 \) | sort -z
  )
done

(( notebook_count > 0 )) || { echo "No .ipynb files were found in the requested source trees." >&2; exit 4; }

if command -v numfmt >/dev/null 2>&1; then
  total_human=$(numfmt --to=iec-i --suffix=B "$total_bytes")
else
  total_human="$total_bytes bytes"
fi

echo "Discovered notebooks: $notebook_count" >&2
echo "Total notebook bytes: $total_human" >&2

if (( ${#very_large_files[@]} > 0 )); then
  cat >&2 <<EOF
Warning: ${#very_large_files[@]} notebook(s) exceed 2 GiB. They are still included,
but the Git LFS per-file limit depends on the GitHub plan. The push will fail
rather than silently omit a file if the account limit is exceeded:
EOF
  printf '  %s\n' "${very_large_files[@]}" >&2
fi

if (( dry_run )); then
  echo "Dry run complete. Nothing was copied, staged, committed, or pushed." >&2
  exit 0
fi

mv -- "$manifest_tmp" "$manifest_abs"
trap - EXIT

# The committed .gitattributes rules cover every notebook below the archive
# directory. Renormalization also converts a previously copied archive file to
# an LFS pointer instead of excluding it because it is large.
git add -- .gitattributes
git add --renormalize -- "$dest_rel" 2>/dev/null || true
git add -- "$dest_rel"

attribute_failures=0
while IFS= read -r -d '' archived_file; do
  rel=${archived_file#"$repo_root/"}
  attr=$(git check-attr filter -- "$rel" | awk -F': ' '{print $3}')
  if [[ "$attr" != "lfs" ]]; then
    echo "LFS attribute was not applied to: $rel" >&2
    attribute_failures=$((attribute_failures + 1))
  fi
done < <(find "$dest_abs" -type f -name '*.ipynb' -print0)

(( attribute_failures == 0 )) || {
  echo "Refusing to commit because some archived notebooks are not tracked by Git LFS." >&2
  exit 5
}

lfs_count=$(git lfs ls-files --cached | awk -v prefix="$dest_rel/" 'index($0, prefix) {n++} END {print n+0}')
if (( lfs_count < notebook_count )); then
  echo "Warning: $notebook_count source notebooks were inventoried, but $lfs_count archive paths are currently listed by Git LFS." >&2
  echo "This can occur when overlapping source roots map to existing identical paths; inspect $manifest_rel before pushing." >&2
fi

echo "Copied or updated: $copied_count" >&2
echo "Already identical: $unchanged_count" >&2
echo "Manifest: $manifest_rel" >&2
echo "Git LFS archive entries: $lfs_count" >&2

if (( do_commit )); then
  if git diff --cached --quiet; then
    echo "No staged changes to commit." >&2
  else
    git commit -m "Archive all Sherlock notebooks with Git LFS"
  fi
fi

if (( do_push )); then
  branch=$(git branch --show-current)
  [[ -n "$branch" ]] || { echo "Cannot push from detached HEAD." >&2; exit 6; }
  echo "Uploading Git LFS objects for branch: $branch" >&2
  git lfs push origin "$branch"
  git push origin "$branch"
fi

echo "Notebook sync complete. No file-size selection cutoff was applied." >&2
