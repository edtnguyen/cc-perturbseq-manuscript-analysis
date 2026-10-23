# Uploading all Sherlock notebooks

The Sherlock import must not use a `-size -10M` filter or any equivalent size-based selection rule. Every `.ipynb` file under the explicitly selected source roots is copied and staged.

## Storage policy

- Curated manuscript drivers under `notebooks/` remain in ordinary Git so GitHub can render and diff them.
- Full-fidelity historical/source notebooks imported from Sherlock are stored under `recovery/sherlock_notebooks/`.
- The archive preserves the original absolute Sherlock path below that directory, preventing filename collisions and retaining provenance.
- All archived `.ipynb` files are tracked with Git LFS. Their outputs are not stripped and duplicate-content notebooks at distinct source paths are not discarded.
- File size affects only whether GitHub and the account's Git LFS plan can accept the object. The sync script never silently omits a notebook because it is large.

GitHub blocks ordinary Git blobs larger than 100 MiB, whereas Git LFS stores a small pointer in Git and uploads the underlying object separately. Git LFS per-file limits depend on the GitHub plan. See GitHub's current documentation:

- https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github
- https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-git-large-file-storage

## One-time setup on Sherlock

```bash
cd /path/to/cc-perturbseq-manuscript-analysis
git fetch origin
git switch curation/figure-reproduction-v1
git pull --ff-only origin curation/figure-reproduction-v1
```

Install Git LFS in the environment used for the push. One option is:

```bash
mamba install -c conda-forge git-lfs
git lfs install
```

Confirm that it is available:

```bash
git lfs version
```

## Define the source roots

Create an untracked local file containing one source directory per line:

```bash
cat > config/sherlock_notebook_roots.local.txt <<'EOF'
/oak/stanford/groups/engreitz/Users/tri/CodeBase
/oak/stanford/groups/engreitz/Users/ymo/ipsc_ec
/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec
# Add every other explicitly in-scope source tree here.
EOF
```

The file is ignored by Git because it contains machine-specific paths. Source roots may overlap; the script de-duplicates only the same absolute file encountered through overlapping roots. It does not drop distinct files with identical contents.

## Inventory before copying

```bash
bash scripts/sync_all_sherlock_notebooks.sh \
  --source-list config/sherlock_notebook_roots.local.txt \
  --dry-run
```

The dry run reports the number of notebooks and their total size. There is no 10 MB cutoff.

## Copy, verify, commit, and push everything

```bash
bash scripts/sync_all_sherlock_notebooks.sh \
  --source-list config/sherlock_notebook_roots.local.txt \
  --push
```

The script performs the following steps:

1. Recursively finds every `.ipynb` file below the listed roots, excluding only `.git` and `.ipynb_checkpoints` directories.
2. Copies each notebook byte-for-byte into `recovery/sherlock_notebooks/<absolute Sherlock path>`.
3. Writes `recovery/sherlock_notebooks/MANIFEST.tsv` with source path, repository path, size and SHA-256 hash.
4. Applies the committed Git LFS attributes and stages all copied notebooks.
5. Verifies that every archived notebook receives the `filter=lfs` attribute.
6. Commits the archive, explicitly uploads the branch's LFS objects, and pushes the current branch.

The operation is restartable. Re-running it leaves byte-identical files unchanged, and Git LFS skips objects that were already uploaded.

## Verification after the push

On Sherlock:

```bash
git status --short
git lfs ls-files | wc -l
git lfs fsck
```

Check that the manifest count matches the source inventory:

```bash
python - <<'PY'
from pathlib import Path

manifest = Path('recovery/sherlock_notebooks/MANIFEST.tsv')
rows = max(sum(1 for _ in manifest.open()) - 1, 0)
archive = sum(1 for _ in Path('recovery/sherlock_notebooks').rglob('*.ipynb'))
print({'manifest_notebooks': rows, 'archived_notebooks': archive})
assert rows == archive
PY
```

After a fresh clone, retrieve the full notebook objects with:

```bash
git lfs pull
```

Without Git LFS, a clone contains pointer files rather than the original large notebooks.

## Failure behavior

The script fails rather than skipping a notebook when:

- Git LFS is unavailable;
- a configured source root does not exist;
- an archived notebook is not covered by the LFS attributes;
- the GitHub account's storage, bandwidth or per-file limit rejects an object; or
- the network upload fails.

After fixing the cause, rerun the same command. For an interrupted or incomplete LFS upload, the explicit recovery command is:

```bash
git lfs push --all origin curation/figure-reproduction-v1
git push origin curation/figure-reproduction-v1
```
