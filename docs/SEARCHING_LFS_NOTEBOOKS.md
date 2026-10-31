# Searching the Sherlock notebook archive

The full-fidelity notebooks under `recovery/sherlock_notebooks/` are stored with Git LFS. GitHub's Contents and raw-file APIs return the small LFS pointer rather than the notebook payload, so connected code-search tools cannot inspect those notebooks directly.

The repository therefore includes `scripts/export_sherlock_notebook_text.py`, which creates an ordinary-Git, UTF-8 mirror containing:

- every Markdown cell;
- every code cell;
- textual stream, error, `text/plain`, and `text/markdown` outputs;
- `NOTEBOOK_INDEX.tsv` with cell counts, sizes, hashes, and language metadata;
- `ABSOLUTE_PATH_INDEX.tsv` containing absolute paths mentioned in code, Markdown, or textual output.

The original `.ipynb` files remain unchanged in Git LFS.

## Generate the mirror on Sherlock

```bash
module load system git/2.45.1
hash -r

cd /oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction

git fetch origin
git checkout curation/figure-reproduction-v1
git pull --ff-only origin curation/figure-reproduction-v1

git lfs pull

python3 scripts/export_sherlock_notebook_text.py
```

The default output is:

```text
recovery/sherlock_notebook_text/
```

Review the generated indexes:

```bash
column -ts $'\t' recovery/sherlock_notebook_text/NOTEBOOK_INDEX.tsv | less -S
column -ts $'\t' recovery/sherlock_notebook_text/ABSOLUTE_PATH_INDEX.tsv | less -S
```

Commit and push the searchable mirror:

```bash
git add recovery/sherlock_notebook_text
git commit -m "Add searchable text mirror of Sherlock notebooks"
git push origin curation/figure-reproduction-v1
```

The script applies conservative redaction to obvious token, API-key, and password assignments. Disable this only after deliberate review:

```bash
python3 scripts/export_sherlock_notebook_text.py --no-redaction
```

Textual outputs are included by default because executed notebooks often contain the resolved input paths, table dimensions, selected parameters, and final numerical values needed for figure reproduction. To omit them:

```bash
python3 scripts/export_sherlock_notebook_text.py --no-text-outputs
```
