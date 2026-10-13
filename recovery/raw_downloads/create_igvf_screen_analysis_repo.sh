#!/usr/bin/env bash
# Create a private EngreitzLab GitHub repository containing the existing
# genome-wide and validation-screen analysis code documented in Dulguun's
# offboarding file. This script copies CODE and documentation only; it does
# not copy raw sequencing data or analysis result tables.

set -Eeuo pipefail
IFS=$'\n\t'

ORG="${ORG:-EngreitzLab}"
REPO_NAME="${REPO_NAME:-hiPSC-EC-genomewide-validation-screen-analysis}"
LOCAL_ROOT="${LOCAL_ROOT:-/oak/stanford/groups/engreitz/Users/tri/repos}"
VISIBILITY="${VISIBILITY:-private}"
CREATE_REMOTE=1
FORCE_REBUILD=0
USE_EXISTING_REMOTE=0
ALLOW_EMPTY=0
REVIEWER_GITHUB_USER="${REVIEWER_GITHUB_USER:-}"
GH_ENV="${GH_ENV:-/oak/stanford/groups/engreitz/Users/tri/envs/github-cli}"
MAX_COPY_MB="${MAX_COPY_MB:-20}"

usage() {
  cat <<'EOF'
Usage:
  bash create_igvf_screen_analysis_repo.sh [options]

Options:
  --repo-name NAME        Repository name (default: hiPSC-EC-genomewide-validation-screen-analysis)
  --local-root PATH       Parent directory for the local clone
  --org NAME              GitHub organization (default: EngreitzLab)
  --private               Create a private repository (default)
  --internal              Create an internal repository
  --public                Create a public repository (not recommended for first pass)
  --local-only            Build and commit locally, but do not create/push a GitHub repository
  --force-rebuild         Move an existing local directory to a timestamped backup and rebuild
  --use-existing-remote   Push to an already-existing ORG/REPO instead of refusing
  --allow-empty           Permit remote creation even if no source analysis code was found
  --reviewer USER         Optionally invite a GitHub user with pull access after creation
  -h, --help              Show this help

Environment overrides:
  ORG, REPO_NAME, LOCAL_ROOT, VISIBILITY, REVIEWER_GITHUB_USER,
  GH_ENV, MAX_COPY_MB
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo-name) REPO_NAME="$2"; shift 2 ;;
    --local-root) LOCAL_ROOT="$2"; shift 2 ;;
    --org) ORG="$2"; shift 2 ;;
    --private) VISIBILITY="private"; shift ;;
    --internal) VISIBILITY="internal"; shift ;;
    --public) VISIBILITY="public"; shift ;;
    --local-only) CREATE_REMOTE=0; shift ;;
    --force-rebuild) FORCE_REBUILD=1; shift ;;
    --use-existing-remote) USE_EXISTING_REMOTE=1; shift ;;
    --allow-empty) ALLOW_EMPTY=1; shift ;;
    --reviewer) REVIEWER_GITHUB_USER="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

case "$VISIBILITY" in
  private|internal|public) ;;
  *) echo "VISIBILITY must be private, internal, or public." >&2; exit 2 ;;
esac

LOCAL_REPO="$LOCAL_ROOT/$REPO_NAME"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

command -v python3 >/dev/null 2>&1 || die "python3 is required."
command -v git >/dev/null 2>&1 || die "git is required."

# Install/use GitHub CLI through a small mamba environment if it is absent.
if [[ "$CREATE_REMOTE" -eq 1 ]] && ! command -v gh >/dev/null 2>&1; then
  command -v mamba >/dev/null 2>&1 || die "gh is missing and mamba is unavailable. Install GitHub CLI first."
  if [[ ! -x "$GH_ENV/bin/gh" ]]; then
    log "Installing GitHub CLI in $GH_ENV"
    mamba create -y -p "$GH_ENV" -c conda-forge gh git
  fi
  export PATH="$GH_ENV/bin:$PATH"
  hash -r
fi

if [[ "$CREATE_REMOTE" -eq 1 ]]; then
  gh auth status --hostname github.com >/dev/null 2>&1 || {
    cat >&2 <<'EOF'
GitHub CLI is not authenticated. Run:

  gh auth login --hostname github.com --git-protocol ssh --web

Then rerun this script.
EOF
    exit 1
  }
  GH_LOGIN="$(gh api user --jq .login)"
  log "GitHub CLI authenticated as: $GH_LOGIN"
  gh api "orgs/$ORG" >/dev/null 2>&1 || die "Cannot access GitHub organization $ORG with the current account."
  gh config set git_protocol ssh --host github.com >/dev/null
fi

mkdir -p "$LOCAL_ROOT"
if [[ -e "$LOCAL_REPO" ]]; then
  if [[ "$FORCE_REBUILD" -eq 1 ]]; then
    BACKUP="${LOCAL_REPO}.backup_${TIMESTAMP}"
    log "Moving existing local directory to $BACKUP"
    mv "$LOCAL_REPO" "$BACKUP"
  else
    die "$LOCAL_REPO already exists. Use --force-rebuild or choose another --repo-name."
  fi
fi

log "Creating repository scaffold at $LOCAL_REPO"
mkdir -p \
  "$LOCAL_REPO/docs" \
  "$LOCAL_REPO/config" \
  "$LOCAL_REPO/environment" \
  "$LOCAL_REPO/scripts" \
  "$LOCAL_REPO/genome_wide/notebooks" \
  "$LOCAL_REPO/genome_wide/plotting" \
  "$LOCAL_REPO/genome_wide/discovered" \
  "$LOCAL_REPO/validation/plotting" \
  "$LOCAL_REPO/validation/discovered"

cat > "$LOCAL_REPO/.gitignore" <<'EOF'
# Raw sequencing and large genomics files
*.fastq
*.fastq.gz
*.fq
*.fq.gz
*.bam
*.bai
*.cram
*.crai
*.sam
*.bed.gz
*.gtf.gz
*.tar
*.tar.gz

# Generated results and local data
/data/
/results/
/output/
/outputs/
/figures/
/logs/
/tmp/
*.log

# Credentials and local configuration
.env
.env.*
!.env.example
*secret*
*credential*

# Python / R / notebooks
__pycache__/
*.py[cod]
.ipynb_checkpoints/
.Rhistory
.RData
.Rproj.user/

# OS/editor
.DS_Store
.vscode/
.idea/
EOF

cat > "$LOCAL_REPO/README.md" <<'EOF'
# hiPSC-EC genome-wide and validation CRISPRi screen analysis

This repository is a provenance snapshot of the analysis code used for the
hiPSC-to-endothelial-cell genome-wide CRISPRi differentiation screen and its
focused validation screen.

## Biological comparisons

### Genome-wide screen

- Primary differentiation comparison: Day 3 CDH5-positive versus Day 3
  CDH5-negative fractions.
- Survival comparison: Day 3 total (combined CDH5-positive and CDH5-negative)
  versus Day 0.
- Analysis methods documented in the source handoff: MAGeCK and a U-test
  workflow.

### Validation screen

- Primary comparison: Day 3 CDH5-positive versus Day 3 CDH5-negative fractions.
- Analysis methods documented in the source handoff: MAGeCK and a U-test
  workflow, followed by cross-screen summaries and volcano/correlation plots.

## Repository organization

- `genome_wide/`: copied genome-wide analysis notebooks and plotting scripts.
- `validation/`: copied validation-screen guide-counting, MAGeCK, U-test, and
  plotting code.
- `docs/SOURCE_PATHS.tsv`: source inputs, outputs, and code locations recorded
  in the offboarding document.
- `docs/CODE_SNAPSHOT_MANIFEST.tsv`: exact source-to-repository copy manifest,
  file sizes, and SHA-256 checksums.
- `docs/REPRODUCIBILITY_GAPS.md`: items that still need confirmation before the
  repository can be described as a fully reproducible rerun.
- `config/paths.example.env`: configurable Sherlock/OAK paths.
- `scripts/check_source_files.sh`: verifies documented Sherlock paths.
- `scripts/snapshot_versions.sh`: records currently available software versions.
- `scripts/prepush_safety_check.sh`: rejects obvious raw-data or credential
  files before a push.

## Data policy

Raw FASTQs, count matrices, MAGeCK result tables, and other large or controlled
files are not committed. The analysis code expects those files to remain on
Sherlock/OAK or to be supplied through paths in `config/paths.example.env`.

## Current limitations

This first-pass repository preserves the existing analysis code rather than
rewriting undocumented commands. Exact MAGeCK versions, full command lines,
U-test parameters, sample-group definitions, and the original runtime
environment must be confirmed from the copied code and run logs before claiming
full computational reproducibility.

A known source note states that `NEDD8` is duplicated in some genome-wide result
tables because of inconsistent use of an underscore in the target name.
EOF

cat > "$LOCAL_REPO/docs/WORKFLOW.md" <<'EOF'
# Analysis workflow represented by this repository

## Step 1: Guide-read processing and guide quantification

Inputs are guide-sequencing FASTQs, sample metadata, and the corresponding
library design. The output is a guide-by-sample count table used by MAGeCK and
the U-test workflow.

## Step 2: Genome-wide screen analysis

The genome-wide screen evaluates:

1. Day 3 CDH5-positive versus Day 3 CDH5-negative guide frequencies.
2. Day 3 total versus Day 0 guide frequencies as a survival-related comparison.

The documented methods are MAGeCK and a U-test analysis. Native result files
are kept outside Git and are listed in `docs/SOURCE_PATHS.tsv`.

## Step 3: Validation-screen analysis

The validation screen evaluates Day 3 CDH5-positive versus Day 3 CDH5-negative
fractions with MAGeCK and U-test workflows. Downstream code combines validation
and genome-wide effect estimates for summary plots and hit comparisons.

## Step 4: Reporting

R scripts generate volcano plots, correlation plots, rank plots, and additional
screen summaries. Plot outputs are not committed by default.

## Proposed IGVF workflow mapping

For IGVF metadata, the code can be represented by separate AnalysisStep records
for guide quantification, MAGeCK promoter-level enrichment, U-test
promoter-level enrichment, and optional reporting/cross-screen comparison.
The precise SoftwareVersion and AnalysisStepVersion records should be based on
versions and parameters confirmed from this repository snapshot and the
original run logs.
EOF

cat > "$LOCAL_REPO/docs/REPRODUCIBILITY_GAPS.md" <<'EOF'
# Reproducibility gaps to resolve before final IGVF workflow submission

The offboarding document identifies code and result locations but does not fully
specify every executable command. Confirm the following from copied notebooks,
shell histories, run logs, or the original analysts:

1. Exact MAGeCK version.
2. Exact `mageck count` and `mageck test` commands.
3. Treatment/control sample names and replicate grouping for each comparison.
4. Normalization method and any control-guide list.
5. Filtering thresholds for guides or promoters.
6. Guide-to-promoter aggregation rules, including P1/P2/P1P2 handling.
7. Exact U-test implementation, grouping rule, sidedness, and multiple-testing
   correction.
8. The normalization factor used in the validation summary and how it was
   derived.
9. R and Python package versions used for plots and table assembly.
10. Whether the additional sequencing directories were merged or one was a
    redundant rerun.
11. Resolution of the documented `NEDD8` duplicate-name issue.
12. A definitive primary output for the validation-screen per-promoter MAGeCK
    analysis rather than only a gene-level combined summary.

Do not fill these gaps by guessing. Record confirmed values in a release or an
IGVF AnalysisStepVersion description.
EOF

cat > "$LOCAL_REPO/docs/REFERENCES.md" <<'EOF'
# Internal reference documents

These links were recorded in Dulguun Amgalan's offboarding document and may
require Stanford/Engreitz Lab access.

- Genome-wide screen Results/Methods:
  https://docs.google.com/document/d/1NgqKp3x4moxZVe19AacY-AE8uGzuYr9KlWQ18wTtLjQ/edit
- Genome-wide screen Drive folder:
  https://drive.google.com/drive/folders/1zPQ12mFrq6nJQifHXoeaXxJ9tBJjy9Ph
- MAGeCK pipeline explainer:
  https://docs.google.com/presentation/d/1dM6c_juSq6LS5wA-v3Xyd4dbH3F_3LmN4_c5P8YYXwg/edit
- Validation-screen Results/Methods:
  https://docs.google.com/document/d/1d3IrYeiHzimtv5uKlFN6avjEScjcnQycprJSJZ7qnhs/edit
- Validation-screen experimental protocol:
  https://docs.google.com/document/d/1Y-SYDgI2tM39KwZFbOC8CFm0jOHtFtAEVJlkMvpjpcI/edit
- Genome-wide validation Drive folder:
  https://drive.google.com/drive/folders/1zH5YU3KINArxxl34YdgM_82fp_BbE5VH
EOF

cat > "$LOCAL_REPO/docs/SOURCE_PATHS.tsv" <<'EOF'
screen	category	role	path
genome_wide	input	guide_count_input	/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/Files_used_for_mageck/combined_mageck_input.txt
genome_wide	code	utest_pipeline_notebook	/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/from_chad/mageck/230523_Pipeline.ipynb
genome_wide	output	mageck_d3pos_vs_d3neg	/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/Pos_vs_Neg/Pos_vs_Neg.gene_summary.txt
genome_wide	output	utest_d3pos_vs_d3neg	/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/Pos_vs_Neg/combined_utest_P1P2_adj_test.txt
genome_wide	output	mageck_d3total_vs_d0	/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/D0_vs_Total/D0_vs_EC_total.gene_summary.txt
genome_wide	output	utest_d3total_vs_d0	/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/D0_vs_Total/combined_utest_P1P2_adj_test.txt
genome_wide	documentation	mageck_output_readme	/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/Mageck_output_README.xlsx
genome_wide	code	volcano_plot_script	/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/generatePlots.R
genome_wide	code	mageckflute_plot_script	/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/mageckflute/230920_MAGECK_Flute.R
validation	code	guide_read_processing	/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/CRISPRiFlowFish
validation	code	mageck_analysis	/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/MageckAnalysis
validation	input	mageck_guide_counts	/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/MageckAnalysis/mageck_sgRNA_readcount.txt
validation	code	utest_analysis	/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/UTestAnalysis
validation	input	utest_guide_counts	/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/UTestAnalysis/guideCounts.csv
validation	output	combined_validation_summary	/oak/stanford/groups/engreitz/Users/dulguun/231023_Validation_Screen_Files/d3pos_d3neg_mageckLFC.csv
validation	output	hit_summary	/oak/stanford/groups/engreitz/Users/dulguun/231023_Validation_Screen_Files/ec_screen_hits.csv
validation	code	plotting_script	/oak/stanford/groups/engreitz/Users/dulguun/231023_Validation_Screen_Files/generatePlots.R
EOF

cat > "$LOCAL_REPO/config/paths.example.env" <<'EOF'
# Copy this file to a location outside Git and edit as needed.

export GW_ROOT=/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis
export GW_GUIDE_COUNTS="$GW_ROOT/full_result_files_from_Chad/Files_used_for_mageck/combined_mageck_input.txt"
export GW_MAGECK_POS_NEG="$GW_ROOT/full_result_files_from_Chad/Pos_vs_Neg/Pos_vs_Neg.gene_summary.txt"
export GW_UTEST_POS_NEG="$GW_ROOT/full_result_files_from_Chad/Pos_vs_Neg/combined_utest_P1P2_adj_test.txt"
export GW_MAGECK_D0_TOTAL="$GW_ROOT/full_result_files_from_Chad/D0_vs_Total/D0_vs_EC_total.gene_summary.txt"
export GW_UTEST_D0_TOTAL="$GW_ROOT/full_result_files_from_Chad/D0_vs_Total/combined_utest_P1P2_adj_test.txt"

export VALIDATION_ROOT=/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation
export VALIDATION_MAGECK_COUNTS="$VALIDATION_ROOT/MageckAnalysis/mageck_sgRNA_readcount.txt"
export VALIDATION_UTEST_COUNTS="$VALIDATION_ROOT/UTestAnalysis/guideCounts.csv"
export VALIDATION_SUMMARY=/oak/stanford/groups/engreitz/Users/dulguun/231023_Validation_Screen_Files/d3pos_d3neg_mageckLFC.csv
EOF

cat > "$LOCAL_REPO/environment/environment.template.yml" <<'EOF'
# Provisional environment only. Pin versions after confirming the original run.
name: hipsc-ec-screen-analysis
channels:
  - conda-forge
  - bioconda
dependencies:
  - python
  - jupyterlab
  - pandas
  - numpy
  - scipy
  - statsmodels
  - openpyxl
  - mageck
  - r-base
  - r-tidyverse
  - r-readxl
  - r-ggrepel
EOF

cat > "$LOCAL_REPO/environment/README.md" <<'EOF'
# Environment notes

`environment.template.yml` is a convenience template, not a reconstruction of
the original environment. Run `scripts/snapshot_versions.sh` in the original
Sherlock environment and pin versions only after confirming which environment
was used to generate the submitted results.
EOF

cat > "$LOCAL_REPO/scripts/check_source_files.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TABLE="$ROOT/docs/SOURCE_PATHS.tsv"

printf 'status\tscreen\tcategory\trole\tpath\n'
tail -n +2 "$TABLE" | while IFS=$'\t' read -r screen category role path; do
  if [[ -e "$path" ]]; then
    status=FOUND
  else
    status=MISSING
  fi
  printf '%s\t%s\t%s\t%s\t%s\n' "$status" "$screen" "$category" "$role" "$path"
done
EOF
chmod +x "$LOCAL_REPO/scripts/check_source_files.sh"

cat > "$LOCAL_REPO/scripts/snapshot_versions.sh" <<'EOF'
#!/usr/bin/env bash
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/environment/software_versions_$(date +%Y%m%d_%H%M%S).txt"
{
  echo "Snapshot date: $(date -Is)"
  echo "Host: $(hostname)"
  echo
  for cmd in python python3 R Rscript mageck git; do
    echo "===== $cmd ====="
    if command -v "$cmd" >/dev/null 2>&1; then
      command -v "$cmd"
      case "$cmd" in
        R|Rscript) "$cmd" --version 2>&1 | head -n 4 ;;
        mageck) "$cmd" --version 2>&1 || "$cmd" -v 2>&1 || true ;;
        *) "$cmd" --version 2>&1 | head -n 4 ;;
      esac
    else
      echo "NOT FOUND"
    fi
    echo
  done
  if command -v conda >/dev/null 2>&1; then
    echo "===== conda list ====="
    conda list 2>&1 || true
  elif command -v mamba >/dev/null 2>&1; then
    echo "===== mamba list ====="
    mamba list 2>&1 || true
  fi
} > "$OUT"
echo "$OUT"
EOF
chmod +x "$LOCAL_REPO/scripts/snapshot_versions.sh"

cat > "$LOCAL_REPO/scripts/prepush_safety_check.sh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

bad=0

# Reject large tracked files.
while IFS= read -r -d '' file; do
  size=$(stat -c %s "$file")
  if (( size > 20 * 1024 * 1024 )); then
    echo "LARGE TRACKED FILE: $file ($size bytes)" >&2
    bad=1
  fi
done < <(git ls-files -z)

# Reject common raw-data extensions if tracked.
if git ls-files | grep -Ei '\.(fastq|fq)(\.gz)?$|\.(bam|cram|sam)$|\.tar(\.gz)?$' >/dev/null; then
  echo "Raw sequencing/archive file appears to be tracked:" >&2
  git ls-files | grep -Ei '\.(fastq|fq)(\.gz)?$|\.(bam|cram|sam)$|\.tar(\.gz)?$' >&2
  bad=1
fi

# Look for common credential assignments. This is intentionally conservative.
if git grep -nEI '(IGVF_SECRET_KEY|ANTHROPIC_API_KEY|GITHUB_TOKEN|ghp_[A-Za-z0-9]{20,}|github_pat_)' -- . \
   ':(exclude)scripts/prepush_safety_check.sh' >/tmp/igvf_repo_secret_scan.$$ 2>/dev/null; then
  echo "Possible credential material detected:" >&2
  cat /tmp/igvf_repo_secret_scan.$$ >&2
  bad=1
fi
rm -f /tmp/igvf_repo_secret_scan.$$

if (( bad )); then
  echo "Safety check FAILED." >&2
  exit 1
fi

echo "Safety check passed."
EOF
chmod +x "$LOCAL_REPO/scripts/prepush_safety_check.sh"

log "Copying documented analysis code and creating a checksum manifest"
REPO_DIR="$LOCAL_REPO" MAX_COPY_MB="$MAX_COPY_MB" python3 <<'PY'
from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
from pathlib import Path

repo = Path(os.environ["REPO_DIR"])
max_bytes = int(os.environ.get("MAX_COPY_MB", "20")) * 1024 * 1024
manifest_path = repo / "docs" / "CODE_SNAPSHOT_MANIFEST.tsv"

allowed_suffixes = {
    ".py", ".r", ".R", ".ipynb", ".sh", ".bash", ".zsh", ".pl", ".awk",
    ".yaml", ".yml", ".json", ".toml", ".ini", ".cfg", ".md", ".Rmd",
    ".rmd", ".qmd",
}
text_name_tokens = ("readme", "command", "parameter", "config", "pipeline", "script", "sample")
skip_parts = {".git", "__pycache__", ".ipynb_checkpoints", "node_modules", ".pixi", ".conda"}

exact_files = [
    (
        "genome_wide/notebooks/230523_Pipeline.ipynb",
        [
            "/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/from_chad/mageck/230523_Pipeline.ipynb",
            "/oak/stanford/groups/engreitz/User/dulguun/230916_Genome-wide_Screen_Analysis/from_chad/mageck/230523_Pipeline.ipynb",
        ],
        "core genome-wide U-test/pipeline notebook",
    ),
    (
        "genome_wide/plotting/generatePlots.R",
        ["/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/generatePlots.R"],
        "genome-wide volcano plotting",
    ),
    (
        "genome_wide/plotting/230920_MAGECK_Flute.R",
        ["/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/mageckflute/230920_MAGECK_Flute.R"],
        "genome-wide rank and MAGeCKFlute plots",
    ),
    (
        "validation/plotting/generatePlots.R",
        ["/oak/stanford/groups/engreitz/Users/dulguun/231023_Validation_Screen_Files/generatePlots.R"],
        "validation volcano and correlation plots",
    ),
    (
        "docs/source/Mageck_output_README.xlsx",
        ["/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/Mageck_output_README.xlsx"],
        "MAGeCK output column documentation",
    ),
]

scan_dirs = [
    (
        "genome_wide/discovered/mageck",
        [
            "/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/from_chad/mageck",
            "/oak/stanford/groups/engreitz/User/dulguun/230916_Genome-wide_Screen_Analysis/from_chad/mageck",
        ],
        "genome-wide MAGeCK/U-test analysis directory",
    ),
    (
        "genome_wide/discovered/files_used_for_mageck",
        [
            "/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/Files_used_for_mageck",
            "/oak/stanford/groups/engreitz/User/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/Files_used_for_mageck",
        ],
        "code/configuration adjacent to genome-wide MAGeCK inputs",
    ),
    (
        "validation/discovered/guide_read_processing",
        ["/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/CRISPRiFlowFish"],
        "validation FASTQ/guide-count processing",
    ),
    (
        "validation/discovered/mageck",
        ["/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/MageckAnalysis"],
        "validation MAGeCK analysis",
    ),
    (
        "validation/discovered/utest",
        ["/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/UTestAnalysis"],
        "validation U-test analysis",
    ),
    (
        "validation/discovered/summary_and_plots",
        ["/oak/stanford/groups/engreitz/Users/dulguun/231023_Validation_Screen_Files"],
        "validation summary/plotting code",
    ),
]

rows: list[dict[str, str | int]] = []
seen_destinations: set[Path] = set()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def copy_notebook_without_outputs(src: Path, dst: Path) -> None:
    try:
        notebook = json.loads(src.read_text(encoding="utf-8"))
        for cell in notebook.get("cells", []):
            if cell.get("cell_type") == "code":
                cell["outputs"] = []
                cell["execution_count"] = None
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(json.dumps(notebook, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    except Exception:
        # Preserve the file if it is not standard JSON, but still avoid modifying the source.
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def copy_one(src: Path, rel_dest: Path, role: str, discovery: str) -> None:
    dst = repo / rel_dest
    if dst in seen_destinations:
        return
    seen_destinations.add(dst)

    if not src.exists():
        rows.append({
            "status": "MISSING",
            "source_path": str(src),
            "repository_path": str(rel_dest),
            "role": role,
            "discovery": discovery,
            "size_bytes": "",
            "sha256": "",
            "note": "source path not found",
        })
        return
    if not src.is_file():
        rows.append({
            "status": "SKIPPED",
            "source_path": str(src),
            "repository_path": str(rel_dest),
            "role": role,
            "discovery": discovery,
            "size_bytes": "",
            "sha256": "",
            "note": "not a regular file",
        })
        return
    size = src.stat().st_size
    if size > max_bytes:
        rows.append({
            "status": "SKIPPED_LARGE",
            "source_path": str(src),
            "repository_path": str(rel_dest),
            "role": role,
            "discovery": discovery,
            "size_bytes": size,
            "sha256": "",
            "note": f"larger than {max_bytes} bytes",
        })
        return

    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.suffix.lower() == ".ipynb":
        copy_notebook_without_outputs(src, dst)
        note = "copied with notebook outputs/execution counts removed"
    else:
        shutil.copy2(src, dst)
        note = "copied unchanged"
    rows.append({
        "status": "COPIED",
        "source_path": str(src),
        "repository_path": str(rel_dest),
        "role": role,
        "discovery": discovery,
        "size_bytes": dst.stat().st_size,
        "sha256": sha256(dst),
        "note": note,
    })


for rel, candidates, role in exact_files:
    existing = next((Path(p) for p in candidates if Path(p).is_file()), None)
    if existing is None:
        # Record the preferred Users path as missing.
        copy_one(Path(candidates[0]), Path(rel), role, "exact")
    else:
        copy_one(existing, Path(rel), role, "exact")


def eligible(path: Path) -> bool:
    if any(part in skip_parts for part in path.parts):
        return False
    if path.suffix in allowed_suffixes or path.suffix.lower() in {s.lower() for s in allowed_suffixes}:
        return True
    if path.name in {"Makefile", "Snakefile"}:
        return True
    if path.suffix.lower() == ".txt" and any(token in path.name.lower() for token in text_name_tokens):
        return True
    return False


for dest_root, candidates, role in scan_dirs:
    source_root = next((Path(p) for p in candidates if Path(p).is_dir()), None)
    if source_root is None:
        rows.append({
            "status": "MISSING_DIRECTORY",
            "source_path": candidates[0],
            "repository_path": dest_root,
            "role": role,
            "discovery": "recursive_allowlist",
            "size_bytes": "",
            "sha256": "",
            "note": "no candidate directory found",
        })
        continue
    for src in sorted(source_root.rglob("*")):
        if not src.is_file() or not eligible(src):
            continue
        relative = src.relative_to(source_root)
        copy_one(src, Path(dest_root) / relative, role, "recursive_allowlist")

manifest_path.parent.mkdir(parents=True, exist_ok=True)
fieldnames = [
    "status", "source_path", "repository_path", "role", "discovery",
    "size_bytes", "sha256", "note",
]
with manifest_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

copied = sum(row["status"] == "COPIED" for row in rows)
missing = sum(str(row["status"]).startswith("MISSING") for row in rows)
print(f"Manifest: {manifest_path}")
print(f"Copied files: {copied}")
print(f"Missing source entries/directories: {missing}")
PY

COPIED_COUNT="$(awk -F'\t' 'NR>1 && $1=="COPIED" {n++} END {print n+0}' "$LOCAL_REPO/docs/CODE_SNAPSHOT_MANIFEST.tsv")"
MISSING_COUNT="$(awk -F'\t' 'NR>1 && $1 ~ /^MISSING/ {n++} END {print n+0}' "$LOCAL_REPO/docs/CODE_SNAPSHOT_MANIFEST.tsv")"
log "Code files copied: $COPIED_COUNT; missing source entries/directories: $MISSING_COUNT"

if [[ "$COPIED_COUNT" -eq 0 && "$ALLOW_EMPTY" -ne 1 ]]; then
  die "No analysis code was found. Review docs/CODE_SNAPSHOT_MANIFEST.tsv or rerun with --allow-empty."
fi

# Initialize and commit.
cd "$LOCAL_REPO"
git init -b main >/dev/null

GIT_NAME="$(git config --global user.name || true)"
GIT_EMAIL="$(git config --global user.email || true)"
if [[ -z "$GIT_NAME" ]]; then
  read -r -p "Git commit name: " GIT_NAME
fi
if [[ -z "$GIT_EMAIL" ]]; then
  read -r -p "Git commit email: " GIT_EMAIL
fi
[[ -n "$GIT_NAME" ]] || die "Git commit name is required."
[[ -n "$GIT_EMAIL" ]] || die "Git commit email is required."
git config user.name "$GIT_NAME"
git config user.email "$GIT_EMAIL"

# Run the safety scan before the first commit.
git add .
"$LOCAL_REPO/scripts/prepush_safety_check.sh"
git commit -m "Initial analysis-code provenance snapshot" >/dev/null

log "Local repository created and committed: $LOCAL_REPO"
log "Review the copy manifest: $LOCAL_REPO/docs/CODE_SNAPSHOT_MANIFEST.tsv"

if [[ "$CREATE_REMOTE" -eq 0 ]]; then
  log "Local-only mode selected. No GitHub repository was created."
  exit 0
fi

REMOTE_SLUG="$ORG/$REPO_NAME"
if gh repo view "$REMOTE_SLUG" >/dev/null 2>&1; then
  if [[ "$USE_EXISTING_REMOTE" -ne 1 ]]; then
    die "GitHub repository $REMOTE_SLUG already exists. Use --use-existing-remote or choose another name."
  fi
  log "Using existing GitHub repository: $REMOTE_SLUG"
  git remote add origin "git@github.com:${REMOTE_SLUG}.git"
  git push -u origin main
else
  log "Creating $VISIBILITY GitHub repository: $REMOTE_SLUG"
  gh repo create "$REMOTE_SLUG" \
    "--$VISIBILITY" \
    --description "Analysis code for genome-wide and validation CRISPRi screens of hiPSC-to-endothelial differentiation" \
    --source "$LOCAL_REPO" \
    --remote origin \
    --push
fi

if [[ -n "$REVIEWER_GITHUB_USER" ]]; then
  log "Inviting $REVIEWER_GITHUB_USER with pull access"
  gh api \
    --method PUT \
    -H "Accept: application/vnd.github+json" \
    "repos/$REMOTE_SLUG/collaborators/$REVIEWER_GITHUB_USER" \
    -f permission='pull' >/dev/null
fi

REPO_URL="$(gh repo view "$REMOTE_SLUG" --json url --jq .url)"
cat <<EOF

Repository created successfully:
  $REPO_URL

Next review steps:
  1. Open docs/CODE_SNAPSHOT_MANIFEST.tsv and resolve missing code paths.
  2. Review docs/REPRODUCIBILITY_GAPS.md with the original analysts.
  3. Run scripts/snapshot_versions.sh in the original analysis environment.
  4. Send the repository URL above to Ian.
EOF
