# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/validation_screen_guide_level_EC_efficiency_edgeR.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/validation_screen_guide_level_EC_efficiency_edgeR.ipynb`
- SHA-256: `4931d88d25f68d253fe35ef5b318d9f8ba4b8975e8b6ecdd5737eb59c6519392`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# Validation-screen guide-level EC differentiation analysis with edgeR — fixed

This notebook is a corrected version of the guide-level validation-screen workflow. It fixes the failure:

```text
RuntimeError: No included POS/NEG samples were identified
```

The sample sheet is now created using:

1. exact user-supplied overrides when provided;
2. robust parsing of sample names;
3. the documented validation-screen PCR-library layout when the column names are cryptic.

The documented layout is:

- biological replicate 1: 7 D0, 7 CDH5-positive, and 7 CDH5-negative PCR libraries;
- biological replicate 2: 7 D0, 7 CDH5-positive, and 6 CDH5-negative PCR libraries.

PCR libraries are treated as **technical replicates** and summed before edgeR. The model is fitted to four biological samples using:

```text
~ biological_replicate + fraction
```

with `NEG` as the reference and `fractionPOS` as the tested coefficient.

## Code cell 1

```python
from pathlib import Path
import json
import re
import shutil
import subprocess
import warnings

import numpy as np
import pandas as pd
from IPython.display import display, Markdown

warnings.filterwarnings("ignore")
pd.set_option("display.max_columns", 200)
pd.set_option("display.max_rows", 200)

class cfg:
    # Documented validation-screen count tables on Sherlock.
    RAW_COUNT_CANDIDATES = [
        Path(
            "/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/"
            "MageckAnalysis/mageck_sgRNA_readcount.txt"
        ),
        Path(
            "/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/"
            "UTestAnalysis/guideCounts.csv"
        ),
    ]

    OUTPUT_DIR = Path(
        "/oak/stanford/groups/engreitz/Users/edtng/GenomeWideValidation/"
        "guide_level_EC_efficiency_edgeR_FIXED"
    )

    # Set to a CSV/TSV only when you already have a checked sample sheet.
    # Required columns: sample, fraction, biological_replicate, include.
    SAMPLE_SHEET_PATH = None

    # Exact overrides take priority over every automatic rule.
    # Example:
    # MANUAL_SAMPLE_MAP = {
    #     "exact_count_column_name": ("POS", "B1", 1),
    # }
    MANUAL_SAMPLE_MAP = {}

    # This fallback prevents the original POS/NEG parsing error when the count
    # matrix contains the documented 41 PCR-library columns with cryptic names.
    USE_DOCUMENTED_41_COLUMN_ORDER_FALLBACK = True

    # edgeR settings.
    MIN_TOTAL_COUNT = 10
    R_SCRIPT_BIN = "Rscript"

    # Control guides used to center guide-level effects.
    CONTROL_REGEX = (
        r"(?i)(non.?target|safe.?target|safe[_ -]?targeting|ntc|"
        r"negative.?control|control)"
    )

    # Whole-population validation-screen FACS efficiencies. These are used only
    # for optional approximate conversion of relative odds to EC percentage.
    FACS_EC_FRACTION_BY_REPLICATE = {"B1": 0.844, "B2": 0.721}

cfg.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
print("Output directory:", cfg.OUTPUT_DIR)

```

### Text output (stream)

```text
Output directory: /oak/stanford/groups/engreitz/Users/edtng/GenomeWideValidation/guide_level_EC_efficiency_edgeR_FIXED

```

## Markdown cell 2

## 1. Load the raw guide-count table

The loader accepts either the MAGeCK-formatted count table or the U-test guide-count table. It identifies guide/gene annotation columns and retains numeric sample columns.

## Code cell 3

```python
def read_table_flexible(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path)
    else:
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception:
            df = pd.read_csv(path, sep=None, engine="python")
    df.columns = [str(c).strip() for c in df.columns]
    return df


def choose_count_path(candidates):
    existing = [p for p in candidates if p.exists()]
    if not existing:
        raise FileNotFoundError(
            "None of the documented validation-screen count files exists:\n"
            + "\n".join(str(p) for p in candidates)
        )
    return existing[0]


def identify_annotation_and_sample_columns(df: pd.DataFrame):
    lower_to_original = {str(c).lower(): c for c in df.columns}
    guide_names = ["sgrna", "guide", "guide_id", "id"]
    gene_names = ["gene", "target", "target_gene"]

    guide_col = next(
        (lower_to_original[x] for x in guide_names if x in lower_to_original),
        None,
    )
    gene_col = next(
        (lower_to_original[x] for x in gene_names if x in lower_to_original),
        None,
    )
    if guide_col is None:
        raise ValueError(
            "Could not identify a guide-ID column. Columns are:\n"
            + "\n".join(df.columns)
        )

    excluded = {guide_col}
    if gene_col is not None:
        excluded.add(gene_col)

    annotation_tokens = [
        "sequence", "spacer", "chrom", "strand", "start", "end",
        "symbol", "annotation",
    ]
    sample_cols = []
    for col in df.columns:
        if col in excluded:
            continue
        if any(token in col.lower() for token in annotation_tokens):
            continue
        numeric = pd.to_numeric(df[col], errors="coerce")
        if numeric.notna().mean() >= 0.95:
            sample_cols.append(col)

    if not sample_cols:
        raise ValueError("No numeric sample columns were detected.")
    return guide_col, gene_col, sample_cols


RAW_COUNTS_PATH = choose_count_path(cfg.RAW_COUNT_CANDIDATES)
raw_table = read_table_flexible(RAW_COUNTS_PATH)
GUIDE_COL, GENE_COL, SAMPLE_COLUMNS = identify_annotation_and_sample_columns(raw_table)

raw_counts = (
    raw_table[SAMPLE_COLUMNS]
    .apply(pd.to_numeric, errors="coerce")
    .fillna(0)
)
if (raw_counts < 0).any().any():
    raise ValueError("Negative values were found in the count matrix.")

print("Count file:", RAW_COUNTS_PATH)
print("Guide column:", GUIDE_COL)
print("Gene column:", GENE_COL)
print("Guide rows:", f"{len(raw_table):,}")
print("Numeric sample columns:", len(SAMPLE_COLUMNS))
print("\nSample columns in file order:")
for i, sample in enumerate(SAMPLE_COLUMNS, start=1):
    print(f"{i:2d}. {sample}")

display(
    raw_table[
        [c for c in [GUIDE_COL, GENE_COL] if c is not None]
        + SAMPLE_COLUMNS[:5]
    ].head()
)

```

### Text output (stream)

```text
Count file: /oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/MageckAnalysis/mageck_sgRNA_readcount.txt
Guide column: sgRNA
Gene column: Gene
Guide rows: 1,926
Numeric sample columns: 8

Sample columns in file order:
 1. D3PR1
 2. D3PR2
 3. D0R1
 4. D0R2
 5. D3NR1
 6. D3NR2
 7. D3TOTAL1
 8. D3TOTAL2

```

### Text output (text/plain)

```text
                      sgRNA   Gene  D3PR1  D3PR2  D0R1  D0R2  D3NR1
0  ACAA1_+_38178488.23-P1P2  ACAA1   1880   6101  1248   763   1225
1  ACAA1_+_38178517.23-P1P2  ACAA1   1946   7184  1008   768   1141
2  ACAA1_+_38178559.23-P1P2  ACAA1   1462   4459   447   934   1873
3  ACAA1_+_38178570.23-P1P2  ACAA1   1827   4704   577  1136   1918
4  ACAA1_+_38178620.23-P1P2  ACAA1   1579   5376   535  1036   1648
```

## Markdown cell 4

## 2. Construct and validate the sample sheet

The original notebook stopped because no columns were assigned to `POS` or `NEG`. This implementation recognizes common naming styles and, for a 41-column matrix, can use the documented order:

```text
B1: 7 D0, 7 POS, 7 NEG
B2: 7 D0, 7 POS, 6 NEG
```

Inspect the displayed table before proceeding. If the file order does not match this design, enter exact mappings in `cfg.MANUAL_SAMPLE_MAP`.

## Code cell 5

```python
def canonical_sample_name(name: str) -> str:
    s = str(name).strip().lower()
    # Preserve biological plus/minus semantics as text before punctuation removal.
    s = s.replace("cdh5+", "cdh5_pos").replace("cdh5-", "cdh5_neg")
    s = s.replace("cd144+", "cd144_pos").replace("cd144-", "cd144_neg")
    s = re.sub(r"[^a-z0-9]+", "_", s)
    return re.sub(r"_+", "_", s).strip("_")


def infer_fraction(name: str) -> str:
    s = canonical_sample_name(name)
    compact = s.replace("_", "")

    if (
        re.search(r"(^|_)d0($|_)", s)
        or "baseline" in s
        or re.search(r"(^|_)input($|_)", s)
    ):
        return "D0"

    neg_patterns = [
        r"(^|_)d3_?neg($|_)", r"(^|_)neg($|_)", r"negative",
        r"cdh5_?neg", r"cd144_?neg", r"non_?ec", r"nonec",
    ]
    pos_patterns = [
        r"(^|_)d3_?pos($|_)", r"(^|_)pos($|_)", r"positive",
        r"cdh5_?pos", r"cd144_?pos", r"ec_?pos",
    ]

    if any(re.search(p, s) for p in neg_patterns) or any(
        x in compact for x in ["d3neg", "cdh5neg", "cd144neg"]
    ):
        return "NEG"
    if any(re.search(p, s) for p in pos_patterns) or any(
        x in compact for x in ["d3pos", "cdh5pos", "cd144pos"]
    ):
        return "POS"
    return ""


def infer_biorep(name: str) -> str:
    s = canonical_sample_name(name)
    patterns = {
        "B1": [
            r"biorep_?1", r"bio_?rep_?1", r"bio_?1",
            r"(^|_)br_?1($|_)", r"(^|_)b_?1($|_)",
        ],
        "B2": [
            r"biorep_?2", r"bio_?rep_?2", r"bio_?2",
            r"(^|_)br_?2($|_)", r"(^|_)b_?2($|_)",
        ],
    }
    for rep, pats in patterns.items():
        if any(re.search(p, s) for p in pats):
            return rep
    return ""


def infer_technical_rep(name: str):
    s = canonical_sample_name(name)
    for pat in [
        r"(?:pcr|tech|technical)_?(\d+)",
        r"(?:tube|lib|library)_?(\d+)",
    ]:
        match = re.search(pat, s)
        if match:
            return int(match.group(1))
    return np.nan


def apply_documented_41_column_order(sheet: pd.DataFrame) -> pd.DataFrame:
    if len(sheet) != 41:
        raise ValueError(
            "The documented-order fallback is valid only for exactly 41 sample columns."
        )

    assignments = []
    assignments += [("D0", "B1", i) for i in range(1, 8)]
    assignments += [("POS", "B1", i) for i in range(1, 8)]
    assignments += [("NEG", "B1", i) for i in range(1, 8)]
    assignments += [("D0", "B2", i) for i in range(1, 8)]
    assignments += [("POS", "B2", i) for i in range(1, 8)]
    assignments += [("NEG", "B2", i) for i in range(1, 7)]

    out = sheet.copy()
    for idx, (fraction, biorep, techrep) in enumerate(assignments):
        # The documented order is used only to fill unresolved fields; explicit
        # name inference and manual mappings are never overwritten.
        if out.loc[idx, "fraction"] == "":
            out.loc[idx, "fraction"] = fraction
        if out.loc[idx, "biological_replicate"] == "":
            out.loc[idx, "biological_replicate"] = biorep
        if pd.isna(out.loc[idx, "technical_replicate"]):
            out.loc[idx, "technical_replicate"] = techrep
        if out.loc[idx, "mapping_source"] == "auto":
            out.loc[idx, "mapping_source"] = "documented 41-column order"
    return out


def load_explicit_sample_sheet(path: Path) -> pd.DataFrame:
    sheet = pd.read_csv(path, sep=None, engine="python")
    required = {"sample", "fraction", "biological_replicate", "include"}
    missing = required - set(sheet.columns)
    if missing:
        raise ValueError(
            "Explicit sample sheet is missing columns: " + ", ".join(sorted(missing))
        )
    if "technical_replicate" not in sheet.columns:
        sheet["technical_replicate"] = np.nan
    sheet["mapping_source"] = "explicit sample sheet"
    return sheet


def build_sample_sheet(sample_columns):
    if cfg.SAMPLE_SHEET_PATH is not None:
        source = Path(cfg.SAMPLE_SHEET_PATH)
        if not source.exists():
            raise FileNotFoundError(source)
        sheet = load_explicit_sample_sheet(source)
    else:
        sheet = pd.DataFrame({"sample": list(map(str, sample_columns))})
        sheet["fraction"] = sheet["sample"].map(infer_fraction)
        sheet["biological_replicate"] = sheet["sample"].map(infer_biorep)
        sheet["technical_replicate"] = sheet["sample"].map(infer_technical_rep)
        sheet["mapping_source"] = "auto"

    for sample, (fraction, biorep, techrep) in cfg.MANUAL_SAMPLE_MAP.items():
        if sample not in set(sheet["sample"]):
            raise KeyError(
                f"MANUAL_SAMPLE_MAP contains an unknown count column: {sample}"
            )
        mask = sheet["sample"].eq(sample)
        sheet.loc[
            mask,
            ["fraction", "biological_replicate", "technical_replicate"],
        ] = [str(fraction).upper(), str(biorep).upper(), techrep]
        sheet.loc[mask, "mapping_source"] = "manual override"

    sheet["sample"] = sheet["sample"].astype(str).str.strip()
    sheet["fraction"] = (
        sheet["fraction"].fillna("").astype(str).str.strip().str.upper()
    )
    sheet["biological_replicate"] = (
        sheet["biological_replicate"]
        .fillna("").astype(str).str.strip().str.upper()
    )

    unresolved = sheet["fraction"].eq("") | sheet["biological_replicate"].eq("")
    used_order_fallback = False
    if (
        unresolved.any()
        and len(sheet) == 41
        and cfg.USE_DOCUMENTED_41_COLUMN_ORDER_FALLBACK
    ):
        sheet = apply_documented_41_column_order(sheet)
        used_order_fallback = True

    sheet["include"] = sheet["fraction"].isin(["POS", "NEG"])

    expected = list(map(str, sample_columns))
    if set(sheet["sample"]) != set(expected):
        raise ValueError(
            "The sample sheet and raw-count sample columns are not an exact set match."
        )
    sheet = sheet.set_index("sample").loc[expected].reset_index()

    invalid_fraction = ~sheet["fraction"].isin(["D0", "POS", "NEG"])
    missing_biorep = (
        sheet["include"]
        & ~sheet["biological_replicate"].isin(["B1", "B2"])
    )

    if invalid_fraction.any() or missing_biorep.any():
        edit_path = cfg.OUTPUT_DIR / "sample_sheet_to_edit.csv"
        sheet.to_csv(edit_path, index=False)
        display(sheet)
        display(
            sheet.loc[
                invalid_fraction | missing_biorep,
                ["sample", "fraction", "biological_replicate"],
            ]
        )
        raise RuntimeError(
            "Some sample columns remain unresolved. Add exact entries to "
            "cfg.MANUAL_SAMPLE_MAP or edit the generated sheet and set "
            "cfg.SAMPLE_SHEET_PATH. Template:\n"
            f"{edit_path}"
        )

    counts = (
        sheet.groupby(["biological_replicate", "fraction"])
        .size()
        .unstack(fill_value=0)
        .reindex(
            index=["B1", "B2"],
            columns=["D0", "POS", "NEG"],
            fill_value=0,
        )
    )
    return sheet, counts, used_order_fallback


sample_sheet, sample_group_counts, used_order_fallback = build_sample_sheet(
    SAMPLE_COLUMNS
)
sample_sheet.to_csv(cfg.OUTPUT_DIR / "sample_sheet_used.csv", index=False)

display(sample_sheet)
display(Markdown("### Technical-library counts by replicate and fraction"))
display(sample_group_counts)

if used_order_fallback:
    print(
        "WARNING: the documented 41-column order fallback was used. "
        "Confirm the displayed sample sheet before interpreting results."
    )

# Validate the four-sample D3 design.
d3_sheet = sample_sheet.loc[sample_sheet["include"]].copy()
if set(d3_sheet["fraction"]) != {"POS", "NEG"}:
    raise RuntimeError("Both POS and NEG fractions are required.")
for rep in ["B1", "B2"]:
    present = set(
        d3_sheet.loc[
            d3_sheet["biological_replicate"].eq(rep), "fraction"
        ]
    )
    if present != {"POS", "NEG"}:
        raise RuntimeError(
            f"{rep} does not contain both POS and NEG samples: {present}"
        )

print("Sample sheet validated successfully.")

```

### Text output (text/plain)

```text
     sample fraction biological_replicate  technical_replicate mapping_source  \
0     D3PR1                                                NaN           auto   
1     D3PR2                                                NaN           auto   
2      D0R1                                                NaN           auto   
3      D0R2                                                NaN           auto   
4     D3NR1                                                NaN           auto   
5     D3NR2                                                NaN           auto   
6  D3TOTAL1                                                NaN           auto   
7  D3TOTAL2                                                NaN           auto   

   include  
0    False  
1    False  
2    False  
3    False  
4    False  
5    False  
6    False  
7    False  
```

### Text output (text/plain)

```text
     sample fraction biological_replicate
0     D3PR1                              
1     D3PR2                              
2      D0R1                              
3      D0R2                              
4     D3NR1                              
5     D3NR2                              
6  D3TOTAL1                              
7  D3TOTAL2                              
```

### Text output (error)

```text
[0;31m---------------------------------------------------------------------------[0m
[0;31mRuntimeError[0m                              Traceback (most recent call last)
Cell [0;32mIn[3], line 203[0m
[1;32m    190[0m     counts [38;5;241m=[39m (
[1;32m    191[0m         sheet[38;5;241m.[39mgroupby([[38;5;124m"[39m[38;5;124mbiological_replicate[39m[38;5;124m"[39m, [38;5;124m"[39m[38;5;124mfraction[39m[38;5;124m"[39m])
[1;32m    192[0m         [38;5;241m.[39msize()
[0;32m   (...)[0m
[1;32m    198[0m         )
[1;32m    199[0m     )
[1;32m    200[0m     [38;5;28;01mreturn[39;00m sheet, counts, used_order_fallback
[0;32m--> 203[0m sample_sheet, sample_group_counts, used_order_fallback [38;5;241m=[39m [43mbuild_sample_sheet[49m[43m([49m
[1;32m    204[0m [43m    [49m[43mSAMPLE_COLUMNS[49m
[1;32m    205[0m [43m)[49m
[1;32m    206[0m sample_sheet[38;5;241m.[39mto_csv(cfg[38;5;241m.[39mOUTPUT_DIR [38;5;241m/[39m [38;5;124m"[39m[38;5;124msample_sheet_used.csv[39m[38;5;124m"[39m, index[38;5;241m=[39m[38;5;28;01mFalse[39;00m)
[1;32m    208[0m display(sample_sheet)

Cell [0;32mIn[3], line 183[0m, in [0;36mbuild_sample_sheet[0;34m(sample_columns)[0m
[1;32m    176[0m     display(sheet)
[1;32m    177[0m     display(
[1;32m    178[0m         sheet[38;5;241m.[39mloc[
[1;32m    179[0m             invalid_fraction [38;5;241m|[39m missing_biorep,
[1;32m    180[0m             [[38;5;124m"[39m[38;5;124msample[39m[38;5;124m"[39m, [38;5;124m"[39m[38;5;124mfraction[39m[38;5;124m"[39m, [38;5;124m"[39m[38;5;124mbiological_replicate[39m[38;5;124m"[39m],
[1;32m    181[0m         ]
[1;32m    182[0m     )
[0;32m--> 183[0m     [38;5;28;01mraise[39;00m [38;5;167;01mRuntimeError[39;00m(
[1;32m    184[0m         [38;5;124m"[39m[38;5;124mSome sample columns remain unresolved. Add exact entries to [39m[38;5;124m"[39m
[1;32m    185[0m         [38;5;124m"[39m[38;5;124mcfg.MANUAL_SAMPLE_MAP or edit the generated sheet and set [39m[38;5;124m"[39m
[1;32m    186[0m         [38;5;124m"[39m[38;5;124mcfg.SAMPLE_SHEET_PATH. Template:[39m[38;5;130;01m\n[39;00m[38;5;124m"[39m
[1;32m    187[0m         [38;5;124mf[39m[38;5;124m"[39m[38;5;132;01m{[39;00medit_path[38;5;132;01m}[39;00m[38;5;124m"[39m
[1;32m    188[0m     )
[1;32m    190[0m counts [38;5;241m=[39m (
[1;32m    191[0m     sheet[38;5;241m.[39mgroupby([[38;5;124m"[39m[38;5;124mbiological_replicate[39m[38;5;124m"[39m, [38;5;124m"[39m[38;5;124mfraction[39m[38;5;124m"[39m])
[1;32m    192[0m     [38;5;241m.[39msize()
[0;32m   (...)[0m
[1;32m    198[0m     )
[1;32m    199[0m )
[1;32m    200[0m [38;5;28;01mreturn[39;00m sheet, counts, used_order_fallback

[0;31mRuntimeError[0m: Some sample columns remain unresolved. Add exact entries to cfg.MANUAL_SAMPLE_MAP or edit the generated sheet and set cfg.SAMPLE_SHEET_PATH. Template:
/oak/stanford/groups/engreitz/Users/edtng/GenomeWideValidation/guide_level_EC_efficiency_edgeR_FIXED/sample_sheet_to_edit.csv
```

## Markdown cell 6

## 3. Sum PCR technical replicates

PCR products from the same biological replicate and fraction are technical replicates, so they are summed. D0 libraries are retained in the sample sheet for provenance but are not included in the POS-versus-NEG edgeR model.

## Code cell 7

```python
annotation_cols = [GUIDE_COL] + ([GENE_COL] if GENE_COL is not None else [])
aggregated = raw_table[annotation_cols].copy()

for (biorep, fraction), rows in d3_sheet.groupby(
    ["biological_replicate", "fraction"]
):
    columns = rows["sample"].tolist()
    aggregated[f"{biorep}_{fraction}"] = (
        raw_counts[columns].sum(axis=1).astype(np.int64)
    )

aggregated_sample_columns = [
    f"{rep}_{fraction}"
    for rep in ["B1", "B2"]
    for fraction in ["NEG", "POS"]
]
missing = [c for c in aggregated_sample_columns if c not in aggregated.columns]
if missing:
    raise RuntimeError(
        "Missing aggregated D3 count columns: " + ", ".join(missing)
    )

edgeR_counts = (
    aggregated[[GUIDE_COL] + aggregated_sample_columns]
    .set_index(GUIDE_COL)
)
edgeR_sample_metadata = pd.DataFrame(
    {
        "sample": aggregated_sample_columns,
        "biological_replicate": [x.split("_")[0] for x in aggregated_sample_columns],
        "fraction": [x.split("_")[1] for x in aggregated_sample_columns],
    }
)
library_qc = pd.DataFrame(
    {
        "sample": aggregated_sample_columns,
        "library_size": edgeR_counts.sum(axis=0).to_numpy(),
        "nonzero_guides": (edgeR_counts > 0).sum(axis=0).to_numpy(),
    }
)

display(edgeR_sample_metadata)
display(library_qc)
display(aggregated.head())

aggregated.to_csv(
    cfg.OUTPUT_DIR / "technical_replicates_summed_counts.csv", index=False
)
edgeR_counts.to_csv(cfg.OUTPUT_DIR / "edgeR_guide_count_matrix.csv")
edgeR_sample_metadata.to_csv(
    cfg.OUTPUT_DIR / "edgeR_sample_metadata.csv", index=False
)
library_qc.to_csv(cfg.OUTPUT_DIR / "edgeR_library_QC.csv", index=False)

```

## Markdown cell 8

## 4. Identify negative-control guides

Safe-targeting and non-targeting guides are used to center guide-level effects.

## Code cell 9

```python
guide_text = aggregated[GUIDE_COL].astype(str)
if GENE_COL is not None:
    guide_text = guide_text + " " + aggregated[GENE_COL].astype(str)

aggregated["is_negative_control"] = guide_text.str.contains(
    cfg.CONTROL_REGEX, regex=True, na=False
)

print("Negative-control guides detected:", int(aggregated["is_negative_control"].sum()))
print("Targeting/other guides:", int((~aggregated["is_negative_control"]).sum()))
display(
    aggregated.loc[
        aggregated["is_negative_control"],
        annotation_cols + ["is_negative_control"],
    ].head(20)
)

if aggregated["is_negative_control"].sum() < 20:
    raise RuntimeError(
        "Too few negative-control guides were detected. Inspect the guide/gene "
        "labels and update cfg.CONTROL_REGEX."
    )

```

## Markdown cell 10

## 5. Run edgeR

The quasi-likelihood model tests `fractionPOS` while accounting for biological replicate. Positive log2 fold-change indicates enrichment in the CDH5-positive fraction; negative values indicate depletion from the EC fraction.

## Code cell 11

```python
rscript = shutil.which(cfg.R_SCRIPT_BIN)
if rscript is None:
    raise RuntimeError(
        "Rscript was not found. Load an R environment containing edgeR on "
        "Sherlock, or set cfg.R_SCRIPT_BIN to the full Rscript path."
    )

r_code = r'''
suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
counts_path <- args[[1]]
metadata_path <- args[[2]]
output_path <- args[[3]]
min_total_count <- as.numeric(args[[4]])

counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample
metadata <- metadata[colnames(counts), , drop = FALSE]

metadata$biological_replicate <- factor(metadata$biological_replicate)
metadata$fraction <- factor(metadata$fraction, levels = c("NEG", "POS"))

design <- model.matrix(~ biological_replicate + fraction, data = metadata)
if (!("fractionPOS" %in% colnames(design))) {
    stop("fractionPOS is absent from the design matrix.")
}

y <- DGEList(counts = counts)
keep <- rowSums(counts) >= min_total_count & filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
test <- glmQLFTest(fit, coef = "fractionPOS")
result <- topTags(test, n = Inf, sort.by = "none")$table
result$guide_id <- rownames(result)
result$tested_by_edgeR <- TRUE

all_guides <- data.frame(guide_id = rownames(counts), stringsAsFactors = FALSE)
output <- merge(all_guides, result, by = "guide_id", all.x = TRUE, sort = FALSE)
output$tested_by_edgeR[is.na(output$tested_by_edgeR)] <- FALSE
write.csv(output, output_path, row.names = FALSE)
'''

r_script_path = cfg.OUTPUT_DIR / "run_edgeR_guide_level.R"
r_script_path.write_text(r_code)
edgeR_result_path = cfg.OUTPUT_DIR / "edgeR_guide_level_POS_vs_NEG.csv"

command = [
    rscript,
    str(r_script_path),
    str(cfg.OUTPUT_DIR / "edgeR_guide_count_matrix.csv"),
    str(cfg.OUTPUT_DIR / "edgeR_sample_metadata.csv"),
    str(edgeR_result_path),
    str(cfg.MIN_TOTAL_COUNT),
]
print("Running:", " ".join(command))
completed = subprocess.run(command, capture_output=True, text=True)
print(completed.stdout)
if completed.returncode != 0:
    print(completed.stderr)
    raise RuntimeError("edgeR failed; see the R error above.")

edgeR_results = pd.read_csv(edgeR_result_path)
print("Guides tested by edgeR:", int(edgeR_results["tested_by_edgeR"].sum()))
display(edgeR_results.sort_values("FDR").head(20))

```

## Markdown cell 12

## 6. Center guide effects on negative controls

The primary output is a control-centered relative CDH5-positive odds index:

\[
100\times 2^{\mathrm{logFC}_g-\mathrm{median}(\mathrm{logFC}_{NTC})}.
\]

A value of 100% denotes control-like odds. This is not the directly observed fraction of EC cells. For interpretation only, the notebook also supplies approximate EC percentages calibrated to the whole-population FACS efficiencies.

## Code cell 13

```python
guide_annotation = aggregated[
    annotation_cols + ["is_negative_control"]
].rename(columns={GUIDE_COL: "guide_id"})

results = guide_annotation.merge(edgeR_results, on="guide_id", how="left")
for column in ["logFC", "PValue", "FDR"]:
    results[column] = pd.to_numeric(results[column], errors="coerce")

control_tested = (
    results["is_negative_control"]
    & results["tested_by_edgeR"].fillna(False)
    & results["logFC"].notna()
)
if control_tested.sum() < 20:
    raise RuntimeError(
        "Too few tested control guides are available to center the effects."
    )

control_median_logFC = float(results.loc[control_tested, "logFC"].median())
results["log2_relative_CDH5_pos_odds"] = (
    results["logFC"] - control_median_logFC
)
results["relative_CDH5_pos_odds_percent"] = (
    100 * np.power(2.0, results["log2_relative_CDH5_pos_odds"])
)


def odds_shift_to_probability(log2_relative_odds, baseline_fraction):
    baseline_odds = baseline_fraction / (1 - baseline_fraction)
    guide_odds = baseline_odds * np.power(2.0, log2_relative_odds)
    return guide_odds / (1 + guide_odds)


for rep, baseline in cfg.FACS_EC_FRACTION_BY_REPLICATE.items():
    results[f"estimated_EC_percent_calibrated_to_{rep}"] = (
        100
        * odds_shift_to_probability(
            results["log2_relative_CDH5_pos_odds"], baseline
        )
    )

calibration_columns = [
    f"estimated_EC_percent_calibrated_to_{rep}"
    for rep in cfg.FACS_EC_FRACTION_BY_REPLICATE
]
results["estimated_EC_percent_mean_FACS_calibration"] = results[
    calibration_columns
].mean(axis=1)

results["total_raw_count"] = (
    edgeR_counts.sum(axis=1).reindex(results["guide_id"]).to_numpy()
)

print(f"Median control-guide edgeR logFC: {control_median_logFC:.4f}")
display(results.sort_values("FDR").head(25))
results.to_csv(
    cfg.OUTPUT_DIR / "guide_level_EC_efficiency_edgeR_results.csv",
    index=False,
)

```

## Markdown cell 14

## 7. Summarize guide concordance by target

This guide-level summary is complementary to, not a replacement for, the manuscript’s original target-level Mann–Whitney analysis.

## Code cell 15

```python
if GENE_COL is None:
    print("No gene column was detected; target-level summary was skipped.")
    target_summary = pd.DataFrame()
else:
    targeting = results.loc[~results["is_negative_control"]].copy()
    targeting[GENE_COL] = targeting[GENE_COL].astype(str)

    target_summary = (
        targeting.groupby(GENE_COL, as_index=False)
        .agg(
            n_guides=("guide_id", "size"),
            n_guides_tested=("tested_by_edgeR", "sum"),
            mean_edgeR_log2FC=("logFC", "mean"),
            median_edgeR_log2FC=("logFC", "median"),
            mean_log2_relative_CDH5_pos_odds=(
                "log2_relative_CDH5_pos_odds", "mean"
            ),
            median_relative_CDH5_pos_odds_percent=(
                "relative_CDH5_pos_odds_percent", "median"
            ),
            mean_estimated_EC_percent=(
                "estimated_EC_percent_mean_FACS_calibration", "mean"
            ),
            minimum_guide_FDR=("FDR", "min"),
            n_guides_FDR_lt_0_05=(
                "FDR",
                lambda x: int((pd.to_numeric(x, errors="coerce") < 0.05).sum()),
            ),
            n_guides_depleted_from_POS=(
                "log2_relative_CDH5_pos_odds",
                lambda x: int((pd.to_numeric(x, errors="coerce") < 0).sum()),
            ),
            n_guides_enriched_in_POS=(
                "log2_relative_CDH5_pos_odds",
                lambda x: int((pd.to_numeric(x, errors="coerce") > 0).sum()),
            ),
        )
    )
    target_summary["guide_direction_concordance"] = (
        target_summary[
            ["n_guides_depleted_from_POS", "n_guides_enriched_in_POS"]
        ].max(axis=1)
        / target_summary["n_guides_tested"].replace(0, np.nan)
    )
    target_summary.to_csv(
        cfg.OUTPUT_DIR / "target_level_summary_from_guide_edgeR.csv",
        index=False,
    )
    display(
        target_summary.sort_values(
            ["n_guides_FDR_lt_0_05", "minimum_guide_FDR"],
            ascending=[False, True],
        ).head(30)
    )

```

## Markdown cell 16

## 8. Diagnostic plots

## Code cell 17

```python
import matplotlib.pyplot as plt

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

fig, ax = plt.subplots(figsize=(5.8, 3.8))
ax.bar(library_qc["sample"], library_qc["library_size"])
ax.set_ylabel("Total guide reads")
ax.set_title(
    "Library sizes after PCR-replicate aggregation", fontweight="bold"
)
ax.tick_params(axis="x", rotation=35)
fig.tight_layout()
for extension in ["png", "pdf", "svg"]:
    fig.savefig(
        cfg.OUTPUT_DIR / f"library_sizes_aggregated.{extension}",
        dpi=300 if extension == "png" else None,
        bbox_inches="tight",
    )
plt.show()

plot_df = results.loc[
    results["tested_by_edgeR"].fillna(False) & results["FDR"].notna()
].copy()
plot_df["neglog10FDR"] = -np.log10(plot_df["FDR"].clip(lower=1e-300))
plot_df["significant"] = plot_df["FDR"] < 0.05

fig, ax = plt.subplots(figsize=(6.2, 5.0))
for significant, size, alpha in [(False, 10, 0.35), (True, 18, 0.85)]:
    subset = plot_df.loc[plot_df["significant"].eq(significant)]
    ax.scatter(
        subset["log2_relative_CDH5_pos_odds"],
        subset["neglog10FDR"],
        s=size,
        alpha=alpha,
        label="FDR < 0.05" if significant else "Not significant",
    )
ax.axvline(0, linestyle="--", linewidth=0.8, color="gray")
ax.axhline(-np.log10(0.05), linestyle="--", linewidth=0.8, color="gray")
ax.set_xlabel("Control-centered log2 CDH5-positive odds")
ax.set_ylabel("-log10(FDR)")
ax.set_title("Guide-level POS-versus-NEG effects", fontweight="bold")
ax.legend(frameon=False)
fig.tight_layout()
for extension in ["png", "pdf", "svg"]:
    fig.savefig(
        cfg.OUTPUT_DIR / f"guide_level_edgeR_volcano.{extension}",
        dpi=300 if extension == "png" else None,
        bbox_inches="tight",
    )
plt.show()

fig, ax = plt.subplots(figsize=(6.2, 4.0))
ax.hist(
    results.loc[
        results["is_negative_control"]
        & results["log2_relative_CDH5_pos_odds"].notna(),
        "log2_relative_CDH5_pos_odds",
    ],
    bins=30,
    alpha=0.65,
    label="Negative controls",
)
ax.hist(
    results.loc[
        ~results["is_negative_control"]
        & results["log2_relative_CDH5_pos_odds"].notna(),
        "log2_relative_CDH5_pos_odds",
    ],
    bins=60,
    alpha=0.45,
    label="Targeting guides",
)
ax.axvline(0, linestyle="--", linewidth=0.8, color="black")
ax.set_xlabel("Control-centered log2 CDH5-positive odds")
ax.set_ylabel("Guides")
ax.set_title("Guide-effect distribution", fontweight="bold")
ax.legend(frameon=False)
fig.tight_layout()
for extension in ["png", "pdf", "svg"]:
    fig.savefig(
        cfg.OUTPUT_DIR / f"guide_effect_distribution.{extension}",
        dpi=300 if extension == "png" else None,
        bbox_inches="tight",
    )
plt.show()

```

## Markdown cell 18

## 9. Output manifest

## Code cell 19

```python
manifest = []
for path in sorted(cfg.OUTPUT_DIR.glob("*")):
    if path.is_file():
        manifest.append({"file": path.name, "bytes": path.stat().st_size})
manifest_df = pd.DataFrame(manifest)
manifest_df.to_csv(cfg.OUTPUT_DIR / "output_manifest.csv", index=False)
display(manifest_df)

analysis_summary = {
    "raw_counts_path": str(RAW_COUNTS_PATH),
    "n_guides": int(len(raw_table)),
    "n_raw_PCR_library_columns": int(len(SAMPLE_COLUMNS)),
    "documented_order_fallback_used": bool(used_order_fallback),
    "aggregated_samples": aggregated_sample_columns,
    "n_negative_control_guides": int(results["is_negative_control"].sum()),
    "median_control_logFC": float(control_median_logFC),
    "FACS_EC_fraction_B1": cfg.FACS_EC_FRACTION_BY_REPLICATE["B1"],
    "FACS_EC_fraction_B2": cfg.FACS_EC_FRACTION_BY_REPLICATE["B2"],
}
(cfg.OUTPUT_DIR / "analysis_summary.json").write_text(
    json.dumps(analysis_summary, indent=2)
)
print(json.dumps(analysis_summary, indent=2))

```

