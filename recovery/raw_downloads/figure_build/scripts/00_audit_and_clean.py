from __future__ import annotations

from pathlib import Path
import json
import math
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source_data"
CLEAN = ROOT / "cleaned_data"
AUDIT = ROOT / "data_audit"

SCREEN_SOURCE = SOURCE / "BatchA_screen_hits_244_derived_source.csv"
ALIAS_SOURCE = SOURCE / "gene_aliases_applied_from_guide_metadata.csv"

CLEAN.mkdir(parents=True, exist_ok=True)
AUDIT.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Load the locally available Batch A hit-level source.
# -----------------------------
df = pd.read_csv(SCREEN_SOURCE)
required = {
    "gene_mapped", "screen_log2FC", "screen_p", "screen_fdr",
    "original_gene", "original_target", "n_screen_rows"
}
missing = sorted(required - set(df.columns))
if missing:
    raise ValueError(f"Missing required columns in {SCREEN_SOURCE.name}: {missing}")

# Clean gene-level table while preserving source identifiers.
clean = pd.DataFrame({
    "gene": df["gene_mapped"].astype(str),
    "original_gene": df["original_gene"].astype(str),
    "original_target_ids": df["original_target"].astype(str),
    "n_source_target_records": pd.to_numeric(df["n_screen_rows"], errors="coerce").fillna(1).astype(int),
    "u_test_effect_log2_CDH5pos_vs_CDH5neg": pd.to_numeric(df["screen_log2FC"], errors="coerce"),
    "u_test_p_value": pd.to_numeric(df["screen_p"], errors="coerce"),
    "u_test_fdr": pd.to_numeric(df["screen_fdr"], errors="coerce"),
})
clean["u_test_hit_fdr_lt_0_05"] = clean["u_test_fdr"] < 0.05
clean["final_union_hit"] = True
clean["screen_direction"] = np.where(
    clean["u_test_effect_log2_CDH5pos_vs_CDH5neg"] < 0,
    "Required for EC differentiation",
    "Restricts EC differentiation",
)
clean["direction_basis"] = (
    "sign of available U-test log2 enrichment in CDH5+ versus CDH5-; "
    "final union-table direction requires confirmation"
)
clean["display_gene"] = clean["gene"]
clean["gene_upper"] = clean["gene"].str.upper()

# Deterministic ranking for the limited hit-only display.
clean["rank_by_u_test_effect_ascending"] = (
    clean["u_test_effect_log2_CDH5pos_vs_CDH5neg"].rank(method="min", ascending=True).astype(int)
)
clean["minus_log10_u_test_fdr"] = -np.log10(clean["u_test_fdr"].clip(lower=1e-300))

# Manuscript-defined known EC regulators, because no author-approved full list was found.
known = {"CDH5", "KDR", "ETV2"}
clean["limited_manuscript_known_regulator"] = clean["gene"].isin(known)

clean = clean.sort_values("gene").reset_index(drop=True)
clean.to_csv(CLEAN / "gw_screen_hits_244_gene_level.tsv", sep="\t", index=False)

# Target-level audit: preserve both NEDD8 target identifiers rather than losing the duplicate.
target_rows = []
for _, row in clean.iterrows():
    target_ids = [x for x in str(row["original_target_ids"]).split("|") if x]
    if not target_ids:
        target_ids = [""]
    for tid in target_ids:
        target_rows.append({
            "gene": row["gene"],
            "original_target_id": tid,
            "source_original_gene": row["original_gene"],
            "u_test_effect_log2_CDH5pos_vs_CDH5neg": row["u_test_effect_log2_CDH5pos_vs_CDH5neg"],
            "u_test_p_value": row["u_test_p_value"],
            "u_test_fdr": row["u_test_fdr"],
            "screen_direction": row["screen_direction"],
            "final_union_hit": True,
        })
target_audit = pd.DataFrame(target_rows)
target_audit.to_csv(CLEAN / "gw_screen_hits_target_level_audit.tsv", sep="\t", index=False)

# Alias mapping used in the prior cleaning step.
alias = pd.read_csv(ALIAS_SOURCE)
alias.to_csv(CLEAN / "gene_alias_mapping.tsv", sep="\t", index=False)

# -----------------------------
# Author-review templates
# -----------------------------
annotation_template = pd.DataFrame({
    "gene": clean["gene"],
    "original_target_id": clean["original_target_ids"],
    "screen_direction": clean["screen_direction"],
    "primary_category": "",
    "secondary_category": "",
    "known_EC_regulator": np.where(clean["gene"].isin(known), "TRUE", ""),
    "predefined_positive_control": np.where(clean["gene"].isin(known), "REVIEW", ""),
    "broad_essential_or_viability": "",
    "must_label": np.where(clean["gene"].isin(known), "TRUE", ""),
    "label_priority": np.where(clean["gene"].isin(known), "1", ""),
    "annotation_status": np.where(clean["gene"].isin(known), "limited manuscript-defined", "needs author review"),
    "notes": np.where(
        clean["gene"].isin(known),
        "Known EC regulator explicitly named in the attached manuscript/brief; verify final approved list.",
        "",
    ),
})
annotation_template.to_csv(CLEAN / "screen_regulators_244_annotation_template.tsv", sep="\t", index=False)

known_template = pd.DataFrame({
    "gene": clean["gene"],
    "known_regulator_status": np.where(clean["gene"].isin(known), "manuscript-defined known EC regulator", ""),
    "must_label": np.where(clean["gene"].isin(known), "TRUE", ""),
    "reason_for_label": np.where(
        clean["gene"].isin(known),
        "Explicitly identified as a recovered known EC regulator in the manuscript/figure redesign brief.",
        "",
    ),
    "label_priority": np.where(clean["gene"].isin(known), "1", ""),
    "author_decision": "",
    "author_notes": "",
})
known_template.to_csv(CLEAN / "known_regulator_review_template.tsv", sep="\t", index=False)

# -----------------------------
# Data audit and file-to-panel map
# -----------------------------
# File-library attachments are included from the attachment audit even when not locally mounted.
audit_rows = [
    {
        "file_name": "CC-Perturb-seq EC Differentiation - Text V2 (1).pdf",
        "file_type": "PDF manuscript",
        "internal_sheets_files": "manuscript text, current figures and legends",
        "rows": "NA",
        "columns": "NA",
        "unit_of_observation": "manuscript pages",
        "important_columns": "screen workflow timing; phenotype definition; current figure legend",
        "candidate_figure_panel": "Fig1a; design/context for Fig1b-e",
        "data_issues": "Parsed through attachment index; raw PDF not locally mounted in build runtime.",
        "sufficient_for_panel": "Yes for conceptual workflow; no for quantitative per-target plotting",
    },
    {
        "file_name": "Figure Redesign Decisions and Precedent Brief.pdf",
        "file_type": "PDF design brief",
        "internal_sheets_files": "executive summary, design decisions, precedent table, layout options",
        "rows": "NA",
        "columns": "NA",
        "unit_of_observation": "document pages",
        "important_columns": "244-regulator count; balanced Figure 1 layout; panel roles",
        "candidate_figure_panel": "All Figure 1 panels and assembly",
        "data_issues": "None for design authority.",
        "sufficient_for_panel": "Yes for design/narrative",
    },
    {
        "file_name": "gw_245_hits.txt",
        "file_type": "tab-delimited text in attachment library",
        "internal_sheets_files": "single table",
        "rows": 245,
        "columns": 5,
        "unit_of_observation": "screen target identifier",
        "important_columns": "Target; u_test_hit; mageck_hit; pos_neg; Gene",
        "candidate_figure_panel": "Fig1b; Fig1c; templates for Fig1d",
        "data_issues": "NEDD8 occurs under two target identifiers; attachment not locally mounted for exact row-level merge.",
        "sufficient_for_panel": "Partially; confirms 245 target rows/244 genes and hit-call schema",
    },
    {
        "file_name": SCREEN_SOURCE.name,
        "file_type": "CSV derived from prior Batch A analysis",
        "internal_sheets_files": "single table",
        "rows": len(df),
        "columns": len(df.columns),
        "unit_of_observation": "unique biological gene",
        "important_columns": "gene_mapped; original_target; screen_log2FC; screen_p; screen_fdr",
        "candidate_figure_panel": "temporary Fig1b; Fig1c; Fig1d review template",
        "data_issues": "Contains 244 final hits only, not the full tested genome-wide universe; MaGECK per-target statistics absent.",
        "sufficient_for_panel": "Yes for hit-only rank and known-regulator recovery; no for final Fig1b volcano",
    },
    {
        "file_name": ALIAS_SOURCE.name,
        "file_type": "CSV",
        "internal_sheets_files": "single mapping table",
        "rows": len(alias),
        "columns": len(alias.columns),
        "unit_of_observation": "gene alias mapping",
        "important_columns": "; ".join(alias.columns),
        "candidate_figure_panel": "identifier harmonization",
        "data_issues": "RBM14 resolved previously to RBM14-RBM4 from guide metadata; author review still recommended.",
        "sufficient_for_panel": "Yes for current cleaned source",
    },
    {
        "file_name": "Fig 2 - CC-Perturb-seq V6.pdf",
        "file_type": "PDF existing figure artwork in attachment library",
        "internal_sheets_files": "current screen schematic, hit plot, category summary, concordance, program atlas",
        "rows": "NA",
        "columns": "NA",
        "unit_of_observation": "figure panels",
        "important_columns": "current Figure 2 screen artwork and counts",
        "candidate_figure_panel": "visual reference for Fig1a-b; categories require final annotation provenance",
        "data_issues": "Existing category panel appears to mix screen and expanded-library bookkeeping; not used as final source.",
        "sufficient_for_panel": "Reference only",
    },
    {
        "file_name": "Pos_vs_Neg.zip",
        "file_type": "ZIP",
        "internal_sheets_files": "not found in accessible Batch A attachments",
        "rows": "NA",
        "columns": "NA",
        "unit_of_observation": "unknown",
        "important_columns": "expected screen direction support",
        "candidate_figure_panel": "Fig1b-d",
        "data_issues": "Missing/not accessible.",
        "sufficient_for_panel": "No",
    },
    {
        "file_name": "validation screen results.zip",
        "file_type": "ZIP",
        "internal_sheets_files": "not found in accessible Batch A attachments",
        "rows": "NA",
        "columns": "NA",
        "unit_of_observation": "expected target/gene validation effects and replicates",
        "important_columns": "genome-wide effect; validation effect; replicate 1; replicate 2; p/q values",
        "candidate_figure_panel": "Fig1e",
        "data_issues": "Missing/not accessible.",
        "sufficient_for_panel": "No",
    },
    {
        "file_name": "final_244_regulator_annotation.tsv",
        "file_type": "expected table",
        "internal_sheets_files": "not found",
        "rows": "expected 244",
        "columns": "expected category/known-control fields",
        "unit_of_observation": "unique biological gene",
        "important_columns": "primary_category; secondary_category; known_EC_regulator; must_label",
        "candidate_figure_panel": "Fig1d; author-approved labels for Fig1b-c",
        "data_issues": "Missing; review template generated.",
        "sufficient_for_panel": "No",
    },
    {
        "file_name": "author_approved_known_regulator_list.tsv",
        "file_type": "expected table",
        "internal_sheets_files": "not found",
        "rows": "unknown",
        "columns": "expected gene/status/must_label/reason/priority",
        "unit_of_observation": "gene",
        "important_columns": "known_regulator_status; must_label; label_priority",
        "candidate_figure_panel": "Fig1b-c",
        "data_issues": "Missing; limited manuscript-defined set used and review template generated.",
        "sufficient_for_panel": "No for final labels; partial for CDH5/KDR/ETV2",
    },
]

audit = pd.DataFrame(audit_rows)
audit.to_csv(AUDIT / "data_audit.tsv", sep="\t", index=False)

file_panel = audit[["file_name", "candidate_figure_panel", "sufficient_for_panel", "data_issues"]].copy()
file_panel.to_csv(AUDIT / "file_to_panel_map.tsv", sep="\t", index=False)

# Transformation log
transforms = pd.DataFrame([
    {
        "step": 1,
        "input": SCREEN_SOURCE.name,
        "output": "gw_screen_hits_244_gene_level.tsv",
        "operation": "Selected screen columns; preserved original gene/target IDs; derived U-test hit status and biological direction from U-test effect sign.",
        "caveat": "Available source contains final hits only and does not carry full MaGECK statistics.",
    },
    {
        "step": 2,
        "input": "original_target joined identifiers",
        "output": "gw_screen_hits_target_level_audit.tsv",
        "operation": "Split pipe-delimited target identifiers to preserve target-level source rows.",
        "caveat": "NEDD8 remains two target records in target-level audit but one gene in main count.",
    },
    {
        "step": 3,
        "input": ALIAS_SOURCE.name,
        "output": "gene_alias_mapping.tsv",
        "operation": "Copied previously documented guide-metadata alias mapping.",
        "caveat": "RBM14/RBM14-RBM4 relationship retained as documented and flagged for review.",
    },
    {
        "step": 4,
        "input": "gw_screen_hits_244_gene_level.tsv",
        "output": "screen_regulators_244_annotation_template.tsv",
        "operation": "Created one-row-per-gene author review template; populated only source-supported fields.",
        "caveat": "Functional categories intentionally blank.",
    },
    {
        "step": 5,
        "input": "gw_screen_hits_244_gene_level.tsv",
        "output": "known_regulator_review_template.tsv",
        "operation": "Created review template; pre-populated CDH5, KDR and ETV2 as manuscript-defined known regulators.",
        "caveat": "No author-approved comprehensive label list was attached.",
    },
])
transforms.to_csv(AUDIT / "transformation_log.tsv", sep="\t", index=False)

# Summary metrics for README and legends.
summary = {
    "n_gene_level_hits": int(clean["gene"].nunique()),
    "n_target_level_records_reconstructed": int(len(target_audit)),
    "duplicate_genes_target_level": sorted(target_audit.loc[target_audit.duplicated("gene", keep=False), "gene"].unique().tolist()),
    "n_u_test_fdr_lt_0_05": int(clean["u_test_hit_fdr_lt_0_05"].sum()),
    "n_negative_u_test_effect": int((clean["u_test_effect_log2_CDH5pos_vs_CDH5neg"] < 0).sum()),
    "n_positive_u_test_effect": int((clean["u_test_effect_log2_CDH5pos_vs_CDH5neg"] > 0).sum()),
    "known_regulators_available": sorted(clean.loc[clean["limited_manuscript_known_regulator"], "gene"].tolist()),
}
(ROOT / "data_audit" / "summary_metrics.json").write_text(json.dumps(summary, indent=2))

# Unresolved issues and consolidated request.
unresolved = f"""# Unresolved issues — Figure 1 build

## Quantitative screen universe

- The locally available quantitative source contains **{summary['n_gene_level_hits']} final screen-derived hit genes only**. It does not contain the full genome-wide tested set required for the final Figure 1b.
- A file-library attachment named `gw_245_hits.txt` confirms 245 target records, 244 unique genes and columns `Target`, `u_test_hit`, `mageck_hit`, `pos_neg`, and `Gene`, but the full row-level file was not locally mounted in this build runtime.
- Per-target MaGECK effect/rank statistics and MaGECK FDR were not accessible. Candidate B is therefore a placeholder.

## Direction-count discrepancy

- In the locally available U-test effect table, {summary['n_negative_u_test_effect']} genes have negative effects and {summary['n_positive_u_test_effect']} have positive effects.
- Existing draft artwork in the attachment index reports 238 significant decreases and 6 significant increases. These quantities may differ because of hit-calling method, target/gene unit, or final union direction rules. The build does **not** silently reconcile them.

## Tested-set count discrepancy

- The redesign brief describes approximately 20,528 protein-coding genes.
- Existing draft artwork reports 20,411 non-significant + 238 decreases + 6 increases = 20,655 plotted target records.
- This likely reflects a difference between unique genes and promoter/target identifiers; final Figure 1b must state the unit explicitly.

## Gene aliases

- NEDD8 is one biological gene with two source target identifiers and is counted once in main displays.
- The available alias table maps RBM14 to RBM14-RBM4. Because the instructions flag RBM14/RBM14-P2/RBM14-RBM4/RBM4 as potentially ambiguous, final author confirmation is still recommended.

## Missing Figure 1 inputs

1. Full genome-wide U-test table for all tested targets/genes.
2. Full MaGECK table with a coherent per-target effect/rank statistic and FDR.
3. Author-approved known-regulator/label list.
4. Completed 244-regulator functional annotation table.
5. Secondary validation-screen effects and replicate-level data.
"""
(ROOT / "unresolved_issues.md").write_text(unresolved)

missing_request = """# Consolidated missing-input request for final Figure 1

The current build is complete for every panel supported by the accessible Batch A files. To replace the temporary placeholders and finalize Figure 1, provide the following files if they are not already available in a non-expired attachment:

## 1. Full genome-wide U-test results
Required columns:

- `target_id`
- `gene`
- `u_test_effect` or `log2FC_CDH5pos_vs_CDH5neg`
- `u_test_p_value`
- `u_test_adjusted_p_value` or `u_test_fdr`
- optional `n_CDH5pos`, `n_CDH5neg`

The table must include the complete tested set, not only hits.

## 2. Full MaGECK results
Required columns:

- `target_id`
- `gene`
- coherent MaGECK effect/rank statistic
- `mageck_p_value`
- `mageck_fdr`
- final MaGECK hit flag

## 3. Completed 244-regulator annotation
Use the generated file `cleaned_data/screen_regulators_244_annotation_template.tsv` and complete at minimum:

- `primary_category`
- `secondary_category`
- `known_EC_regulator`
- `predefined_positive_control`
- `broad_essential_or_viability`
- `must_label`
- `label_priority`
- `notes`

## 4. Author-approved known-regulator labels
Use `cleaned_data/known_regulator_review_template.tsv`, or provide:

- `gene`
- `known_regulator_status`
- `must_label`
- `reason_for_label`
- `label_priority`

## 5. Secondary validation-screen results
Required columns:

- `gene` and/or `target_id`
- `genomewide_screen_effect`
- `validation_screen_effect`
- `validation_p_value`
- `validation_q_value`
- `replicate1_effect`
- `replicate2_effect`
- target class (`screen_hit`, `random_nonhit`, `additional_curated`, `control`)

Earlier ZIP attachments may have expired; re-upload only the files not already accessible in the current Batch A set.
"""
(ROOT / "data_audit" / "consolidated_missing_input_request.md").write_text(missing_request)

print(json.dumps(summary, indent=2))
