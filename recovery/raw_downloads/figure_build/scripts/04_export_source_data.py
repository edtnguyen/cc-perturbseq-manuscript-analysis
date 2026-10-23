from __future__ import annotations

from pathlib import Path
import json
import platform
import sys
import pandas as pd
import numpy as np
import matplotlib

ROOT = Path(__file__).resolve().parents[1]
CLEAN = ROOT / "cleaned_data"
SOURCE = ROOT / "source_data"
LEGENDS = ROOT / "legends"
AUDIT = ROOT / "data_audit"

screen = pd.read_csv(CLEAN / "gw_screen_hits_244_gene_level.tsv", sep="\t")
summary = json.loads((AUDIT / "summary_metrics.json").read_text())

# Fig1a source timeline.
fig1a = pd.DataFrame([
    {"order": 1, "time": "day -4", "event": "infect doxycycline-inducible CRISPRi hiPSCs with genome-wide lentiviral gRNA pool", "source_status": "manuscript-supported"},
    {"order": 2, "time": "day -2", "event": "induce dCas9-KRAB expression with doxycycline", "source_status": "manuscript-supported"},
    {"order": 3, "time": "D0", "event": "hiPSC / differentiation start", "source_status": "manuscript-supported"},
    {"order": 4, "time": "D1", "event": "primitive streak", "source_status": "manuscript-supported"},
    {"order": 5, "time": "D2", "event": "lateral mesoderm", "source_status": "manuscript-supported"},
    {"order": 6, "time": "D3", "event": "endothelial differentiation; MACS separation of CDH5+ and CDH5- populations", "source_status": "manuscript-supported"},
    {"order": 7, "time": "D3", "event": "gRNA sequencing and regulator identification", "source_status": "manuscript-supported"},
])
fig1a.to_csv(SOURCE / "Fig1a_source_data.tsv", sep="\t", index=False)

# Fig1b source data: exact plotted temporary rank panel.
fig1b = screen.sort_values("u_test_effect_log2_CDH5pos_vs_CDH5neg").reset_index(drop=True).copy()
fig1b["display_rank"] = np.arange(1, len(fig1b) + 1)
fig1b["display_label"] = np.where(fig1b["gene"].isin(["CDH5", "KDR", "ETV2"]), fig1b["gene"], "")
fig1b["panel_status"] = "temporary hit-only representation; full tested set missing"
fig1b.to_csv(SOURCE / "Fig1b_source_data.tsv", sep="\t", index=False)

# Fig1c source data.
fig1c = screen[screen["gene"].isin(["CDH5", "KDR", "ETV2"])].copy()
fig1c["known_regulator_source"] = "limited manuscript-defined set"
fig1c["author_approval_status"] = "pending"
fig1c.to_csv(SOURCE / "Fig1c_source_data.tsv", sep="\t", index=False)

# Fig1d placeholder source metadata; actual template is in cleaned_data.
fig1d = pd.DataFrame([
    {
        "panel_status": "temporary placeholder",
        "n_unique_screen_regulators": summary["n_gene_level_hits"],
        "required_input": "completed screen_regulators_244_annotation_template.tsv",
        "preferred_plot": "horizontal diverging stacked bars by functional category and biological direction",
        "categories_plotted": "none; background bars are schematic only",
    }
])
fig1d.to_csv(SOURCE / "Fig1d_source_data.tsv", sep="\t", index=False)

# Fig1e placeholder source metadata.
fig1e = pd.DataFrame([
    {
        "panel_status": "temporary placeholder",
        "required_columns": "gene/target_id; genomewide_screen_effect; validation_screen_effect; validation_p_value; validation_q_value; replicate1_effect; replicate2_effect; target_class",
        "preferred_main_plot": "genome-wide effect versus secondary-validation effect",
        "preferred_inset": "replicate 1 versus replicate 2 effect",
    }
])
fig1e.to_csv(SOURCE / "Fig1e_source_data.tsv", sep="\t", index=False)

# Candidate comparison.
candidate_comparison = pd.DataFrame([
    {
        "candidate": "A",
        "representation": "U-test volcano",
        "statistical_coherence": "coherent U-test effect and U-test FDR",
        "available_input": "only 244 final hit genes",
        "limitation": "does not show the complete genome-wide tested set",
        "main_text_decision": "not selected as final; exported for review",
    },
    {
        "candidate": "B",
        "representation": "MaGECK result plot",
        "statistical_coherence": "would use MaGECK effect/rank and MaGECK FDR",
        "available_input": "missing",
        "limitation": "per-target MaGECK statistics not accessible",
        "main_text_decision": "placeholder only",
    },
    {
        "candidate": "C",
        "representation": "ranked U-test effect",
        "statistical_coherence": "single coherent U-test effect quantity",
        "available_input": "244 final hit genes",
        "limitation": "hit-only; cannot display non-hits or final hit-calling methods",
        "main_text_decision": "selected as temporary Figure 1b because it is readable and avoids mixed-method statistics",
    },
])
candidate_comparison.to_csv(SOURCE / "Fig1b_candidate_comparison.tsv", sep="\t", index=False)

legend = f"""# Draft Figure 1 legend

**Figure 1 | A genome-wide CRISPRi screen identifies regulators of endothelial-cell differentiation.**

**a, Genome-wide CRISPRi screen design.** Doxycycline-inducible CRISPRi hiPSCs were infected with a genome-wide lentiviral gRNA pool four days before differentiation. dCas9–KRAB expression was induced with doxycycline two days before differentiation and maintained throughout differentiation. Cells progressed through hiPSC (D0), primitive streak (D1), lateral mesoderm (D2) and endothelial differentiation (D3). D3 cells were separated by MACS into CDH5-positive EC and CDH5-negative non-EC populations, followed by gRNA sequencing and regulator identification.

**b, Temporary hit-only ranking of screen-derived regulators.** The {summary['n_gene_level_hits']} unique screen-derived hit genes are ranked by the available Mann–Whitney U-test effect, defined as log2 gRNA enrichment in CDH5-positive versus CDH5-negative cells. Each point represents one unique biological gene. Blue denotes a negative available U-test effect, interpreted as knockdown shifting toward the CDH5-negative/non-EC population and therefore a gene required for EC differentiation. Vermillion denotes a positive available U-test effect, interpreted as knockdown shifting toward the CDH5-positive/EC population and therefore a gene that restricts EC differentiation. KDR, CDH5 and ETV2 are labeled as the limited manuscript-defined known-regulator set. This panel is explicitly temporary because the available quantitative source contains only the 244 final hits; the full tested-target universe and per-target MaGECK results are required for the final statistically coherent genome-wide panel. In the available U-test table, {summary['n_negative_u_test_effect']} genes have negative effects and {summary['n_positive_u_test_effect']} have positive effects; this direction count remains to be reconciled against the final union hit table.

**c, Recovery of known EC regulators.** U-test effects for KDR, CDH5 and ETV2, the known EC regulators explicitly identified in the attached manuscript and design brief, are shown together with their rank among the 244 final hit genes. The set is deliberately limited pending author approval of the complete known-regulator and labeling list.

**d, Functional organization of the screen-derived regulators.** Temporary placeholder for an author-reviewed functional-category summary of the 244 unique screen-derived regulators. The faint diverging bars are schematic and do not encode data. A one-row-per-gene review file, `screen_regulators_244_annotation_template.tsv`, was generated; unsupported category assignments were left blank.

**e, Secondary validation and reproducibility.** Temporary placeholder for comparison of genome-wide and secondary-validation screen effects, with a biological-replicate concordance inset. Validation and replicate-level result tables were not present in the accessible Batch A files.

The main regulator count is at the unique biological-gene level. NEDD8 is represented by two original target identifiers in the target-level audit but is counted once in the 244-gene display.
"""
(LEGENDS / "Figure1_legend_draft.md").write_text(legend)

software = {
    "python": sys.version.split()[0],
    "platform": platform.platform(),
    "pandas": pd.__version__,
    "numpy": np.__version__,
    "matplotlib": matplotlib.__version__,
}

readme = f"""# Figure build README — Figure 1 Batch A

## Build status

Figure 1 has been assembled immediately from the accessible Batch A files. Panels supported by the available data were completed; unsupported panels contain clearly labeled temporary placeholders.

### Completed or partially completed panels

- **Fig1a:** complete conceptual screen workflow using manuscript-supported timing and D3 CDH5 MACS sorting.
- **Fig1b:** temporary hit-only ranked U-test-effect panel; candidate A–C comparison exported.
- **Fig1c:** completed limited manuscript-defined known-regulator recovery panel for CDH5, KDR and ETV2.
- **Fig1d:** placeholder plus a populated 244-row author-review template.
- **Fig1e:** placeholder because secondary-validation and replicate data were not accessible.

## Panel-to-input mapping

- Fig1a: `CC-Perturb-seq EC Differentiation - Text V2 (1).pdf`; `Figure Redesign Decisions and Precedent Brief.pdf`.
- Fig1b: `source_data/BatchA_screen_hits_244_derived_source.csv`; attachment-index schema from `gw_245_hits.txt`.
- Fig1c: same screen-hit source plus the manuscript-defined genes CDH5, KDR and ETV2.
- Fig1d: `cleaned_data/screen_regulators_244_annotation_template.tsv`; final category fields pending author completion.
- Fig1e: secondary validation result ZIP/table not accessible; placeholder only.

See `data_audit/data_audit.tsv` and `data_audit/file_to_panel_map.tsv` for the full audit.

## Cleaning and deduplication rules

- Main count uses unique biological genes.
- The locally available table contains {summary['n_gene_level_hits']} unique genes.
- Pipe-separated original target identifiers were expanded into `gw_screen_hits_target_level_audit.tsv`.
- NEDD8 has two source target identifiers and is retained as two target-level rows but one gene-level row.
- Original gene and target IDs are preserved.
- Alias decisions are in `cleaned_data/gene_alias_mapping.tsv`.
- RBM14→RBM14-RBM4 is preserved from the available guide-derived alias table and remains flagged for author confirmation.

## Statistical quantities

- Temporary Fig1b and Fig1c use the Mann–Whitney U-test log2 enrichment in CDH5-positive versus CDH5-negative populations.
- Available U-test adjusted P values are stored as `u_test_fdr`.
- The final 244-hit set is a union of U-test and MaGECK calls, but exact MaGECK per-target statistics were not locally accessible.
- Candidate A is coherent for the U-test but contains hits only.
- Candidate B is a placeholder.
- Candidate C was selected as the temporary main panel because it uses one coherent effect quantity and avoids combining an effect from one method with significance from another.

## Hit-direction definition

- Negative available U-test effect: knockdown shifts toward CDH5-negative/non-EC; displayed as **Required for EC differentiation**.
- Positive available U-test effect: knockdown shifts toward CDH5-positive/EC; displayed as **Restricts EC differentiation**.
- Final direction labels must be verified against the complete union-table documentation because current draft artwork reports a different 238/6 split than the 236/8 sign split in the available U-test table.

## Category provenance

No final author-reviewed 244-gene functional annotation was attached. Categories remain blank in the review template; the placeholder bars in Fig1d are explicitly schematic.

## Color system

- Required for EC differentiation: blue `#0072B2`.
- Restricts EC differentiation: vermillion `#D55E00`.
- Temporary placeholders: muted rose `#CC6677`.
- Background/non-data structure: neutral gray.

## Output formats

Each individual panel is exported as SVG, vector PDF and transparent PNG at 600 dpi. The assembled figure is exported as SVG, vector PDF and white-background PNG at 600 dpi. A PowerPoint version is exported with each panel inserted as a separate SVG object.

## Software versions

```json
{json.dumps(software, indent=2)}
```

## Regeneration commands

From the `figure_build/scripts` directory:

```bash
python 00_audit_and_clean.py
python 01_build_figure1.py
python 03_assemble_figures.py
python 04_export_source_data.py
node 05_build_figure1_pptx.js
```

## Placeholders remaining

- Final full-universe Figure 1b.
- Final functional categories in Figure 1d.
- Final secondary validation/reproducibility panel in Figure 1e.

See `unresolved_issues.md` and `data_audit/consolidated_missing_input_request.md`.
"""
(ROOT / "README_figure_build.md").write_text(readme)

print("Exported Figure 1 source data, README and draft legend.")
