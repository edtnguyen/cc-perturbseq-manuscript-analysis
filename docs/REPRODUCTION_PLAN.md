# CC-Perturb-seq manuscript figure reproduction plan

## Scope

This branch organizes the code and provenance needed to reproduce the current main and supplementary figures in **Working CC-Perturb-seq EC Differentiation - Text V3**.

The raw recovery sources remain unchanged:

- `cc-perturbseq-manuscript-analysis:import/manuscript-scripts-20261013-121611`
- `programDE:sherlock-sync-20261013-070509`

Do not delete, rename, or rewrite files on those source branches. Canonical scripts will be copied or refactored onto this curation branch only after their panel provenance is established.

## Ownership rule

A panel is not assumed to be Tri-generated merely because it appears in the manuscript. Each panel is classified as one of:

1. **Tri / CC-Perturb-seq analysis** — recover and reproduce from existing code.
2. **Mixed analysis** — combine Tri-generated analysis with collaborator-generated inputs or subpanels.
3. **Collaborator-generated** — request the original source file, data, plotting code, and statistical output; do not reconstruct from a raster manuscript image unless the source is unavailable.
4. **Manual schematic or assembly** — preserve the editable source artwork and separately version any code-generated components.

The initial classifications in `PANEL_REPRODUCTION_MATRIX.csv` are provisional and inferred from the current manuscript legends, Methods, author contributions, and recovered filenames.

## Reproduction states

- `UNTRIAGED`: no candidate source inspected.
- `CANDIDATES_FOUND`: one or more possible scripts identified.
- `CANONICAL_SELECTED`: one source version selected and justified.
- `INPUTS_MAPPED`: every required input has a documented path or acquisition source.
- `RUNNABLE`: script/notebook runs from a clean entry point.
- `OUTPUT_MATCHED`: regenerated panel agrees with the current manuscript panel in data, labels, and visual encoding.
- `BLOCKED_COLLABORATOR`: requires a collaborator-owned source asset or dataset.
- `BLOCKED_INPUT`: required local/Sherlock input has not been located.
- `SUPERSEDED`: historical source retained but not used.

## Execution order

### Phase 1 — provenance inventory

For every panel in Figures 1–5 and Supplementary Figures S1–S13:

- identify likely author or analysis group;
- identify candidate scripts/notebooks in both source repositories;
- identify raw and intermediate inputs;
- identify the expected panel output;
- record unresolved statistical or labeling decisions.

### Phase 2 — canonical-source selection

For each candidate family:

- compare code and parameters rather than selecting by filename or modification date alone;
- check that statistical definitions agree with the current Methods;
- determine whether the output corresponds to the current panel rather than an older draft;
- retain provenance back to the original source path and commit.

### Phase 3 — first reproducible panel set

Prioritize code-derived CC-Perturb-seq panels that are both central and likely recoverable:

- Fig. 2a, 2b, 2d, 2e, 2f
- Fig. 3b, 3c, 3f
- Fig. 4c, 4d
- Fig. S7, S8, S9, S12

These panels exercise the principal shared components: cNMF usage, programDE/CRT results, guide-level EC fractions, carrier-seeded neighborhood enrichment, and Hippo-modulation single-cell projection.

### Phase 4 — collaborator-source collection

Create a request list for panels based on flow cytometry, qPCR, immunostaining, mechanical perturbation, zebrafish experiments, genome-wide-screen analyses, or multiome/ChromBPNet analyses that were generated primarily by collaborators.

For each requested panel, obtain where applicable:

- source data or the minimal plotted table;
- original plotting code;
- editable vector output (`.svg`, `.pdf`, `.ai`, or equivalent);
- exact statistical-test output;
- sample/replicate labels and ordering;
- software and package versions.

### Phase 5 — clean repository layout

Canonical code will be organized under:

```text
analysis/
  01_genome_wide_screen/
  02_cc_perturbseq_qc/
  03_program_effects/
  04_cell_state_and_neighborhoods/
  05_hippo_small_molecule/
  06_multiome_and_motifs/
figures/
  main/figure_01 ... figure_05/
  supplementary/figure_S01 ... figure_S13/
src/cc_perturbseq/
config/
docs/
```

Large data and generated results must not be committed. Machine-specific paths belong in an ignored local configuration file; the repository will track only `config/paths.example.yaml`.

## Output requirements

- Save editable vector output as PDF and/or SVG.
- Keep text editable in Illustrator when technically possible.
- Rasterize only dense point layers when necessary for file size.
- Record random seeds and statistical-test parameters.
- Generate a run receipt listing script commit, input paths/checksums where feasible, output paths, runtime, and success/failure state.
- Preserve panel-level data tables used to plot statistics when they are small and suitable for public release.

## Manuscript inconsistencies

The current V3 document contains tracked-change remnants and unresolved alternatives in several places. Figure reconstruction must not silently choose between conflicting counts, labels, or tests. Such cases are recorded as manuscript decisions in the panel matrix and resolved explicitly before a panel is marked `OUTPUT_MATCHED`.
