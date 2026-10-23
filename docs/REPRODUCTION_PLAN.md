# CC-Perturb-seq manuscript figure reproduction plan

## Authoritative scope

This branch now targets only the figures requested for remaking from **Working CC-Perturb-seq EC Differentiation - Text V3**:

- Main figures: **Figure 1, Figure 2, Figure 3b-d, and Figure 4d-e**.
- Supplementary figures: **Figures S1-S9**.

Panels outside this scope are retained in the historical inventory but are not part of the current execution plan.

The raw recovery sources remain unchanged:

- `cc-perturbseq-manuscript-analysis:import/manuscript-scripts-20261013-121611`
- `programDE:sherlock-sync-20261013-070509`

Do not delete, rename, or rewrite files on those source branches. Canonical scripts will be copied or refactored onto this curation branch only after their panel provenance is established.

## Upstream-data rule

The genome-wide and focused validation screens were processed by other authors. Their result tables are treated as **fixed upstream inputs**. The present task is to reproduce the manuscript panels from those result tables, not to rerun FASTQ processing, guide counting, MAGeCK, or the original U-test by default.

For every upstream screen input we will:

1. record the exact Sherlock path and file checksum;
2. retain the original column names and analysis provenance;
3. audit any gene-name remapping or deduplication before plotting;
4. verify that the plotted effect, significance field, and direction agree with the current manuscript legend;
5. avoid silently recomputing or substituting statistics.

Rerunning an upstream screen analysis is required only if a discrepancy cannot be resolved from the supplied result tables and original analysis scripts.

## Ownership rule

A panel is not assumed to be Tri-generated merely because it appears in the manuscript. Each panel is classified as one of:

1. **Tri / CC-Perturb-seq analysis** — recover and reproduce from existing code.
2. **Mixed analysis** — combine Tri-generated analysis with collaborator-generated inputs or subpanels.
3. **Collaborator-generated upstream result** — use the supplied processed table or request the original source asset; do not reverse-engineer values from a raster panel.
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

## Work packages and execution order

### WP1 — Screen-result panels

Panels:

- Fig. 1b
- Fig. S1c-d
- Fig. S7a

Inputs are the processed genome-wide and validation-screen result tables documented in `SHERLOCK_INPUT_MANIFEST.md`. These panels can be reproduced without rerunning the screen pipelines.

Fig. S1a-b are handled separately because they require flow-cytometry source material rather than the screen result tables.

### WP2 — Core CC-Perturb-seq cell-state and perturbation panels

Panels:

- Fig. 1c-d
- Fig. 2a-b and 2d-f
- Fig. 3b-c
- Fig. 4d
- Fig. S3
- Fig. S4
- Fig. S7b
- Fig. S8b-g
- Fig. S9

Shared dependencies include the final day-specific AnnData objects, SCEPTRE guide assignments, K = 50 cNMF usages and spectra, programDE/CRT outputs, guide-level EC-percentage results, the integrated UMAP, and carrier-seeded neighborhood outputs.

### WP3 — Program annotation and multiome motif panels

Panels:

- Fig. 3d
- Fig. 4e
- Fig. S5
- Fig. S6

These panels combine CC-Perturb-seq perturbation effects and cNMF loadings with scE2G/ChromBPNet/TF-MoDISco/Fi-NeMo motif results. The exact final multiome input set must be confirmed before these panels are marked runnable.

### WP4 — CC-Perturb-seq development and benchmarking

Panels:

- Fig. S2b-n

Subgroups:

- S2b-h: circular capture versus CROP-seq benchmarking in hiPSCs and TeloHAECs.
- S2i-n: Parse Evercode v3 compatibility experiment in hiPSCs and K562 cells.

The canonical source notebooks and their Sherlock input paths must be selected separately for these two experiments.

### WP5 — Schematics and source-asset panels

Panels:

- Fig. 1a
- Fig. 2c
- Fig. S1a-b
- Fig. S2a
- Fig. S8a

The schematics require editable source artwork or a deliberate redraw. Fig. S1a-b require original FCS files, FlowJo workspaces, or collaborator-exported vector plots with documented gates.

## Canonical-source selection

For each candidate family:

- compare code and parameters rather than selecting by filename or modification date alone;
- check that statistical definitions agree with the current Methods;
- determine whether the output corresponds to the current panel rather than an older draft;
- retain provenance back to the original source path and commit;
- separate data computation from panel styling and assembly.

## Clean repository layout

Canonical code will be organized under:

```text
analysis/
  01_screen_results/
  02_cc_perturbseq_qc/
  03_program_effects/
  04_cell_state_and_neighborhoods/
  05_multiome_and_motifs/
  06_method_benchmarking/
figures/
  main/figure_01/
  main/figure_02/
  main/figure_03/
  main/figure_04/
  supplementary/figure_S01/ ... figure_S09/
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
- Preserve panel-level source-data tables when they are small and suitable for public release.

## Manuscript and legacy-code checks

The current manuscript contains tracked-change remnants and unresolved alternatives. Figure reconstruction must not silently choose between conflicting counts, labels, correlations, or tests.

Known screen-panel checks include:

- directional MAGeCK FDR versus a single combined FDR;
- gene/promoter identifier cleanup and the duplicated `NEDD8` entry;
- Pearson correlation versus linear-model R-squared in Fig. S1c-d;
- the exact pair of genome-wide and validation-screen effect columns used for Fig. S1d;
- whether significance colors represent genome-wide hits, validation hits, or their intersection.

Such cases are recorded as explicit decisions before a panel is marked `OUTPUT_MATCHED`.
