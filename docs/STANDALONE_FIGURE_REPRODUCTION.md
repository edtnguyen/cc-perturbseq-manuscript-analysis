# Re-scripted figure notebooks — source audit

The seven retained figure notebooks are executable analysis entry points. None
opens, hashes for admission, patches, executes, or copies plots from an old
notebook. Recovery archives are provenance records only and may be absent.

| Entry point | Actual computation | Audit action |
| --- | --- | --- |
| Figure_01 | Local screen, EC fraction, and UMAP producers | Audited; accepted code/style unchanged |
| Figure_02 | Local effect matrices, statistics, feature UMAPs and carrier enrichment | Audited; accepted code/style unchanged |
| Figure_03 | Local Hippo program annotation, neighborhood and time-course producers | Audited; accepted code/style unchanged |
| Figure_04 | Local early-time-course, UMAP and annotation producers | Audited; accepted code/style unchanged |
| Figure_S02 | Local RNA/guide processing plus original plotnine panel code | Re-scripted from the later Parse-v3 paper-panel notebook; no old-notebook hash dependency |
| Figure_S04 | Recompute Q-Q coordinates, target-program joins, correlations, plots from numeric files | Removed old notebook executor and figure copier; ported CRT backend and plotting functions |
| Figure_S06 | Load six scientific inputs, select features, calculate shared scales and render logos/heatmaps | Removed old notebook executor; ported original selection and drawing functions |

## S2

Scientific sources are the recovered `K562_hiPSC_ParseV3_plotnine_figure_reproduction.ipynb`
and the underlying Parse processing. The previous curated plotting module was
an approximate seaborn redraw; it is now replaced by the later source's
plotnine functions. Paper expected-pool/specificity/top-guide panels use the
source palettes, layers, facets, axes and themes. S2l retains the no-guide
control through the source function's existing `include_no_guide` option.
Transcript yield/saturation and guide mixing use the source cells 21 and 27.
The mixing marginals are exported separately, as in that source.

The raw-GEX route labels whole Leiden clusters from their lineage-marker scores.
A per-cell zero-marker override from the separate low-memory export method had
been mixed into that route; it is removed. These are not interchangeable
annotation methods. The complete final metrics/labels are saved as
`Figure_S02_all_cell_metrics.csv.gz`. Supply this file to the GEO exporter to
freeze the accepted labels rather than reclustering during release preparation.
The paper source can also use older cached labels; no equivalence to such a
cache is claimed without comparing its cell IDs and labels.

## S4

Default: read `all_days_QQ_points.csv` and `replicate_beta_manifest.csv` with the
referenced numeric beta matrices. `CCPERTURBSEQ_S4_RESULTS` selects this input
model directory. Existing PDF/SVG/PNG files are never inputs. Q-Q coordinates,
real-target pair matching, Pearson/Spearman summaries, and plots are regenerated.

Set `RECOMPUTE_CRT=True` in the notebook to run the ported backend from the
original H5MU/usage/guide inputs. It retains B=10000 for the full-day CRT,
B_REP=127 for replicate beta extraction, full covariates, and original NTC
labels. This optional raw-analysis mode requires the scientific programDE
software and its data, but not the old notebook or its generated Python script.
No hidden long CRT run occurs in the default plotting mode.

## S6

`figure_s06_annotations.py` contains the original helper, feature-selection,
shared-scale, vector-logo and rendering functions. The notebook itself performs
the source's load/selection/export/render steps. The active min10 enrichment TSV,
current cluster map and standard compiled MoDISco H5 remain selected. The cNMF
workbook, program annotation workbook and consensus spectra are scientific
inputs; the previous notebook and its SHA are not execution gates.

## Validation and limits

`STANDALONE_FIGURE_SOURCE_AUDIT.json` records original source identities, cell
numbers and normalized Python AST hashes for the ported functions. Tests check
those functions against their recorded source hashes, not against an installed
copy of the old notebook.

The end-to-end tests create a separate checkout containing only current code,
configuration and notebooks. It contains no `recovery/` directory, no previous
figure files and no private scientific input files. Fresh Python processes
execute every curated S2/S4/S6 code cell against synthetic numerical inputs.
The tests verify source tables, correlations, vector logos and image exports.
GEO H5MU roundtrip tests are rerun against the same S2/Hippo producers.

These tests establish executable independence and source-port consistency,
not numerical agreement with the real Oak data or visual identity to the
assembled manuscript. That final comparison requires the real-data rerun.
Small historical cell-count differences remain diagnostic, not fatal.

The Hippo rebuild notebook optionally runs the *current* Figure 3/4 notebooks
as orchestration. It does not execute the archived production notebooks; the
rebuild/projection mathematics already live in ordinary Python modules.

The GEO preparation/testing notebooks likewise use repository functions and
submitted H5MU files, not old notebooks. Re-export the mixed-cell release after
accepting the updated S2 outputs, since its stored labels/source tables refer
to the producer version used at export. Hippo scientific calculations are
unchanged.

A rendering-only compatibility guard disables invisible automatic log minor
breaks in S2l. The original theme hides these grids; major ticks, limits,
values, layers, and source functions are unchanged. This prevents a mizani
IndexError for a narrow range containing only one major tick.

The rendered S6 previews exposed a DPI-dependent placement defect inherited
from the source's fixed pixel-coordinate logo rotation. The port expresses the
same 90-degree rotation in axes coordinates with physical-aspect correction.
This retains the motif matrices, glyphs, colors, orientation and intended
placement while keeping logos in position for PDF/SVG/600-dpi PNG exports.
The source and adapted function hashes are recorded separately, and an actual
logomaker glyph test verifies identical physical bounds at 72, 100 and 600 dpi.
