# Standalone figure validation — 2026-10-30

Tested scientific-code commit: `bdb20be054fb65656e5fc30720a1cdbd6f07a5e2`.

Validation run: https://github.com/edtnguyen/cc-perturbseq-manuscript-analysis/actions/runs/41045607082

The complete repository suite passed before the candidate was pushed. An earlier run of the same code also reported `77 passed, 36 warnings`; its tests and main-figure preflight passed, but publication failed because the Actions token could not modify workflow files. Publication was then separated from workflow registration; no analysis logic was changed to resolve the permission issue.

## Executed checks

- JSON/nbformat and code compilation for curated notebooks.
- All seven figure entry points audited for archived-notebook execution, dynamic source execution, and copying old figures.
- Entire S2, S4 and S6 notebook code-cell sequences executed in separate Python processes with synthetic numerical inputs. Their temporary checkouts contained no `recovery/` directory, no old notebooks, and no previous figure outputs. File-access guards rejected Oak, scratch and recovery accesses.
- S2 used synthetic expression MTX, metadata, guide counts/reference and all four saturation reports through the actual Scanpy classification and panel producers.
- S4 regenerated Q-Q coordinates, matched real-target/program effects, correlations and plots from numerical p-values and beta matrices. The optional full-scale CRT computation was NOT run on the actual study data.
- S6 read synthetic workbook, consensus, enrichment, cluster-map and HDF5 CWM files, selected features, calculated shared scales, and produced all day panels and the combined panel.
- S2 narrow-log-range rendering and S6 motif-logo physical placement at 72/100/600 dpi were exercised.
- Both GEO H5MU roundtrips passed after removal of original fixture inputs.
- Main Figure 1–4 static preflight passed; their accepted notebooks and scientific producer code were not modified by this repair.

Rendered synthetic previews were inspected. They are validation fixtures, not manuscript results or GEO submission data.

## Limits

No actual Oak study matrix or full-study CRT run was executed here. These checks establish source-port consistency, executable independence and tested synthetic behavior, not numerical or pixel-for-pixel agreement with the assembled manuscript. Run the revised supplementary notebooks on the real inputs and review the outputs before accepting the release. Historical cell totals remain advisory rather than fixed acceptance gates.
