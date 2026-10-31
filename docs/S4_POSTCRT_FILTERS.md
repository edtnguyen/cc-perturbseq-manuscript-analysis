# S4c final post-CRT selection

Source: executed Sherlock version of
`iPSC_EC_QQ_full_covariates_with_replicate_correlation.ipynb`, cells 29-31;
SHA-256 6e80625abd9b5e71d45da6dd19e888898318d39a159a7ac35600fc1d85489e87.
The shorter initial raw-download notebook lacks these later cells. The prior
port selected the earlier unfiltered scatter rather than the final panel.

The standalone module `figure_s04_postcrt.py` now applies the final rules:
positive guide assignments; unique target carriers present in the ORIGINAL
cNMF usage index; >=20 such cells in each replicate; pooled full-day BH q<=0.10.
BH includes all finite real target-program tests within day BEFORE pair/count
filters. NTC and non-finite effects are excluded. Beta values are not refitted
or trimmed. There is no automatic top-N fallback or requirement for equal signs.
The original selected-scatter function, formatting, dimensions and PDF/SVG/PNG
exports are ported. Source-function AST hashes are separately recorded/tested.

S4b remains empirical. The archived final S4c cell preferred skew-calibrated
pooled p-values, falling back to empirical when absent; the actual type/path is
recorded in S4c_pooled_pvalue_input_audit.csv. This choice is explicit rather
than claiming empirical and skew selections are equivalent.

Completed full-day and replicate results are reusable with RECOMPUTE_CRT=False.
Set INPUT_MODEL_ROOT (or CCPERTURBSEQ_S4_RESULTS) to the completed model directory,
including k_50.dt_0_2. Counts are reconstructed from recorded H5MU metadata,
usage IDs and assignment CSVs, or loaded from a same-model source diagnostic
when the original inputs are unavailable. TARGET_COUNTS_CSV provides an explicit
numerical override. Missing counts raise an error, never skip the cell filter.
No original notebook is a runtime dependency and no previous figures are copied.

Use pooled_selected_* figures. Old unfiltered files are not deleted but are
excluded from new receipts. S4c_filter_audit.csv and S4c_filter_counts.csv record
inclusion/exclusion and counts by day; historical totals are not exact gates.

Tests cover the prepared filtering rules and original-function equivalence,
the full numerical-only notebook, actual H5MU metadata reconstruction, and the
recomputation notebook with a deterministic programDE stub. The separate real
programDE tests remain in the full suite. Synthetic tests do not establish
study-level numerical or visual agreement with actual Oak data.

Source-helper fingerprints use sorted compact JSON of AST node names and fields,
excluding empty type_params (added in Python 3.12). They were independently
computed from code cell 31 in the executed source text, Git blob
333850b851ddd13586f5635a3adefb33ad7e151e, and compared to the port before publication.
They do not depend on Python-version-specific ast.dump formatting.
