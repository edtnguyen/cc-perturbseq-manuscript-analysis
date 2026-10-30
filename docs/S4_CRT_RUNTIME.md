# S4 CRT runtime correction

The source program CRT notebook's executed command used
`/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python` and
`--guide-long-csv-template`. The separate `ipsc_ec_sceptre.ipynb` ran under
`/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/bin/python` with
rpy2/R. S4 consumes its saved guide assignments; it does not rerun R-SCEPTRE.

The re-scripted notebook passed the invalid `--guide-long-template` spelling.
It now uses a parser-validated builder and the original canonical flag. The
backend accepts the shorter spelling as a compatibility alias; other invalid
flags still fail. Explicit argument lists exclude Jupyter connection arguments.

`CRT_PYTHON` defaults to the notebook's `sys.executable`, matching the source.
Set it in the notebook, or set `CCPERTURBSEQ_S4_PYTHON`, to the sc-dl Python path
above to run CRT independently of the plotting kernel. Both full-day CRT and
replicate beta extraction use this same Python in a fresh process. There is no
old-notebook execution. Missing dependencies, wrong programDE imports and
missing scientific input files are checked before CRT. Logs and the runtime
report are written to `recomputed_crt/s4_crt.log` and `s4_crt_runtime.json`.

The original calculation and plotting functions are preserved. Defaults remain
B=10000, B_REP=127, minimum 200 cells per replicate, and the original covariates
and NTC labels. K, density threshold, overrides and labels now propagate to both
stages. The overwrite setting also regenerates replicate beta files instead of
silently reusing old matrices. A worker error stops subsequent plotting.

Tests cover the actual notebook command, both flag spellings, Jupyter argv,
interpreter selection, error propagation and parameters. The integration test
executes the entire notebook with RECOMPUTE_CRT=True using small H5MU/CSV/TSV
files, four full-day and eight replicate calls, and checks generated tables and
figures. Only the programDE numerical engine is a deterministic test stub; this
is an orchestration test, not validation of study CRT statistics. Production
imports the real user programDE checkout. No actual Oak study run is claimed.
