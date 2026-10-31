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

## Output contract repair and resume

The previous test stub accepted a combination rejected by real programDE:
`calibrate_skew_normal=False, return_raw_pvals=True`. The new shared call boundary
sets both supplementary output flags to the calibration flag. For replicate
betas all three are False. In no-skew full-day mode, empirical p-values come
from pvals_df and no skew output is labeled or exported. Skew-enabled full-day
behavior is unchanged. B_REP=127 and the statistical/plotting methods are unchanged.

The guide loader recognizes intended_target_name, ignores missing target values
rather than replacing control labels with 'nan', and rejects conflicting guide
annotations or invalid explicit column selections. The worker audits mappings
before computation. It retains the original S4 cNMF inputs; it does not replace
that model with a different run merely because another file is newer.

Set CRT_STAGE='replicates' to avoid repeating a completed full-day CRT after a
replicate-only failure. Resume checks prior arguments, each completed day record,
required output files, input timestamps, and the mapping audit. Legacy v43
fallback runs are reusable only when all metadata mappings equal prior inferred
names. Otherwise rerun both stages with CRT_STAGE='all'. These checks verify
recorded provenance, not mathematical correctness of existing study results.

Regression tests now call the actual programDE implementation pinned at
11cfa3f5d9c8d400c15cd15e1e550fc11fe0f4f9 on small synthetic inputs. The invalid
combination must fail there, the corrected call must return finite beta and
empirical outputs, and beta results must match the direct API call. Existing
full-notebook integration stubs also enforce the real flag constraint.
No execution on the actual Oak datasets is claimed.

### Preserve the source's target grouping

The real v43 audit found differences between guide-prefix names and metadata
labels, including gene aliases and per-guide SAFE target identifiers. Therefore
recognizing intended_target_name must not silently change the reproduction.
GUIDE_MAPPING_POLICY='source_names' is now explicit in the notebook and applied
to BOTH stages. Metadata is read and audited but original target names, control
groups, and gene-name-dependent random seeds are preserved. The optional
metadata policy is an explicit alternative, keeps controls canonical, and is
not used to resume a source-names full-day export. Actual v43 metadata and a
SAFE_TARGETING example are exercised in regression tests.
