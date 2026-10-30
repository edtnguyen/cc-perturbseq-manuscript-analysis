from pathlib import Path
import shutil
import nbformat
r=Path.cwd()
p=r/'src/cc_perturbseq/crt_qq_backend.py'
s=p.read_text()
old='parser = argparse.ArgumentParser(\n        description='
assert old in s
s=s.replace(old, 'parser = argparse.ArgumentParser(\n        prog="cc_perturbseq.crt_qq_backend",\n        allow_abbrev=False,\n        description=',1)
old='parser.add_argument("--guide-long-csv-template", default=DEFAULT_GUIDE_LONG_TEMPLATE)'
assert old in s
s=s.replace(old, '''parser.add_argument(
        "--guide-long-csv-template", "--guide-long-template",
        dest="guide_long_csv_template", default=DEFAULT_GUIDE_LONG_TEMPLATE,
        help="SCEPTRE-assigned guide CSV template. The shorter flag is a compatibility alias.",
    )''',1)
p.write_text(s)
p=r/'src/cc_perturbseq/figure_s04_compute.py';s=p.read_text()
a=s.index('def compute_replicate_betas('); b=s.index('    from pathlib import Path',a)
s=s[:a]+'''def compute_replicate_betas(
    *, programde_root, data_dir, cnmf_dir, guide_long_template,
    guide_metadata_tsv, output_dir, days=("d0", "d1", "d2", "d3"), n_jobs=4,
    k=50, density_threshold="0_2", ntc_labels=("SAFE", "non-targeting", "NTC"),
    overwrite=False, cnmf_file=None, usage_file_template=None,
    gex_modality=None, guide_modality=None,
    guide_metadata_guide_col=None, guide_metadata_target_col=None,
):
    PROGRAMDE_ROOT, DATA_DIR, CNMF_DIR = str(programde_root), str(data_dir), str(cnmf_dir)
    GUIDE_LONG_TEMPLATE = str(guide_long_template)
    GUIDE_METADATA_TSV = str(guide_metadata_tsv) if guide_metadata_tsv is not None else None
    OUT_DIR, DAYS, N_JOBS = str(output_dir), list(days), int(n_jobs)
    K, DENSITY_THRESHOLD, NTC_LABELS = int(k), str(density_threshold), list(ntc_labels)
'''+s[b:]
s=s.replace('SKIP_EXISTING_REP_BETAS = True','SKIP_EXISTING_REP_BETAS = not overwrite',1)
s=s.replace('cnmf_file=None,\n        usage_file_template=None,','cnmf_file=cnmf_file,\n        usage_file_template=usage_file_template,',1)
s=s.replace('gex_modality=None,\n        guide_modality=None,','gex_modality=gex_modality,\n        guide_modality=guide_modality,',1)
s=s.replace('guide_metadata_guide_col=None,\n        guide_metadata_target_col=None,','guide_metadata_guide_col=guide_metadata_guide_col,\n        guide_metadata_target_col=guide_metadata_target_col,',1)
s=s.replace('and beta_path.exists():','and beta_path.is_file() and beta_path.stat().st_size > 0:',1)
p.write_text(s)
p=r/'notebooks/supplementary/Figure_S04.ipynb';nb=nbformat.read(p,as_version=4)
s=nb.cells[1].source
s=s.replace('from cc_perturbseq.figure_s04_compute import compute_replicate_betas','from cc_perturbseq.figure_s04_runtime import build_crt_argv, launch_crt_recomputation')
s=s.replace('RECOMPUTE_CRT = False', '''RECOMPUTE_CRT = False
# Original S4 CRT used sc-dl; R-SCEPTRE guide assignment is a separate upstream step.
# Explicit override: /oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python
CRT_PYTHON = os.environ.get("CCPERTURBSEQ_S4_PYTHON", sys.executable)
print("Notebook Python:", sys.executable)
print("CRT worker Python:", CRT_PYTHON)''')
nb.cells[1].source=s
nb.cells[2].source='''## Numerical inputs and runtime

The archived program CRT/replicate notebook ran with
`/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python`.
The separate `ipsc_ec_sceptre.ipynb` used
`/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/bin/python`
and `rpy2`/R to assign guides. This S4 notebook **reads those saved assignments**;
it does not run R-SCEPTRE or require switching its plotting kernel to R.

With `RECOMPUTE_CRT=True`, the code below validates the actual CLI arguments,
then starts the repository's Python CRT worker using `CRT_PYTHON`. The worker
checks its dependencies, programDE import location, and all requested day inputs
before any CRT. Both full-day CRT and replicate beta estimation run in that same
interpreter. No archived notebook is executed. The original calculation and
plotting functions remain in repository modules.

The default `CRT_PYTHON=sys.executable` follows the original notebook. Set it to
the explicit sc-dl Python path above to keep plotting in another kernel. A failed
worker stops the notebook; its traceback is retained in `recomputed_crt/s4_crt.log`.
The runtime/import/input report is `recomputed_crt/s4_crt_runtime.json`.
'''
s=nb.cells[3].source;end=s.index('\nqq_csv =')
nb.cells[3].source='''if RECOMPUTE_CRT:
    argv = build_crt_argv(
        programde_root=PROGRAMDE_ROOT, data_dir=DATA_DIR, cnmf_dir=CNMF_DIR,
        guide_long_template=GUIDE_LONG_TEMPLATE, guide_metadata_tsv=GUIDE_METADATA_TSV,
        output_dir=RECOMPUTED_ROOT, days=DAYS, k=K,
        density_threshold=DENSITY_THRESHOLD, B=B, n_jobs=N_JOBS,
        ntc_labels=NTC_LABELS, overwrite=True,
    )
    # Uses the source parser's --guide-long-csv-template spelling. It is parsed
    # before a child process is started; Jupyter's own argv is never consumed.
    runtime_report = launch_crt_recomputation(argv, python_executable=CRT_PYTHON)
    display(runtime_report)
    INPUT_MODEL_ROOT = RECOMPUTED_ROOT / qqb.model_tag(K, qqb.normalize_density_token(DENSITY_THRESHOLD))
'''+s[end:]
nbformat.validate(nb)
for i,c in enumerate(nb.cells):
    if c.cell_type=='code':compile(c.source,f'S4 cell {i}','exec')
nbformat.write(nb,p)
shutil.copy2('/tmp/s4_runtime.py',r/'src/cc_perturbseq/figure_s04_runtime.py')
shutil.copy2('/tmp/test_s4_crt_runtime.py',r/'tests/test_s4_crt_runtime.py')
p=r/'config/geo-reproduction-requirements.txt'
s=p.read_text();assert 'muon==' not in s
p.write_text(s+'\n# Original S4 raw-CRT H5MU reader; included so recomputation tests run.\nmuon==0.1.7\n')
(r/'docs/S4_CRT_RUNTIME.md').write_text('''# S4 CRT runtime correction

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
''')
for p in (r/'src').rglob('*.py'):compile(p.read_text(),str(p),'exec')
print('S4 runtime correction applied and compiled.')
