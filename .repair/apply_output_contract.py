from pathlib import Path
import ast
import shutil
import nbformat

root=Path.cwd()
def replace(text, old, new, count=1):
    observed=text.count(old)
    if observed != count:
        raise RuntimeError(f'Expected {count} occurrences, found {observed}: {old[:100]!r}')
    return text.replace(old,new,count)

shutil.copy2('/tmp/s4_output_contract.py',root/'src/cc_perturbseq/figure_s04_contract.py')
shutil.copy2('/tmp/test_s4_output_contract.py',root/'tests/test_s4_output_contract.py')
p=root/'src/cc_perturbseq/figure_s04_compute.py';s=p.read_text()
s=replace(s,'from cc_perturbseq import crt_qq_backend as qqb',
          'from cc_perturbseq import crt_qq_backend as qqb\n    from cc_perturbseq.figure_s04_contract import run_program_crt')
s=replace(s,'''out = programde.run_all_genes_union_crt(
                inputs=inputs,
                B=B_REP,
                n_jobs=N_JOBS,
                calibrate_skew_normal=False,
                return_raw_pvals=True,
                return_skew_normal=False,
            )''','''# betas_df is always returned. Raw/skew supplement flags must both
            # be False when skew calibration is disabled in real programDE.
            out = run_program_crt(
                programde, inputs=inputs, B=B_REP, n_jobs=N_JOBS,
                calibrate_skew_normal=False,
            )''')
p.write_text(s)
p=root/'src/cc_perturbseq/crt_qq_backend.py';s=p.read_text()
s=replace(s,'        "grna_target",\n    ]','        "grna_target",\n        "intended_target_name",\n    ]')
s=replace(s,'''    raw_mapping = dict(
        zip(
            table[resolved_guide_col].astype(str),
            table[resolved_target_col].astype(str),
        )
    )''','''    if guide_col and resolved_guide_col != guide_col:
        raise KeyError(f"Explicit guide column {guide_col!r} was not found.")
    if target_col and resolved_target_col != target_col:
        raise KeyError(f"Explicit target column {target_col!r} was not found.")
    # Do not overwrite inferred SAFE/non-targeting/NTC labels with string 'nan'.
    selected = table[[resolved_guide_col, resolved_target_col]].dropna().astype(str)
    selected[resolved_guide_col] = selected[resolved_guide_col].str.strip()
    selected[resolved_target_col] = selected[resolved_target_col].str.strip()
    missing = {"", "nan", "none", "na", "n/a", "null"}
    selected = selected.loc[
        ~selected[resolved_guide_col].str.lower().isin(missing)
        & ~selected[resolved_target_col].str.lower().isin(missing)
    ]
    conflicts = selected.groupby(resolved_guide_col)[resolved_target_col].nunique()
    if (conflicts > 1).any():
        raise ValueError(f"Conflicting target annotations for guides: {list(conflicts[conflicts > 1].index[:10])}")
    raw_mapping = dict(zip(selected[resolved_guide_col], selected[resolved_target_col]))''')
s=replace(s,'except Exception as exc:  # keep fallback mapping usable\n                warnings.warn(',
          'except Exception as exc:  # retain legacy fallback only when no columns were explicit\n                if guide_col or target_col:\n                    raise\n                warnings.warn(')
s=replace(s,'''    out = programde.run_all_genes_union_crt(
        inputs=inputs,
        B=args.B,
        n_jobs=args.n_jobs,
        calibrate_skew_normal=not args.no_skew,
        return_raw_pvals=True,
        return_skew_normal=not args.no_skew,
    )

    pvals_empirical = out["pvals_raw_df"]
    pvals_skew = out.get("pvals_df")
    if pvals_skew is None:
        pvals_skew = out.get("pvals_skew_df")''','''    from .figure_s04_contract import run_program_crt, extract_pvalue_outputs
    out = run_program_crt(
        programde, inputs=inputs, B=args.B, n_jobs=args.n_jobs,
        calibrate_skew_normal=not args.no_skew,
    )
    pvals_empirical, pvals_skew = extract_pvalue_outputs(
        out, calibrate_skew_normal=not args.no_skew,
    )''')
s=s.replace('source="run_all_genes_union_crt pvals_raw_df",',
            'source="run_all_genes_union_crt " + ("pvals_df" if args.no_skew else "pvals_raw_df"),')
s=s.replace('source="observed target rows from pvals_raw_df",',
            'source="observed target rows from " + ("pvals_df" if args.no_skew else "pvals_raw_df"),')
p.write_text(s)
p=root/'src/cc_perturbseq/figure_s04_runtime.py';s=p.read_text()
s=replace(s,'from . import crt_qq_backend as qqb',
          'from . import crt_qq_backend as qqb\nfrom .figure_s04_contract import audit_guide_metadata, validate_completed_full_day')
s=replace(s,'"sceptre_R_required": False, "inputs": records,',
          '"sceptre_R_required": False, "inputs": records,\n            "guide_metadata_audit": audit_guide_metadata(args),')
s=replace(s,'def run_worker(argv, *, preflight_only=False):',
          'def run_worker(argv, *, preflight_only=False, stage="all"):')
s=replace(s,'    args = qqb.parse_args(list(argv))',
          '    if stage not in {"all", "replicates"}:\n        raise ValueError("stage must be all or replicates")\n    args = qqb.parse_args(list(argv))')
s=replace(s,'    report = check_crt_runtime(args)',
          '    report = check_crt_runtime(args)\n    resume = validate_completed_full_day(args, report["guide_metadata_audit"]) if stage == "replicates" else None')
s=replace(s,'"B": args.B, "B_REP": 127})',
          '"B": args.B, "B_REP": 127, "stage": stage, "full_day_resume": resume})')
s=replace(s,'''        print("S4 full-day program CRT starting", flush=True)
        code = qqb.main(list(argv))
        if code != 0:
            raise RuntimeError(f"Full-day CRT returned exit status {code}")''','''        if stage == "all":
            print("S4 full-day program CRT starting", flush=True)
            code = qqb.main(list(argv))
            if code != 0:
                raise RuntimeError(f"Full-day CRT returned exit status {code}")
        else:
            print("Reusing verified completed full-day export; running replicates only.", flush=True)''')
s=replace(s,'def launch_crt_recomputation(argv, *, python_executable=None, preflight_only=False):',
          'def launch_crt_recomputation(argv, *, python_executable=None, preflight_only=False, stage="all"):')
s=replace(s,'    argv = list(map(str, argv))',
          '    if stage not in {"all", "replicates"}:\n        raise ValueError("stage must be all or replicates")\n    argv = list(map(str, argv))')
s=replace(s,'"check" if preflight_only else "run", "--", *argv]',
          '("check-" if preflight_only else "") + ("replicates" if stage == "replicates" else "run"), "--", *argv]')
s=replace(s,'''    if len(argv) < 2 or argv[0] not in {"check", "run"} or argv[1] != "--":
        raise ValueError("Usage: python -m cc_perturbseq.figure_s04_runtime {check|run} -- BACKEND_OPTIONS")
    run_worker(argv[2:], preflight_only=(argv[0] == "check"))''','''    modes = {"check": (True, "all"), "check-run": (True, "all"),
             "run": (False, "all"), "replicates": (False, "replicates"),
             "check-replicates": (True, "replicates")}
    if len(argv) < 2 or argv[0] not in modes or argv[1] != "--":
        raise ValueError("Usage: python -m cc_perturbseq.figure_s04_runtime {check|run|replicates} -- BACKEND_OPTIONS")
    preflight_only, stage = modes[argv[0]]
    run_worker(argv[2:], preflight_only=preflight_only, stage=stage)''')
p.write_text(s)
p=root/'notebooks/supplementary/Figure_S04.ipynb';nb=nbformat.read(p,as_version=4)
for cell in nb.cells:
    if cell.cell_type!='code':continue
    if 'RECOMPUTE_CRT = False' in cell.source:
        cell.source=replace(cell.source,'RECOMPUTE_CRT = False',
            'RECOMPUTE_CRT = False\n# Set to "replicates" to resume a completed full-day export after a replicate failure.\n# Its prior runtime/configuration/mapping records are verified before any output is changed.\nCRT_STAGE = "all"')
    if 'runtime_report = launch_crt_recomputation' in cell.source:
        cell.source=replace(cell.source,'launch_crt_recomputation(argv, python_executable=CRT_PYTHON)',
                            'launch_crt_recomputation(argv, python_executable=CRT_PYTHON, stage=CRT_STAGE)')
    compile(cell.source,'S4 code','exec')
nb.cells.insert(3,nbformat.v4.new_markdown_cell('''### Replicate-output flags and resuming

Replicate beta extraction uses `calibrate_skew_normal=False`,
`return_raw_pvals=False`, and `return_skew_normal=False`. `betas_df` is still
returned by programDE. The full-day skew-enabled call retains its original flags.
The optional no-skew path reads empirical results from `pvals_df`, not `pvals_raw_df`.

After a full-day export completed but replicate estimation failed, set
`RECOMPUTE_CRT=True` and `CRT_STAGE="replicates"`. The same output directory,
scientific settings, inputs, and guide mapping must match the completed run.
Missing/mismatched records stop before modifying outputs; use `CRT_STAGE="all"`
for a genuinely new or changed analysis. No completed p-values are silently
combined with a different guide mapping.

The v43 metadata schema (`guide_id`, `intended_target_name`) is now recognized.
Blank control targets retain their inferred control labels. The runtime report
records any target mapping differences from the previous name-based fallback.
'''))
nbformat.validate(nb);nbformat.write(nb,p)
p=root/'tests/test_s4_crt_runtime.py';s=p.read_text()
s=replace(s,'ns=dict(RECOMPUTE_CRT=True, build_crt_argv=runtime.build_crt_argv,',
          'ns=dict(RECOMPUTE_CRT=True, CRT_STAGE="all", build_crt_argv=runtime.build_crt_argv,')
s=replace(s,'def run_all_genes_union_crt(*,inputs,B,**kwargs):',
'''def run_all_genes_union_crt(*,inputs,B,calibrate_skew_normal=False,
                            return_raw_pvals=False,return_skew_normal=False,**kwargs):
    if (return_raw_pvals or return_skew_normal) and not calibrate_skew_normal:
        raise ValueError("raw/skew supplemental outputs require calibration")''')
s=replace(s,"    return {'betas_df':beta,'pvals_raw_df':p,'pvals_df':p}",
'''    result={'betas_df':beta,'pvals_df':p}
    if return_raw_pvals: result['pvals_raw_df']=p
    if return_skew_normal: result['pvals_skew_df']=p
    return result''')
p.write_text(s)
p=root/'docs/S4_CRT_RUNTIME.md';p.write_text(p.read_text()+'''
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
''')
for p in (root/'src/cc_perturbseq').glob('*.py'):compile(p.read_text(),str(p),'exec')
print('S4 output contract repair applied and notebook validated.')
