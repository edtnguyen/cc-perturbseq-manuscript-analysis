"""Apply the reviewed S4c notebook patch and strengthen integration fixtures."""
from pathlib import Path
import ast
import json
import nbformat

ROOT = Path.cwd()
nb_path = ROOT / 'notebooks/supplementary/Figure_S04.ipynb'
nb = nbformat.read(nb_path, as_version=4)
settings = next(c for c in nb.cells if c.cell_type == 'code' and 'RECOMPUTE_CRT = ' in c.source and 'INPUT_MODEL_ROOT = ' in c.source)
assert 'POOLED_Q_CUTOFF = ' not in settings.source
settings.source += '''
# S4c uses the FINAL executed source selection, not the earlier unfiltered panel.
POOLED_PVALUE_TYPE = "skew"
POOLED_Q_CUTOFF = 0.10
MIN_TARGET_CELLS_PER_REP = 20
TARGET_COUNTS_CSV = None  # optional explicit source target-count CSV
'''
plot = next(c for c in nb.cells if c.get('id') == 'd809a541')
plot.source = '''# Executed source cells 29-31: post-CRT filtering, not a new CRT.
from cc_perturbseq.figure_s04_postcrt import (
    load_or_rebuild_target_counts, read_pooled_results,
    select_effect_pairs, render_selected_effects,
)
PLOT_OUT_DIR = PANEL_C_DIR
unfiltered_effects = build_effect_pair_table(manifest)
unfiltered_path = PANEL_C_DIR / "S4c_unfiltered_real_effect_pairs.csv"
unfiltered_effects.to_csv(unfiltered_path, index=False)

target_counts = load_or_rebuild_target_counts(
    INPUT_MODEL_ROOT, DAYS, PANEL_C_DIR,
    counts_csv=TARGET_COUNTS_CSV,
    mapping_policy=GUIDE_MAPPING_POLICY,
)
pooled, pooled_inputs = read_pooled_results(
    INPUT_MODEL_ROOT, DAYS, preferred=POOLED_PVALUE_TYPE,
)
pooled.to_csv(PANEL_C_DIR / "S4c_pooled_full_day_pvalues_and_qvalues.csv", index=False)
pooled_inputs.to_csv(PANEL_C_DIR / "S4c_pooled_pvalue_input_audit.csv", index=False)
selected, all_pairs_audit, filter_summary = select_effect_pairs(
    unfiltered_effects, pooled, target_counts,
    q_cutoff=POOLED_Q_CUTOFF, min_target_cells=MIN_TARGET_CELLS_PER_REP,
)
all_pairs_audit.to_csv(PANEL_C_DIR / "S4c_filter_audit.csv", index=False)
filter_summary.to_csv(PANEL_C_DIR / "S4c_filter_counts.csv", index=False)
display(pooled_inputs)
display(filter_summary)
# Never replace an empty selection with an unfiltered or top-N panel.
effects, summary, combined_plot = render_selected_effects(
    selected, PANEL_C_DIR,
    q_cutoff=POOLED_Q_CUTOFF, min_target_cells=MIN_TARGET_CELLS_PER_REP,
)
display(summary)
display(combined_plot)
selection_tag = (
    f"q_q{str(POOLED_Q_CUTOFF).replace('.', '_')}_top500_mincells{MIN_TARGET_CELLS_PER_REP}"
)
# Source filename contains top500 even though selection mode is q; no fallback.
selection_settings_path = PANEL_C_DIR / "S4c_selection_settings.json"
selection_settings_path.write_text(json.dumps({
    "source_executed_notebook_cells": [29, 30, 31],
    "selection_mode": "q", "pooled_q_cutoff": POOLED_Q_CUTOFF,
    "min_target_cells_per_replicate": MIN_TARGET_CELLS_PER_REP,
    "pooled_pvalue_requested": POOLED_PVALUE_TYPE,
    "pooled_pvalue_inputs": pooled_inputs.to_dict("records"),
    "counts_definition": "unique target carriers with original cNMF usage, each replicate",
    "bh_family": "all finite real target-program tests within day, before pair/count filters",
    "no_crt_rerun_required_for_selection": True,
    "historical_counts_are_diagnostic_only": True,
}, indent=2) + "\\n")
selected_output_files = [
    unfiltered_path,
    PANEL_C_DIR / "DIAGNOSTIC_target_cell_counts_by_replicate.csv",
    PANEL_C_DIR / "target_cell_count_provenance.json",
    PANEL_C_DIR / "S4c_pooled_full_day_pvalues_and_qvalues.csv",
    PANEL_C_DIR / "S4c_pooled_pvalue_input_audit.csv",
    PANEL_C_DIR / "S4c_filter_audit.csv",
    PANEL_C_DIR / "S4c_filter_counts.csv",
    selection_settings_path,
    PANEL_C_DIR / f"pooled_selected_replicate_effect_pairs_{selection_tag}.csv",
    PANEL_C_DIR / f"pooled_selected_replicate_effect_correlation_summary_{selection_tag}.csv",
]
stems = [f"pooled_selected_biological_replicate_DE_program_effect_correlation_all_days_{selection_tag}"]
stems += [f"{day}_pooled_selected_biological_replicate_DE_program_effect_correlation_{selection_tag}"
          for day in effects["day"].drop_duplicates()]
selected_output_files += [PANEL_C_DIR / f"{stem}.{suffix}"
                          for stem in stems for suffix in ("pdf", "svg", "png")]
print("Selected S4c outputs:", PANEL_C_DIR)
print("Use pooled_selected_* figures, not unfiltered figures from an older run.")
'''
receipt = next(c for c in nb.cells if c.cell_type == 'code' and 'write_receipt(' in c.source and '"recomputed_crt"' in c.source)
receipt.source = '''for panel in ("b", "c"):
    if panel == "b":
        files = [Path(p) for day_outputs in all_saved_files.values() for p in day_outputs.values()]
        files.append(PANEL_B_DIR / "all_days_QQ_points.csv")
    else:
        files = selected_output_files
    if not any(p.suffix == ".svg" and p.is_file() for p in files):
        raise RuntimeError(f"Panel {panel}: no SVG output from this execution")
    write_receipt(
        ctx, panel, status="RESCRIPTED_NUMERICAL_INPUTS_RENDERED",
        inputs={"qq_points": str(qq_csv), "replicate_manifest": str(manifest_path)},
        outputs={p.name: str(p) for p in files if p.is_file()},
        parameters={
            "runtime_source_notebook_required": False,
            "recomputed_crt": RECOMPUTE_CRT,
            "source_style_changed": False,
            "panels": ["b", "c"],
            "s4c_selection": "executed source cells 29-31; pooled q and target-cell support",
            "selection_settings": str(selection_settings_path),
        },
    )
'''
for c in nb.cells:
    if c.cell_type == 'markdown' and c.source.startswith('## S4c'):
        c.source = '''## S4c — final post-CRT selection from executed source cells 29–31

The earlier unfiltered scatterplot was not the final source panel. This cell
joins actual target-carrier counts with original cNMF usage in each replicate,
then selects pooled full-day BH q <= 0.10 and >=20 carriers in BOTH replicates.
BH is computed separately within each day over all real target-program tests,
before any pair/count filtering. There is no beta trimming or top-N fallback.

The archived S4c selection preferred skew-calibrated pooled p-values, with an
explicitly recorded empirical fallback when skew files were absent. S4b remains
empirical. The input audit records the actual p-value files used. The source
scatter style (point size 1.15, alpha 0.82) and beta values are preserved.

For a completed run, leave RECOMPUTE_CRT=False and set INPUT_MODEL_ROOT to its
model directory. Reconstructing counts from recorded inputs is not a CRT rerun.
Missing counts never silently disable the >=20 filter. Historical exact counts
and correlations are diagnostic, not acceptance gates.
'''
    if c.cell_type == 'code':
        c.execution_count = None
        c.outputs = []
        compile(c.source, f'S4:{c.get("id")}', 'exec')
nb.metadata.setdefault('cc_perturbseq', {})['post_crt_selection_source_cells'] = [29,30,31]
nbformat.validate(nb)
nbformat.write(nb, nb_path)

hashes = {
  'normalize_program_name': '91e449c50564145d7607fa60379bf7e002f4c79b85652e71e0b0a6eb5a31edad',
  'bh_qvalue': 'bed04dc1beb561aee1d512cdf7b1118e122889b0f9ffc0f413ddb7eb3e117e23',
  'read_pvalue_long_or_wide': '33d99bf5ab800440201468d706258e63e257deaf5db16ea65fe4b91fc310a5c0',
  'summarize_replicate_correlations': 'd6972753085005d4841c5484d65c883a8427b4cdd69b41fcd219ad02ea747ca7',
  'format_r': 'cebca50ac8513d666bbb879bfe500c566a009347bbc318909e5b5ba03a11a3e5',
  'save_plotnine_vector_outputs': '5b4abbeadb0fe3cbccbe3c44f346c1abc8559309a0b4d24ca02fc851296e3e9f',
  'make_replicate_scatter': '7a76fce9d41bf22db9790c04354e68f4b441bd86a0b4026814a47803913085d2',
}
(ROOT/'docs/S4_POSTCRT_SOURCE_FUNCTION_HASHES.json').write_text(json.dumps(hashes,indent=2)+'\n')
(ROOT/'docs/S4_POSTCRT_FILTERS.md').write_text('''# S4c final post-CRT selection

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
''')

# Upgrade the existing entire-notebook numerical fixture. It must now supply
# ALL inputs of the final source panel, not silently omit the filters.
p = ROOT/'tests/test_standalone_figure_ports.py'
s = p.read_text()
old = "    output=root/'outputs'\n    _run_notebook(root,'Figure_S04.ipynb',{'outputs':{'root':str(output)}},"
assert s.count(old) == 1
new = '''    counts = []
    for day in ['d0','d1','d2','d3']:
        day_dir = model/day
        day_dir.mkdir(parents=True, exist_ok=True)
        pvalues = pd.DataFrame(
            [[0.001,0.002,0.8], [0.0001,0.02,0.8], [0.001,0.3,0.8], [1e-100]*3],
            index=['A','B','C','NTC'], columns=['P0','P1','P2'])
        pvalues.to_csv(day_dir/f'{day}_target_pvalues_skew_wide.tsv', sep='\\t')
        for sample in ['r1','r2']:
            for target in ['A','B','C']:
                counts.append({'day':day, 'biological_sample':sample, 'target_name':target,
                    'n_target_cells_with_cnmf_usage':19 if (sample,target)==('r2','B') else 20})
    diagnostic = model/'biological_replicate_effects/plots/DIAGNOSTIC_target_cell_counts_by_replicate.csv'
    diagnostic.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(counts).to_csv(diagnostic, index=False)
    output=root/'outputs'
    _run_notebook(root,'Figure_S04.ipynb',{'outputs':{'root':str(output)}},'''
s = s.replace(old, new)
old = "    pairs=pd.read_csv(c/'biological_replicate_DE_program_effect_pairs.csv')\n    assert len(pairs)==36\n    assert 'NTC' not in set(pairs.target_name)"
assert old in s
s = s.replace(old, '''    unfiltered = pd.read_csv(c/'S4c_unfiltered_real_effect_pairs.csv')
    pairs = pd.read_csv(next(c.glob('pooled_selected_replicate_effect_pairs_*.csv')))
    assert len(unfiltered)==36 and len(pairs)==12
    assert set(zip(pairs.target_name,pairs.program_name))=={('A','P0'),('A','P1'),('C','P0')}
    assert set(pairs.target_name)=={'A','C'}
    assert pairs.pooled_qvalue.le(0.10).all()
    assert pairs.n_target_cells_x.ge(20).all() and pairs.n_target_cells_y.ge(20).all()
    filters=pd.read_csv(c/'S4c_filter_counts.csv')
    assert filters.after_target_cell_filter.tolist()==[6]*4
    assert filters.selected_pairs.tolist()==[3]*4''')
p.write_text(s)

# Make the existing recomputation integration test verify the real H5MU
# metadata-to-count route, not only the resulting number of SVG files.
p=ROOT/'tests/test_s4_crt_runtime.py'
s=p.read_text()
s+='''
    panel_c=output/'Figure_S04/panel_c'
    provenance=json.loads((panel_c/'target_cell_count_provenance.json').read_text())
    assert provenance['mode']=='reconstructed_from_recorded_inputs'
    counts=pd.read_csv(panel_c/'DIAGNOSTIC_target_cell_counts_by_replicate.csv')
    assert counts.n_target_cells_with_cnmf_usage.eq(100).all()
    assert len(counts)==16
    selected=pd.read_csv(next(panel_c.glob('pooled_selected_replicate_effect_pairs_*.csv')))
    assert selected.pooled_qvalue.le(0.10).all()
    assert selected.n_target_cells_x.ge(20).all() and selected.n_target_cells_y.ge(20).all()
    assert set(selected.day)=={'d0','d1','d2','d3'}
'''
p.write_text(s)
for p in (ROOT/'src').rglob('*.py'):
    compile(p.read_text(),str(p),'exec')
for p in (ROOT/'tests').glob('*.py'):
    compile(p.read_text(),str(p),'exec')
print('S4 post-CRT notebook patch applied; integration fixtures require the final selection inputs.')
