# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure_1/Figure1_CC_Perturbseq_plot_input_package_paper_panels/CC_Figure1_reprocess_and_plot.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure_1/Figure1_CC_Perturbseq_plot_input_package_paper_panels/CC_Figure1_reprocess_and_plot.ipynb`
- SHA-256: `920b638c310c713d92d93af8a77956651fe594df38f322f1dbe7a96f6def05e1`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# Figure 1 CC-Perturb-seq capture benchmarking

This notebook reprocesses the Figure 1 capture-comparison datasets on Sherlock and writes an uploadable plot-input bundle. Send back the generated `figure1_plot_input_bundle.zip` if you want me to make the final plots from the Sherlock-processed data.

The current package includes the original Figure 1 panels plus paper-inspired QC/supplemental candidates: guide assignment classes, per-guide recovery curves, cell throughput, transcriptome QC summaries, guide saturation, primer optimization, guide-read source fractions, replicate correlations, bulk-vs-single-cell guide recovery, and optional target-knockdown QC.

The upload bundle contains derived TSV/TSV.GZ tables only, not raw MatrixMarket, h5ad, h5mu, or FASTQ files. It also can generate Illustrator-editable SVG panels locally with `plotnine`.


## Markdown cell 1

## Remaining inputs to fill before final plotting

The teloHAEC h5ad/kite paths and the iPSC pilot guide annotation paths are filled in from the paths you provided. Remaining optional inputs:

1. **Transcriptome downsampling outputs** for circular capture and CROP-seq in hiPSCs, if panel 1d must be a true GEX saturation curve rather than the transcriptome UMI/QC fallback.
2. **K562 or other additional-cell-type paths** if panel 1e should include more than teloHAEC and the iPSC pilot.
3. **Optional QC tables** for primer optimization, oligo-dT versus guide-RT primer read fractions, guide-count replicate correlations, bulk guide counts, and target-knockdown QC. These are configured under `cfg["optional_panels"]`.

Note: the path you provided for `CROP_GEX1.h5ad` points to one file inside the teloHAEC folder. The config uses its parent directory, `/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/tony_results`, as `root_dir`.


## Code cell 2

```python
from pathlib import Path
import json
import pandas as pd

# Make the helper module importable when the notebook is run from its own directory.
import sys
sys.path.insert(0, str(Path.cwd()))

from figure1_cc_perturbseq_pipeline import (
    default_config,
    check_config_paths,
    run_reprocessing,
    export_plot_input_bundle,
    load_outputs_from_dir,
    generate_figure1_panels,
    save_plot,
    plot_guide_umis_distribution,
    plot_assignment_metrics,
    plot_transcriptome_saturation,
    plot_transcriptome_qc_distribution,
    plot_other_celltype_guide_umis,
    plot_assignment_classes_stacked,
    plot_guide_detection_curve,
    plot_cell_throughput,
    plot_primer_optimization,
    plot_read_origin_fraction,
    plot_replicate_correlation_heatmap,
    plot_bulk_singlecell_correlation,
    plot_target_knockdown,
    write_sceptre_inputs_from_parse,
    SCEPTRE_R_TEMPLATE,
)

pd.set_option("display.max_colwidth", 160)

```

## Markdown cell 3

## 1. Configure datasets

The hiPSC single-guide, teloHAEC, iPSC pilot, teloHAEC guide-saturation, and iPSC pilot guide-annotation paths are prefilled. Edit only if files move or if you want to disable a dataset.


## Code cell 4

```python
cfg = default_config()

# Recommended: write all outputs into a dedicated project folder on Sherlock.
cfg["out_dir"] = "/oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/figure1_cc_capture"

# --- Dataset 1: single-guide infection in iPSC ---
cfg["datasets"]["ipsc_single_cc"]["enabled"] = True
cfg["datasets"]["ipsc_single_crop"]["enabled"] = True

# --- Dataset 2: teloHAEC circular capture vs CROP-seq ---
# You provided a path to CROP_GEX1.h5ad; the analysis needs the parent directory containing all h5ad files.
telohaec_root = "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/tony_results"
for key in ["telohaec_cc", "telohaec_crop"]:
    cfg["datasets"][key]["root_dir"] = telohaec_root
    cfg["datasets"][key]["enabled"] = True

# --- Dataset 3: iPSC ~10-gene pilot ---
cfg["datasets"]["ipsc_pilot_cc"]["enabled"] = True
cfg["datasets"]["ipsc_pilot_cc"]["guide_target_table"] = (
    "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure_1/eTN3_gRNA_group_table.csv"
)
cfg["datasets"]["ipsc_pilot_cc"]["guide_group_pairs"] = (
    "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure_1/eTN3_gRNA_group_pairs.csv"
)

# --- teloHAEC guide saturation curves ---
# These are the GUIDE_S1.*_kb_kite and CROP_GUIDE1_S4.*_kb_kite folders from splitseq_analysis.ipynb.
cfg["saturation"]["enabled"] = True
cfg["saturation"]["read_counts"] = [1000000, 3000000, 10000000, 30000000]
cfg["saturation"]["curves"] = [
    {
        "dataset_id": "telohaec_capture_comparison",
        "sample_id": "teloHAEC_circular_capture_guide_saturation",
        "method": "Circular capture",
        "cell_type": "teloHAEC",
        "modality": "guide",
        "pattern": f"{telohaec_root}/GUIDE_S1.{{reads}}_kb_kite/counts_unfiltered_modified/adata.h5ad",
        "cells_from_sample_id": "teloHAEC_circular_capture",
    },
    {
        "dataset_id": "telohaec_capture_comparison",
        "sample_id": "teloHAEC_CROP_seq_guide_saturation",
        "method": "CROP-seq",
        "cell_type": "teloHAEC",
        "modality": "guide",
        "pattern": f"{telohaec_root}/CROP_GUIDE1_S4.{{reads}}_kb_kite/counts_unfiltered_modified/adata.h5ad",
        "cells_from_sample_id": "teloHAEC_CROP_seq",
    },
]

# --- Upload bundle settings ---
# The full panel-specific derived tables are exported; the sampled preview tables are bounded by this setting.
cfg["upload_bundle"]["max_cells_per_sample"] = 50000
cfg["upload_bundle"]["include_full_cell_metrics"] = False
cfg["upload_bundle"]["histogram_bins"] = 80

Path(cfg["out_dir"]).mkdir(parents=True, exist_ok=True)
with open(Path(cfg["out_dir"]) / "figure1_config_used.json", "w") as handle:
    json.dump(cfg, handle, indent=2)

cfg


```

## Markdown cell 5

## 2. Check paths before reprocessing

Run this cell first on Sherlock. Disabled datasets are shown as disabled; missing enabled paths need to be fixed before final outputs can be generated.

## Code cell 6

```python
path_check = check_config_paths(cfg)
path_check
```

## Markdown cell 7

## 3. Reprocess count matrices into common metrics

This step standardizes the Parse-style and h5ad/h5mu inputs into shared per-cell, per-guide, per-sample, and saturation tables.

## Code cell 8

```python
outputs = run_reprocessing(cfg)
outputs["path_check"] = path_check

print("Tables written to:", cfg["out_dir"])
for name, df in outputs.items():
    print(f"{name}: {df.shape}")

outputs["sample_summary"]
```

## Markdown cell 9

## 4. Create uploadable plot-input bundle

This is the main handoff file to send back for plotting. It contains derived per-cell metrics, exact summary tables, guide-level detection metrics, saturation rows if available, and diagnostics. It does not contain raw matrices.

## Code cell 10

```python
bundle_zip = export_plot_input_bundle(
    outputs,
    cfg,
    # Default is 50,000 cells per panel/sample stratum. Lower this if the ZIP is too large to upload.
    max_rows_per_panel_stratum=cfg.get("upload_bundle", {}).get("max_cells_per_sample", 50000),
)

print("Upload this file back to ChatGPT:")
print(bundle_zip)

# Optional: inspect the manifest table sizes.
manifest_path = Path(cfg["out_dir"]) / "plot_input_bundle" / "manifest.json"
with open(manifest_path) as handle:
    manifest = json.load(handle)
manifest["table_shapes"]


```

## Markdown cell 11

### Optional: verify the upload bundle can be loaded

## Code cell 12

```python
# Optional verification
# reloaded = load_plot_input_bundle(bundle_zip)
# {name: df.shape for name, df in reloaded.items()}
```

## Markdown cell 13

### Command-line equivalent on Sherlock

This route skips local SVG generation and just creates the plot-input ZIP.

## Code cell 14

```python
# python figure1_cc_perturbseq_pipeline.py --config-json figure1_config_template.json --check-paths
# python figure1_cc_perturbseq_pipeline.py --config-json figure1_config_template.json --skip-plots
# Optional smaller upload: add --max-rows-per-panel-stratum 10000
# Optional later plotting from the bundle: python figure1_cc_perturbseq_pipeline.py --from-bundle /path/to/figure1_plot_input_bundle.zip

```

## Markdown cell 15

## 5. Generate Figure 1 SVG panels locally, optional

The SVG export keeps text editable in Illustrator. If true transcriptome saturation files are absent, the notebook writes `fig1d_transcriptome_umis_hiPSC.svg` as a fallback QC panel.

## Code cell 16

```python
manifest = generate_figure1_panels(outputs, Path(cfg["out_dir"]) / "svg")
manifest
```

## Markdown cell 17

## 6. Optional manual panel regeneration

Use these cells if you want to tune panel dimensions, titles, or which metric is plotted without rerunning the reprocessing step.

## Code cell 18

```python
cell_metrics = outputs["cell_metrics"]
summary = outputs["sample_summary"]
saturation = outputs["saturation"]
svg_dir = Path(cfg["out_dir"]) / "svg"
svg_dir.mkdir(parents=True, exist_ok=True)

# Panel b: guide UMIs/cell in hiPSC single-guide comparison.
if cell_metrics.empty:
    print("No cell_metrics table yet. Fix enabled paths and rerun Section 3.")
else:
    p1b = plot_guide_umis_distribution(
        cell_metrics,
        dataset_id="ipsc_single_guide",
        cell_type="hiPSC",
        y_col="guide_umis_top",
        title="hiPSC guide capture",
    )
    save_plot(p1b, svg_dir / "fig1b_guide_umis_hiPSC.svg", width=2.25, height=2.0)
    p1b
```

## Code cell 19

```python
# Panel c: assignment and guide detection metrics.
if summary.empty:
    print("No sample_summary table yet. Fix enabled paths and rerun Section 3.")
else:
    p1c = plot_assignment_metrics(summary, dataset_id="ipsc_single_guide", title="Guide assignment")
    save_plot(p1c, svg_dir / "fig1c_assignment_metrics.svg", width=3.0, height=2.2)
    p1c
```

## Code cell 20

```python
# Panel d: true transcriptome saturation if supplied; otherwise fallback transcriptome UMI distribution.
if cell_metrics.empty and saturation.empty:
    print("No cell_metrics or saturation table yet. Fix enabled paths and rerun Section 3.")
elif not saturation.empty and (saturation["modality"].str.lower().isin(["gex", "transcriptome", "mrna"]).any()):
    p1d = plot_transcriptome_saturation(saturation, dataset_id="ipsc_single_guide", title="Transcriptome saturation")
    save_plot(p1d, svg_dir / "fig1d_transcriptome_saturation_hiPSC.svg", width=2.5, height=2.0)
    p1d
else:
    p1d = plot_transcriptome_qc_distribution(cell_metrics, dataset_id="ipsc_single_guide", title="Transcriptome complexity")
    save_plot(p1d, svg_dir / "fig1d_transcriptome_umis_hiPSC.svg", width=2.25, height=2.0)
    p1d
```

## Code cell 21

```python
# Panel e: guide capture in additional cell types/datasets.
if cell_metrics.empty or (cell_metrics["dataset_id"] != "ipsc_single_guide").sum() == 0:
    print("No additional-cell-type data yet. Enable teloHAEC and/or confirm iPSC pilot paths.")
else:
    p1e = plot_other_celltype_guide_umis(
        cell_metrics,
        exclude_dataset_id="ipsc_single_guide",
        y_col="guide_umis_top",
        title="Guide capture across cell types",
    )
    save_plot(p1e, svg_dir / "fig1e_guide_umis_other_celltypes.svg", width=3.0, height=2.2)
    p1e
```

## Markdown cell 22

## 7. Optional SCEPTRE export for the iPSC pilot

This section writes MatrixMarket files and metadata for SCEPTRE from the Parse-style iPSC pilot. It is not required for the Figure 1 capture metrics unless you want SCEPTRE-based assignment or knockdown results.

## Code cell 23

```python
pilot = cfg["datasets"]["ipsc_pilot_cc"]
sceptre_out = Path(cfg["out_dir"]) / "sceptre_ipsc_pilot"

guide_target_table = pilot.get("guide_target_table")
guide_group_pairs = pilot.get("guide_group_pairs")  # available for downstream checks; not required by write_sceptre_inputs_from_parse()

print("Guide target table:", guide_target_table)
print("Guide group pairs:", guide_group_pairs)

# Uncomment to export SCEPTRE inputs.
# write_sceptre_inputs_from_parse(
#     gex_dir=pilot["gex_dir"],
#     guide_dir=pilot["guide_dir"],
#     out_dir=sceptre_out,
#     guide_target_table=guide_target_table,
#     guide_id_col="grna_id",
#     target_col="grna_target",
#     qc_filter=cfg["qc_filter_parse"],
# )
# print("SCEPTRE input written to", sceptre_out)


```

## Code cell 24

```python
# Write an R script for running SCEPTRE after the export cell above has been run.
r_script = SCEPTRE_R_TEMPLATE.format(data_dir=str(sceptre_out), p_mito_threshold=0.22)
Path(sceptre_out).mkdir(parents=True, exist_ok=True)
with open(sceptre_out / "run_sceptre.R", "w") as handle:
    handle.write(r_script)
print(sceptre_out / "run_sceptre.R")
```

