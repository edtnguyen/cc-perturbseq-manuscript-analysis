# Motif-panel input manifest

The retained motif panels are Figure 3d, Figure 4e, and Figure S6. They share
the following motif-specific inputs.

## Selected files

| Role | Selected path | Selection basis |
|---|---|---|
| Day/program motif enrichment | `/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec/analysis/2.tf_enrichment_in_programs/enhancer_enrichment_results_min10.tsv` | Active non-archive result from 2026-10-21; the manuscript panels use the min-10-carrier analysis. |
| Motif-name to cluster mapping | `/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec/analysis/0.0.motif_id_annotation/tables/cluster_to_tf_motif_id_mapping.tsv` | Active non-archive mapping table from 2026-08-26. |
| Motif contribution-weight matrices | `/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec/results/compendium/modisco_compiled/modisco_compiled.h5` | Standard compiled compendium referenced by the recovered Figure 3/4/S6 notebooks. |

The later
`results/compendium/modisco_compiled_profile/modisco_compiled.h5` file is a
separate profile compendium. It is not treated as a renamed replacement for
the standard compiled H5.

## Deliberately excluded alternatives

The figure code rejects the following as silent substitutions:

- `analysis/archive_v1/2.tf_enrichment_in_programs/enhancer_enrichment_results.tsv`
- `analysis/archive_v1/2.tf_enrichment_in_programs_union_links/enhancer_enrichment_results.tsv`
- `analysis/archive_v1/2.tf_enrichment_in_programs_pos_only_2026-10-20/enhancer_enrichment_results.tsv`
- `analysis/archive_v1/2.tf_enrichment_in_programs_zero_filled_2026-10-21/enhancer_enrichment_results.tsv`
- `analysis/archive_v1/2.tf_enrichment_in_programs_carriers_min5_2026-10-21/enhancer_enrichment_results.tsv`
- archived sensitivity files `enhancer_enrichment_results_min1.tsv` and `enhancer_enrichment_results_min20.tsv`
- `analysis/archive_v1/tables/cluster_to_tf_motif_id_mapping.tsv`
- `results/compendium/modisco_compiled_profile/modisco_compiled.h5`

## Additional non-motif inputs for the composite panels

The complete regulator/loading/motif panels also use:

- `/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50.spectra.k_50.dt_0_2.consensus.txt`
- `/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/data/cNMF_50_0_2_032526_6co.xlsx`
- `/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/data/Programs_annotations.xlsx` for Figure S6.

The canonical paths are enforced in `src/cc_perturbseq/motif_sources.py`.
Figure 3d and Figure 4e audit all three configured motif paths. Figure S6
patches these exact paths into the frozen recovered source notebook before
execution.
