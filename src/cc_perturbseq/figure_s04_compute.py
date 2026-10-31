"""Explicit opt-in CRT recomputation, ported from source cell 19."""
def compute_replicate_betas(
    *, programde_root, data_dir, cnmf_dir, guide_long_template,
    guide_metadata_tsv, output_dir, days=("d0", "d1", "d2", "d3"), n_jobs=4,
    k=50, density_threshold="0_2", ntc_labels=("SAFE", "non-targeting", "NTC"),
    overwrite=False, cnmf_file=None, usage_file_template=None,
    gex_modality=None, guide_modality=None,
    guide_metadata_guide_col=None, guide_metadata_target_col=None, guide_mapping_policy="source_names",
):
    PROGRAMDE_ROOT, DATA_DIR, CNMF_DIR = str(programde_root), str(data_dir), str(cnmf_dir)
    GUIDE_LONG_TEMPLATE = str(guide_long_template)
    GUIDE_METADATA_TSV = str(guide_metadata_tsv) if guide_metadata_tsv is not None else None
    OUT_DIR, DAYS, N_JOBS = str(output_dir), list(days), int(n_jobs)
    K, DENSITY_THRESHOLD, NTC_LABELS = int(k), str(density_threshold), list(ntc_labels)
    from pathlib import Path
    from types import SimpleNamespace
    import gc
    import importlib
    import json
    import time

    import muon as mu
    import numpy as np
    import pandas as pd
    from IPython.display import display

    from cc_perturbseq import crt_qq_backend as qqb
    from cc_perturbseq.figure_s04_contract import run_program_crt

    # For the replicate-correlation plot, beta estimates are the relevant output.
    # Keep B_REP small for a QC run; use a larger value only if you also need replicate p-values.
    B_REP = 127
    MIN_CELLS_PER_REPLICATE = 200
    SKIP_EXISTING_REP_BETAS = not overwrite

    # Output paths.
    density_token = qqb.normalize_density_token(DENSITY_THRESHOLD)
    MODEL_ROOT = Path(OUT_DIR) / qqb.model_tag(K, density_token)
    REPL_ROOT = MODEL_ROOT / "biological_replicate_effects"
    REPL_BETA_DIR = REPL_ROOT / "betas_by_replicate"
    REPL_BETA_DIR.mkdir(parents=True, exist_ok=True)

    programde = qqb.load_programde(PROGRAMDE_ROOT)

    args_for_helpers = SimpleNamespace(
        cnmf_file=cnmf_file,
        usage_file_template=usage_file_template,
        cnmf_dir=CNMF_DIR,
        k=K,
        gex_modality=gex_modality,
        guide_modality=guide_modality,
        ntc_label=NTC_LABELS,
        guide_metadata_tsv=GUIDE_METADATA_TSV,
        guide_metadata_guide_col=guide_metadata_guide_col,
        guide_metadata_target_col=guide_metadata_target_col,
        require_covariates=True,
    )

    replicate_rows = []

    for day in DAYS:
        print(f"\n=== {day}: loading full day object ===")

        h5mu_path = qqb.h5mu_path_for_day(Path(DATA_DIR), day)
        guide_csv = qqb.guide_long_path_for_day(GUIDE_LONG_TEMPLATE, day)

        mdata = mu.read(h5mu_path)

        gex_modality = qqb.resolve_modality(
            mdata,
            args_for_helpers.gex_modality,
            ["GEX", "gex", "rna", "RNA"],
            "GEX/RNA",
        )
        guide_modality = qqb.resolve_modality(
            mdata,
            args_for_helpers.guide_modality,
            ["guide", "Guide", "gRNA", "grna"],
            "guide",
        )

        adata = mdata[gex_modality].copy()

        usage_candidates = qqb.usage_candidates_from_args(
            args_for_helpers,
            day,
            density_token,
        )
        usage_path, matched_usage_cells = qqb.choose_usage_file(
            usage_candidates,
            pd.Index(adata.obs_names),
            day,
        )

        qqb.load_and_attach_usage(adata, usage_path)

        guide_names, guide2gene = qqb.load_and_attach_guides(
            adata,
            mdata,
            guide_modality,
            guide_csv,
            NTC_LABELS,
            args_for_helpers.guide_metadata_tsv,
            args_for_helpers.guide_metadata_guide_col,
            args_for_helpers.guide_metadata_target_col,
            mapping_policy=guide_mapping_policy,
        )

        if "biological_sample" not in adata.obs:
            raise KeyError(f"{day}: adata.obs does not contain 'biological_sample'.")

        biological_samples = sorted(adata.obs["biological_sample"].astype(str).unique())
        print(f"{day}: biological samples = {biological_samples}")

        for biological_sample in biological_samples:
            sample_safe = qqb.safe_name(biological_sample)
            day_rep_dir = REPL_BETA_DIR / day
            day_rep_dir.mkdir(parents=True, exist_ok=True)

            beta_path = day_rep_dir / f"{day}__{sample_safe}__target_betas_wide.tsv"
            meta_path = day_rep_dir / f"{day}__{sample_safe}__run_metadata.json"

            if SKIP_EXISTING_REP_BETAS and beta_path.is_file() and beta_path.stat().st_size > 0:
                print(f"  {day} / {biological_sample}: using existing {beta_path}")
                replicate_rows.append({
                    "day": day,
                    "biological_sample": biological_sample,
                    "status": "existing",
                    "n_cells": None,
                    "beta_path": str(beta_path),
                    "metadata_path": str(meta_path),
                })
                continue

            mask = adata.obs["biological_sample"].astype(str).eq(biological_sample).to_numpy()
            n_cells = int(mask.sum())

            if n_cells < MIN_CELLS_PER_REPLICATE:
                print(f"  {day} / {biological_sample}: skipped, only {n_cells:,} cells")
                replicate_rows.append({
                    "day": day,
                    "biological_sample": biological_sample,
                    "status": "skipped_too_few_cells",
                    "n_cells": n_cells,
                    "beta_path": None,
                    "metadata_path": None,
                })
                continue

            print(f"  {day} / {biological_sample}: running CRT on {n_cells:,} cells")

            adata_rep = adata[mask].copy()

            # biological_sample is constant within this subset and will be dropped by attach_covariates().
            covariates = qqb.attach_covariates(
                adata_rep,
                require_covariates=True,
            )

            inputs = qqb.prepare_inputs(programde, adata_rep)

            t0 = time.time()
            # betas_df is always returned. Raw/skew supplement flags must both
            # be False when skew calibration is disabled in real programDE.
            out = run_program_crt(
                programde, inputs=inputs, B=B_REP, n_jobs=N_JOBS,
                calibrate_skew_normal=False,
            )

            betas = out.get("betas_df")
            if betas is None:
                raise RuntimeError(f"{day} / {biological_sample}: no betas_df returned.")

            qqb.save_wide_frame(betas, beta_path)

            meta = {
                "day": day,
                "biological_sample": biological_sample,
                "n_cells": n_cells,
                "B_REP": B_REP,
                "usage_path": str(usage_path),
                "matched_usage_cells": int(matched_usage_cells),
                "h5mu_path": str(h5mu_path),
                "guide_csv": str(guide_csv),
                "n_guides": int(len(guide_names)),
                "n_targets": int(betas.shape[0]),
                "n_programs": int(betas.shape[1]),
                "covariate_columns": list(covariates.columns),
                "covariate_sources": adata_rep.uns.get("covariate_sources", {}),
                "runtime_seconds": round(time.time() - t0, 2),
                "beta_path": str(beta_path),
            }
            meta_path.write_text(json.dumps(meta, indent=2))

            replicate_rows.append({
                "day": day,
                "biological_sample": biological_sample,
                "status": "written",
                "n_cells": n_cells,
                "beta_path": str(beta_path),
                "metadata_path": str(meta_path),
            })

            del adata_rep, inputs, out, betas
            gc.collect()

        del mdata, adata
        gc.collect()

    replicate_manifest = pd.DataFrame(replicate_rows)
    replicate_manifest_path = REPL_ROOT / "replicate_beta_manifest.csv"
    replicate_manifest.to_csv(replicate_manifest_path, index=False)

    print("\nSaved manifest:", replicate_manifest_path)
    display(replicate_manifest)
    return replicate_manifest
