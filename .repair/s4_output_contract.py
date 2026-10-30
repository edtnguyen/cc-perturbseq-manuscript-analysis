"""programDE output contract and checks for reusing a completed full-day S4 export."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from . import crt_qq_backend as qqb


def run_program_crt(programde, *, inputs, B, n_jobs, calibrate_skew_normal):
    """Request raw/skew supplements only when programDE's skew fit is enabled.

    With calibration disabled, empirical p-values are already in pvals_df and
    betas_df remains available. Enabling a skew fit is NOT needed to obtain beta.
    """
    calibrated = bool(calibrate_skew_normal)
    return programde.run_all_genes_union_crt(
        inputs=inputs, B=B, n_jobs=n_jobs,
        calibrate_skew_normal=calibrated,
        return_raw_pvals=calibrated,
        return_skew_normal=calibrated,
    )


def extract_pvalue_outputs(out, *, calibrate_skew_normal):
    key = "pvals_raw_df" if calibrate_skew_normal else "pvals_df"
    if key not in out or out[key] is None:
        raise RuntimeError(f"programDE did not return required empirical output {key!r}.")
    skew = None
    if calibrate_skew_normal:
        skew = out.get("pvals_df")
        if skew is None:
            skew = out.get("pvals_skew_df")
    return out[key], skew


def audit_guide_metadata(args):
    """Compare the now-recognized metadata against the previous name inference.

    This audits every metadata guide, a conservative superset of guides used on
    an individual day. Missing control targets retain their inferred labels.
    """
    if not args.guide_metadata_tsv:
        return {"mode": "explicit_name_inference", "legacy_name_difference_count": 0}
    path = Path(args.guide_metadata_tsv).expanduser().resolve()
    table = qqb.read_table_auto(path)
    columns = {str(c).lower().replace(" ", "_"): c for c in table.columns}
    candidates = [args.guide_metadata_guide_col, "guide_name", "guide_id", "guide", "grna", "grna_name", "grna_id", "sgrna", "protospacer"]
    column = next((columns[str(c).lower().replace(" ", "_")] for c in candidates
                   if c and str(c).lower().replace(" ", "_") in columns), None)
    if column is None:
        raise KeyError(f"Cannot resolve guide column in {path}: {list(table.columns)}")
    guides = table[column].dropna().astype(str).tolist()
    before = qqb.infer_guide2gene_from_names(guides, args.ntc_label)
    mapped = qqb.load_guide2gene_metadata(path, guides, args.guide_metadata_guide_col,
                                        args.guide_metadata_target_col)
    after = dict(before)
    after.update(mapped)
    changes = [{"guide": g, "previous_inferred_target": before[g], "metadata_target": after[g]}
               for g in sorted(after) if after[g] != before[g]]
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "guide_column": str(column), "n_guides": len(after),
            "n_metadata_mapped": len(mapped),
            "n_name_inferred_including_controls": len(after) - len(mapped),
            "legacy_name_difference_count": len(changes),
            "legacy_name_differences": changes,
            "effective_mapping_sha256": hashlib.sha256(
                json.dumps(after, sort_keys=True).encode()).hexdigest()}


def validate_completed_full_day(args, current_mapping_audit):
    """Resume only the same completed export, never a different model or mapping.

    Old runtime receipts lacked mapping fingerprints. They are accepted only
    for the exact legacy v43 schema that produced the user's fallback warning,
    and only when the corrected mapping equals previous name inference.
    """
    root = Path(args.out_dir).expanduser().resolve()
    report_path = root / "s4_crt_runtime.json"
    if not report_path.is_file():
        raise ValueError("Cannot verify completed full-day inputs: runtime receipt missing. Use CRT_STAGE='all'.")
    previous = json.loads(report_path.read_text())
    if "argv" not in previous:
        raise ValueError("Completed full-day receipt lacks arguments. Use CRT_STAGE='all'.")
    old = qqb.parse_args(previous["argv"])
    ignored = {"n_jobs", "overwrite", "skip_if_exists", "no_plot", "save_original_qq", "dry_run"}
    differences = {k: [getattr(old, k, None), v] for k, v in vars(args).items()
                   if k not in ignored and getattr(old, k, None) != v}
    if differences:
        raise ValueError(f"Full-day scientific inputs/settings changed; use CRT_STAGE='all': {differences}")
    old_audit = previous.get("guide_metadata_audit")
    if old_audit is not None:
        if old_audit != current_mapping_audit:
            raise ValueError("Guide metadata changed since full-day CRT; use CRT_STAGE='all'.")
    elif args.guide_metadata_tsv:
        columns = set(qqb.read_table_auto(Path(args.guide_metadata_tsv)).columns)
        legacy_target_aliases = {"target_name", "target", "gene", "gene_symbol", "symbol", "grna_target"}
        known_legacy_fallback = (old.guide_metadata_target_col is None and
                                "intended_target_name" in columns and
                                not columns.intersection(legacy_target_aliases))
        if not known_legacy_fallback or current_mapping_audit["legacy_name_difference_count"]:
            raise ValueError("Cannot establish unchanged guide mapping for the old full-day export. Use CRT_STAGE='all'; inspect guide_metadata_audit.")
    model = root / qqb.model_tag(args.k, qqb.normalize_density_token(args.density_threshold))
    required = [model / "all_days_QQ_points.csv", model / "QQ_export_manifest.csv"]
    records = []
    for day in args.days:
        directory = model / day
        config_path = directory / f"{day}_run_config.json"
        if not config_path.is_file():
            raise ValueError(f"Completed full-day configuration missing: {config_path}. Use CRT_STAGE='all'.")
        config = json.loads(config_path.read_text())
        expected = {"day": day, "k": args.k, "B": args.B,
                    "density_threshold_token": qqb.normalize_density_token(args.density_threshold),
                    "ntc_labels": args.ntc_label}
        for key, value in expected.items():
            if config.get(key) != value:
                raise ValueError(f"Cannot resume {day}: {key} changed; use CRT_STAGE='all'.")
        for key in ("h5mu_path", "usage_path", "guide_long_csv"):
            p = Path(config[key]).expanduser()
            if not p.is_file() or p.stat().st_mtime_ns > config_path.stat().st_mtime_ns:
                raise ValueError(f"Input missing or newer than full-day completion: {p}. Use CRT_STAGE='all'.")
        if args.guide_metadata_tsv:
            p = Path(args.guide_metadata_tsv).expanduser()
            if p.stat().st_mtime_ns > config_path.stat().st_mtime_ns:
                raise ValueError(f"Metadata newer than full-day completion: {p}. Use CRT_STAGE='all'.")
        for field in ("combined_QQ_points", "target_empirical_wide", "ntc_empirical_wide"):
            required.append(directory / config["pvalue_outputs"][field])
        records.append({"day": day, "configuration": str(config_path)})
    for p in required:
        if not p.is_file() or p.stat().st_size == 0:
            raise ValueError(f"Completed full-day output missing/empty: {p}. Use CRT_STAGE='all'.")
    manifest = qqb.pd.read_csv(model / "QQ_export_manifest.csv")
    completed = set(manifest.loc[manifest["status"].eq("completed"), "day"].astype(str))
    if not set(args.days).issubset(completed):
        raise ValueError("Full-day manifest does not mark all requested days completed. Use CRT_STAGE='all'.")
    return {"status": "verified_completed_full_day", "model_root": str(model), "days": records,
            "input_verification": "matching recorded parameters/paths and unmodified timestamps; metadata content/mapping audit"}
