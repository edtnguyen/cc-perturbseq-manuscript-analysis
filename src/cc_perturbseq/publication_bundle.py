"""Maintainer-only collection of exact figure companions into repository data.

Institutional source paths are read ONLY by this opt-in preparation tool, not
by any figure notebook. Scientific matrices are neither fitted nor modified.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import shutil
import tempfile

from .igvf_download import IgvfClient, digest_file

DAYS = ("d0", "d1", "d2", "d3")
COUNTS_MAP = {"D01": "D0R1", "D02": "D0R2", "D3neg1": "D3NR1", "D3neg2": "D3NR2", "D3pos1": "D3PR1", "D3pos2": "D3PR2"}


def _valid(path):
    return Path(path).is_file() and Path(path).stat().st_size > 0


def _read_obs_index(path, *, modality=None):
    import h5py
    import pandas as pd
    with h5py.File(path, "r") as handle:
        group = handle[f"mod/{modality}/obs"] if modality else handle["obs"]
        name = group.attrs.get("_index", "_index")
        if isinstance(name, bytes):
            name = name.decode()
        ids = group[name].asstr()[()]
    index = pd.Index(ids.astype(str))
    if not index.is_unique:
        raise ValueError(f"Source cell ordering contains duplicate identifiers: {path}")
    return index


def validation_bridge(original_path, portal_path):
    """Keep original labels and QC-only fields; all six plotted counts come from IGVF.

Mapping requires identical numeric count vectors and an identical analysis-gene
label (allowing the documented NEDD8 alias), then an exact normalized guide ID
when counts are not unique. Ambiguous identities are rejected, not guessed.
"""
    import numpy as np
    import pandas as pd
    from .figure1_panel_c import _read_count_table, _infer_column
    original = _read_count_table(original_path)
    deposited = pd.read_csv(portal_path, sep="\t")
    guide = _infer_column(original, ["sgRNA", "guide", "guide_id", "grna"])
    gene = _infer_column(original, ["Gene", "target", "target_gene", "gene_symbol"])
    required = list(COUNTS_MAP.values()) + ["D3TOTAL1", "D3TOTAL2"]
    if set(required) - set(original):
        raise ValueError("Original validation counts lack the QC totals required for the companion table.")
    if len(original) != len(deposited) or deposited.guide_id.duplicated().any():
        raise ValueError("Validation source/deposit guide counts or identities differ.")
    def normalized_gene(v):
        text = str(v).strip().upper().replace("_", "-")
        if "NON-TARGET" in text:
            return "NON-TARGETING"
        if "SAFE" in text:
            return "SAFE"
        return "NEDD8" if text == "NEDD8-MDP1" else text
    def signature(label, counts):
        a = np.asarray(counts, dtype=float)
        if not np.isfinite(a).all() or (a < 0).any():
            raise ValueError("Invalid validation count vector.")
        return normalized_gene(label), tuple(a)
    lookup = defaultdict(list)
    for _, row in deposited.iterrows():
        lookup[signature(row["gene"], row[list(COUNTS_MAP)])].append(str(row["guide_id"]))
    assigned = []
    for _, row in original.iterrows():
        options = lookup[signature(row[gene], row[list(COUNTS_MAP.values())])]
        exact = [x for x in options if x.strip() == str(row[guide]).strip()]
        chosen = exact[0] if len(exact) == 1 else options[0] if len(options) == 1 else None
        if chosen is None:
            raise ValueError(f"Cannot uniquely reconcile validation guide {row[guide]!r} with deposited counts. No mapping was guessed.")
        assigned.append(chosen)
        options.remove(chosen)
    result = pd.DataFrame({"portal_guide_id": assigned, "sgRNA": original[guide].astype(str),
        "Gene": original[gene].astype(str), "D3TOTAL1": original["D3TOTAL1"], "D3TOTAL2": original["D3TOTAL2"]})
    aligned = deposited.set_index("guide_id").loc[assigned]
    for portal, local in COUNTS_MAP.items():
        np.testing.assert_array_equal(aligned[portal].to_numpy(float), original[local].to_numpy(float), err_msg=local)
    return result


def _beta_sources(model_root):
    import pandas as pd
    model_root = Path(model_root)
    path = model_root / "biological_replicate_effects/replicate_beta_manifest.csv"
    table = pd.read_csv(path)
    table = table.loc[table["status"].isin(["written", "existing"]) & table.beta_path.notna()].copy()
    if table.empty or set(table.day.astype(str)) != set(DAYS):
        raise ValueError("S4 replicate manifest must contain usable results for all four days.")
    results = []
    for i, row in table.iterrows():
        source = Path(str(row.beta_path))
        choices = [source] if source.is_absolute() else [path.parent / source, model_root / source]
        choices.append(path.parent / "betas_by_replicate" / str(row.day) / source.name)
        found = next((p for p in choices if _valid(p)), None)
        if found is None:
            raise FileNotFoundError(f"Missing S4 numerical beta file for {row.day}: {source.name}")
        relative = Path("betas_by_replicate") / str(row.day) / source.name
        table.loc[i, "beta_path"] = relative.as_posix()
        results.append((str(row.day), relative, found))
    # Other saved paths are audit fields, not runtime inputs. Do not publish
    # unusable absolute paths in the relocated numerical manifest.
    for col in table:
        if col != "beta_path" and ("path" in str(col).lower() or "dir" in str(col).lower()):
            table[col] = table[col].map(lambda x: Path(str(x)).name if pd.notna(x) else x)
    return table, results



def _geo_sources(config, *, hippo_h5mu=None, mixed_h5mu=None):
    configured = config.get("geo_plot_sources", {})
    sources = {
        "hippo_timecourse": hippo_h5mu or configured.get("hippo_timecourse"),
        "ipsc_k562_mix": mixed_h5mu or configured.get("ipsc_k562_mix"),
    }
    missing = []
    for dataset, source in sources.items():
        path = Path(source).expanduser().resolve() if source else None
        if path is None or not _valid(path):
            missing.append(f"{dataset}: {path or 'not configured'}")
        else:
            sources[dataset] = path
    if missing:
        raise FileNotFoundError("Submitted objects needed to prepare plotting tables:\n  " + "\n  ".join(missing)
            + "\nUse --hippo-h5mu and --mixed-h5mu to select the exact submitted integrated files. "
            "The original Split-pipe experiment directories are not needed.")
    return sources


def _build_geo_plots(sources, parent):
    from .geo_plot_data import build_plot_bundle
    reports = {}
    for dataset, source in sources.items():
        print(f"Generating {dataset} plotting inputs from {source.name}", flush=True)
        manifest = build_plot_bundle(source, Path(parent) / dataset, dataset)
        reports[dataset] = {"source_sha256": manifest["source"]["sha256"],
                            "n_tables": len(manifest["tables"])}
    return reports


def collect_geo_plots_only(root, source_config, *, hippo_h5mu=None, mixed_h5mu=None, apply=False, overwrite=False):
    """Add/update just the two independent plot bundles, preserving other data."""
    root = Path(root).expanduser().resolve()
    config = json.loads(Path(source_config).read_text())
    sources = _geo_sources(config, hippo_h5mu=hippo_h5mu, mixed_h5mu=mixed_h5mu)
    for dataset, source in sources.items():
        print(f"FOUND {dataset} submitted object: {source}")
    if not apply:
        return {"status": "dry_run", "datasets": list(sources), "raw_objects_uploaded": False}
    destination = root / "data/reproduction/geo_plots"
    if destination.exists() and not overwrite:
        raise FileExistsError("geo_plots already exists. Use --overwrite to back it up and replace it.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".geo-plots-", dir=destination.parent))
    try:
        reports = _build_geo_plots(sources, staging)
        if destination.exists():
            backup = root / "data" / ("geo_plots.backup-" + datetime.now().strftime("%Y%m%d-%H%M%S-%f"))
            destination.rename(backup)
            print("Previous plot bundles retained:", backup)
        staging.rename(destination)
        print("Plotting data saved to", destination, "; no Git changes were staged or pushed.")
        return {"status": "collected", "datasets": reports, "raw_objects_uploaded": False}
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def collect(root, source_config, *, apply=False, overwrite=False, s4_results=None, target_counts=None, hippo_h5mu=None, mixed_h5mu=None):
    import pandas as pd
    root = Path(root).expanduser().resolve()
    config = json.loads(Path(source_config).read_text())
    geo_sources = _geo_sources(config, hippo_h5mu=hippo_h5mu, mixed_h5mu=mixed_h5mu)
    model_root = Path(s4_results or config["s4_results_source"]).expanduser().resolve()
    target_counts = target_counts or config.get("s4_target_counts_source")
    sources = [(row["key"], Path(row["source"]).expanduser()) for row in config["files"]]
    sources += [("expression cell order", Path(config["integrated_cell_order_source"])),
                ("cNMF cell order", Path(config["cnmf_cell_order_source"])),
                ("validation compatibility", Path(config["validation_guide_counts_source"])),
                ("S4 Q–Q", model_root / "all_days_QQ_points.csv"),
                ("S4 replicate manifest", model_root / "biological_replicate_effects/replicate_beta_manifest.csv")]
    problems = []
    for name, path in sources:
        present = _valid(path)
        print(f"{'FOUND' if present else 'MISSING'} {name}: {path}")
        if not present:
            problems.append(f"{name}: {path}")
    pvalue_files = []
    for day in DAYS:
        found = [model_root / day / f"{day}_target_pvalues_{kind}_{shape}.tsv"
                 for kind in ("skew", "empirical") for shape in ("long", "wide")]
        found = [p for p in found if _valid(p)]
        if not found:
            problems.append(f"S4 pooled p-values for {day} under {model_root}")
        pvalue_files.extend(found)
    if problems:
        raise FileNotFoundError("Collection preflight found missing inputs:\n  " + "\n  ".join(problems))
    table, beta_sources = _beta_sources(model_root)
    if not apply:
        print(f"Dry run: {len(config['files'])} exact companion files, two cell-order tables, validation adapter, and S4 numerical data would be collected. Use --apply to copy. Both GEO plotting bundles will be regenerated; integrated H5MU objects and raw CRT inputs are not copied.")
        return {"status": "dry_run", "s4_model": str(model_root)}
    destination = root / "data/reproduction"
    if destination.exists() and any(p.name != "README.md" for p in destination.iterdir()) and not overwrite:
        raise FileExistsError("data/reproduction already exists. Verify it or rerun with --overwrite; the existing bundle will be backed up.")
    (root / "data").mkdir(exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".publication-", dir=root / "data"))
    files = []
    def stage_path(relative):
        rel = Path(relative)
        if rel.is_absolute() or ".." in rel.parts or not rel.is_relative_to("data/reproduction"):
            raise ValueError(f"Invalid publication destination {relative}")
        p = staging / rel.relative_to("data/reproduction")
        p.parent.mkdir(parents=True, exist_ok=True)
        return p
    def record(key, p, **extra):
        if not _valid(p):
            raise ValueError(f"Empty publication file: {key}")
        if p.stat().st_size >= 95 * 1024 * 1024:
            raise ValueError(f"{key} exceeds the 95 MiB publication guard. Use compression or configure Git LFS before publishing; no large file was staged to Git.")
        files.append({"key": key, "path": "data/reproduction/" + p.relative_to(staging).as_posix(),
            "size_bytes": p.stat().st_size, "sha256": digest_file(p, "sha256"), **extra})
    def copy(key, source, relative, *, compress=False):
        out = stage_path(relative)
        if compress:
            with Path(source).open("rb") as inp, out.open("wb") as raw:
                with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as target:
                    shutil.copyfileobj(inp, target, 4 * 1024 * 1024)
            with gzip.open(out, "rb") as stream:
                import hashlib
                check = hashlib.sha256()
                for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
                    check.update(block)
            if check.hexdigest() != digest_file(source, "sha256"):
                raise ValueError("Compressed copy is not lossless.")
        else:
            shutil.copy2(source, out)
            if digest_file(out, "sha256") != digest_file(source, "sha256"):
                raise ValueError("Copied file checksum differs from source.")
        record(key, out, source_basename=Path(source).name, source_sha256=digest_file(source, "sha256"), compression="gzip" if compress else "unchanged")
        return out
    try:
        if (destination / "README.md").is_file():
            shutil.copy2(destination / "README.md", staging / "README.md")
        for row in config["files"]:
            copy(row["key"], row["source"], row["destination"], compress=row.get("gzip", False))
        for kind, source, modality in (("expression", config["integrated_cell_order_source"], None), ("cnmf", config["cnmf_cell_order_source"], "cNMF")):
            ids = _read_obs_index(source, modality=modality)
            path = stage_path(f"data/reproduction/programs/{kind}_cell_order.tsv.gz")
            pd.DataFrame({"cell_id": ids}).to_csv(path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
            record("cell_order." + kind, path, n_cells=len(ids), source_basename=Path(source).name)
        # A small public deposit is used to prove the column/guide bridge.
        portal_manifest = json.loads((root / "config/public_data_manifest.json").read_text())
        validation = IgvfClient(allow_prompt=True).download(portal_manifest["files"]["IGVFFI5157STFU"], root / "data/igvf")
        bridge = validation_bridge(config["validation_guide_counts_source"], validation)
        path = stage_path("data/reproduction/screen/validation_guide_compatibility.tsv")
        bridge.to_csv(path, sep="\t", index=False)
        record("validation.compatibility", path, count_source="IGVFFI5157STFU", all_six_count_columns_verified_equal=True)
        prefix = "data/reproduction/s4/k_50.dt_0_2/"
        copy("s4.qq", model_root / "all_days_QQ_points.csv", prefix + "all_days_QQ_points.csv")
        out = stage_path(prefix + "biological_replicate_effects/replicate_beta_manifest.csv")
        table.to_csv(out, index=False)
        record("s4.manifest", out)
        for i, (day, relative, source) in enumerate(beta_sources):
            copy(f"s4.beta.{day}.{i}", source, prefix + "biological_replicate_effects/" + relative.as_posix())
        for p in pvalue_files:
            copy("s4.pvalues." + p.name, p, prefix + p.relative_to(model_root).as_posix())
        from .figure_s04_postcrt import load_or_rebuild_target_counts
        counts_out = stage_path(prefix + "biological_replicate_effects/plots/DIAGNOSTIC_target_cell_counts_by_replicate.csv")
        # Prefer the same-model saved numerical counts; rebuild on the source
        # machine only when absent. No CRT estimation is performed.
        saved = model_root / "biological_replicate_effects/plots/DIAGNOSTIC_target_cell_counts_by_replicate.csv"
        selected_counts = Path(target_counts) if target_counts else saved if _valid(saved) else None
        counts = load_or_rebuild_target_counts(model_root, DAYS, counts_out.parent, counts_csv=selected_counts)
        counts.to_csv(counts_out, index=False)
        record("s4.counts", counts_out)
        # Do not ship path-based raw-recomputation receipts as runtime inputs.
        for extra in counts_out.parent.glob("*.json"):
            extra.unlink()
        for extra in counts_out.parent.glob("DIAGNOSTIC_replicate_cnmf_usage_overlap.csv"):
            extra.unlink()
        # The two plotting bundles have their own manifests/checksums, so they
        # can subsequently be refreshed with --geo-plots-only without invalidating
        # the independent core companion manifest.
        geo_reports = _build_geo_plots(geo_sources, staging / "geo_plots")
        result = {"schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
            "files": files, "geo_datasets_included": False,
            "geo_plot_bundles": geo_reports,
            "cnmf_tsv_policy": "Exact original representations retained; no unverified substitution from H5MU loadings or usage.",
            "s4_policy": "Completed numerical results; >=20 target cells and pooled q<=0.10 selection unchanged; no raw CRT run."}
        (staging / "bundle_manifest.json").write_text(json.dumps(result, indent=2) + "\n")
        if destination.exists():
            backup = root / "data" / ("reproduction.backup-" + datetime.now().strftime("%Y%m%d-%H%M%S-%f"))
            destination.rename(backup)
            print("Previous bundle retained:", backup)
        staging.rename(destination)
        print(f"Collected {len(files)} verified files into {destination}. No Git changes were staged or pushed.")
        return result
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--source-config", type=Path)
    parser.add_argument("--s4-results", type=Path, help="Exact completed S4 model directory, not its plots subdirectory")
    parser.add_argument("--target-counts", type=Path, help="Same-model numerical target-cell support table")
    parser.add_argument("--apply", action="store_true", help="Copy and verify data; otherwise perform a dry run")
    parser.add_argument("--overwrite", action="store_true", help="Back up and replace an existing collected bundle")
    parser.add_argument("--hippo-h5mu", type=Path, help="Exact submitted Hippo integrated H5MU")
    parser.add_argument("--mixed-h5mu", type=Path, help="Exact submitted hiPSC-K562 integrated H5MU")
    parser.add_argument("--geo-plots-only", action="store_true", help="Collect only the two precomputed plotting bundles, leaving other companion data unchanged")
    args = parser.parse_args(argv)
    if args.geo_plots_only:
        return collect_geo_plots_only(args.repo, args.source_config or args.repo / "config/publication_source_files.json",
            hippo_h5mu=args.hippo_h5mu, mixed_h5mu=args.mixed_h5mu, apply=args.apply, overwrite=args.overwrite)
    return collect(args.repo, args.source_config or args.repo / "config/publication_source_files.json", apply=args.apply,
            overwrite=args.overwrite, s4_results=args.s4_results, target_counts=args.target_counts,
            hippo_h5mu=args.hippo_h5mu, mixed_h5mu=args.mixed_h5mu)


if __name__ == "__main__":
    main()
