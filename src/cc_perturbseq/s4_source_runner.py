"""Run the exact prior Figure S4b/c CRT notebook and collect its outputs."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import os
from pathlib import Path
import shutil
from typing import Iterable

from .source_notebook_runner import git_blob_sha1_file, sha256_file


SOURCE_ENVIRONMENT_VARIABLE = "CCPERTURBSEQ_FIGURE_S04_SOURCE_NOTEBOOK"
REUSE_ENVIRONMENT_VARIABLE = "CCPERTURBSEQ_FIGURE_S04_REUSE_SOURCE_OUTPUTS"

OAK_SOURCE_NOTEBOOK = Path(
    "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/"
    "iPSC_EC_QQ_full_covariates_with_replicate_correlation_bundle/"
    "iPSC_EC_QQ_full_covariates_with_replicate_correlation.ipynb"
)
OAK_SOURCE_SHA256 = (
    "6e80625abd9b5e71d45da6dd19e888898318d39a159a7ac35600fc1d85489e87"
)
ARCHIVE_SOURCE_RELATIVE = Path(
    "recovery/raw_downloads/"
    "iPSC_EC_QQ_full_covariates_with_replicate_correlation_bundle/"
    "iPSC_EC_QQ_full_covariates_with_replicate_correlation.ipynb"
)
ARCHIVE_SOURCE_GIT_BLOB_SHA1 = "04e5c97707b847e09a89f7ca7b4ca2c83bd3b015"

MODEL_OUTPUT_RELATIVE = Path(
    "results/ipsc_ec_qq_exports_full_covariates/k_50.dt_0_2"
)


@dataclass(frozen=True)
class S4SourceNotebook:
    path: str
    source_kind: str
    verification_algorithm: str
    expected_hash: str
    observed_hash: str


def _verify_source(path: Path) -> S4SourceNotebook:
    path = path.expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Figure S4 source notebook is missing or empty: {path}")

    observed_sha256 = sha256_file(path)
    if observed_sha256 == OAK_SOURCE_SHA256:
        return S4SourceNotebook(
            path=str(path),
            source_kind="executed_sherlock_source",
            verification_algorithm="sha256",
            expected_hash=OAK_SOURCE_SHA256,
            observed_hash=observed_sha256,
        )

    observed_git_sha1 = git_blob_sha1_file(path)
    if observed_git_sha1 == ARCHIVE_SOURCE_GIT_BLOB_SHA1:
        return S4SourceNotebook(
            path=str(path),
            source_kind="archived_code_source",
            verification_algorithm="git_blob_sha1",
            expected_hash=ARCHIVE_SOURCE_GIT_BLOB_SHA1,
            observed_hash=observed_git_sha1,
        )

    raise ValueError(
        "Figure S4 source notebook does not match either frozen prior version. "
        f"Path: {path}; SHA-256: {observed_sha256}; Git blob SHA-1: {observed_git_sha1}."
    )


def resolve_s4_source_notebook(
    repo_root: str | Path,
    explicit_path: str | Path | None = None,
) -> S4SourceNotebook:
    """Resolve only the exact previous notebook; never substitute a generic plotter."""
    repo_root = Path(repo_root).expanduser().resolve()
    explicit = explicit_path or os.environ.get(SOURCE_ENVIRONMENT_VARIABLE)
    if explicit:
        return _verify_source(Path(explicit))

    candidates = [OAK_SOURCE_NOTEBOOK, repo_root / ARCHIVE_SOURCE_RELATIVE]
    failures: list[str] = []
    for candidate in candidates:
        if not candidate.is_file() or candidate.stat().st_size == 0:
            failures.append(f"missing: {candidate}")
            continue
        try:
            return _verify_source(candidate)
        except Exception as exc:
            failures.append(f"rejected: {candidate}: {exc}")

    raise FileNotFoundError(
        "The exact prior Figure S4b/c source notebook could not be resolved. "
        f"Set {SOURCE_ENVIRONMENT_VARIABLE} to the frozen notebook. Checked: "
        + " | ".join(failures)
    )


def source_model_output_root(source: S4SourceNotebook | str | Path) -> Path:
    path = Path(source.path if isinstance(source, S4SourceNotebook) else source)
    return path.expanduser().resolve().parent / MODEL_OUTPUT_RELATIVE


def _copy_files(paths: Iterable[Path], destination: Path) -> dict[str, str]:
    destination.mkdir(parents=True, exist_ok=True)
    copied: dict[str, str] = {}
    for path in sorted({Path(p).resolve() for p in paths}):
        if not path.is_file() or path.stat().st_size == 0:
            continue
        target = destination / path.name
        shutil.copy2(path, target)
        key = target.stem
        suffix = target.suffix.lstrip(".") or "file"
        candidate = f"{key}_{suffix}"
        index = 2
        while candidate in copied:
            candidate = f"{key}_{suffix}_{index}"
            index += 1
        copied[candidate] = str(target)
    return copied


def collect_s4_source_outputs(
    source: S4SourceNotebook,
    *,
    panel_b_directory: str | Path,
    panel_c_directory: str | Path,
) -> dict[str, dict[str, str]]:
    """Copy the exact source notebook's S4b/c files into the figure output tree."""
    model_root = source_model_output_root(source)
    if not model_root.is_dir():
        raise FileNotFoundError(
            "The source notebook did not produce its expected model-output directory: "
            f"{model_root}"
        )

    panel_b_source = model_root / "illustrator_editable_full_covariates"
    b_paths = list(panel_b_source.glob("d[0-3]_QQ_full_covariates.*"))
    for name in ("all_days_QQ_points.csv", "QQ_export_manifest.csv"):
        b_paths.append(model_root / name)
    for day in ("d0", "d1", "d2", "d3"):
        b_paths.extend((model_root / day).glob(f"{day}_*QQ_points.csv"))
        b_paths.append(model_root / day / f"{day}_run_config.json")

    replicate_root = model_root / "biological_replicate_effects"
    panel_c_source = replicate_root / "plots"
    c_paths = list(
        panel_c_source.glob("*biological_replicate_DE_program_effect_correlation*")
    )
    c_paths.extend(
        [
            panel_c_source / "biological_replicate_DE_program_effect_pairs.csv",
            panel_c_source / "biological_replicate_DE_program_effect_correlation_summary.csv",
            replicate_root / "replicate_beta_manifest.csv",
        ]
    )

    copied_b = _copy_files(b_paths, Path(panel_b_directory).expanduser().resolve())
    copied_c = _copy_files(c_paths, Path(panel_c_directory).expanduser().resolve())

    b_vector = [path for path in copied_b.values() if Path(path).suffix in {".pdf", ".svg", ".eps"}]
    c_vector = [path for path in copied_c.values() if Path(path).suffix in {".pdf", ".svg", ".eps"}]
    if not b_vector:
        raise RuntimeError(
            f"No prior-notebook vector output was collected for Figure S4b from {panel_b_source}."
        )
    if not c_vector:
        raise RuntimeError(
            f"No prior-notebook vector output was collected for Figure S4c from {panel_c_source}."
        )
    if not any(Path(path).name == "all_days_QQ_points.csv" for path in copied_b.values()):
        raise RuntimeError("Figure S4b source-data table all_days_QQ_points.csv was not collected.")
    if not any("effect_pairs" in Path(path).name for path in copied_c.values()):
        raise RuntimeError("Figure S4c matched target-program effect-pair table was not collected.")

    return {
        "panel_b": copied_b,
        "panel_c": copied_c,
        "source": {
            **{key: str(value) for key, value in asdict(source).items()},
            "model_output_root": str(model_root),
        },
    }


__all__ = [
    "ARCHIVE_SOURCE_GIT_BLOB_SHA1",
    "ARCHIVE_SOURCE_RELATIVE",
    "MODEL_OUTPUT_RELATIVE",
    "OAK_SOURCE_NOTEBOOK",
    "OAK_SOURCE_SHA256",
    "REUSE_ENVIRONMENT_VARIABLE",
    "SOURCE_ENVIRONMENT_VARIABLE",
    "S4SourceNotebook",
    "collect_s4_source_outputs",
    "resolve_s4_source_notebook",
    "source_model_output_root",
]
