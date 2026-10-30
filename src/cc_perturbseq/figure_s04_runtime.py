"""S4 CRT command/runtime boundary; no archived notebook or R execution.

The archived S4 run used sc-dl Python. R-SCEPTRE in sceptre_env produced the
upstream guide-assignment CSVs; programDE's Python src.sceptre consumes them.
"""
from __future__ import annotations

from collections import deque
import importlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from . import crt_qq_backend as qqb
from .figure_s04_contract import audit_guide_metadata, validate_completed_full_day


def build_crt_argv(*, programde_root, data_dir, cnmf_dir, guide_long_template,
                   guide_metadata_tsv, output_dir, days=("d0", "d1", "d2", "d3"),
                   k=50, density_threshold="0_2", B=10000, n_jobs=4,
                   ntc_labels=("SAFE", "non-targeting", "NTC"), overwrite=True):
    """Build and parse the original S4 command before starting any computation."""
    argv = [
        "--programde-root", str(programde_root), "--data-dir", str(data_dir),
        "--cnmf-dir", str(cnmf_dir), "--guide-long-csv-template", str(guide_long_template),
        "--out-dir", str(output_dir), "--days", *map(str, days), "--k", str(k),
        "--density-threshold", str(density_threshold), "--B", str(B),
        "--n-jobs", str(n_jobs), "--require-covariates", "--no-plot",
        "--overwrite" if overwrite else "--skip-if-exists",
    ]
    if guide_metadata_tsv is not None:
        argv += ["--guide-metadata-tsv", str(guide_metadata_tsv)]
    for label in ntc_labels:
        argv += ["--ntc-label", str(label)]
    qqb.parse_args(argv)
    return argv


def check_crt_runtime(args):
    """Check the selected interpreter and source inputs, without running CRT."""
    versions = {}
    for name in ("numpy", "pandas", "scipy", "anndata", "mudata", "muon"):
        try:
            module = importlib.import_module(name)
        except ImportError as exc:
            raise RuntimeError(
                f"CRT interpreter {sys.executable} cannot import {name}. "
                "Use the original sc-dl Python via CRT_PYTHON; R-SCEPTRE is not "
                "run by S4. No CRT has started."
            ) from exc
        versions[name] = getattr(module, "__version__", "unknown")

    programde_root = Path(args.programde_root).expanduser().resolve()
    if not (programde_root / "src" / "sceptre").is_dir():
        raise FileNotFoundError(f"No programDE Python src/sceptre package: {programde_root}")
    qqb.load_programde(str(programde_root))
    package = importlib.import_module("src.sceptre")
    package_file = getattr(package, "__file__", None)
    if not package_file or not Path(package_file).resolve().is_relative_to(programde_root):
        raise RuntimeError(
            f"Wrong programDE import: {package_file}; expected beneath {programde_root}. "
            "The CRT worker must import the requested checkout. No CRT has started."
        )

    records, missing = [], []
    if args.guide_metadata_tsv:
        metadata = Path(args.guide_metadata_tsv).expanduser()
        if not metadata.is_file() or metadata.stat().st_size == 0:
            missing.append(f"guide metadata: {metadata}")
    density = qqb.normalize_density_token(args.density_threshold)
    for day in args.days:
        h5mu = qqb.h5mu_path_for_day(Path(args.data_dir).expanduser(), day)
        guides = qqb.guide_long_path_for_day(args.guide_long_csv_template, day)
        for label, path in (("H5MU", h5mu), ("SCEPTRE guide assignments", guides)):
            if not path.is_file() or path.stat().st_size == 0:
                missing.append(f"{day} {label}: {path}")
        candidates = qqb.usage_candidates_from_args(args, day, density)
        usable = [p for p in candidates if p.is_file() and p.stat().st_size > 0]
        if not usable:
            missing.append(f"{day} cNMF usage: no nonempty candidate under {args.cnmf_dir}")
        records.append({"day": day, "h5mu": str(h5mu), "guide_assignments": str(guides),
                        "usage_candidates": [str(p) for p in usable]})
    if missing:
        raise FileNotFoundError(
            "S4 CRT input preflight failed before computation:\n  " + "\n  ".join(missing)
            + "\nS4 consumes saved SCEPTRE assignments; it does not regenerate them or silently "
            "substitute another assignment file. Correct GUIDE_LONG_TEMPLATE for your actual exports."
        )
    return {"python_executable": sys.executable, "python_version": sys.version,
            "programde_module": str(Path(package_file).resolve()), "versions": versions,
            "sceptre_R_required": False, "inputs": records,
            "guide_metadata_audit": audit_guide_metadata(args),
            "usage_cell_alignment": "checked by original backend after reading each day"}


def run_worker(argv, *, preflight_only=False, stage="all"):
    """Run both original Python stages in the same interpreter after preflight."""
    if stage not in {"all", "replicates"}:
        raise ValueError("stage must be all or replicates")
    args = qqb.parse_args(list(argv))
    if args.dry_run:
        raise ValueError("Use preflight_only=True for this worker, not backend --dry-run.")
    report = check_crt_runtime(args)
    resume = validate_completed_full_day(args, report["guide_metadata_audit"]) if stage == "replicates" else None
    out = Path(args.out_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    report_path = out / "s4_crt_runtime.json"
    report.update({"status": "preflight_passed", "argv": list(argv),
                   "B": args.B, "B_REP": 127, "stage": stage, "full_day_resume": resume})
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    if preflight_only:
        return report
    try:
        if stage == "all":
            print("S4 full-day program CRT starting", flush=True)
            code = qqb.main(list(argv))
            if code != 0:
                raise RuntimeError(f"Full-day CRT returned exit status {code}")
        else:
            print("Reusing verified completed full-day export; running replicates only.", flush=True)
        from .figure_s04_compute import compute_replicate_betas
        print("S4 replicate beta estimation starting", flush=True)
        compute_replicate_betas(
            programde_root=args.programde_root, data_dir=args.data_dir, cnmf_dir=args.cnmf_dir,
            guide_long_template=args.guide_long_csv_template,
            guide_metadata_tsv=args.guide_metadata_tsv, output_dir=args.out_dir,
            days=args.days, n_jobs=args.n_jobs, k=args.k,
            density_threshold=qqb.normalize_density_token(args.density_threshold),
            ntc_labels=args.ntc_label, overwrite=args.overwrite,
            cnmf_file=args.cnmf_file, usage_file_template=args.usage_file_template,
            gex_modality=args.gex_modality, guide_modality=args.guide_modality,
            guide_metadata_guide_col=args.guide_metadata_guide_col,
            guide_metadata_target_col=args.guide_metadata_target_col,
        )
    except BaseException as exc:
        report.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        report_path.write_text(json.dumps(report, indent=2) + "\n")
        raise
    report["status"] = "completed"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    return report


def launch_crt_recomputation(argv, *, python_executable=None, preflight_only=False, stage="all"):
    """Stream a fresh Python worker's output; never consume Jupyter's sys.argv."""
    if stage not in {"all", "replicates"}:
        raise ValueError("stage must be all or replicates")
    argv = list(map(str, argv))
    args = qqb.parse_args(argv)
    requested = str(python_executable or sys.executable)
    executable = shutil.which(os.path.expanduser(requested))
    if executable is None:
        raise FileNotFoundError(f"CRT Python is missing or not executable: {requested}")
    # Do not resolve a venv's Python symlink: that can change its environment.
    executable = os.path.abspath(executable)
    source_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(source_root) + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env["PATH"] = str(Path(executable).parent) + os.pathsep + env.get("PATH", "")
    env["PYTHONUNBUFFERED"] = "1"
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env.setdefault(key, "1")
    command = [executable, "-u", "-m", "cc_perturbseq.figure_s04_runtime",
               ("check-" if preflight_only else "") + ("replicates" if stage == "replicates" else "run"), "--", *argv]
    out = Path(args.out_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    log_path = out / "s4_crt.log"
    print("CRT command:", shlex.join(command), flush=True)
    print("CRT log:", log_path, flush=True)
    tail = deque(maxlen=40)
    with log_path.open("w") as log:
        log.write(shlex.join(command) + "\n")
        with subprocess.Popen(command, env=env, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True, bufsize=1) as worker:
            try:
                for line in worker.stdout:
                    log.write(line); log.flush()
                    print(line, end="", flush=True)
                    tail.append(line)
                code = worker.wait()
            except BaseException:
                worker.terminate()
                raise
    if code != 0:
        raise RuntimeError(
            f"S4 CRT failed under {executable} (exit {code}); see {log_path}.\n"
            + "".join(tail)
        )
    return json.loads((out / "s4_crt_runtime.json").read_text())


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    modes = {"check": (True, "all"), "check-run": (True, "all"),
             "run": (False, "all"), "replicates": (False, "replicates"),
             "check-replicates": (True, "replicates")}
    if len(argv) < 2 or argv[0] not in modes or argv[1] != "--":
        raise ValueError("Usage: python -m cc_perturbseq.figure_s04_runtime {check|run|replicates} -- BACKEND_OPTIONS")
    preflight_only, stage = modes[argv[0]]
    run_worker(argv[2:], preflight_only=preflight_only, stage=stage)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
