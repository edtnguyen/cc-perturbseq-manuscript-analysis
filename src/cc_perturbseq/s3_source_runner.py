"""Capture S3b/c from matching cells in an existing QC notebook."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import re

from .source_notebook_runner import _contains_group, _load_notebook, sha256_file


@dataclass(frozen=True)
class SourceCellMatch:
    panel: str
    cell_index: int
    cell_id: str | None
    score: int
    matched_tokens: tuple[str, ...]


def _window(notebook, index, radius=2):
    return "\n".join(
        notebook.cells[i].source
        for i in range(max(0, index-radius), min(len(notebook.cells), index+radius+1))
        if notebook.cells[i].cell_type in {"code", "markdown"}
    ).lower()


def identify_s3_bc_cells(notebook):
    matches = {"b": [], "c": []}
    for index, cell in enumerate(notebook.cells):
        if cell.cell_type != "code":
            continue
        source, window = cell.source.lower(), _window(notebook, index)
        b = []
        if any(t in source for t in ("scatter(", "geom_point(", "sns.scatterplot")):
            if _contains_group(window, ("total_counts", "total_umis", "ncount_rna")): b.append("total counts")
            if _contains_group(window, ("n_genes_by_counts", "detected_genes", "nfeature_rna")): b.append("detected genes")
            if _contains_group(window, ("round1", "round_1", "round 1", "bci", "barcode index", "first-round barcode")): b.append("round-1 barcode")
        if len(b) == 3:
            matches["b"].append(SourceCellMatch("b", index, cell.get("id"), 3, tuple(b)))

        c = []
        if any(t in source for t in ("violinplot", "geom_violin", "violin(")):
            if _contains_group(window, ("total_counts", "total_umis", "ncount_rna")): c.append("total counts")
            if _contains_group(window, ("n_genes_by_counts", "detected_genes", "nfeature_rna")): c.append("detected genes")
            if _contains_group(window, ("pct_counts_mt", "mitochondrial_fraction")): c.append("mitochondrial fraction")
            if _contains_group(window, ("pct_counts_ribo", "ribosomal_fraction")): c.append("ribosomal fraction")
            if _contains_group(window, ("median_abs_deviation", "median absolute deviation", "mad", "qc_pass", "before", "after")): c.append("QC threshold/stage")
        if len(c) == 5:
            matches["c"].append(SourceCellMatch("c", index, cell.get("id"), 5, tuple(c)))
    missing = [panel for panel, rows in matches.items() if not rows]
    if missing:
        raise ValueError(
            f"No identifiable original plotting cells for S3 panels {missing}. "
            "No generic replacement will be generated."
        )
    return matches


_SETUP = r'''from pathlib import Path as _CCPath
import json as _cc_json
import matplotlib.pyplot as _cc_plt
_CC_ROOT = _CCPath({output_dir!r}); _CC_ROOT.mkdir(parents=True, exist_ok=True)
_CC_COUNTS = {{"b": 0, "c": 0}}; _CC_FILES = []; _CC_SEEN = set()
_CC_SHOW = _cc_plt.show
def _cc_tracking_show(*args, **kwargs):
    result = _CC_SHOW(*args, **kwargs)
    _CC_SEEN.update(id(_cc_plt.figure(n)) for n in _cc_plt.get_fignums())
    return result
_cc_plt.show = _cc_tracking_show
def _cc_capture_s3(panel):
    figures = [_cc_plt.figure(n) for n in _cc_plt.get_fignums() if id(_cc_plt.figure(n)) not in _CC_SEEN]
    for figure in figures:
        _CC_COUNTS[panel] += 1
        stem = f"Figure_S03{{panel}}_source_{{_CC_COUNTS[panel]:02d}}"
        for suffix, kwargs in (("pdf", {{}}), ("svg", {{}}), ("png", {{"dpi": 300}})):
            path = _CC_ROOT / f"{{stem}}.{{suffix}}"
            figure.savefig(path, bbox_inches="tight", facecolor="white", transparent=False, **kwargs)
            _CC_FILES.append(str(path))
        _CC_SEEN.add(id(figure))
'''


def _inject(source, panel):
    call = f"_cc_capture_s3({panel!r})"
    pattern = re.compile(r"(?m)^(?P<i>[ \t]*)(?:plt|pyplot)\.show\(\s*\)")
    if pattern.search(source):
        return pattern.sub(lambda m: f"{m.group('i')}{call}\n{m.group(0)}", source)
    return source.rstrip() + "\n\n" + call + "\n"


def prepare_s3_source_reproduction(source_path, prepared_path, *, capture_output_dir):
    import nbformat

    source_path = Path(source_path).expanduser().resolve()
    notebook = _load_notebook(source_path)
    matches = identify_s3_bc_cells(notebook)
    setup = nbformat.v4.new_code_cell(_SETUP.format(output_dir=str(Path(capture_output_dir).resolve())))
    setup.metadata["tags"] = ["cc-perturbseq-source-capture"]
    notebook.cells.insert(0, setup)
    for panel, rows in matches.items():
        for row in rows:
            cell = notebook.cells[row.cell_index + 1]
            cell.source = _inject(cell.source, panel)
            cell.metadata.setdefault("tags", []).append(f"figure-s03{panel}-source-cell")
    notebook.cells.append(nbformat.v4.new_code_cell(
        "_cc_manifest = _CC_ROOT / 'Figure_S03_source_capture_manifest.json'\n"
        "_cc_manifest.write_text(_cc_json.dumps({'captured_files': _CC_FILES, "
        "'capture_counts': _CC_COUNTS}, indent=2) + '\\n')\nprint(_cc_manifest)"
    ))
    metadata = {
        "source_notebook": str(source_path), "source_sha256": sha256_file(source_path),
        "panels": ["b", "c"],
        "source_cell_matches": {p: [asdict(row) for row in rows] for p, rows in matches.items()},
    }
    notebook.metadata.setdefault("cc_perturbseq", {}).update(metadata)
    prepared_path = Path(prepared_path).expanduser().resolve()
    prepared_path.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(notebook, prepared_path)
    return notebook, {
        "source_notebook": str(source_path), "source_sha256": metadata["source_sha256"],
        "prepared_notebook": str(prepared_path), "matches": metadata["source_cell_matches"],
    }


__all__ = ["SourceCellMatch", "identify_s3_bc_cells", "prepare_s3_source_reproduction"]
