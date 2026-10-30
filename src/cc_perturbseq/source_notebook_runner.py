"""Run existing source notebooks without replacing their scientific logic."""
from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import re
from typing import Mapping, Sequence


def sha256_file(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_notebook(path: str | Path):
    import nbformat

    path = Path(path).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Source notebook is missing or empty: {path}")
    return nbformat.read(path, as_version=4)


_LITERAL = r"(?:(?:[rRuUbBfF]{0,2})\"(?:\\.|[^\"\\])*\"|(?:[rRuUbBfF]{0,2})'(?:\\.|[^'\\])*')"


def replace_path_assignment(notebook, variable, value, *, expected_replacements=1):
    pattern = re.compile(
        rf"(?ms)^(?P<i>[ \t]*){re.escape(variable)}[ \t]*=[ \t]*"
        rf"Path\([ \t\r\n]*(?:{_LITERAL}[ \t\r\n]*)+\)"
    )
    replacement = json.dumps(str(Path(value).expanduser()))
    count = 0
    for cell in notebook.cells:
        if cell.cell_type != "code":
            continue

        def repl(match):
            nonlocal count
            count += 1
            return f"{match.group('i')}{variable} = Path({replacement})"

        cell.source = pattern.sub(repl, cell.source)
    if count != expected_replacements:
        raise ValueError(
            f"Expected {expected_replacements} assignment(s) for {variable}; found {count}."
        )
    return count


def patch_source_notebook(
    source_path: str | Path,
    output_path: str | Path,
    *,
    path_replacements: Mapping[str, str | Path],
    expected_sha256: str | None = None,
):
    import nbformat

    source_path = Path(source_path).expanduser().resolve()
    observed = sha256_file(source_path)
    if expected_sha256 and observed != expected_sha256:
        raise ValueError(
            f"Source SHA-256 mismatch: expected {expected_sha256}; observed {observed}."
        )
    notebook = _load_notebook(source_path)
    counts = {
        name: replace_path_assignment(notebook, name, path)
        for name, path in path_replacements.items()
    }
    notebook.metadata.setdefault("cc_perturbseq", {}).update(
        {
            "source_notebook": str(source_path),
            "source_sha256": observed,
            "path_only_replacements": {
                name: str(Path(path).expanduser())
                for name, path in path_replacements.items()
            },
        }
    )
    output_path = Path(output_path).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(notebook, output_path)
    return notebook, {
        "source_notebook": str(source_path),
        "source_sha256": observed,
        "patched_notebook": str(output_path),
        "replacement_counts": counts,
    }


def execute_notebook(notebook_or_path, output_path, *, working_directory, timeout=14_400):
    import nbformat
    from nbclient import NotebookClient

    notebook = (
        _load_notebook(notebook_or_path)
        if isinstance(notebook_or_path, (str, Path))
        else notebook_or_path
    )
    output_path = Path(output_path).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    client = NotebookClient(
        notebook,
        timeout=timeout,
        kernel_name=notebook.metadata.get("kernelspec", {}).get("name", "python3"),
        resources={"metadata": {"path": str(Path(working_directory).resolve())}},
        allow_errors=False,
        record_timing=True,
    )
    executed = client.execute()
    nbformat.write(executed, output_path)
    return executed


def _notebook_text(path: Path) -> str:
    try:
        notebook = _load_notebook(path)
    except Exception:
        return ""
    return "\n".join(
        cell.source for cell in notebook.cells
        if cell.cell_type in {"code", "markdown"}
    ).lower()


def _contains_group(text: str, alternatives: Sequence[str]) -> bool:
    return any(token.lower() in text for token in alternatives)


def discover_notebook(
    *, explicit_path, search_roots, required_token_groups,
    preferred_tokens=(), excluded_tokens=(), max_notebooks=20_000,
):
    if explicit_path:
        path = Path(explicit_path).expanduser().resolve()
        text = _notebook_text(path)
        missing = [tuple(g) for g in required_token_groups if not _contains_group(text, g)]
        if missing:
            raise ValueError(f"Explicit source notebook lacks token groups: {missing}")
        return path, [{"path": str(path), "score": sum(t.lower() in text for t in preferred_tokens), "explicit": True}]

    skip = {".git", ".ipynb_checkpoints", "__pycache__", "archive", "archives", "cache", "env", "envs", "outputs", "results", "tmp"}
    rows = []
    inspected = 0
    for raw_root in search_roots:
        root = Path(raw_root).expanduser()
        if not root.is_dir():
            continue
        for directory, dirs, files in os.walk(root):
            dirs[:] = [d for d in dirs if d.lower() not in skip and not d.startswith(".")]
            for filename in files:
                if not filename.endswith(".ipynb") or filename.startswith("."):
                    continue
                inspected += 1
                if inspected > max_notebooks:
                    raise RuntimeError("Notebook search limit reached; set the source environment variable.")
                path = Path(directory) / filename
                text = _notebook_text(path)
                if not text or any(t.lower() in text for t in excluded_tokens):
                    continue
                if not all(_contains_group(text, group) for group in required_token_groups):
                    continue
                matched = [t for t in preferred_tokens if t.lower() in text]
                rows.append({"path": str(path.resolve()), "score": len(matched), "matched_preferred_tokens": matched, "explicit": False})
    if not rows:
        raise FileNotFoundError("No existing source notebook matched; set the explicit source-notebook environment variable.")
    rows.sort(key=lambda row: (-row["score"], row["path"]))
    best = [row for row in rows if row["score"] == rows[0]["score"]]
    if len(best) != 1:
        raise RuntimeError("Source-notebook discovery is ambiguous; set the explicit source-notebook environment variable: " + json.dumps(best, indent=2))
    return Path(best[0]["path"]), rows


from .s3_source_runner import SourceCellMatch, identify_s3_bc_cells, prepare_s3_source_reproduction  # noqa: E402

__all__ = [
    "SourceCellMatch", "discover_notebook", "execute_notebook",
    "identify_s3_bc_cells", "patch_source_notebook",
    "prepare_s3_source_reproduction", "replace_path_assignment", "sha256_file",
]
