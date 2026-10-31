#!/usr/bin/env python3
"""Export Git-LFS notebook archives to an ordinary-Git searchable text mirror.

The full notebooks remain byte-for-byte under recovery/sherlock_notebooks and
are stored with Git LFS. This script reads the materialized local notebooks and
writes small UTF-8 Markdown files containing code, Markdown, and textual output.
It also builds notebook and absolute-path indexes for repository search.

Run this only after `git lfs pull` has materialized the notebook objects.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Iterable

LFS_PREFIX = "version https://git-lfs.github.com/spec/v1"
ABSOLUTE_PATH_RE = re.compile(
    r"(?P<path>/(?:oak|home|scratch|tmp|var|data|mnt|Users)/[^\s\"'`<>|;,)]*)"
)
SECRET_PATTERNS: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(r"(?i)(github[_-]?token\s*[=:]\s*)[^\s\"']+"), r"\1<REDACTED>"),
    (re.compile(r"(?i)(api[_-]?key\s*[=:]\s*)[^\s\"']+"), r"\1<REDACTED>"),
    (re.compile(r"(?i)(password\s*[=:]\s*)[^\s\"']+"), r"\1<REDACTED>"),
    (re.compile(r"gh[pousr]_[A-Za-z0-9_]{20,}"), "<REDACTED_GITHUB_TOKEN>"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("recovery/sherlock_notebooks"),
        help="Materialized notebook archive root.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("recovery/sherlock_notebook_text"),
        help="Ordinary-Git text mirror destination.",
    )
    parser.add_argument(
        "--no-text-outputs",
        action="store_true",
        help="Omit textual cell outputs; code and Markdown are always retained.",
    )
    parser.add_argument(
        "--no-redaction",
        action="store_true",
        help="Disable conservative token/password redaction in the text mirror.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_lfs_pointer(path: Path) -> bool:
    with path.open("rb") as handle:
        prefix = handle.read(256).decode("utf-8", errors="replace")
    return prefix.startswith(LFS_PREFIX)


def redact(text: str, enabled: bool) -> str:
    if not enabled:
        return text
    for pattern, replacement in SECRET_PATTERNS:
        text = pattern.sub(replacement, text)
    return text


def source_text(cell: dict[str, Any]) -> str:
    source = cell.get("source", "")
    return "".join(source) if isinstance(source, list) else str(source)


def textual_outputs(cell: dict[str, Any]) -> Iterable[tuple[str, str]]:
    for output in cell.get("outputs", []):
        output_type = output.get("output_type", "output")
        if output_type == "stream":
            text = output.get("text", "")
            yield "stream", "".join(text) if isinstance(text, list) else str(text)
            continue
        if output_type == "error":
            traceback = output.get("traceback", [])
            yield "error", "\n".join(str(item) for item in traceback)
            continue
        data = output.get("data", {})
        for mime in ("text/plain", "text/markdown"):
            if mime in data:
                value = data[mime]
                yield mime, "".join(value) if isinstance(value, list) else str(value)


def kernel_language(payload: dict[str, Any]) -> str:
    metadata = payload.get("metadata", {})
    language = metadata.get("language_info", {}).get("name")
    if language:
        return str(language)
    return str(metadata.get("kernelspec", {}).get("language", "text"))


def extract_paths(text: str) -> list[str]:
    paths: list[str] = []
    for match in ABSOLUTE_PATH_RE.finditer(text):
        value = match.group("path").rstrip(".:")
        if value not in paths:
            paths.append(value)
    return paths


def export_notebook(
    notebook: Path,
    source_root: Path,
    output_root: Path,
    *,
    include_outputs: bool,
    redact_enabled: bool,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    payload = json.loads(notebook.read_text(encoding="utf-8"))
    relative = notebook.relative_to(source_root)
    output_path = output_root / relative.with_suffix(".ipynb.md")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    language = kernel_language(payload)
    rendered: list[str] = [
        f"# Notebook mirror: `{relative.as_posix()}`\n",
        "\n",
        f"- Source archive path: `{notebook.as_posix()}`\n",
        f"- SHA-256: `{sha256(notebook)}`\n",
        f"- Kernel language: `{language}`\n",
        "- Full-fidelity notebook storage: Git LFS\n",
        "\n",
    ]
    path_rows: list[dict[str, Any]] = []
    code_cells = 0
    markdown_cells = 0
    output_blocks = 0

    for index, cell in enumerate(payload.get("cells", [])):
        cell_type = str(cell.get("cell_type", "unknown"))
        text = redact(source_text(cell), redact_enabled)
        if cell_type == "code":
            code_cells += 1
            rendered.extend(
                [f"## Code cell {index}\n\n", f"```{language}\n", text, "\n```\n\n"]
            )
        elif cell_type == "markdown":
            markdown_cells += 1
            rendered.extend([f"## Markdown cell {index}\n\n", text, "\n\n"])
        else:
            rendered.extend([f"## {cell_type.title()} cell {index}\n\n", text, "\n\n"])

        for absolute_path in extract_paths(text):
            path_rows.append(
                {
                    "notebook": relative.as_posix(),
                    "cell_index": index,
                    "cell_type": cell_type,
                    "absolute_path": absolute_path,
                }
            )

        if include_outputs and cell_type == "code":
            for label, output_text in textual_outputs(cell):
                output_text = redact(output_text, redact_enabled)
                if not output_text.strip():
                    continue
                output_blocks += 1
                rendered.extend(
                    [
                        f"### Text output ({label})\n\n",
                        "```text\n",
                        output_text,
                        "\n```\n\n",
                    ]
                )
                for absolute_path in extract_paths(output_text):
                    path_rows.append(
                        {
                            "notebook": relative.as_posix(),
                            "cell_index": index,
                            "cell_type": f"output:{label}",
                            "absolute_path": absolute_path,
                        }
                    )

    output_path.write_text("".join(rendered), encoding="utf-8")
    summary = {
        "notebook": relative.as_posix(),
        "mirror": output_path.relative_to(output_root.parent).as_posix(),
        "size_bytes": notebook.stat().st_size,
        "sha256": sha256(notebook),
        "kernel_language": language,
        "cells_total": len(payload.get("cells", [])),
        "code_cells": code_cells,
        "markdown_cells": markdown_cells,
        "text_output_blocks": output_blocks,
        "absolute_path_occurrences": len(path_rows),
    }
    return summary, path_rows


def write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    source_root = args.source.resolve()
    output_root = args.output.resolve()
    if not source_root.is_dir():
        raise SystemExit(f"Notebook archive not found: {source_root}")

    notebooks = sorted(source_root.rglob("*.ipynb"))
    if not notebooks:
        raise SystemExit(f"No notebooks found below {source_root}")

    pointers = [path for path in notebooks if is_lfs_pointer(path)]
    if pointers:
        examples = "\n".join(f"  - {path}" for path in pointers[:10])
        raise SystemExit(
            "Git LFS objects are not materialized. Run `git lfs pull` first.\n"
            f"Pointer files detected ({len(pointers)}):\n{examples}"
        )

    output_root.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, Any]] = []
    path_rows: list[dict[str, Any]] = []
    for notebook in notebooks:
        summary, notebook_paths = export_notebook(
            notebook,
            source_root,
            output_root,
            include_outputs=not args.no_text_outputs,
            redact_enabled=not args.no_redaction,
        )
        summaries.append(summary)
        path_rows.extend(notebook_paths)

    write_tsv(
        output_root / "NOTEBOOK_INDEX.tsv",
        summaries,
        [
            "notebook",
            "mirror",
            "size_bytes",
            "sha256",
            "kernel_language",
            "cells_total",
            "code_cells",
            "markdown_cells",
            "text_output_blocks",
            "absolute_path_occurrences",
        ],
    )
    write_tsv(
        output_root / "ABSOLUTE_PATH_INDEX.tsv",
        path_rows,
        ["notebook", "cell_index", "cell_type", "absolute_path"],
    )

    unique_paths = len({row["absolute_path"] for row in path_rows})
    print(
        json.dumps(
            {
                "notebooks_exported": len(summaries),
                "text_mirror_root": str(output_root),
                "path_occurrences": len(path_rows),
                "unique_absolute_paths": unique_paths,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
