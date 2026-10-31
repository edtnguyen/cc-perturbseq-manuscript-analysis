#!/usr/bin/env python3
"""Insert or refresh a repository-path bootstrap cell in curated figure notebooks."""

from __future__ import annotations

import json
from pathlib import Path


NOTEBOOK_DIRECTORIES = (
    Path("notebooks/main"),
    Path("notebooks/supplementary"),
)
BOOTSTRAP_TAG = "repo-path-bootstrap"
BOOTSTRAP_MARKER = "# Repository-path bootstrap for local support modules."
BOOTSTRAP_SOURCE = [
    f"{BOOTSTRAP_MARKER}\n",
    "from pathlib import Path\n",
    "import sys\n",
    "\n",
    "\n",
    "def _find_repo_root(start=None):\n",
    "    start = Path(start or Path.cwd()).resolve()\n",
    "    for candidate in (start, *start.parents):\n",
    "        if (candidate / \"src\" / \"cc_perturbseq\").is_dir():\n",
    "            return candidate\n",
    "    raise FileNotFoundError(\n",
    "        \"Could not locate the repository root containing \"\n",
    "        \"src/cc_perturbseq. Run this notebook from within the \"\n",
    "        \"cc-perturbseq-manuscript-analysis checkout.\"\n",
    "    )\n",
    "\n",
    "\n",
    "REPO_ROOT = _find_repo_root()\n",
    "for _support_path in (REPO_ROOT, REPO_ROOT / \"src\"):\n",
    "    _support_path_text = str(_support_path)\n",
    "    if _support_path_text not in sys.path:\n",
    "        sys.path.insert(0, _support_path_text)\n",
    "\n",
    "print(\n",
    "    \"Repository support paths added:\",\n",
    "    REPO_ROOT,\n",
    "    \"and\",\n",
    "    REPO_ROOT / \"src\",\n",
    ")\n",
]


def bootstrap_cell() -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "id": "repo-path-bootstrap",
        "metadata": {"tags": [BOOTSTRAP_TAG]},
        "outputs": [],
        "source": BOOTSTRAP_SOURCE,
    }


def is_bootstrap(cell: dict) -> bool:
    tags = cell.get("metadata", {}).get("tags", [])
    source = "".join(cell.get("source", []))
    return BOOTSTRAP_TAG in tags or BOOTSTRAP_MARKER in source


def insertion_index(cells: list[dict]) -> int:
    """Place the bootstrap after leading title/description cells and before code."""
    for index, cell in enumerate(cells):
        if cell.get("cell_type") == "code":
            return index
    return len(cells)


def validate_notebook(payload: dict, path: Path) -> None:
    if payload.get("nbformat") != 4 or not isinstance(payload.get("cells"), list):
        raise ValueError(f"{path}: unsupported or invalid notebook structure")
    for index, cell in enumerate(payload["cells"]):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        try:
            compile(source, f"{path}:cell-{index}", "exec")
        except SyntaxError as exc:
            raise SyntaxError(f"{path}: code cell {index} failed compilation: {exc}") from exc


def serialize(payload: dict, original_text: str) -> str:
    # Preserve the repository's existing compact-versus-pretty convention.
    if "\n" not in original_text.strip():
        text = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    else:
        text = json.dumps(payload, ensure_ascii=False, indent=1)
    return text + ("\n" if original_text.endswith("\n") else "")


def update_notebook(path: Path) -> bool:
    original_text = path.read_text(encoding="utf-8")
    payload = json.loads(original_text)
    cells = payload["cells"]
    matching = [index for index, cell in enumerate(cells) if is_bootstrap(cell)]

    if matching:
        first = matching[0]
        cells[first] = bootstrap_cell()
        for index in reversed(matching[1:]):
            del cells[index]
    else:
        cells.insert(insertion_index(cells), bootstrap_cell())

    validate_notebook(payload, path)
    updated_text = serialize(payload, original_text)
    if updated_text == original_text:
        return False
    path.write_text(updated_text, encoding="utf-8")
    return True


def main() -> None:
    paths = sorted(
        path
        for directory in NOTEBOOK_DIRECTORIES
        for path in directory.glob("*.ipynb")
    )
    if not paths:
        raise FileNotFoundError("No curated figure notebooks were found.")

    changed = []
    for path in paths:
        if update_notebook(path):
            changed.append(path)

    print(f"Validated {len(paths)} curated notebooks.")
    print(f"Updated {len(changed)} notebooks.")
    for path in changed:
        print(path)


if __name__ == "__main__":
    main()
