from __future__ import annotations

from collections import OrderedDict
from collections.abc import Sequence
from pathlib import Path
import hashlib
from typing import Any

import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from .program_io import significance_stars

FIGURE2_TARGET_GROUPS: "OrderedDict[str, list[str]]" = OrderedDict(
    [
        ("BMP/TGFβ pathway", ["BMPR1A", "SMAD1", "SMAD5"]),
        ("Mesoderm regulators", ["MIXL1", "EOMES"]),
        ("LM + EC regulators", ["KDR", "KCTD10", "PRKCD", "ETV2", "PLCG1"]),
        ("Hippo pathway", ["TAOK1", "AMOTL2", "LATS2", "NF2"]),
    ]
)

FIGURE2_PROGRAMS: list[dict[str, Any]] = [
    {"day": "D1", "program": 13, "label": "P13 Posterior\nprimitive streak", "column_group": "Day 1"},
    {"day": "D1", "program": 30, "label": "P30 Primitive streak,\nBMP-repressed", "column_group": "Day 1"},
    {"day": "D1", "program": 48, "label": "P48 Primitive streak", "column_group": "Day 1"},
    {"day": "D2", "program": 15, "label": "P15 Vascular\nprogenitor", "column_group": "Day 2"},
    {"day": "D2", "program": 14, "label": "P14 Lateral\nmesoderm", "column_group": "Day 2"},
    {"day": "D2", "program": 17, "label": "P17 Lateral mesoderm,\nYAP/TAZ-induced", "column_group": "Day 2"},
    {"day": "D3", "program": 2, "label": "P2 Artery\nendothelium 1", "column_group": "Day 3 — EC"},
    {"day": "D3", "program": 3, "label": "P3 Endothelial\nfunction 1", "column_group": "Day 3 — EC"},
    {"day": "D3", "program": 8, "label": "P8 Artery\nendothelium 2", "column_group": "Day 3 — EC"},
    {"day": "D3", "program": 11, "label": "P11 Endothelial\nfunction 2", "column_group": "Day 3 — EC"},
    {"day": "D3", "program": 33, "label": "P33 Anterior-like\nmesoderm", "column_group": "Non-EC"},
    {"day": "D3", "program": 34, "label": "P34 Mesenchyme,\nYAP/TAZ-induced", "column_group": "Non-EC"},
    {"day": "D3", "program": 27, "label": "P27 Mesenchyme,\nECM", "column_group": "Non-EC"},
]

FIGURE2D_TARGETS = ["KDR", "PLCG1", "PRKCD", "ETV2"]

FIGURE2_DIVERGING = LinearSegmentedColormap.from_list(
    "figure2_diverging",
    [
        "#213460", "#2A62A5", "#71A6CA", "#CDE2ED", "#FAFBFA",
        "#F2DAD2", "#E3AA97", "#CA6E5D", "#9A2A30", "#621422",
    ],
)

def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def bh_adjust(pvalues: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1 or not np.isfinite(p).all():
        raise ValueError("BH adjustment requires one-dimensional finite P values.")
    order = np.argsort(p, kind="stable")
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.empty_like(p)
    adjusted[order] = np.clip(np.minimum.accumulate(ranked[::-1])[::-1], 0, 1)
    return adjusted
