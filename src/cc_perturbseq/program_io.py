from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
import re
from typing import Any

import numpy as np
import pandas as pd


def config_value(config: Mapping[str, Any], dotted_key: str, *, required: bool = True) -> Any:
    node: Any = config
    for part in dotted_key.split("."):
        if not isinstance(node, Mapping) or part not in node:
            if required:
                raise KeyError(f"Missing configuration key: {dotted_key}")
            return None
        node = node[part]
    if node in (None, "") and required:
        raise ValueError(f"Empty configuration value: {dotted_key}")
    return node


def parse_program(value: Any) -> int | None:
    text = str(value)
    for pattern in (
        r"\bP\s*[_-]?\s*(\d+)\b",
        r"(?:cnmf|program)\s*[_-]?\s*(\d+)\b",
        r"^\s*(\d+)\s*$",
    ):
        match = re.search(pattern, text, flags=re.I)
        if match:
            return int(match.group(1))
    return None


def clean_day(value: Any) -> str | None:
    text = str(value).strip().lower().replace("day", "d").replace(" ", "")
    match = re.search(r"(?:^|[^0-9])d?([0-3])(?:$|[^0-9])", text)
    return f"D{match.group(1)}" if match else None


def significance_stars(q: float | int | None) -> str:
    if q is None:
        return ""
    try:
        finite = np.isfinite(float(q))
    except (TypeError, ValueError):
        return ""
    if not finite:
        return ""
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


def infer_column(
    df: pd.DataFrame,
    candidates: Sequence[str],
    *,
    required: bool = True,
) -> str | None:
    exact = {str(c).strip().lower(): str(c) for c in df.columns}
    for candidate in candidates:
        if candidate.lower() in exact:
            return exact[candidate.lower()]
    for column in df.columns:
        if any(candidate.lower() in str(column).lower() for candidate in candidates):
            return str(column)
    if required:
        raise ValueError(
            f"Could not infer a column from {list(candidates)}. "
            f"Available: {list(df.columns)}"
        )
    return None


def read_table(path: str | Path, **kwargs: Any) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, **kwargs)
    if suffix in {".parquet", ".pq"}:
        return pd.read_parquet(path, **kwargs)
    if suffix in {".pkl", ".pickle"}:
        return pd.read_pickle(path, **kwargs)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", **kwargs)
    return pd.read_csv(path, **kwargs)


def infer_day_from_path(path: Path) -> str | None:
    for part in reversed(path.parts):
        day = clean_day(part)
        if day:
            return day
    return None


def standardize_program_effects(
    effect_table: pd.DataFrame,
    qvalue_table: pd.DataFrame | None = None,
    *,
    inferred_day: str | None = None,
) -> pd.DataFrame:
    table = effect_table.copy().dropna(axis=0, how="all").dropna(axis=1, how="all")
    target_col = infer_column(
        table,
        ["target", "target_name", "gene", "gene_symbol", "perturbation"],
        required=False,
    )
    program_col = infer_column(
        table,
        ["program", "program_id", "program_name", "cnmf_program", "response", "component"],
        required=False,
    )
    effect_col = infer_column(
        table,
        ["effect", "obs_coef", "approx_log2fc", "log2fc", "coefficient", "coef", "beta", "estimate"],
        required=False,
    )
    if target_col and program_col and effect_col:
        day_col = infer_column(
            table, ["day", "timepoint", "sample", "diff_day"], required=False
        )
        q_col = infer_column(
            table,
            ["q_value", "qvalue", "qval", "fdr", "adj_pval", "padj", "adjusted_p_value"],
            required=False,
        )
        p_col = infer_column(
            table, ["p_value", "pvalue", "pval"], required=False
        )
        out = pd.DataFrame(
            {
                "target": table[target_col].astype(str).str.strip(),
                "program": table[program_col].map(parse_program),
                "effect": pd.to_numeric(table[effect_col], errors="coerce"),
                "q_value": (
                    pd.to_numeric(table[q_col], errors="coerce") if q_col else np.nan
                ),
                "p_value": (
                    pd.to_numeric(table[p_col], errors="coerce") if p_col else np.nan
                ),
                "day": table[day_col].map(clean_day) if day_col else inferred_day,
            }
        )
    else:
        target_col = target_col or str(table.columns[0])
        records: list[dict[str, Any]] = []
        indexed = table.set_index(target_col)
        for column in indexed.columns:
            program = parse_program(column)
            if program is None:
                continue
            day = clean_day(column) or inferred_day
            for target, value in indexed[column].items():
                records.append(
                    {
                        "target": str(target),
                        "program": program,
                        "day": day,
                        "effect": value,
                    }
                )
        out = pd.DataFrame(records)
        if out.empty:
            raise ValueError(
                "Unrecognized program-effect table; expected long columns or "
                "wide program-labelled columns."
            )
        out["effect"] = pd.to_numeric(out["effect"], errors="coerce")
        out["q_value"] = np.nan
        out["p_value"] = np.nan

    out = out.dropna(subset=["target", "program", "effect"]).copy()
    out["program"] = out["program"].astype(int)
    out["day"] = out["day"].map(clean_day)
    if qvalue_table is not None:
        q = standardize_program_effects(
            qvalue_table, inferred_day=inferred_day
        ).rename(columns={"effect": "q_from_table"})
        out = out.merge(
            q[["target", "day", "program", "q_from_table"]],
            on=["target", "day", "program"],
            how="left",
        )
        out["q_value"] = out["q_value"].where(
            out["q_value"].notna(), out.pop("q_from_table")
        )
    return out.sort_values(
        ["day", "program", "target"], na_position="last"
    ).reset_index(drop=True)


def load_program_effects(
    effect_path: str | Path,
    qvalue_path: str | Path | None = None,
    *,
    filename: str | None = None,
) -> pd.DataFrame:
    effect_path = Path(effect_path)
    qvalue_path = Path(qvalue_path) if qvalue_path else None
    if effect_path.is_dir():
        files = sorted(effect_path.rglob(filename or "*.csv"))
        if not files:
            raise FileNotFoundError(f"No matching effect files under {effect_path}")
        pieces = [
            standardize_program_effects(
                read_table(path), inferred_day=infer_day_from_path(path)
            )
            for path in files
        ]
        return pd.concat(pieces, ignore_index=True).drop_duplicates(
            ["target", "day", "program"], keep="last"
        )
    q_table = (
        read_table(qvalue_path)
        if qvalue_path is not None and qvalue_path.exists()
        else None
    )
    return standardize_program_effects(
        read_table(effect_path),
        q_table,
        inferred_day=infer_day_from_path(effect_path),
    )


def effect_matrix(effects: pd.DataFrame) -> pd.DataFrame:
    data = effects.copy()
    data["feature"] = data.apply(
        lambda row: (
            f"{row.day}_P{int(row.program)}"
            if pd.notna(row.day)
            else f"P{int(row.program)}"
        ),
        axis=1,
    )
    return data.pivot_table(
        index="target", columns="feature", values="effect", aggfunc="mean"
    )


def load_spectra(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep=r"\s+|\t|,", engine="python", index_col=0)
    row_ids = [parse_program(value) for value in table.index]
    col_ids = [parse_program(value) for value in table.columns]
    if sum(value is not None for value in col_ids) > sum(
        value is not None for value in row_ids
    ):
        table = table.T
        row_ids = [parse_program(value) for value in table.index]
    table.index = [
        program if program is not None else value
        for program, value in zip(row_ids, table.index)
    ]
    table = table.loc[
        [value for value in table.index if isinstance(value, (int, np.integer))]
    ]
    table.index = table.index.astype(int)
    return table.apply(pd.to_numeric, errors="coerce")


def standardize_motif_enrichment(
    table: pd.DataFrame,
    *,
    inferred_day: str | None = None,
) -> pd.DataFrame:
    day_col = infer_column(table, ["day", "timepoint"], required=False)
    program_col = infer_column(
        table, ["program", "program_id", "cnmf"], required=False
    )
    motif_col = infer_column(
        table, ["motif", "tf_motif_id", "tf", "cluster"], required=False
    )
    q_col = infer_column(
        table, ["q_value", "qvalue", "fdr", "padj", "adjusted"], required=False
    )
    score_col = infer_column(
        table, ["score", "effect", "statistic", "u_stat", "signed"], required=False
    )
    if not (program_col and motif_col and q_col):
        raise ValueError(
            "Motif enrichment must be long with program, motif and q/FDR columns."
        )
    return pd.DataFrame(
        {
            "day": (
                table[day_col].map(clean_day) if day_col else inferred_day
            ),
            "program": table[program_col].map(parse_program),
            "motif": table[motif_col].astype(str),
            "q_value": pd.to_numeric(table[q_col], errors="coerce"),
            "score": (
                pd.to_numeric(table[score_col], errors="coerce")
                if score_col
                else 1.0
            ),
        }
    ).dropna(subset=["program", "motif", "q_value"])


def load_motif_enrichment(
    path: str | Path | None = None,
    *,
    fdr_template: str | None = None,
    score_template: str | None = None,
    days: Sequence[str] = ("D0", "D1", "D2", "D3"),
) -> pd.DataFrame:
    if path is not None and Path(path).exists():
        return standardize_motif_enrichment(
            read_table(path), inferred_day=infer_day_from_path(Path(path))
        )
    if not fdr_template:
        raise FileNotFoundError(
            "No readable motif table or FDR template configured."
        )
    records: list[dict[str, Any]] = []
    for day in days:
        day_number = int(day[1:])
        fdr_path = Path(fdr_template.format(day=day_number, DAY=day))
        if not fdr_path.exists():
            continue
        fdr = read_table(fdr_path, index_col=0)
        scores = None
        if score_template:
            score_path = Path(score_template.format(day=day_number, DAY=day))
            scores = read_table(score_path, index_col=0) if score_path.exists() else None
        for row in fdr.index:
            program = parse_program(row)
            if program is None:
                continue
            for motif in fdr.columns:
                q_value = pd.to_numeric(
                    pd.Series([fdr.loc[row, motif]]), errors="coerce"
                ).iloc[0]
                score = (
                    1.0
                    if scores is None
                    or row not in scores.index
                    or motif not in scores.columns
                    else pd.to_numeric(
                        pd.Series([scores.loc[row, motif]]), errors="coerce"
                    ).iloc[0]
                )
                records.append(
                    {
                        "day": day,
                        "program": program,
                        "motif": str(motif),
                        "q_value": q_value,
                        "score": score,
                    }
                )
    if not records:
        raise FileNotFoundError("No motif matrices resolved from template.")
    return pd.DataFrame(records)


def load_program_annotations(
    path: str | Path | None,
) -> tuple[dict[int, str], dict[str, list[int]]]:
    if path is None or not Path(path).exists():
        return {}, {}
    table = read_table(path)
    program_col = infer_column(table, ["program_id", "program", "cnmf"])
    label_col = infer_column(
        table, ["program_label", "annotation", "name", "label"], required=False
    )
    day_col = infer_column(
        table, ["day", "stage", "timepoint"], required=False
    )
    labels: dict[int, str] = {}
    by_day: dict[str, list[int]] = {}
    for _, row in table.iterrows():
        program = parse_program(row[program_col])
        day = clean_day(row[day_col]) if day_col else None
        if program is None:
            continue
        labels[program] = str(row[label_col]) if label_col else f"P{program}"
        if day:
            by_day.setdefault(day, []).append(program)
    return labels, by_day


def load_usage(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".h5ad":
        import anndata as ad

        adata = ad.read_h5ad(path)
        if "cnmf_usage" not in adata.obsm:
            raise KeyError("AnnData lacks obsm['cnmf_usage']")
        matrix = np.asarray(adata.obsm["cnmf_usage"])
        names = adata.uns.get("cnmf_program_names")
        if names is None or len(names) != matrix.shape[1]:
            names = [f"P{i + 1}" for i in range(matrix.shape[1])]
        return pd.DataFrame(
            matrix,
            index=adata.obs_names,
            columns=[str(value) for value in names],
        )

    table = read_table(path)
    program_columns = [
        column for column in table.columns if parse_program(column) is not None
    ]
    if program_columns:
        usage = table[program_columns].apply(pd.to_numeric, errors="coerce")
        usage.columns = [f"P{parse_program(column)}" for column in program_columns]
    else:
        usage = table.select_dtypes(include=[np.number])
    usage = usage.dropna(axis=1, how="all")
    if usage.empty:
        raise ValueError(f"No numeric program-usage columns in {path}")
    return usage


def save_source_tables(
    tables: Mapping[str, pd.DataFrame],
    output_dir: str | Path,
    prefix: str,
) -> dict[str, str]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, str] = {}
    for name, table in tables.items():
        path = output_dir / f"{prefix}_{name}.tsv"
        table.to_csv(path, sep="\t")
        outputs[name] = str(path)
    return outputs


def select_top_targets(
    effects: pd.DataFrame,
    day: str,
    programs: Sequence[int],
    *,
    n_per_program: int = 4,
    q_threshold: float = 0.05,
) -> list[str]:
    selected: list[str] = []
    day = clean_day(day) or day
    data = effects[
        (effects.day == day) & effects.program.isin(programs)
    ].copy()
    for program in programs:
        part = data[data.program == int(program)].copy()
        significant = part[part.q_value < q_threshold]
        if not significant.empty:
            part = significant
        part["magnitude"] = part.effect.abs()
        for target in part.nlargest(n_per_program, "magnitude").target:
            if target not in selected:
                selected.append(str(target))
    return selected
