"""S4c post-CRT selection, ported from executed source cells 29-31.

The initial raw-download notebook ends before these cells. No old notebook is
loaded or executed at runtime. Exact original plotting/BH/summary helpers below.
Scientific settings: pooled q <= 0.10; >=20 unique target carriers whose IDs
occur in the original cNMF usage file in EACH biological replicate; prefer
skew-calibrated pooled p-values as in the source. S4b remains empirical.
"""
from pathlib import Path
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib as mpl
from plotnine import (
    ggplot, aes, geom_abline, geom_hline, geom_vline, geom_point, facet_wrap,
    scale_x_continuous, scale_y_continuous, coord_equal, labs, theme_bw, theme,
    element_text, element_line, element_rect, element_blank,
)
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
NTC_KEY_SET = {'safe', 'non-targeting', 'ntc'}
POINT_COLOR = '#17365D'
GRID_COLOR = '#E8EDF3'
ZERO_LINE_COLOR = '#D5DCE3'
DIAGONAL_COLOR = '#AEB7C2'
BORDER_COLOR = '#182635'

def normalize_program_name(x):
    s = str(x).strip()
    s = pd.Series([s]).str.replace(
        r"^(program|Program|usage|Usage|cnmf|cNMF)[_\s\-]*",
        "",
        regex=True,
    ).iloc[0]
    s = pd.Series([s]).str.replace(r"\.0$", "", regex=True).iloc[0]
    return s


def bh_qvalue(pvalues):
    p = pd.to_numeric(pd.Series(pvalues), errors="coerce").to_numpy(dtype=float)
    q = np.full(p.shape, np.nan, dtype=float)

    finite = np.isfinite(p)
    if finite.sum() == 0:
        return q

    pv = np.clip(p[finite], 0, 1)
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]

    q_finite = np.empty_like(pv)
    q_finite[order] = np.clip(ranked, 0, 1)
    q[finite] = q_finite
    return q


def read_pvalue_long_or_wide(path, day):
    path = Path(path)
    sep = "\t" if path.suffix.lower() in [".tsv", ".txt"] else ","

    header = pd.read_csv(path, sep=sep, nrows=0).columns.tolist()
    lower = {c.lower(): c for c in header}

    target_col = None
    for candidate in ["target_name", "target", "gene", "gene_name"]:
        if candidate in lower:
            target_col = lower[candidate]
            break

    program_col = None
    for candidate in ["program_name", "program", "cnmf_program", "component"]:
        if candidate in lower:
            program_col = lower[candidate]
            break

    pvalue_col = None
    for candidate in ["pvalue", "p_value", "pval", "p_val"]:
        if candidate in lower:
            pvalue_col = lower[candidate]
            break

    # Long format.
    if target_col is not None and program_col is not None and pvalue_col is not None:
        df = pd.read_csv(path, sep=sep)
        df = df[[target_col, program_col, pvalue_col]].copy()
        df = df.rename(
            columns={
                target_col: "target_name",
                program_col: "program_name",
                pvalue_col: "pooled_pvalue",
            }
        )

    # Wide target x program matrix.
    else:
        df_wide = pd.read_csv(path, sep=sep, index_col=0)
        df_wide.index = df_wide.index.astype(str)
        df_wide.columns = df_wide.columns.astype(str)
        df_wide.index.name = "target_name"

        df = (
            df_wide
            .reset_index()
            .melt(
                id_vars="target_name",
                var_name="program_name",
                value_name="pooled_pvalue",
            )
        )

    df["day"] = str(day)
    df["target_name"] = df["target_name"].astype(str).str.strip()
    df["program_name"] = df["program_name"].map(normalize_program_name)
    df["pooled_pvalue"] = pd.to_numeric(df["pooled_pvalue"], errors="coerce")
    df = df[np.isfinite(df["pooled_pvalue"])].copy()
    df["pooled_pvalue"] = df["pooled_pvalue"].clip(1e-300, 1.0)
    df["pooled_pvalue_source"] = path.name

    df = df[
        ~df["target_name"].str.strip().str.lower().isin(NTC_KEY_SET)
    ].copy()

    return df


def summarize_replicate_correlations(df):
    rows = []

    for (day, pair), g in df.groupby(["day", "replicate_pair"], sort=False):
        signable = g[
            np.isfinite(g["beta_x"])
            & np.isfinite(g["beta_y"])
            & (g["beta_x"] != 0)
            & (g["beta_y"] != 0)
        ]

        sign_agreement = np.nan
        if len(signable) > 0:
            sign_agreement = (
                np.sign(signable["beta_x"].to_numpy())
                == np.sign(signable["beta_y"].to_numpy())
            ).mean()

        rows.append({
            "day": day,
            "replicate_pair": pair,
            "n_target_program_pairs": int(len(g)),
            "n_targets": int(g["target_name"].nunique()),
            "n_programs": int(g["program_name"].nunique()),
            "pearson_r": g["beta_x"].corr(g["beta_y"], method="pearson"),
            "spearman_r": g["beta_x"].corr(g["beta_y"], method="spearman"),
            "sign_agreement": sign_agreement,
            "min_pooled_pvalue": g["pooled_pvalue"].min(),
            "median_pooled_pvalue": g["pooled_pvalue"].median(),
            "max_pooled_pvalue": g["pooled_pvalue"].max(),
            "min_pooled_qvalue": g["pooled_qvalue"].min(),
            "median_pooled_qvalue": g["pooled_qvalue"].median(),
            "max_pooled_qvalue": g["pooled_qvalue"].max(),
        })

    return pd.DataFrame(rows)


def format_r(x):
    if pd.isna(x) or not np.isfinite(x):
        return "NA"
    return f"{x:.2f}"


def save_plotnine_vector_outputs(p, output_stem, width, height):
    output_stem = Path(output_stem)
    output_stem.parent.mkdir(parents=True, exist_ok=True)

    # PDF/SVG are Illustrator-editable. PNG is just a preview.
    p.save(output_stem.with_suffix(".pdf"), width=width, height=height, units="in", verbose=False)
    p.save(output_stem.with_suffix(".svg"), width=width, height=height, units="in", verbose=False)
    p.save(output_stem.with_suffix(".png"), width=width, height=height, units="in", dpi=300, verbose=False)


def make_replicate_scatter(plot_df, title, subtitle, width, height, ncol=2):
    axis_limit = float(
        np.nanmax(np.abs(plot_df[["beta_x", "beta_y"]].to_numpy(dtype=float)))
    )
    axis_limit = axis_limit * 1.05 if np.isfinite(axis_limit) and axis_limit > 0 else 1.0

    return (
        ggplot(plot_df, aes("beta_x", "beta_y"))
        + geom_hline(yintercept=0, color=ZERO_LINE_COLOR, size=0.35)
        + geom_vline(xintercept=0, color=ZERO_LINE_COLOR, size=0.35)
        + geom_abline(
            intercept=0,
            slope=1,
            linetype="dashed",
            color=DIAGONAL_COLOR,
            size=0.75,
        )
        + geom_point(
            size=1.15,
            alpha=0.82,
            color=POINT_COLOR,
            stroke=0,
        )
        + facet_wrap("~facet_label", ncol=ncol)
        + scale_x_continuous(limits=(-axis_limit, axis_limit))
        + scale_y_continuous(limits=(-axis_limit, axis_limit))
        + coord_equal()
        + labs(
            title=title,
            subtitle=subtitle,
            x="DE program beta effect in biological replicate 1",
            y="DE program beta effect in biological replicate 2",
        )
        + theme_bw(base_size=10)
        + theme(
            figure_size=(width, height),
            plot_title=element_text(size=13, weight="bold", ha="center"),
            plot_subtitle=element_text(size=8.5, ha="center"),
            strip_text=element_text(size=8.5),
            axis_title=element_text(size=10),
            axis_text=element_text(size=8),
            panel_grid_major=element_line(color=GRID_COLOR, size=0.45),
            panel_grid_minor=element_blank(),
            panel_border=element_rect(color=BORDER_COLOR, size=0.65),
        )
    )


def _nonempty(path):
    path = Path(path).expanduser()
    return path.is_file() and path.stat().st_size > 0


def _record(path):
    path = Path(path).expanduser().resolve()
    return {'path': str(path), 'size_bytes': path.stat().st_size,
            'mtime_ns': path.stat().st_mtime_ns}


def read_pooled_results(model_root, days, preferred='skew'):
    """Original day-level BH family, before joining/filtering replicate effects.

    Original exact file priority is retained (long then wide; preferred p-type,
    then the other type). No arbitrary glob fallback or NTC pseudotarget files.
    """
    if preferred not in {'skew', 'empirical'}:
        raise ValueError('preferred must be skew or empirical')
    frames, records = [], []
    for day in days:
        order = list(dict.fromkeys([preferred, 'skew', 'empirical']))
        choices = [(Path(model_root) / day / f'{day}_target_pvalues_{kind}_{shape}.tsv', kind)
                   for kind in order for shape in ('long', 'wide')]
        chosen = next(((p, kind) for p, kind in choices if _nonempty(p)), None)
        if chosen is None:
            raise FileNotFoundError(f'No pooled target p-value table for {day}. Checked: '
                                    + ', '.join(str(p) for p, _ in choices))
        path, kind = chosen
        if kind != preferred:
            warnings.warn(f'{day}: requested {preferred} pooled p-values absent; using {kind}, '
                          'as allowed by the source. This fallback is recorded.')
        frame = read_pvalue_long_or_wide(path, day)
        if frame.duplicated(['day', 'target_name', 'program_name']).any():
            raise ValueError(f'Duplicate target/program IDs after normalization: {path}')
        # All finite REAL target x program tests for this day; not just pairs
        # with enough cells, pairs present in both replicates, or significant pairs.
        frame['pooled_qvalue'] = bh_qvalue(frame['pooled_pvalue'])
        frames.append(frame)
        records.append({'day': day, 'requested_pvalue_type': preferred,
                        'actual_pvalue_type': kind, 'n_tests_for_bh': len(frame), **_record(path)})
    return pd.concat(frames, ignore_index=True), pd.DataFrame(records)


def count_target_cells(obs, guide_long, usage_cells, guide2gene, day):
    """Cell-29 counting semantics: positive assignments, union within target,
    actual raw-usage index membership (not the subsequently zero-filled usage).
    Returns all carriers and carriers with cNMF usage; no CRT is performed.
    """
    obs = obs.copy()
    if 'biological_sample' not in obs:
        raise KeyError(f'{day}: biological_sample missing from expression observations')
    obs.index = obs.index.astype(str)
    if not obs.index.is_unique:
        raise ValueError(f'{day}: duplicated expression cell IDs')
    cells = pd.DataFrame({'cell_id': obs.index,
                          'biological_sample': obs['biological_sample'].astype(str).to_numpy()})
    cells['has_cnmf_usage'] = cells.cell_id.isin(set(map(str, usage_cells)))
    guide_col = next((k for k in ('grna_target', 'guide', 'guide_id') if k in guide_long), None)
    if 'cell_id' not in guide_long or guide_col is None:
        raise KeyError('Assignments require cell_id and grna_target/guide/guide_id')
    guides = guide_long.copy()
    guides['cell_id'] = guides.cell_id.astype(str)
    guides[guide_col] = guides[guide_col].astype(str)
    if 'value' in guides:
        guides = guides.loc[pd.to_numeric(guides.value, errors='coerce').fillna(0) > 0].copy()
    guides['target_name'] = guides[guide_col].map(guide2gene)
    missing = guides.target_name.isna()
    guides.loc[missing, 'target_name'] = guides.loc[missing, guide_col].str.split('_', n=1).str[0]
    guides = guides.loc[~guides.target_name.str.strip().str.lower().isin(NTC_KEY_SET)].copy()
    guides = guides.merge(cells, on='cell_id', how='inner', validate='many_to_one')
    keys = ['biological_sample', 'target_name']
    unique = guides.drop_duplicates([*keys, 'cell_id'])
    all_counts = unique.groupby(keys, dropna=False).size().rename('n_target_cells_all')
    usage_counts = unique.loc[unique.has_cnmf_usage].groupby(keys, dropna=False).size().rename('n_target_cells_with_cnmf_usage')
    counts = pd.concat([all_counts, usage_counts], axis=1).fillna(0).astype(int).reset_index()
    counts.insert(0, 'day', str(day))
    overlap = cells.groupby('biological_sample', dropna=False).agg(
        n_cells=('cell_id', 'size'), n_cells_with_cnmf_usage=('has_cnmf_usage', 'sum')).reset_index()
    overlap['pct_cnmf_usage'] = 100 * overlap.n_cells_with_cnmf_usage / overlap.n_cells
    overlap.insert(0, 'day', str(day))
    return counts, overlap


def _read_h5mu_metadata(path, gex_modality=None, guide_modality=None):
    """Read obs and guide names only, avoiding loading the full RNA matrix."""
    import h5py
    try:
        from anndata.io import read_elem
    except ImportError:
        from anndata._io.specs import read_elem
    with h5py.File(path, 'r') as handle:
        mods = handle['mod']
        def choose(explicit, candidates, label):
            if explicit:
                if explicit not in mods:
                    raise KeyError(f'{path}: requested {label} modality {explicit!r} absent')
                return explicit
            selected = next((key for key in candidates if key in mods), None)
            if selected is None:
                raise KeyError(f'{path}: no {label} modality; found {list(mods.keys())}')
            return selected
        rna = choose(gex_modality, ['GEX', 'gex', 'rna', 'RNA'], 'RNA')
        guides = choose(guide_modality, ['guide', 'Guide', 'gRNA', 'grna'], 'guide')
        obs = read_elem(mods[rna]['obs'])
        var = read_elem(mods[guides]['var'])
    return obs, list(var.index.astype(str))


def load_or_rebuild_target_counts(model_root, days, output_dir, counts_csv=None,
                                  mapping_policy='source_names'):
    """Use an explicit numerical count table, or rebuild source cell-29 counts.

    Without an explicit table, prefer rebuilding from per-day CRT provenance.
    An existing same-model source diagnostic is usable when raw inputs are gone.
    Count filtering is NEVER silently disabled when the counts are unavailable.
    """
    model_root, output_dir = Path(model_root), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    saved = model_root / 'biological_replicate_effects/plots/DIAGNOSTIC_target_cell_counts_by_replicate.csv'
    configs, can_rebuild = [], True
    for day in days:
        path = model_root / day / f'{day}_run_config.json'
        if not _nonempty(path):
            can_rebuild = False
            break
        conf = json.loads(path.read_text())
        if not all(k in conf and _nonempty(conf[k]) for k in ('h5mu_path', 'usage_path', 'guide_long_csv')):
            can_rebuild = False
            break
        conf['day'] = day
        configs.append(conf)
    if counts_csv is not None or (not can_rebuild and _nonempty(saved)):
        path = Path(counts_csv) if counts_csv is not None else saved
        if not _nonempty(path):
            raise FileNotFoundError(f'Target-cell count table missing or empty: {path}')
        table = pd.read_csv(path)
        receipt = {'mode': 'existing_numerical_counts', **_record(path)}
    elif can_rebuild:
        from . import crt_qq_backend as qqb
        frames, overlaps, records = [], [], []
        for conf in configs:
            policy = conf.get('guide_mapping_policy', mapping_policy)
            if policy != mapping_policy:
                raise ValueError('Count mapping policy differs from the recorded CRT. '
                                 'Use the same mapping policy, not a new gene alias mapping.')
            if policy != 'source_names':
                raise ValueError('Automatic reconstruction currently supports the source_names '
                                 'reproduction policy only. Supply a count table explicitly for a different analysis.')
            obs, guide_names = _read_h5mu_metadata(conf['h5mu_path'], conf.get('gex_modality'), conf.get('guide_modality'))
            mapping = qqb.infer_guide2gene_from_names(guide_names, conf.get('ntc_labels', ['SAFE','non-targeting','NTC']))
            usage_index = pd.read_csv(conf['usage_path'], sep='\t', usecols=[0], dtype=str).iloc[:, 0]
            assignments = pd.read_csv(conf['guide_long_csv'])
            frame, overlap = count_target_cells(obs, assignments, usage_index, mapping, conf['day'])
            frames.append(frame); overlaps.append(overlap)
            records.append({'day': conf['day'], 'mapping_policy': policy,
                            'inputs': [_record(conf[k]) for k in ('h5mu_path','usage_path','guide_long_csv')]})
        table = pd.concat(frames, ignore_index=True)
        pd.concat(overlaps, ignore_index=True).to_csv(output_dir/'DIAGNOSTIC_replicate_cnmf_usage_overlap.csv', index=False)
        receipt = {'mode': 'reconstructed_from_recorded_inputs', 'days': records}
    else:
        raise FileNotFoundError(
            'S4c requires target-carrier counts WITH ORIGINAL cNMF usage. Neither '
            'readable per-day recorded inputs nor the same-model source diagnostic '
            f'were found. Set TARGET_COUNTS_CSV to an explicit numerical count table. Expected: {saved}. '
            'The >=20 filter has NOT been silently skipped. CRT need not be rerun to generate counts.')
    required = ['day', 'biological_sample', 'target_name', 'n_target_cells_with_cnmf_usage']
    missing = set(required)-set(table)
    if missing:
        raise KeyError(f'Target count table lacks {sorted(missing)}')
    for col in required[:3]:
        table[col] = table[col].astype(str)
    if table.duplicated(required[:3]).any():
        raise ValueError('Duplicate day/replicate/target rows in target count table')
    absent = set(map(str,days))-set(table.day)
    if absent:
        raise ValueError(f'Target-count table does not cover days {sorted(absent)}')
    values = pd.to_numeric(table.n_target_cells_with_cnmf_usage, errors='coerce')
    if not (np.isfinite(values) & values.ge(0) & values.eq(np.floor(values))).all():
        raise ValueError('Target-cell counts must be finite nonnegative integers')
    table['n_target_cells_with_cnmf_usage'] = values.astype(int)
    table.to_csv(output_dir/'DIAGNOSTIC_target_cell_counts_by_replicate.csv', index=False)
    (output_dir/'target_cell_count_provenance.json').write_text(json.dumps(receipt, indent=2)+'\n')
    return table


def select_effect_pairs(effects, pooled, target_counts, *, q_cutoff=0.10, min_target_cells=20):
    """Source's active SELECTION_MODE='q'. No beta trimming or top-N fallback.

    BH q-values must already be computed over ALL real tests within each day.
    Selection changes only the compared effect pairs, never their beta values.
    """
    if not 0 <= q_cutoff <= 1 or min_target_cells < 0:
        raise ValueError('Invalid q cutoff or target-cell threshold')
    frame = effects.copy().drop(columns=['facet_label'], errors='ignore')
    for c in ['day','replicate_x','replicate_y','replicate_pair','target_name','program_name']:
        frame[c] = frame[c].astype(str).str.strip()
    frame['program_name'] = frame.program_name.map(normalize_program_name)
    frame = frame.loc[~frame.target_name.str.lower().isin(NTC_KEY_SET)
                      & np.isfinite(frame.beta_x) & np.isfinite(frame.beta_y)].copy()
    keys = ['day','replicate_pair','target_name','program_name']
    if frame.duplicated(keys).any():
        raise ValueError('Duplicate normalized effect-pair IDs')
    initial = frame.groupby(['day','replicate_pair'], sort=False).size().rename('matched_finite_real_pairs')
    if pooled.duplicated(['day','target_name','program_name']).any():
        raise ValueError('Duplicate normalized pooled p-value IDs')
    merged = frame.merge(pooled, on=['day','target_name','program_name'], how='left', validate='many_to_one')
    for side in ['x','y']:
        columns = ['day','biological_sample','target_name','n_target_cells_with_cnmf_usage']
        count = target_counts[columns].copy()
        for col in columns[:3]: count[col] = count[col].astype(str)
        if count.duplicated(columns[:3]).any(): raise ValueError('Duplicate count keys')
        count = count.rename(columns={'biological_sample': f'replicate_{side}',
                                     'n_target_cells_with_cnmf_usage': f'n_target_cells_{side}'})
        merged = merged.drop(columns=[f'n_target_cells_{side}'], errors='ignore').merge(
            count, on=['day',f'replicate_{side}','target_name'], how='left', validate='many_to_one')
        merged[f'n_target_cells_{side}'] = pd.to_numeric(merged[f'n_target_cells_{side}'], errors='coerce').fillna(0).astype(int)
    merged['has_pooled_pvalue'] = np.isfinite(merged.pooled_pvalue)
    merged['passes_target_cell_filter'] = (merged.n_target_cells_x >= min_target_cells) & (merged.n_target_cells_y >= min_target_cells)
    merged['passes_pooled_q_filter'] = merged.pooled_qvalue.le(q_cutoff)
    merged['selected_for_S4c'] = merged.has_pooled_pvalue & merged.passes_target_cell_filter & merged.passes_pooled_q_filter
    if (~merged.has_pooled_pvalue).any():
        warnings.warn(f'{(~merged.has_pooled_pvalue).sum()} pairs lack a matching pooled p-value; see filter audit.')
    selected = merged.loc[merged.selected_for_S4c].copy()
    selected['selection_rule'] = f'pooled q <= {q_cutoff}'
    report = initial.to_frame()
    for name, mask in [('with_pooled_pvalue', merged.has_pooled_pvalue),
                       ('after_target_cell_filter', merged.has_pooled_pvalue & merged.passes_target_cell_filter),
                       ('selected_pairs', merged.selected_for_S4c)]:
        report[name] = merged.loc[mask].groupby(['day','replicate_pair']).size().reindex(report.index, fill_value=0)
    # Return even when empty so callers can persist the exclusion report first.
    return selected, merged, report.reset_index()


def render_selected_effects(selected, output_dir, *, q_cutoff=0.10, min_target_cells=20):
    """Use cell-31 scatter helper, physical dimensions, titles and filenames."""
    if selected.empty:
        raise ValueError('No effects pass the source S4c filters. Inspect S4c_filter_audit.csv; '
                         'do not silently switch to an unfiltered or top-N panel.')
    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    tag = f"q_q{str(q_cutoff).replace('.', '_')}_top500_mincells{min_target_cells}"
    summary = summarize_replicate_correlations(selected)
    summary['facet_label'] = summary.apply(lambda r: (
        f"{r['day']}: {r['replicate_pair']}\nPearson r={format_r(r['pearson_r'])}; "
        f"Spearman rho={format_r(r['spearman_r'])}; n={int(r['n_target_program_pairs']):,}"), axis=1)
    plotted = selected.drop(columns=['facet_label'], errors='ignore').merge(
        summary[['day','replicate_pair','facet_label']], on=['day','replicate_pair'], how='left', validate='many_to_one')
    plotted.to_csv(out/f'pooled_selected_replicate_effect_pairs_{tag}.csv', index=False)
    summary.to_csv(out/f'pooled_selected_replicate_effect_correlation_summary_{tag}.csv', index=False)
    subtitle = f'pooled full-day q <= {q_cutoff}; target cells >= {min_target_cells} in both replicates'
    nfacets = int(plotted.facet_label.nunique()); ncol = 2 if nfacets > 1 else 1
    width = 8.5 if ncol == 2 else 5.2
    height = float(max(4.6, 3.25*np.ceil(nfacets/ncol)))
    combined = make_replicate_scatter(plotted, 'Biological replicate correlation of DE program effects', subtitle, width, height, ncol)
    save_plotnine_vector_outputs(combined, out/f'pooled_selected_biological_replicate_DE_program_effect_correlation_all_days_{tag}', width, height)
    for day, group in plotted.groupby('day', sort=False):
        width, height = 5.2, float(max(4.4, 3.6*group.facet_label.nunique()))
        plot = make_replicate_scatter(group, f'{day}: biological replicate correlation of DE program effects', subtitle, width, height, 1)
        save_plotnine_vector_outputs(plot, out/f'{day}_pooled_selected_biological_replicate_DE_program_effect_correlation_{tag}', width, height)
    return plotted, summary, combined
