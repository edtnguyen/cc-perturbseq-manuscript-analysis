from __future__ import annotations

import csv
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plotnine import (
    aes,
    annotate,
    coord_fixed,
    element_blank,
    element_rect,
    element_text,
    geom_segment,
    geom_text,
    geom_tile,
    ggplot,
    scale_fill_cmap,
    theme,
    theme_void,
)

INPUT_ASSOC = Path('/mnt/data/cNMF_50_0_2_032526_6co(3).xlsx')
ASSOC_SHEET = 'Perturbation Association 1'
INPUT_EC = Path('/mnt/data/d3_EC_percent_by_guide_with_individual_NT_excl_Hippo(3).csv')
OUTDIR = Path('/mnt/data/gene_program_EC_heatmap_plotnine_revised')
OUTDIR.mkdir(parents=True, exist_ok=True)

CMAP_NAME = 'RdBu_r'
COLOR_LIMIT = 3.0
PERMUTATIONS = 100_000
PERMUTATION_SEED = 20_260_806

GENE_GROUPS: list[tuple[str, list[str]]] = [
    ('BMP/TGFb\npathway', ['BMPR1A', 'SMAD1', 'SMAD5', 'ACVR2A', 'SMAD4']),
    ('EC+ LM\nregulators', ['KDR', 'KCTD10', 'PRKCD', 'ETV2', 'PLCG1']),
    ('Hippo\npathway', ['TAOK1', 'AMOTL2', 'LATS2', 'NF2']),
    ('MIXL1/EOMES\ncluster', ['MIXL1', 'EOMES', 'ACVR1B', 'SMAD3']),
]

# Tighter within-block spacing than the previous version.
PROGRAMS: list[dict[str, object]] = [
    {'day': 'D1', 'program': 13, 'label': 'P13 Posterior\nprimitive streak', 'x': 8.34},
    {'day': 'D1', 'program': 30, 'label': 'P30 Primitive streak,\nBMP-repressed', 'x': 9.32},
    {'day': 'D1', 'program': 48, 'label': 'P48 Primitive streak', 'x': 10.30},
    {'day': 'D2', 'program': 15, 'label': 'P15 Vascular\nprogenitor', 'x': 13.35},
    {'day': 'D2', 'program': 14, 'label': 'P14 Lateral\nmesoderm', 'x': 14.33},
    {'day': 'D2', 'program': 17, 'label': 'P17 Lateral mesoderm,\nYAP/TAZ-induced', 'x': 15.31},
    {'day': 'D3', 'program': 2, 'label': 'P2 Artery\nendothelium 1', 'x': 18.55},
    {'day': 'D3', 'program': 3, 'label': 'P3 Endothelial\nfunction 1', 'x': 19.53},
    {'day': 'D3', 'program': 8, 'label': 'P8 Artery\nendothelium 2', 'x': 20.51},
    {'day': 'D3', 'program': 11, 'label': 'P11 Endothelial\nfunction 2', 'x': 21.49},
    {'day': 'D3', 'program': 33, 'label': 'P33 Anterior-like\nmesoderm', 'x': 22.47},
    {'day': 'D3', 'program': 34, 'label': 'P34 Mesenchyme,\nYAP/TAZ-induced', 'x': 23.45},
    {'day': 'D3', 'program': 27, 'label': 'P27 Mesenchyme,\nECM', 'x': 24.43},
]

EC_X = 28.55
TILE_WIDTH = 0.93
TILE_HEIGHT = 0.93

NS_MAIN = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
NS_REL = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
NS_PKG = 'http://schemas.openxmlformats.org/package/2006/relationships'


def _column_index(cell_ref: str) -> int:
    match = re.match(r'([A-Z]+)', cell_ref)
    if match is None:
        raise ValueError(f'Invalid Excel cell reference: {cell_ref}')
    index = 0
    for char in match.group(1):
        index = index * 26 + ord(char) - 64
    return index - 1


def _sheet_xml_path(zf: zipfile.ZipFile, sheet_name: str) -> str:
    workbook = ET.fromstring(zf.read('xl/workbook.xml'))
    relationships = ET.fromstring(zf.read('xl/_rels/workbook.xml.rels'))
    rel_map = {
        rel.attrib['Id']: rel.attrib['Target']
        for rel in relationships.findall(f'{{{NS_PKG}}}Relationship')
    }
    sheets = workbook.find(f'{{{NS_MAIN}}}sheets')
    if sheets is None:
        raise ValueError('Workbook has no sheets element')
    for sheet in sheets:
        if sheet.attrib.get('name') == sheet_name:
            rel_id = sheet.attrib[f'{{{NS_REL}}}id']
            target = rel_map[rel_id]
            return target if target.startswith('xl/') else f'xl/{target}'
    raise KeyError(f'Sheet not found: {sheet_name}')


def _load_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    if 'xl/sharedStrings.xml' not in zf.namelist():
        return []
    root = ET.fromstring(zf.read('xl/sharedStrings.xml'))
    return [
        ''.join(text.text or '' for text in item.iter(f'{{{NS_MAIN}}}t'))
        for item in root.findall(f'{{{NS_MAIN}}}si')
    ]


def read_association_rows(path: Path, sheet_name: str) -> list[dict[str, str]]:
    with zipfile.ZipFile(path) as zf:
        shared = _load_shared_strings(zf)
        root = ET.fromstring(zf.read(_sheet_xml_path(zf, sheet_name)))

    parsed_rows: list[list[str | None]] = []
    for row in root.iter(f'{{{NS_MAIN}}}row'):
        values: list[str | None] = [None] * 9
        for cell in row.findall(f'{{{NS_MAIN}}}c'):
            col = _column_index(cell.attrib['r'])
            if col >= len(values):
                continue
            cell_type = cell.attrib.get('t')
            raw_value = cell.find(f'{{{NS_MAIN}}}v')
            inline = cell.find(f'{{{NS_MAIN}}}is')
            if cell_type == 's' and raw_value is not None:
                value = shared[int(raw_value.text)]
            elif cell_type == 'inlineStr' and inline is not None:
                value = ''.join(text.text or '' for text in inline.iter(f'{{{NS_MAIN}}}t'))
            elif raw_value is not None:
                value = raw_value.text
            else:
                value = None
            values[col] = value
        parsed_rows.append(values)

    if not parsed_rows:
        raise ValueError(f'No rows found in {sheet_name}')
    header = [str(value or '') for value in parsed_rows[0]]
    expected = [
        '', 'target_name', 'program_name', 'specificity_scores',
        'beta_hat', 'p-value', 'adj_pval', 'approx_log2FC', 'Sample',
    ]
    if header != expected:
        raise ValueError(f'Unexpected association-sheet header: {header}')

    result: list[dict[str, str]] = []
    for values in parsed_rows[1:]:
        result.append({header[i]: '' if values[i] is None else str(values[i]) for i in range(len(header))})
    return result


def build_gene_y_positions() -> tuple[dict[str, float], dict[str, float]]:
    y_by_gene: dict[str, float] = {}
    y_by_group: dict[str, float] = {}
    y = 18.8
    group_gap = 0.80
    for group_index, (group_name, genes) in enumerate(GENE_GROUPS):
        ys = []
        for gene in genes:
            y_by_gene[gene] = y
            ys.append(y)
            y -= 1.0
        y_by_group[group_name] = float(np.mean(ys))
        if group_index < len(GENE_GROUPS) - 1:
            y -= group_gap
    return y_by_gene, y_by_group


def significance_stars(fdr: float) -> str:
    if fdr < 0.001:
        return '***'
    if fdr < 0.01:
        return '**'
    if fdr < 0.05:
        return '*'
    return ''


def benjamini_hochberg(p_values: dict[str, float]) -> dict[str, float]:
    names = list(p_values)
    p = np.asarray([p_values[name] for name in names], dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    m = len(ranked)
    adjusted_sorted = np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1]
    adjusted = np.empty_like(adjusted_sorted)
    adjusted[order] = np.minimum(adjusted_sorted, 1.0)
    return dict(zip(names, adjusted))


def build_program_tables(
    association_rows: list[dict[str, str]],
    y_by_gene: dict[str, float],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected_genes = {gene for _, genes in GENE_GROUPS for gene in genes}
    selected_keys = {(str(program['day']), int(program['program'])) for program in PROGRAMS}
    lookup: dict[tuple[str, str, int], tuple[float, float]] = {}

    for row in association_rows:
        target = row['target_name']
        if target not in selected_genes:
            continue
        sample = row['Sample']
        program = int(float(row['program_name']))
        if (sample, program) not in selected_keys:
            continue
        lookup[(target, sample, program)] = (float(row['approx_log2FC']), float(row['adj_pval']))

    tile_rows: list[dict[str, object]] = []
    star_rows: list[dict[str, object]] = []
    for _, genes in GENE_GROUPS:
        for gene in genes:
            for program in PROGRAMS:
                key = (gene, str(program['day']), int(program['program']))
                if key not in lookup:
                    raise KeyError(f'Missing association result: {key}')
                effect, adj_pval = lookup[key]
                x = float(program['x'])
                y = y_by_gene[gene]
                stars = significance_stars(adj_pval)
                tile_rows.append({
                    'gene': gene,
                    'day': str(program['day']),
                    'program': int(program['program']),
                    'program_label': str(program['label']).replace('\n', ' '),
                    'x': x,
                    'y': y,
                    'approx_log2FC': effect,
                    'display_value': float(np.clip(effect, -COLOR_LIMIT, COLOR_LIMIT)),
                    'adj_pval': adj_pval,
                    'significance': stars,
                })
                if stars:
                    star_rows.append({
                        'x': x,
                        'y': y,
                        'label': stars,
                        'text_color': 'white' if abs(effect) >= 1.45 else 'black',
                    })

    return pd.DataFrame(tile_rows), pd.DataFrame(star_rows)


def guide_level_permutation_pvalues(
    targeting_values: dict[str, np.ndarray],
    nt_values: np.ndarray,
    n_resamples: int,
    seed: int,
) -> dict[str, float]:
    """Two-sided randomization test on the unweighted mean guide EC percentage."""
    rng = np.random.default_rng(seed)
    p_values: dict[str, float] = {}
    nt_n = len(nt_values)

    guide_counts = sorted({len(values) for values in targeting_values.values()})
    for m in guide_counts:
        genes = [gene for gene, values in targeting_values.items() if len(values) == m]
        total_n = m + nt_n

        # Generate random subsets once for all genes with the same guide count.
        random_scores = rng.random((n_resamples, total_n), dtype=np.float32)
        subset_indices = np.argpartition(random_scores, m - 1, axis=1)[:, :m]
        del random_scores

        for gene in genes:
            gene_values = targeting_values[gene]
            combined = np.concatenate([gene_values, nt_values])
            total_sum = float(combined.sum())
            observed = float(gene_values.mean() - nt_values.mean())
            selected_sum = combined[subset_indices].sum(axis=1)
            permuted = selected_sum / m - (total_sum - selected_sum) / nt_n
            exceedances = int(np.count_nonzero(np.abs(permuted) >= abs(observed)))
            p_values[gene] = (exceedances + 1.0) / (n_resamples + 1.0)

        del subset_indices

    return p_values


def read_ec_summary(
    path: Path,
    selected_genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, float, int]:
    targeting: dict[str, list[dict[str, float | str]]] = {}
    individual_nt: list[float] = []

    with path.open(newline='') as handle:
        reader = csv.DictReader(handle)
        required = {
            'row_type', 'guide_id', 'target', 'n_D3_cells',
            'n_D3_EC_cells', 'n_D3_nonEC_cells', 'pct_EC_D3',
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f'EC CSV is missing required columns: {sorted(required)}')
        for row in reader:
            row_type = row['row_type']
            if row_type == 'individual_NT_guide_excl_Hippo':
                if row['pct_EC_D3'] not in ('', None):
                    individual_nt.append(float(row['pct_EC_D3']))
            elif row_type == 'targeting_guide' and row['pct_EC_D3'] not in ('', None):
                gene = row['target']
                targeting.setdefault(gene, []).append({
                    'guide_id': row['guide_id'],
                    'pct_EC_D3': float(row['pct_EC_D3']),
                    'n_D3_cells': float(row['n_D3_cells']),
                    'n_D3_EC_cells': float(row['n_D3_EC_cells']),
                    'n_D3_nonEC_cells': float(row['n_D3_nonEC_cells']),
                })

    if not individual_nt:
        raise ValueError('No usable individual NT-guide rows found in EC CSV')
    nt_values = np.asarray(individual_nt, dtype=float)
    nt_mean_pct = float(nt_values.mean())

    targeting_values = {
        gene: np.asarray([float(row['pct_EC_D3']) for row in guides], dtype=float)
        for gene, guides in targeting.items()
    }
    p_values = guide_level_permutation_pvalues(
        targeting_values,
        nt_values,
        n_resamples=PERMUTATIONS,
        seed=PERMUTATION_SEED,
    )
    q_values = benjamini_hochberg(p_values)

    all_rows: list[dict[str, object]] = []
    for gene in sorted(targeting):
        guides = targeting[gene]
        guide_pct = targeting_values[gene]
        n_cells = np.asarray([float(row['n_D3_cells']) for row in guides], dtype=float)
        n_ec = np.asarray([float(row['n_D3_EC_cells']) for row in guides], dtype=float)
        gene_mean_pct = float(guide_pct.mean())
        pooled_pct = float(100.0 * n_ec.sum() / n_cells.sum()) if n_cells.sum() > 0 else float('nan')
        log2_ratio = float(np.log2(gene_mean_pct / nt_mean_pct))
        q = float(q_values[gene])
        all_rows.append({
            'gene': gene,
            'n_guides': len(guides),
            'mean_guide_pct_EC_D3': gene_mean_pct,
            'median_guide_pct_EC_D3': float(np.median(guide_pct)),
            'sd_guide_pct_EC_D3': float(np.std(guide_pct, ddof=1)) if len(guide_pct) > 1 else float('nan'),
            'pooled_cell_pct_EC_D3': pooled_pct,
            'mean_individual_NT_guide_pct_EC_D3': nt_mean_pct,
            'mean_pct_difference_vs_NT': float(gene_mean_pct - nt_mean_pct),
            'log2_gene_mean_over_NT_mean': log2_ratio,
            'permutation_pvalue': float(p_values[gene]),
            'permutation_fdr_bh_298_genes': q,
            'significance': significance_stars(q),
        })

    all_summary = pd.DataFrame(all_rows)
    selected_summary = (
        all_summary.set_index('gene')
        .loc[selected_genes]
        .reset_index()
    )
    return selected_summary, all_summary, nt_mean_pct, len(nt_values)


def colorbar_df() -> pd.DataFrame:
    values = np.linspace(-COLOR_LIMIT, COLOR_LIMIT, 241)
    x0, x1 = 31.8, 37.1
    step = (x1 - x0) / len(values)
    return pd.DataFrame({
        'x': x0 + (np.arange(len(values)) + 0.5) * step,
        'y': np.full(len(values), 14.55),
        'display_value': values,
    })


def make_plot(
    program_tiles: pd.DataFrame,
    program_star_rows: pd.DataFrame,
    ec_summary: pd.DataFrame,
    y_by_gene: dict[str, float],
    y_by_group: dict[str, float],
):
    ec_lookup = ec_summary.set_index('gene').to_dict('index')
    ec_rows = []
    ec_star_rows = []
    for _, genes in GENE_GROUPS:
        for gene in genes:
            row = ec_lookup[gene]
            value = float(row['log2_gene_mean_over_NT_mean'])
            stars = str(row['significance'])
            ec_rows.append({
                'gene': gene,
                'x': EC_X,
                'y': y_by_gene[gene],
                'ec_log2_ratio': value,
                'display_value': float(np.clip(value, -COLOR_LIMIT, COLOR_LIMIT)),
                'permutation_pvalue': float(row['permutation_pvalue']),
                'permutation_fdr': float(row['permutation_fdr_bh_298_genes']),
                'significance': stars,
            })
            if stars:
                ec_star_rows.append({
                    'x': EC_X,
                    'y': y_by_gene[gene],
                    'label': stars,
                    'text_color': 'white' if abs(value) >= 1.45 else 'black',
                })

    ec_tiles = pd.DataFrame(ec_rows)
    ec_stars = pd.DataFrame(ec_star_rows)
    all_stars = pd.concat([program_star_rows, ec_stars], ignore_index=True)
    white_stars = all_stars.loc[all_stars['text_color'].eq('white')]
    black_stars = all_stars.loc[all_stars['text_color'].eq('black')]
    cb = colorbar_df()

    p = (
        ggplot()
        + geom_tile(
            program_tiles,
            aes('x', 'y', fill='display_value'),
            inherit_aes=False,
            width=TILE_WIDTH,
            height=TILE_HEIGHT,
            color='white',
            size=0.40,
        )
        + geom_tile(
            ec_tiles,
            aes('x', 'y', fill='display_value'),
            inherit_aes=False,
            width=TILE_WIDTH,
            height=TILE_HEIGHT,
            color='white',
            size=0.40,
        )
        + geom_tile(
            cb,
            aes('x', 'y', fill='display_value'),
            inherit_aes=False,
            width=(37.1 - 31.8) / len(cb) + 0.005,
            height=0.70,
            color=None,
        )
        + scale_fill_cmap(cmap_name=CMAP_NAME, limits=(-COLOR_LIMIT, COLOR_LIMIT), guide=None)
        + geom_text(
            white_stars,
            aes('x', 'y', label='label'),
            inherit_aes=False,
            color='white',
            family='Arimo',
            size=7.0,
            fontweight='bold',
            va='center',
        )
        + geom_text(
            black_stars,
            aes('x', 'y', label='label'),
            inherit_aes=False,
            color='black',
            family='Arimo',
            size=7.0,
            fontweight='bold',
            va='center',
        )
        + coord_fixed(ratio=0.78, xlim=(0.0, 43.0), ylim=(-7.2, 21.8), expand=False)
        + theme_void()
        + theme(
            figure_size=(13.6, 7.0),
            plot_background=element_rect(fill='white', color='white'),
            panel_background=element_rect(fill='white', color='white'),
            text=element_text(family='Arimo', color='black'),
            axis_title=element_blank(),
            axis_text=element_blank(),
            axis_ticks=element_blank(),
            axis_line=element_blank(),
            plot_margin=0,
        )
    )

    p = (
        p
        + annotate('text', x=5.9, y=20.55, label='b', family='Arimo', fontweight='bold', size=27, ha='center', va='center')
        + annotate('text', x=16.9, y=20.90, label='Effects on gene programs', family='Arimo', fontweight='bold', size=16.5, ha='center', va='center')
        + annotate('text', x=28.55, y=20.82, label='Effects on EC\ndifferentiation', family='Arimo', fontweight='bold', size=16.2, ha='center', va='center')
        + annotate('text', x=9.32, y=20.0, label='Day 1', family='Arimo', fontweight='bold', size=12.0, ha='center', va='center')
        + annotate('text', x=14.33, y=20.0, label='Day 2', family='Arimo', fontweight='bold', size=12.0, ha='center', va='center')
        + annotate('text', x=21.49, y=20.0, label='Day 3', family='Arimo', fontweight='bold', size=12.0, ha='center', va='center')
        + geom_segment(pd.DataFrame([
            {'x': 7.84, 'xend': 10.79, 'y': 19.55, 'yend': 19.55},
            {'x': 12.85, 'xend': 15.81, 'y': 19.55, 'yend': 19.55},
            {'x': 18.05, 'xend': 24.93, 'y': 19.55, 'yend': 19.55},
            {'x': 28.04, 'xend': 29.06, 'y': 19.55, 'yend': 19.55},
        ]), aes(x='x', xend='xend', y='y', yend='yend'), inherit_aes=False, color='#5A5A5A', size=0.55)
    )

    for group_name, genes in GENE_GROUPS:
        p = p + annotate(
            'text', x=0.35, y=y_by_group[group_name], label=group_name,
            family='Arimo', fontweight='bold', size=10.0, ha='left', va='center',
        )
        for gene in genes:
            p = p + annotate(
                'text', x=6.75, y=y_by_gene[gene], label=gene,
                family='Arimo', size=12.0, ha='right', va='center',
            )

    for program in PROGRAMS:
        p = p + annotate(
            'text', x=float(program['x']) + 0.20, y=-1.35, label=str(program['label']),
            family='Arimo', size=9.0, angle=55, ha='right', va='top',
        )
    p = p + annotate(
        'text', x=28.77, y=-1.35, label='D3 EC\nlog2 ratio vs NT',
        family='Arimo', size=10.0, angle=55, ha='right', va='top',
    )

    p = (
        p
        + annotate('text', x=31.8, y=16.35, label='Effect / EC ratio', family='Arimo', fontweight='bold', size=11.4, ha='left', va='center')
        + annotate('text', x=31.8, y=15.75, label='log2 scale', family='Arimo', size=10.7, ha='left', va='center')
        + annotate('text', x=31.8, y=13.85, label='-3', family='Arimo', size=8.8, ha='center', va='center')
        + annotate('text', x=34.45, y=13.85, label='0', family='Arimo', size=8.8, ha='center', va='center')
        + annotate('text', x=37.1, y=13.85, label='+3', family='Arimo', size=8.8, ha='center', va='center')
        + annotate('text', x=31.8, y=12.15, label='Significance', family='Arimo', fontweight='bold', size=11.4, ha='left', va='center')
        + annotate('text', x=31.8, y=11.35, label='* FDR < 0.05\n** FDR < 0.01\n*** FDR < 0.001', family='Arimo', size=10.3, ha='left', va='top')
        + annotate('text', x=31.8, y=7.95, label='D3 EC ratio', family='Arimo', fontweight='bold', size=11.0, ha='left', va='center')
        + annotate('text', x=31.8, y=7.20, label='log2(gene mean / NT mean)', family='Arimo', size=10.0, ha='left', va='center')
        + annotate('text', x=31.8, y=6.55, label='Mean across guides per gene', family='Arimo', size=9.2, ha='left', va='center')
        + annotate('text', x=31.8, y=5.85, label='EC stars: guide-level permutation', family='Arimo', size=9.2, ha='left', va='center')
        + annotate('text', x=31.8, y=5.20, label='BH FDR across 298 genes', family='Arimo', size=9.2, ha='left', va='center')
    )
    return p, ec_tiles


def write_csv_outputs(
    program_tiles: pd.DataFrame,
    selected_ec_summary: pd.DataFrame,
    all_ec_summary: pd.DataFrame,
    nt_mean_pct: float,
    n_nt_guides: int,
) -> None:
    program_tiles.to_csv(OUTDIR / 'selected_gene_program_effects.csv', index=False)
    selected_ec_summary.to_csv(OUTDIR / 'selected_gene_D3_EC_summary.csv', index=False)
    all_ec_summary.to_csv(OUTDIR / 'all_298_gene_D3_EC_permutation_results.csv', index=False)
    with (OUTDIR / 'METHOD.txt').open('w') as handle:
        handle.write(
            'Program-effect source: Perturbation Association 1 sheet.\n'
            'Program-effect measure: approx_log2FC.\n'
            'Program significance: adj_pval, displayed as * <0.05, ** <0.01, *** <0.001.\n'
            'Color scale: RdBu_r, clipped to [-3, +3].\n'
            'EC source: d3_EC_percent_by_guide_with_individual_NT_excl_Hippo(3).csv.\n'
            'Gene-level EC percentage: unweighted arithmetic mean of pct_EC_D3 across available targeting guides.\n'
            'For all displayed genes, six targeting guides were available.\n'
            'NT baseline: arithmetic mean of usable individual_NT_guide_excl_Hippo percentages.\n'
            f'Usable individual NT guides: {n_nt_guides}.\n'
            f'NT mean percentage used: {nt_mean_pct:.8f}.\n'
            'EC heatmap value: log2(gene mean guide percentage / NT mean guide percentage).\n'
            f'EC significance: two-sided guide-level permutation test on the unweighted mean EC percentage ({PERMUTATIONS} randomizations; seed {PERMUTATION_SEED}).\n'
            'The targeting-guide labels were permuted against the usable individual NT-guide percentages.\n'
            'EC p-values were adjusted by the Benjamini-Hochberg method across all 298 tested genes.\n'
            'Fisher exact test was not used because pooling cell counts would treat cells within the same guide as independent and would not represent guide-to-guide variability.\n'
        )


def main() -> None:
    mpl.rcParams.update({
        'font.family': 'Arimo',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'svg.fonttype': 'none',
    })

    association_rows = read_association_rows(INPUT_ASSOC, ASSOC_SHEET)
    y_by_gene, y_by_group = build_gene_y_positions()
    program_tiles, program_star_rows = build_program_tables(association_rows, y_by_gene)
    selected_genes = [gene for _, genes in GENE_GROUPS for gene in genes]
    selected_ec_summary, all_ec_summary, nt_mean_pct, n_nt_guides = read_ec_summary(INPUT_EC, selected_genes)
    plot, _ = make_plot(program_tiles, program_star_rows, selected_ec_summary, y_by_gene, y_by_group)

    figure = plot.draw(show=False)
    figure.set_size_inches(13.6, 7.0, forward=True)
    outputs = {
        'pdf': OUTDIR / 'gene_program_and_EC_effects_heatmap_plotnine.pdf',
        'svg': OUTDIR / 'gene_program_and_EC_effects_heatmap_plotnine.svg',
        'png': OUTDIR / 'gene_program_and_EC_effects_heatmap_plotnine.png',
    }
    figure.savefig(outputs['pdf'], format='pdf', dpi=300, facecolor='white', bbox_inches='tight')
    figure.savefig(outputs['svg'], format='svg', dpi=300, facecolor='white', bbox_inches='tight')
    figure.savefig(outputs['png'], format='png', dpi=300, facecolor='white', bbox_inches='tight')
    plt.close(figure)

    write_csv_outputs(program_tiles, selected_ec_summary, all_ec_summary, nt_mean_pct, n_nt_guides)
    print(f'NT individual-guide mean EC percentage: {nt_mean_pct:.6f}; n={n_nt_guides}')
    print(selected_ec_summary[[
        'gene', 'n_guides', 'mean_guide_pct_EC_D3', 'log2_gene_mean_over_NT_mean',
        'permutation_pvalue', 'permutation_fdr_bh_298_genes', 'significance',
    ]].to_string(index=False))
    for path in outputs.values():
        print(path)


if __name__ == '__main__':
    main()
