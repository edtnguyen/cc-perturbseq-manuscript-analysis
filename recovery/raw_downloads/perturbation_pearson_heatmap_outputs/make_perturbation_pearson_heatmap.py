from __future__ import annotations

import csv
import math
import os
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

INPUT_XLSX = Path('/mnt/data/cNMF_50_0_2_032526_6co(3).xlsx')
OUTDIR = Path('/mnt/data/perturbation_pearson_heatmap_outputs')
OUTDIR.mkdir(parents=True, exist_ok=True)
SHEET_NAME = 'Perturbation Association 1'
SELECTED = ['EOMES', 'MIXL1', 'TBXT', 'HAND1', 'TBX3', 'SOX4']
NS_MAIN = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
NS_REL = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
NS_PKG = 'http://schemas.openxmlformats.org/package/2006/relationships'


def _column_index(cell_ref: str) -> int:
    match = re.match(r'([A-Z]+)', cell_ref)
    if match is None:
        raise ValueError(f'Invalid Excel cell reference: {cell_ref}')
    idx = 0
    for char in match.group(1):
        idx = idx * 26 + (ord(char) - 64)
    return idx - 1


def _sheet_xml_path(zf: zipfile.ZipFile, sheet_name: str) -> str:
    workbook = ET.fromstring(zf.read('xl/workbook.xml'))
    relationships = ET.fromstring(zf.read('xl/_rels/workbook.xml.rels'))
    rel_map = {
        rel.attrib['Id']: rel.attrib['Target']
        for rel in relationships.findall(f'{{{NS_PKG}}}Relationship')
    }
    for sheet in workbook.find(f'{{{NS_MAIN}}}sheets'):
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


def read_association_sheet(path: Path, sheet_name: str) -> tuple[list[str], list[list[str | None]]]:
    with zipfile.ZipFile(path) as zf:
        shared = _load_shared_strings(zf)
        sheet_path = _sheet_xml_path(zf, sheet_name)
        root = ET.fromstring(zf.read(sheet_path))

    rows: list[list[str | None]] = []
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
        rows.append(values)

    if not rows:
        raise ValueError(f'No rows found in {sheet_name}')
    header = [str(value) if value is not None else '' for value in rows[0]]
    return header, rows[1:]


def build_effect_matrix(rows: list[list[str | None]]) -> tuple[list[str], list[tuple[str, str]], np.ndarray]:
    sample_order = ['D0', 'D1', 'D2', 'D3']
    program_order = [str(i) for i in range(50)]
    features = [(sample, program) for sample in sample_order for program in program_order]
    feature_index = {feature: idx for idx, feature in enumerate(features)}

    targets = sorted({str(row[1]) for row in rows if row[1] not in (None, 'non-targeting')})
    target_index = {target: idx for idx, target in enumerate(targets)}
    matrix = np.full((len(targets), len(features)), np.nan, dtype=float)

    for row in rows:
        target, program, log2fc, sample = row[1], row[2], row[7], row[8]
        if target in (None, 'non-targeting'):
            continue
        if program is None or log2fc is None or sample is None:
            raise ValueError(f'Incomplete row encountered: {row}')
        matrix[target_index[str(target)], feature_index[(str(sample), str(program))]] = float(log2fc)

    if len(targets) != 298:
        raise ValueError(f'Expected 298 perturbed genes after excluding non-targeting; found {len(targets)}')
    if matrix.shape[1] != 200:
        raise ValueError(f'Expected 200 day-program features; found {matrix.shape[1]}')
    missing = int(np.isnan(matrix).sum())
    if missing:
        raise ValueError(f'Effect matrix contains {missing} missing values')
    absent = [gene for gene in SELECTED if gene not in target_index]
    if absent:
        raise ValueError(f'Selected genes absent from workbook: {absent}')
    return targets, features, matrix


def cluster_correlations(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    correlations = np.corrcoef(matrix)
    correlations = np.clip(correlations, -1.0, 1.0)
    np.fill_diagonal(correlations, 1.0)
    distances = np.clip(1.0 - correlations, 0.0, 2.0)
    np.fill_diagonal(distances, 0.0)
    condensed = squareform(distances, checks=False)
    tree = linkage(condensed, method='average', optimal_ordering=True)
    order = leaves_list(tree)
    return correlations, tree, order


def adjusted_label_positions(natural: dict[str, float], lower: float = 0.155, upper: float = 0.895,
                             min_gap: float = 0.057) -> dict[str, float]:
    items = sorted(natural.items(), key=lambda item: item[1], reverse=True)
    placed: dict[str, float] = {}
    previous = upper + min_gap
    for gene, y in items:
        placed[gene] = min(y, previous - min_gap, upper)
        previous = placed[gene]
    bottom = min(placed.values())
    if bottom < lower:
        shift = lower - bottom
        placed = {gene: y + shift for gene, y in placed.items()}
    top = max(placed.values())
    if top > upper:
        shift = top - upper
        placed = {gene: y - shift for gene, y in placed.items()}
    return placed


def save_full_heatmap(targets: list[str], correlations: np.ndarray, order: np.ndarray) -> dict[str, int]:
    ordered_targets = [targets[idx] for idx in order]
    rank = {gene: idx for idx, gene in enumerate(ordered_targets)}
    ordered = correlations[np.ix_(order, order)]
    n = len(ordered_targets)

    mpl.rcParams.update({
        'font.family': 'Arimo',
        'font.size': 12,
        'axes.titlesize': 24,
        'axes.labelsize': 17,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'svg.fonttype': 'none',
    })

    fig = plt.figure(figsize=(12.0, 8.6), facecolor='white')
    ax = fig.add_axes([0.07, 0.17, 0.64, 0.745])
    image = ax.imshow(
        ordered,
        cmap='RdBu_r',
        vmin=-1,
        vmax=1,
        origin='upper',
        interpolation='nearest',
        aspect='equal',
        rasterized=True,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('298 perturbed genes', pad=15)
    ax.set_xlabel('Pearson r - Effect on program expression', labelpad=14)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_alpha(0.45)

    colors = mpl.rcParams['axes.prop_cycle'].by_key()['color']
    natural_y: dict[str, float] = {}
    x_edge_data = n - 0.45
    for idx, gene in enumerate(SELECTED):
        pos = rank[gene]
        color = colors[idx % len(colors)]
        ax.scatter(
            [pos], [pos], marker='s', s=42,
            facecolors='none', edgecolors=color, linewidths=1.8,
            zorder=5, clip_on=True,
        )
        ax.plot(
            [n - 0.5, n + 2.0], [pos, pos],
            linewidth=1.5, color=color, clip_on=False, solid_capstyle='butt',
        )
        _, y_fig = fig.transFigure.inverted().transform(ax.transData.transform((x_edge_data, pos)))
        natural_y[gene] = float(y_fig)

    label_y = adjusted_label_positions(natural_y)
    x0 = 0.710
    x_elbow = 0.745
    x_text = 0.790
    fig.text(0.805, 0.982, 'Selected genes (cluster order)', ha='center', va='top', fontsize=16)

    for idx, gene in enumerate(SELECTED):
        color = colors[idx % len(colors)]
        y0 = natural_y[gene]
        y1 = label_y[gene]
        line = Line2D(
            [x0, x_elbow, x_text - 0.012],
            [y0, y0, y1],
            transform=fig.transFigure,
            color=color,
            linewidth=1.35,
            solid_capstyle='round',
            zorder=6,
        )
        fig.add_artist(line)
        fig.text(
            x_text, y1, f'{gene}  (rank {rank[gene] + 1})',
            ha='left', va='center', fontsize=14, color=color,
        )

    fig.text(
        0.77, 0.145,
        'Average-linkage clustering\non correlation distance (1 - r)',
        ha='center', va='top', fontsize=11,
    )

    cax = fig.add_axes([0.245, 0.060, 0.29, 0.032])
    cbar = fig.colorbar(image, cax=cax, orientation='horizontal', ticks=[-1, 0, 1])
    cbar.set_label('Pearson r', fontsize=15, labelpad=4)
    cbar.ax.tick_params(labelsize=11, length=3)
    cbar.outline.set_visible(False)

    for extension in ('pdf', 'svg', 'png'):
        output = OUTDIR / f'perturbation_gene_pearson_heatmap_annotated.{extension}'
        kwargs = {'facecolor': 'white'}
        if extension == 'png':
            kwargs['dpi'] = 600
        fig.savefig(output, **kwargs)
    plt.close(fig)
    return {gene: rank[gene] + 1 for gene in SELECTED}


def save_selected_pairwise_heatmap(targets: list[str], correlations: np.ndarray) -> np.ndarray:
    target_index = {gene: idx for idx, gene in enumerate(targets)}
    selected_matrix = correlations[np.ix_([target_index[g] for g in SELECTED], [target_index[g] for g in SELECTED])]

    fig = plt.figure(figsize=(7.8, 7.8), facecolor='white')
    ax = fig.add_axes([0.19, 0.245, 0.67, 0.63])
    image = ax.imshow(
        selected_matrix,
        cmap='RdBu_r',
        vmin=-1,
        vmax=1,
        origin='upper',
        interpolation='nearest',
        aspect='equal',
    )
    ax.set_xticks(np.arange(len(SELECTED)))
    ax.set_yticks(np.arange(len(SELECTED)))
    ax.set_xticklabels(SELECTED, rotation=45, ha='right', rotation_mode='anchor', fontsize=12)
    ax.set_yticklabels(SELECTED, fontsize=12)
    ax.tick_params(length=0)
    ax.set_title('Pairwise correlation among selected perturbations', pad=16, fontsize=20)

    for row in range(len(SELECTED)):
        for col in range(len(SELECTED)):
            value = selected_matrix[row, col]
            text_color = 'white' if abs(value) >= 0.55 else 'black'
            ax.text(col, row, f'{value:.2f}', ha='center', va='center', fontsize=11, color=text_color)

    ax.set_xticks(np.arange(-0.5, len(SELECTED), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(SELECTED), 1), minor=True)
    ax.grid(which='minor', linewidth=1.2, color='white')
    ax.tick_params(which='minor', bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    cax = fig.add_axes([0.24, 0.075, 0.57, 0.026])
    cbar = fig.colorbar(image, cax=cax, orientation='horizontal', ticks=[-1, 0, 1])
    cbar.set_label('Pearson r', fontsize=13, labelpad=2)
    cbar.outline.set_visible(False)

    for extension in ('pdf', 'svg', 'png'):
        output = OUTDIR / f'selected_six_gene_pairwise_correlations.{extension}'
        kwargs = {'facecolor': 'white'}
        if extension == 'png':
            kwargs['dpi'] = 600
        fig.savefig(output, **kwargs)
    plt.close(fig)
    return selected_matrix


def write_tables(targets: list[str], correlations: np.ndarray, order: np.ndarray,
                 ranks: dict[str, int], selected_matrix: np.ndarray) -> None:
    ordered_targets = [targets[idx] for idx in order]
    target_index = {gene: idx for idx, gene in enumerate(targets)}
    rank0 = {gene: idx for idx, gene in enumerate(ordered_targets)}

    with (OUTDIR / 'clustered_gene_order.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['cluster_rank', 'target_name', 'selected_gene'])
        for idx, gene in enumerate(ordered_targets, start=1):
            writer.writerow([idx, gene, gene in SELECTED])

    with (OUTDIR / 'selected_six_gene_pairwise_correlations.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['target_name', *SELECTED])
        for gene, row in zip(SELECTED, selected_matrix):
            writer.writerow([gene, *[f'{value:.8f}' for value in row]])

    with (OUTDIR / 'selected_gene_cluster_summary.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow([
            'target_name', 'cluster_rank', 'rank_percentile',
            'immediate_left_gene', 'immediate_right_gene',
            'nearest_selected_gene', 'nearest_selected_r',
            'top_overall_correlated_gene', 'top_overall_r',
            'five_genes_before', 'five_genes_after',
        ])
        for gene in SELECTED:
            pos = rank0[gene]
            index = target_index[gene]
            left = ordered_targets[pos - 1] if pos > 0 else ''
            right = ordered_targets[pos + 1] if pos + 1 < len(ordered_targets) else ''

            selected_others = [other for other in SELECTED if other != gene]
            nearest_selected = max(selected_others, key=lambda other: correlations[index, target_index[other]])
            nearest_selected_r = correlations[index, target_index[nearest_selected]]

            sorted_indices = np.argsort(correlations[index])[::-1]
            top_idx = next(candidate for candidate in sorted_indices if candidate != index)
            top_gene = targets[top_idx]
            top_r = correlations[index, top_idx]

            before = ordered_targets[max(0, pos - 5):pos]
            after = ordered_targets[pos + 1:min(len(ordered_targets), pos + 6)]
            writer.writerow([
                gene,
                ranks[gene],
                f'{100.0 * (ranks[gene] - 1) / (len(ordered_targets) - 1):.2f}',
                left,
                right,
                nearest_selected,
                f'{nearest_selected_r:.8f}',
                top_gene,
                f'{top_r:.8f}',
                ';'.join(before),
                ';'.join(after),
            ])

    with (OUTDIR / 'METHOD.txt').open('w') as handle:
        handle.write(
            'Input sheet: Perturbation Association 1\n'
            'Effect measure: approx_log2FC\n'
            'Rows included: 298 perturbed genes; non-targeting excluded\n'
            'Features per gene: all 50 cNMF programs at D0, D1, D2, and D3 (200 values)\n'
            'Similarity: pairwise Pearson correlation across the 200-value effect vectors\n'
            'Ordering: average-linkage hierarchical clustering using distance = 1 - Pearson r; optimal leaf ordering enabled\n'
            'No p-value filtering, significance thresholding, clipping, or imputation was applied.\n'
        )


def main() -> None:
    header, rows = read_association_sheet(INPUT_XLSX, SHEET_NAME)
    expected = ['target_name', 'program_name', 'specificity_scores', 'beta_hat', 'p-value', 'adj_pval', 'approx_log2FC', 'Sample']
    if header[1:9] != expected:
        raise ValueError(f'Unexpected header: {header}')

    targets, features, matrix = build_effect_matrix(rows)
    correlations, tree, order = cluster_correlations(matrix)
    ranks = save_full_heatmap(targets, correlations, order)
    selected_matrix = save_selected_pairwise_heatmap(targets, correlations)
    write_tables(targets, correlations, order, ranks, selected_matrix)

    print(f'Created outputs in {OUTDIR}')
    print('Selected ranks:', ranks)
    print('Selected pairwise correlations:')
    for gene, row in zip(SELECTED, selected_matrix):
        print(gene, ' '.join(f'{value:+.3f}' for value in row))


if __name__ == '__main__':
    main()
