from __future__ import annotations

from pathlib import Path
import json
import math
import shutil
import zipfile

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import pearsonr, spearmanr

# Illustrator-friendly vector settings
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'DejaVu Sans'

OUTDIR = Path(__file__).resolve().parent
BASE = OUTDIR.parent
INPUT = OUTDIR / 'validation_vs_genomewide_selectedlabels_halfsize_source.tsv'
OUTDIR.mkdir(parents=True, exist_ok=True)

# Exact requested labels; SPIN1 intentionally excluded.
LABEL_GENES = [
    'BMPR1A', 'SMAD1', 'SMAD5', 'ACVR2A', 'SMAD4',
    'KDR', 'KCTD10', 'PRKCD', 'ETV2', 'PLCG1',
    'TAOK1', 'AMOTL2', 'LATS2', 'NF2',
    'MIXL1', 'EOMES', 'ACVR1B', 'SMAD3',
]

COLORS = {
    'Others': '#C9CED6',
    'Significant in GW only': '#D18F2A',
    'Significant in both': '#7751A3',
}
SIZES = {
    # Marker area is one-quarter of the original, which halves marker diameter.
    'Others': 5.75,
    'Significant in GW only': 9.5,
    'Significant in both': 10.5,
}
ORDER = ['Others', 'Significant in GW only', 'Significant in both']


def fmt_p(p: float) -> str:
    return '<1e-300' if p < 1e-300 else (f'{p:.1e}' if p < 1e-4 else f'{p:.4f}')


def pixels_to_data(ax, dx_px: float, dy_px: float) -> tuple[float, float]:
    inv = ax.transData.inverted()
    p0 = inv.transform((0, 0))
    p1 = inv.transform((dx_px, dy_px))
    return float(p1[0] - p0[0]), float(p1[1] - p0[1])


def repel_texts(ax, texts, points, n_iter: int = 90) -> None:
    """Fast adjustText-like repulsion for the small, fixed label set."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    # Initial offsets away from the centre of the plot.
    xmid = np.mean(ax.get_xlim())
    ymid = np.mean(ax.get_ylim())
    for idx, (text, (x, y)) in enumerate(zip(texts, points)):
        sx = 1 if x >= xmid else -1
        sy = 1 if y >= ymid else -1
        # Alternate vertical direction slightly to break identical starts.
        sy = sy if idx % 2 == 0 else -sy
        dx, dy = pixels_to_data(ax, 10 * sx, 8 * sy)
        text.set_position((x + dx, y + dy))
        text.set_ha('left' if sx > 0 else 'right')

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    for iteration in range(n_iter):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        boxes = [t.get_window_extent(renderer=renderer).expanded(1.08, 1.15) for t in texts]
        shifts = np.zeros((len(texts), 2), dtype=float)
        moved = False

        # Text-text repulsion in display space.
        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                bi, bj = boxes[i], boxes[j]
                if not bi.overlaps(bj):
                    continue
                moved = True
                ci = np.array([(bi.x0 + bi.x1) / 2, (bi.y0 + bi.y1) / 2])
                cj = np.array([(bj.x0 + bj.x1) / 2, (bj.y0 + bj.y1) / 2])
                delta = ci - cj
                if np.linalg.norm(delta) < 1e-6:
                    delta = np.array([1.0, 0.5])
                delta = delta / np.linalg.norm(delta)
                # Slightly stronger vertical motion, which is more useful for gene labels.
                delta = delta * np.array([1.2, 1.8])
                shifts[i] += 1.8 * delta
                shifts[j] -= 1.8 * delta

        # Repel labels from all labelled data points.
        point_disp = [ax.transData.transform(p) for p in points]
        for i, box in enumerate(boxes):
            center = np.array([(box.x0 + box.x1) / 2, (box.y0 + box.y1) / 2])
            for p in point_disp:
                if box.x0 - 2 <= p[0] <= box.x1 + 2 and box.y0 - 2 <= p[1] <= box.y1 + 2:
                    delta = center - p
                    if np.linalg.norm(delta) < 1e-6:
                        delta = np.array([1.0, 1.0])
                    shifts[i] += 1.4 * delta / np.linalg.norm(delta)
                    moved = True

        # Weak pull toward each point to avoid very long leaders.
        for i, (text, point) in enumerate(zip(texts, points)):
            tdisp = ax.transData.transform(text.get_position())
            pdisp = ax.transData.transform(point)
            delta = pdisp - tdisp
            dist = np.linalg.norm(delta)
            if dist > 70:
                shifts[i] += 0.055 * delta

        max_shift = 0.0
        for i, text in enumerate(texts):
            shift = np.clip(shifts[i], -8, 8)
            max_shift = max(max_shift, float(np.linalg.norm(shift)))
            if np.linalg.norm(shift) < 0.03:
                continue
            pos_disp = ax.transData.transform(text.get_position()) + shift
            xnew, ynew = ax.transData.inverted().transform(pos_disp)
            # Keep the label anchor inside the axes.
            xnew = min(max(xnew, xmin + 0.01 * (xmax - xmin)), xmax - 0.01 * (xmax - xmin))
            ynew = min(max(ynew, ymin + 0.01 * (ymax - ymin)), ymax - 0.01 * (ymax - ymin))
            text.set_position((xnew, ynew))

        if not moved or max_shift < 0.08:
            break

    fig.canvas.draw()


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)

    df = pd.read_csv(INPUT, sep='\t')
    required = {
        'gene', 'genome_wide_screen_log2fc', 'validation_screen_log2fc',
        'genome_wide_screen_hit', 'validation_screen_hit'
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f'Missing columns: {sorted(missing)}')

    df = df.dropna(subset=['genome_wide_screen_log2fc', 'validation_screen_log2fc']).copy()
    df['plot_category'] = np.select(
        [
            df['genome_wide_screen_hit'].astype(bool) & df['validation_screen_hit'].astype(bool),
            df['genome_wide_screen_hit'].astype(bool) & ~df['validation_screen_hit'].astype(bool),
        ],
        ['Significant in both', 'Significant in GW only'],
        default='Others',
    )

    # Statistics.
    x = df['genome_wide_screen_log2fc'].to_numpy(float)
    y = df['validation_screen_log2fc'].to_numpy(float)
    pear_r, pear_p = pearsonr(x, y)
    spear_r, spear_p = spearmanr(x, y)
    screen = df['genome_wide_screen_hit'].astype(bool)
    val = df['validation_screen_hit'].astype(bool)
    same_direction = np.sign(x) == np.sign(y)

    n_screen = int(screen.sum())
    n_both = int((screen & val).sum())
    n_same = int((screen & same_direction).sum())

    # Figure.
    fig, ax = plt.subplots(figsize=(7.15, 6.05), dpi=180)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    # Reference lines.
    max_abs = float(max(np.abs(x).max(), np.abs(y).max()))
    lim = math.ceil(max_abs * 1.08 * 10) / 10
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axline((0, 0), slope=1, color='#E3E6EA', lw=1.0, ls=':', zorder=0)
    ax.axhline(0, color='#C8CED5', lw=0.8, ls='--', zorder=0)
    ax.axvline(0, color='#C8CED5', lw=0.8, ls='--', zorder=0)

    # Points, layered from background to highlighted.
    for category in ORDER:
        sub = df[df['plot_category'] == category]
        ax.scatter(
            sub['genome_wide_screen_log2fc'],
            sub['validation_screen_log2fc'],
            s=SIZES[category],
            c=COLORS[category],
            alpha=0.72 if category == 'Others' else 0.95,
            edgecolors='none',
            zorder=1 if category == 'Others' else 2,
        )

    # Linear fit across all assayed targets.
    slope, intercept = np.polyfit(x, y, 1)
    xx = np.linspace(x.min(), x.max(), 200)
    ax.plot(xx, intercept + slope * xx, color='#272B2F', lw=1.45, zorder=3)

    # Highlight and label only the exact requested genes.
    label_df = df[df['gene'].isin(LABEL_GENES)].copy()
    found = set(label_df['gene'])
    missing_labels = [g for g in LABEL_GENES if g not in found]

    texts = []
    points = []
    for gene in LABEL_GENES:
        row = label_df[label_df['gene'] == gene]
        if row.empty:
            continue
        r = row.iloc[0]
        px = float(r['genome_wide_screen_log2fc'])
        py = float(r['validation_screen_log2fc'])
        category = str(r['plot_category'])
        ax.scatter(
            [px], [py], s=14.5, c=COLORS[category],
            edgecolors='#20252A', linewidths=0.55, zorder=5,
        )
        t = ax.text(
            px, py, gene,
            fontsize=9.1, color='#1D2024', va='center', zorder=7,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.78, pad=0.28),
        )
        texts.append(t)
        points.append((px, py))

    repel_texts(ax, texts, points, n_iter=90)

    # Leader lines after final placement.
    for text, point in zip(texts, points):
        tx, ty = text.get_position()
        ax.annotate(
            '', xy=point, xytext=(tx, ty),
            arrowprops=dict(arrowstyle='-', color='#9EA7B1', lw=0.55, alpha=0.9),
            zorder=6,
        )

    # Annotations and styling.
    ax.text(
        0.5, 1.025,
        f'{n_both}/{n_screen} screen hits significant; {n_same}/{n_screen} same direction',
        transform=ax.transAxes, ha='center', va='bottom',
        fontsize=11.0, color='#667180',
    )
    ax.text(
        0.025, 0.97,
        f'n = {len(df)}\nPearson r = {pear_r:.2f}',
        transform=ax.transAxes, ha='left', va='top',
        fontsize=10.3, color='#252A30',
    )

    ax.set_xlabel('Genome-wide U-test log2FC', fontsize=12.5, fontweight='bold')
    ax.set_ylabel('Validation U-test log2FC', fontsize=12.5, fontweight='bold')
    ax.tick_params(axis='both', labelsize=10.3, width=0.8, length=4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#30363B')
    ax.spines['bottom'].set_color('#30363B')
    ax.spines['left'].set_linewidth(0.95)
    ax.spines['bottom'].set_linewidth(0.95)

    legend_handles = [
        Line2D([0], [0], marker='o', linestyle='', markersize=3.75,
               markerfacecolor=COLORS[c], markeredgecolor='none',
               label=f'{c} (n={(df["plot_category"] == c).sum()})')
        for c in ORDER
    ]
    ax.legend(
        handles=legend_handles,
        loc='upper center', bbox_to_anchor=(0.5, -0.13),
        ncol=3, frameon=False, fontsize=9.4,
        handletextpad=0.5, columnspacing=1.5,
    )

    fig.subplots_adjust(left=0.13, right=0.98, top=0.90, bottom=0.20)

    # Exports.
    stem = 'validation_vs_genomewide_selectedlabels_halfsize'
    png = OUTDIR / f'{stem}.png'
    pdf = OUTDIR / f'{stem}.pdf'
    svg = OUTDIR / f'{stem}.svg'
    fig.savefig(png, dpi=500, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf, bbox_inches='tight', facecolor='white')
    fig.savefig(svg, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    source_cols = [
        'gene', 'genome_wide_screen_log2fc', 'validation_screen_log2fc',
        'genome_wide_screen_hit', 'validation_screen_hit', 'plot_category'
    ]
    df[source_cols].to_csv(OUTDIR / f'{stem}_source.tsv', sep='\t', index=False)
    label_df[source_cols].to_csv(OUTDIR / f'{stem}_labels.tsv', sep='\t', index=False)

    stats_df = pd.DataFrame([{
        'n_genes': len(df),
        'pearson_r': pear_r,
        'pearson_p': pear_p,
        'spearman_rho': spear_r,
        'spearman_p': spear_p,
        'n_screen_hits': n_screen,
        'n_significant_in_both': n_both,
        'n_significant_in_gw_only': int((screen & ~val).sum()),
        'n_others': int((~screen).sum()),
        'n_same_direction_screen_hits': n_same,
        'n_requested_labels': len(LABEL_GENES),
        'n_found_labels': len(label_df),
        'missing_requested_labels': ', '.join(missing_labels),
        'SPIN1_excluded': True,
        'label_method': 'iterative overlap repulsion for fixed 18-gene label set',
    }])
    stats_df.to_csv(OUTDIR / f'{stem}_stats.tsv', sep='\t', index=False)

    # Build notes and script copy.
    readme = f"""Validation vs genome-wide screen scatter - half-size points

Input:
- {INPUT}

Label set:
- {', '.join(LABEL_GENES)}
- SPIN1 intentionally excluded

Categories:
- Significant in both
- Significant in GW only
- Others

Counts:
- Total genes: {len(df)}
- Screen hits: {n_screen}
- Significant in both: {n_both}
- Significant in GW only: {int((screen & ~val).sum())}
- Others: {int((~screen).sum())}

Point-size change:
- All scatter-marker diameters are 0.5x the original.
- Matplotlib marker areas are therefore 0.25x the original.
- Outlined markers for labelled genes and legend markers were scaled consistently.

Vector exports:
- SVG text is retained as live text where supported.
- PDF uses embedded TrueType fonts.
"""
    (OUTDIR / 'README.md').write_text(readme)
    script_copy = OUTDIR / 'build_validation_vs_genomewide_selectedlabels_halfsize.py'
    if Path(__file__).resolve() != script_copy.resolve():
        shutil.copy2(Path(__file__), script_copy)

    # SVG text QC.
    svg_text = svg.read_text(errors='ignore')
    pd.DataFrame([{
        'file': svg.name,
        'contains_text_elements': '<text' in svg_text,
        'contains_path_elements': '<path' in svg_text,
    }]).to_csv(OUTDIR / 'svg_text_qc.tsv', sep='\t', index=False)

    zip_path = BASE / 'validation_vs_genomewide_selectedlabels_halfsize_20260807_package.zip'
    with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
        for fp in sorted(OUTDIR.rglob('*')):
            if fp.is_file():
                zf.write(fp, arcname=str(fp.relative_to(BASE)))

    print(json.dumps({
        'output_dir': str(OUTDIR),
        'zip': str(zip_path),
        'png': str(png),
        'pdf': str(pdf),
        'svg': str(svg),
        'n_found_labels': len(label_df),
        'missing_labels': missing_labels,
    }, indent=2))


if __name__ == '__main__':
    main()
