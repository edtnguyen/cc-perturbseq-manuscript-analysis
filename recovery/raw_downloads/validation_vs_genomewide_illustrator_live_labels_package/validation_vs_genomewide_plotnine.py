#!/usr/bin/env python3
"""Create an Illustrator-oriented validation-versus-genome-wide scatter plot.

The gray ``Others`` points are rasterized to reduce vector complexity. Red and
black significant points, selected-gene points, regression/reference lines,
axes, legend, title, and labels remain vector. The 14 selected-gene labels are
exported as standard live Helvetica text in PDF and direct Arial/Helvetica text
elements in SVG so their font size can be changed normally in Adobe Illustrator.
Automatic label placement continues to use Plotnine's ``adjust_text`` workflow.
"""
from __future__ import annotations

import argparse
import io
import logging
import re
import zipfile
from pathlib import Path

import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from plotnine import (
    aes,
    annotate,
    coord_cartesian,
    element_blank,
    element_line,
    element_rect,
    element_text,
    geom_abline,
    geom_hline,
    geom_point,
    geom_smooth,
    geom_text,
    geom_vline,
    ggplot,
    guide_legend,
    guides,
    labs,
    scale_alpha_manual,
    scale_color_manual,
    scale_size_manual,
    scale_x_continuous,
    scale_y_continuous,
    theme,
    theme_classic,
)

BASE_NAME = 'validation_vs_genomewide_illustrator_live_labels'
FONT_FAMILY = 'Arimo'
LABEL_FONT_FAMILY = 'Helvetica'
SELECTED_GENES = [
    'BMPR1A', 'SMAD1', 'SMAD5',
    'MIXL1', 'EOMES',
    'KDR', 'KCTD10', 'PRKCD', 'ETV2', 'PLCG1',
    'TAOK1', 'AMOTL2', 'LATS2', 'NF2',
]
CATEGORY_ORDER = [
    'Significant in GW only',
    'Significant in both',
    'Others',
]
CATEGORY_COLORS = {
    'Significant in GW only': '#B22222',
    'Significant in both': '#111111',
    'Others': '#BFC3C7',
}
CATEGORY_SIZES = {
    'Significant in GW only': 1.90,
    'Significant in both': 1.90,
    'Others': 1.20,
}
CATEGORY_ALPHA = {
    'Significant in GW only': 0.94,
    'Significant in both': 0.84,
    'Others': 0.48,
}


def configure_matplotlib() -> None:
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
    mpl.rcParams['font.family'] = FONT_FAMILY
    mpl.rcParams['font.sans-serif'] = [
        'Arimo', 'Liberation Sans', 'Nimbus Sans', 'DejaVu Sans'
    ]
    mpl.rcParams['figure.dpi'] = 300
    mpl.rcParams['savefig.dpi'] = 300
    mpl.rcParams['svg.fonttype'] = 'none'
    # Standard PDF core Helvetica is imported by Illustrator as live text,
    # unlike a subset CID font that may not expose ordinary size controls.
    mpl.rcParams['pdf.use14corefonts'] = True
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['pdf.compression'] = 9
    mpl.rcParams['text.usetex'] = False
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(
        color=['#111111', '#BFC3C7', '#B22222']
    )


def load_input(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == '.zip':
        with zipfile.ZipFile(path) as zf:
            members = [
                n for n in zf.namelist()
                if n.endswith('ec_screen_hits.csv')
                and '__MACOSX' not in n
                and not Path(n).name.startswith('._')
            ]
            if len(members) != 1:
                raise ValueError(f'Expected one ec_screen_hits.csv; found {members}')
            return pd.read_csv(io.BytesIO(zf.read(members[0])))
    return pd.read_csv(path)


def prepare_data(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float]]:
    required = {
        'gene', 'genome_wide_screen_log2fc', 'validation_screen_log2fc',
        'genome_wide_screen_hit', 'validation_screen_hit',
    }
    missing = required.difference(raw.columns)
    if missing:
        raise ValueError(f'Missing required columns: {sorted(missing)}')

    df = raw.copy()
    df['gene'] = df['gene'].astype(str).str.strip()
    df['x'] = pd.to_numeric(df['genome_wide_screen_log2fc'], errors='coerce')
    df['y'] = pd.to_numeric(df['validation_screen_log2fc'], errors='coerce')
    if df[['x', 'y']].isna().any().any():
        raise ValueError('Non-numeric screen effect values were detected.')

    df['category'] = np.select(
        [
            df['genome_wide_screen_hit'] & ~df['validation_screen_hit'],
            df['genome_wide_screen_hit'] & df['validation_screen_hit'],
        ],
        ['Significant in GW only', 'Significant in both'],
        default='Others',
    )
    df['category'] = pd.Categorical(
        df['category'], categories=CATEGORY_ORDER, ordered=True
    )
    # Draw gray points first and emphasized categories last.
    draw_order = {'Others': 0, 'Significant in GW only': 1, 'Significant in both': 2}
    df['_draw_order'] = df['category'].astype(object).map(draw_order)
    df = df.sort_values('_draw_order').reset_index(drop=True)

    selected = df.loc[df['gene'].isin(SELECTED_GENES)].copy()
    missing_genes = sorted(set(SELECTED_GENES).difference(selected['gene']))
    duplicates = selected['gene'].value_counts()
    duplicates = duplicates[duplicates != 1]
    if missing_genes or not duplicates.empty:
        raise ValueError(
            f'Selected-gene mapping failed. Missing={missing_genes}; '
            f'duplicates={duplicates.to_dict()}'
        )
    selected['requested_order'] = pd.Categorical(
        selected['gene'], categories=SELECTED_GENES, ordered=True
    )
    selected = selected.sort_values('requested_order').reset_index(drop=True)
    # Shared, non-gene-specific initialization; final locations are automatic.
    selected['label_start_x'] = selected['x'] - 0.16
    selected['label_start_y'] = selected['y'] + 0.16

    pearson_r = float(df['x'].corr(df['y']))
    slope, intercept = np.polyfit(df['x'].to_numpy(), df['y'].to_numpy(), 1)
    stats = {
        'n': int(len(df)),
        'pearson_r': pearson_r,
        'slope': float(slope),
        'intercept': float(intercept),
    }
    return df, selected, stats


def build_plot(df: pd.DataFrame, selected: pd.DataFrame, stats: dict[str, float]):
    configure_matplotlib()
    counts = df['category'].value_counts().reindex(CATEGORY_ORDER, fill_value=0)
    legend_labels = {
        category: f'{category} (n={int(counts[category])})'
        for category in CATEGORY_ORDER
    }
    annotation = f"n = {stats['n']}\nPearson r = {stats['pearson_r']:.2f}"
    others = df.loc[df['category'] == 'Others'].copy()
    significant = df.loc[df['category'] != 'Others'].copy()

    plot = (
        ggplot(df, aes('x', 'y'))
        + geom_abline(
            slope=1,
            intercept=0,
            color='#D6DBDF',
            linetype='dotted',
            size=0.55,
        )
        + geom_hline(
            yintercept=0,
            color='#B8C0C7',
            linetype='dashed',
            size=0.60,
        )
        + geom_vline(
            xintercept=0,
            color='#B8C0C7',
            linetype='dashed',
            size=0.60,
        )
        # Rasterize only the gray background points.
        + geom_point(
            data=others,
            mapping=aes('x', 'y'),
            inherit_aes=False,
            color=CATEGORY_COLORS['Others'],
            size=CATEGORY_SIZES['Others'],
            alpha=CATEGORY_ALPHA['Others'],
            stroke=0,
            raster=True,
            show_legend=False,
        )
        # Keep red and black significant points as vector objects.
        + geom_point(
            data=significant,
            mapping=aes('x', 'y', color='category', size='category', alpha='category'),
            inherit_aes=False,
            stroke=0,
            show_legend=True,
        )
        + geom_smooth(
            method='lm',
            se=False,
            color='#202020',
            size=0.85,
            fullrange=False,
        )
        + geom_point(
            data=selected,
            mapping=aes('x', 'y'),
            inherit_aes=False,
            color='#111111',
            size=2.35,
            alpha=0.95,
            stroke=0,
            show_legend=False,
        )
        + geom_text(
            data=selected,
            mapping=aes('label_start_x', 'label_start_y', label='gene'),
            inherit_aes=False,
            family=FONT_FAMILY,
            size=6.55,
            color='#111111',
            ha='center',
            va='center',
            adjust_text={
                'x': df['x'].to_numpy(),
                'y': df['y'].to_numpy(),
                'target_x': selected['x'].to_numpy(),
                'target_y': selected['y'].to_numpy(),
                'avoid_self': True,
                'prevent_crossings': True,
                'force_text': (0.75, 1.00),
                'force_static': (0.18, 0.28),
                'force_pull': (0.002, 0.004),
                'force_explode': (0.45, 0.65),
                'pull_threshold': 14,
                'expand': (1.12, 1.28),
                'max_move': (14, 14),
                'ensure_inside_axes': True,
                'only_move': {
                    'text': 'xy',
                    'static': 'xy',
                    'explode': 'xy',
                    'pull': 'xy',
                },
                'min_arrow_len': 0,
                'iter_lim': 1000,
                'arrowprops': {
                    'arrowstyle': '-',
                    'linewidth': 0.42,
                    'color': '#6F6F6F',
                    'shrinkA': 2.0,
                    'shrinkB': 3.5,
                },
            },
        )
        + annotate(
            'text',
            x=-4.30,
            y=4.23,
            label=annotation,
            ha='left',
            va='top',
            family=FONT_FAMILY,
            size=9.0,
            color='#111111',
        )
        + scale_color_manual(
            values=CATEGORY_COLORS,
            breaks=CATEGORY_ORDER,
            labels=[legend_labels[c] for c in CATEGORY_ORDER],
            name=None,
            drop=False,
        )
        + scale_size_manual(
            values=CATEGORY_SIZES,
            breaks=CATEGORY_ORDER,
            guide=None,
        )
        + scale_alpha_manual(
            values=CATEGORY_ALPHA,
            breaks=CATEGORY_ORDER,
            guide=None,
        )
        + guides(
            color=guide_legend(
                override_aes={'size': 2.7, 'alpha': 1.0, 'stroke': 0}
            )
        )
        + scale_x_continuous(
            breaks=list(range(-4, 5)),
            expand=(0, 0),
        )
        + scale_y_continuous(
            breaks=list(range(-4, 5)),
            expand=(0, 0),
        )
        + coord_cartesian(
            xlim=(-4.55, 4.55),
            ylim=(-4.55, 4.45),
        )
        + labs(
            title=' ',
            x=(
                'Genome-wide screen gRNA enrichment in CDH5$^{+}$ vs CDH5$^{-}$\n'
                '($\\log_{2}$ fold-change)'
            ),
            y=(
                'Validation screen: gRNA enrichment\n'
                '($\\log_{2}$ fold-change)'
            ),
        )
        + theme_classic(base_family=FONT_FAMILY)
        + theme(
            figure_size=(8.0, 5.2),
            plot_background=element_rect(fill='white', color=None),
            panel_background=element_rect(fill='white', color=None),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            panel_border=element_blank(),
            axis_line=element_line(color='#111111', size=0.72),
            axis_ticks=element_line(color='#111111', size=0.72),
            axis_ticks_length=4.2,
            axis_text=element_text(
                family=FONT_FAMILY,
                size=8.9,
                color='#111111',
            ),
            axis_title_x=element_text(
                family=FONT_FAMILY,
                size=11.5,
                color='#111111',
                margin={'t': 9},
            ),
            axis_title_y=element_text(
                family=FONT_FAMILY,
                size=11.5,
                color='#111111',
                margin={'r': 10},
            ),
            legend_position=(0.765, 0.135),
            legend_justification=(0.5, 0.5),
            legend_direction='vertical',
            legend_background=element_blank(),
            legend_key=element_blank(),
            legend_title=element_blank(),
            legend_text=element_text(
                family=FONT_FAMILY,
                size=8.9,
                color='#111111',
            ),
            legend_margin=0,
            legend_box_margin=0,
            plot_title_position='plot',
            plot_title=element_text(
                family=FONT_FAMILY,
                size=15.3,
                weight='bold',
                color='white',
                ha='left',
                margin={'b': 8},
            ),
            plot_margin_left=0.17,
            plot_margin_right=0.02,
            plot_margin_top=0.02,
            plot_margin_bottom=0.02,
        )
    )
    return plot



def make_svg_point_labels_illustrator_friendly(svg_path: str | Path) -> None:
    """Make selected-gene labels direct, unclipped SVG text objects."""
    svg_path = Path(svg_path)
    svg = svg_path.read_text(encoding='utf-8')

    pattern = re.compile(
        r'(?P<indent>\s*)<g id="(?P<id>gene_label_[^"]+)">\s*'
        r'(?P<text><text\b[^>]*>.*?</text>)\s*</g>',
        flags=re.S,
    )

    def replace_label(match: re.Match[str]) -> str:
        text = match.group('text')
        text = re.sub(
            r'<text\b',
            f'<text id="{match.group("id")}" data-role="point-label"',
            text,
            count=1,
        )
        text = re.sub(
            r'font-family:\s*[^;]+',
            'font-family: Arial, Helvetica, sans-serif',
            text,
            count=1,
        )
        text = re.sub(
            r'font-size:\s*[0-9.]+px',
            'font-size: 6.55px',
            text,
            count=1,
        )
        return match.group('indent') + text

    svg, count = pattern.subn(replace_label, svg)
    if count != len(SELECTED_GENES):
        raise RuntimeError(
            f'Expected to unwrap {len(SELECTED_GENES)} point labels; found {count}.'
        )
    svg_path.write_text(svg, encoding='utf-8')


def save_outputs(plot, outdir: str | Path):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {
        ext: outdir / f'{BASE_NAME}.{ext}'
        for ext in ('pdf', 'svg', 'png')
    }
    fig = plot.draw(show=False)
    fig.set_size_inches(8.0, 5.2, forward=True)
    fig.text(
        0.003, 0.965, 'd',
        family=FONT_FAMILY, fontsize=19.5, fontweight='bold',
        ha='left', va='top', color='#111111',
    )
    fig.text(
        0.042, 0.958, 'Secondary screening validates effects',
        family=FONT_FAMILY, fontsize=15.3, fontweight='bold',
        ha='left', va='top', color='#111111',
    )
    # Preserve each selected-gene label as a separate live text artist.
    label_count = 0
    for axis in fig.axes:
        for text in axis.texts:
            if text.get_text() in SELECTED_GENES:
                label_count += 1
                text.set_fontfamily(LABEL_FONT_FAMILY)
                text.set_fontsize(6.55)
                text.set_clip_on(False)
                text.set_clip_path(None)
                text.set_rasterized(False)
                text.set_gid(f'gene_label_{text.get_text()}')
    if label_count != len(SELECTED_GENES):
        raise RuntimeError(
            f'Expected {len(SELECTED_GENES)} point labels; found {label_count}.'
        )

    fig.canvas.draw()
    gene_texts = [
        text
        for axis in fig.axes
        for text in axis.texts
        if text.get_text() in SELECTED_GENES
    ]
    adjusted_positions = pd.DataFrame(
        [
            {
                'gene': text.get_text(),
                'adjusted_label_x': float(text.get_position()[0]),
                'adjusted_label_y': float(text.get_position()[1]),
            }
            for text in gene_texts
        ]
    ).sort_values('gene')
    fig.savefig(paths['pdf'], format='pdf', dpi=300, bbox_inches=None)
    fig.savefig(paths['svg'], format='svg', dpi=300, bbox_inches=None)
    make_svg_point_labels_illustrator_friendly(paths['svg'])
    fig.savefig(paths['png'], format='png', dpi=600, bbox_inches=None)
    return fig, paths, adjusted_positions


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input',
        type=Path,
        default=Path(__file__).resolve().parent / 'source' / 'ec_screen_hits.csv',
    )
    parser.add_argument(
        '--outdir',
        type=Path,
        default=Path(__file__).resolve().parent,
    )
    args = parser.parse_args()

    raw = load_input(args.input)
    df, selected, stats = prepare_data(raw)
    plot = build_plot(df, selected, stats)
    fig, paths, adjusted_positions = save_outputs(plot, args.outdir)
    plt.close(fig)

    df.to_csv(args.outdir / f'{BASE_NAME}_plot_data.csv', index=False)
    selected.to_csv(
        args.outdir / f'{BASE_NAME}_selected_gene_statistics.csv',
        index=False,
    )
    pd.DataFrame([stats]).to_csv(
        args.outdir / f'{BASE_NAME}_fit_statistics.csv',
        index=False,
    )
    adjusted_positions.to_csv(
        args.outdir / f'{BASE_NAME}_adjusted_label_positions.csv',
        index=False,
    )
    (
        df['category']
        .value_counts()
        .reindex(CATEGORY_ORDER, fill_value=0)
        .rename_axis('category')
        .reset_index(name='n')
        .to_csv(args.outdir / f'{BASE_NAME}_category_counts.csv', index=False)
    )
    print(f"n={stats['n']}; Pearson r={stats['pearson_r']:.6f}")
    for key, path in paths.items():
        print(f'{key}: {path}')


if __name__ == '__main__':
    main()
