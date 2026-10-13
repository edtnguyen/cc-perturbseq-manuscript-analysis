#!/usr/bin/env python3
from __future__ import annotations

import io
import math
import shutil
import zipfile
from pathlib import Path

import matplotlib as mpl
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
    geom_hline,
    geom_point,
    geom_segment,
    geom_text,
    geom_vline,
    ggplot,
    guides,
    guide_legend,
    labs,
    scale_color_manual,
    scale_fill_manual,
    scale_shape_manual,
    scale_size_continuous,
    scale_x_continuous,
    scale_y_continuous,
    theme,
    theme_classic,
)

INPUT_ZIP = Path('/mnt/data/IP_GV_selected_program_plots.zip')
OUT = Path('/mnt/data/selected_program_vector_plots')
OUT.mkdir(parents=True, exist_ok=True)

# Illustrator-friendly editable text and vector geometry.
mpl.rcParams.update({
    'font.family': 'Arimo',
    'font.sans-serif': ['Arimo', 'Arial', 'Liberation Sans', 'DejaVu Sans'],
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'pdf.use14corefonts': False,
    'axes.unicode_minus': False,
    'savefig.transparent': False,
})

PROGRAM_ORDER = ['P15', 'P9', 'P6', 'P12']
PROGRAM_COLORS = {
    'P15': '#0072B2',  # blue
    'P9': '#E69F00',   # orange
    'P6': '#009E73',   # green
    'P12': '#CC79A7',  # magenta
}
SAMPLE_COLORS = {'IP_GFP': '#0072B2', 'GV_GFP': '#D55E00'}
METHOD_COLORS = {
    'Matched weighted permutation': '#0072B2',
    'MAGMA expression-conditioned': '#E69F00',
    'MAGMA top-50 competitive': '#009E73',
}
METHOD_SHAPES = {
    'Matched weighted permutation': 'o',
    'MAGMA expression-conditioned': 's',
    'MAGMA top-50 competitive': '^',
}


def read_tsv(name: str) -> pd.DataFrame:
    with zipfile.ZipFile(INPUT_ZIP) as zf:
        with zf.open(name) as handle:
            return pd.read_csv(handle, sep='\t')


integrated = read_tsv('selected_program_integrated_evidence.tsv')
pairs = read_tsv('selected_program_nominal_cross_mouse_pairs.tsv')
drivers = read_tsv('selected_program_top_GWAS_driver_genes.tsv')
annotations = read_tsv('selected_program_annotations.tsv')

# Preserve a clean copy of the exact plotting data in the output bundle.
for filename, frame in {
    'selected_program_integrated_evidence.tsv': integrated,
    'selected_program_nominal_cross_mouse_pairs.tsv': pairs,
    'selected_program_top_GWAS_driver_genes.tsv': drivers,
    'selected_program_annotations.tsv': annotations,
}.items():
    frame.to_csv(OUT / filename, sep='\t', index=False)


def base_theme(width: float, height: float, legend_position='right'):
    return (
        theme_classic(base_family='Arimo', base_size=9)
        + theme(
            figure_size=(width, height),
            plot_title=element_text(size=12, weight='bold', ha='left', margin={'b': 2}),
            plot_subtitle=element_text(size=9.5, ha='left', color='#333333', margin={'b': 8}),
            plot_caption=element_text(size=7.5, ha='left', color='#4D4D4D', margin={'t': 8}),
            axis_title=element_text(size=9.5),
            axis_text=element_text(size=8.5, color='#222222'),
            axis_ticks=element_line(color='#333333', size=0.45),
            axis_line=element_line(color='#333333', size=0.55),
            panel_border=element_rect(color='#333333', fill=None, size=0.55),
            legend_position=legend_position,
            legend_title=element_text(size=8.5, weight='bold'),
            legend_text=element_text(size=8),
            legend_key=element_rect(fill='white', color='white'),
            legend_background=element_rect(fill='white', color='white'),
            plot_background=element_rect(fill='white', color='white'),
            panel_background=element_rect(fill='white', color='white'),
            plot_margin=0.08,
        )
    )


def save_plot(plot, stem: str, width: float, height: float):
    pdf = OUT / f'{stem}_editable.pdf'
    png = OUT / f'{stem}_preview.png'
    plot.save(pdf, width=width, height=height, units='in', dpi=300, verbose=False)
    plot.save(png, width=width, height=height, units='in', dpi=300, verbose=False)
    return pdf, png


def fdr_label(value: float) -> str:
    if value < 0.001:
        return f'{value:.1g}'
    if value < 0.01:
        return f'{value:.3f}'.rstrip('0').rstrip('.')
    return f'{value:.3f}'.rstrip('0').rstrip('.')


def p_label(value: float) -> str:
    if value < 0.01:
        return f'{value:.3g}'
    return f'{value:.3f}'.rstrip('0').rstrip('.')


# -----------------------------------------------------------------------------
# 1. Integrated evidence
# -----------------------------------------------------------------------------
p1_df = integrated.copy()
p1_df['matched_neglog10p'] = -np.log10(p1_df['matched_permutation_p'])
p1_df['program'] = pd.Categorical(p1_df['program'], categories=PROGRAM_ORDER, ordered=True)
p1_df['label'] = p1_df.apply(
    lambda row: (
        f"{row['program']}\n"
        f"{int(row['disease_aligned_nominal_pairs'])} aligned target"
        f"{'s' if int(row['disease_aligned_nominal_pairs']) != 1 else ''}\n"
        f"human GSEA FDR={fdr_label(float(row['human_PAH_GSEA_FDR']))}"
    ),
    axis=1,
)
label_positions_1 = {
    'P15': (-1.46, 1.86),
    'P9': (1.34, 0.50),
    'P6': (1.37, 1.69),
    'P12': (1.62, 0.25),
}
p1_df['label_x'] = p1_df['program'].astype(str).map(lambda x: label_positions_1[x][0])
p1_df['label_y'] = p1_df['program'].astype(str).map(lambda x: label_positions_1[x][1])

p1 = (
    ggplot(p1_df, aes('human_PAH_NES', 'matched_neglog10p'))
    + geom_hline(yintercept=-math.log10(0.05), linetype='dashed', color='#666666', size=0.55)
    + geom_vline(xintercept=0, linetype='dashed', color='#666666', size=0.55)
    + geom_point(
        aes(color='program', size='disease_aligned_nominal_pairs'),
        alpha=0.9,
        stroke=0.7,
    )
    + geom_text(
        aes(x='label_x', y='label_y', label='label'),
        ha='left', va='top', size=8.2, lineheight=1.08, color='#202020',
        inherit_aes=False,
    )
    + scale_color_manual(values=PROGRAM_COLORS, name='Program')
    + scale_size_continuous(
        range=(4.2, 8.6),
        breaks=[1, 2, 3],
        limits=(1, 3),
        name='Disease-aligned\nnominal targets',
    )
    + scale_x_continuous(breaks=np.arange(-1.5, 2.6, 0.5))
    + scale_y_continuous(breaks=np.arange(0, 2.01, 0.25))
    + coord_cartesian(xlim=(-1.70, 2.72), ylim=(0.05, 1.95))
    + labs(
        title='Integrated evidence for selected PAH programs',
        subtitle='Human PAH transcriptional direction, common-variant enrichment, and replicated perturbation support',
        x='Human PAH program enrichment (GSE185479 preranked GSEA NES)',
        y='Loading-weighted PAH GWAS enrichment, -log10(p)',
        caption='Dashed horizontal line: nominal GWAS p=0.05. Negative NES indicates a program reduced in human PAH.',
    )
    + guides(
        color=guide_legend(order=1, override_aes={'size': 4.2}),
        size=guide_legend(order=2),
    )
    + base_theme(8.6, 5.6, legend_position='right')
)
save_plot(p1, '01_integrated_program_evidence', 8.6, 5.6)


# -----------------------------------------------------------------------------
# 2. Nominal in-vivo perturbation effects
# -----------------------------------------------------------------------------
pair_order = pairs['pair_label'].tolist()
pairs2 = pairs.copy()
pairs2['pair_label'] = pd.Categorical(
    pairs2['pair_label'], categories=list(reversed(pair_order)), ordered=True
)
pairs2['direction'] = np.where(pairs2['disease_aligned'], 'toward PAH', 'opposite PAH')
pairs2['p_text'] = pairs2.apply(
    lambda row: f"p={p_label(float(row['signed_stouffer_p']))}; {row['direction']}",
    axis=1,
)
long_rows = []
for _, row in pairs2.iterrows():
    for sample, beta_col in [('IP_GFP', 'disease_aligned_IP_beta'), ('GV_GFP', 'disease_aligned_GV_beta')]:
        long_rows.append({
            'pair_label': row['pair_label'],
            'sample': sample,
            'beta': row[beta_col],
        })
pairs_long = pd.DataFrame(long_rows)
pairs_long['sample'] = pd.Categorical(
    pairs_long['sample'], categories=['IP_GFP', 'GV_GFP'], ordered=True
)
pairs_long['pair_label'] = pd.Categorical(
    pairs_long['pair_label'], categories=list(reversed(pair_order)), ordered=True
)

p2 = (
    ggplot(pairs2, aes(y='pair_label'))
    + geom_vline(xintercept=0, linetype='dashed', color='#666666', size=0.55)
    + geom_segment(
        aes(x='disease_aligned_IP_beta', xend='disease_aligned_GV_beta', yend='pair_label'),
        color='#B8C7D9', size=0.8,
    )
    + geom_point(
        data=pairs_long,
        mapping=aes(x='beta', y='pair_label', color='sample', shape='sample'),
        size=3.5, stroke=0.7,
        inherit_aes=False,
    )
    + geom_text(
        aes(x=2.56, label='p_text'),
        ha='left', va='center', size=7.5, color='#222222',
    )
    + scale_color_manual(values=SAMPLE_COLORS, name='Mouse')
    + scale_shape_manual(values={'IP_GFP': 'o', 'GV_GFP': '^'}, name='Mouse')
    + scale_x_continuous(breaks=np.arange(-1.5, 3.1, 0.5))
    + coord_cartesian(xlim=(-1.50, 3.15))
    + labs(
        title='Nominal replicated in-vivo perturbation effects',
        subtitle='Pairs have the same beta sign in both mice and signed-Stouffer p<0.05',
        x='Disease-aligned perturbation effect\n(positive = moves the program toward the human PAH direction)',
        y=None,
        caption='Exploratory nominal tests; no pair is FDR significant. Dagger denotes a disputed PAH gene-disease relationship.',
    )
    + guides(
        color=guide_legend(order=1, nrow=1),
        shape=guide_legend(order=1, nrow=1),
    )
    + base_theme(9.4, 6.1, legend_position='top')
    + theme(
        axis_text_y=element_text(size=8.5),
        legend_direction='horizontal',
        legend_box='horizontal',
        panel_grid_major_y=element_line(color='#ECECEC', size=0.4),
        axis_title_y=element_blank(),
    )
)
save_plot(p2, '02_nominal_in_vivo_perturbation_effects', 9.4, 6.1)


# -----------------------------------------------------------------------------
# 3. Common-variant genetic evidence
# -----------------------------------------------------------------------------
program_label_map = dict(zip(annotations['program'], annotations['program_label']))
base_y = {'P15': 3.0, 'P9': 2.0, 'P6': 1.0, 'P12': 0.0}
methods = [
    ('Matched weighted permutation', 'matched_permutation_p', 'matched_permutation_BH_FDR', 0.18),
    ('MAGMA expression-conditioned', 'MAGMA_expression_conditioned_p', 'MAGMA_expression_conditioned_BH_FDR', 0.0),
    ('MAGMA top-50 competitive', 'MAGMA_top50_p', 'MAGMA_top50_BH_FDR', -0.18),
]
long3 = []
for _, row in integrated.iterrows():
    for method, p_col, q_col, offset in methods:
        long3.append({
            'program': row['program'],
            'program_label': program_label_map[row['program']],
            'method': method,
            'neglog10p': -math.log10(float(row[p_col])),
            'BH_FDR': float(row[q_col]),
            'y': base_y[row['program']] + offset,
            'significant': float(row[q_col]) < 0.10,
        })
p3_df = pd.DataFrame(long3)
p3_df['method'] = pd.Categorical(
    p3_df['method'], categories=list(METHOD_COLORS), ordered=True
)
p3_sig = p3_df[p3_df['significant']].copy()
p3_sig['star_x'] = p3_sig['neglog10p'] + 0.08

p3 = (
    ggplot(p3_df, aes('neglog10p', 'y'))
    + geom_vline(xintercept=-math.log10(0.05), linetype='dashed', color='#666666', size=0.55)
    + geom_point(aes(color='method', shape='method'), size=3.6, stroke=0.7)
    + geom_text(
        data=p3_sig,
        mapping=aes(x='star_x', y='y', label='"*"'),
        inherit_aes=False,
        size=11, ha='left', va='center', color='#111111',
    )
    + scale_color_manual(values=METHOD_COLORS, name='Method')
    + scale_shape_manual(values=METHOD_SHAPES, name='Method')
    + scale_y_continuous(
        breaks=[3, 2, 1, 0],
        labels=[program_label_map[p] for p in PROGRAM_ORDER],
        limits=(-0.48, 3.48),
    )
    + scale_x_continuous(breaks=np.arange(0, 2.81, 0.5), limits=(-0.03, 2.78))
    + labs(
        title='Human common-variant genetic evidence for selected programs',
        subtitle='Three complementary program-level enrichment tests; * BH FDR < 0.10 across all 30 programs',
        x='Common-variant PAH enrichment, -log10(p)',
        y=None,
        caption='Dashed line: nominal p=0.05. The three methods share input genetics and are not independent replications.',
    )
    + guides(
        color=guide_legend(order=1, nrow=1),
        shape=guide_legend(order=1, nrow=1),
    )
    + base_theme(9.5, 5.2, legend_position='top')
    + theme(
        axis_text_y=element_text(size=8.5),
        axis_title_y=element_blank(),
        legend_direction='horizontal',
        legend_box='horizontal',
    )
)
save_plot(p3, '03_common_variant_genetic_evidence', 9.5, 5.2)


# -----------------------------------------------------------------------------
# 4. Human PAH program validation
# -----------------------------------------------------------------------------
p4_df = integrated.copy()
p4_df['program'] = pd.Categorical(p4_df['program'], categories=PROGRAM_ORDER, ordered=True)
p4_df['gsea_strength'] = -np.log10(p4_df['human_PAH_GSEA_FDR'])
p4_df['label'] = p4_df.apply(
    lambda row: (
        f"{row['program']}\nFDR={fdr_label(float(row['human_PAH_GSEA_FDR']))}"
        + ("\ncomplete donor separation" if bool(row['human_complete_separation']) else "")
    ),
    axis=1,
)
label_positions_4 = {
    'P15': (-1.46, -0.31),
    'P9': (1.57, 0.72),
    'P6': (1.33, 0.43),
    'P12': (1.26, 3.16),
}
p4_df['label_x'] = p4_df['program'].astype(str).map(lambda x: label_positions_4[x][0])
p4_df['label_y'] = p4_df['program'].astype(str).map(lambda x: label_positions_4[x][1])

p4 = (
    ggplot(p4_df, aes('human_PAH_NES', 'human_donor_Hedges_g'))
    + geom_hline(yintercept=0, linetype='dashed', color='#666666', size=0.55)
    + geom_vline(xintercept=0, linetype='dashed', color='#666666', size=0.55)
    + geom_point(aes(color='program', size='gsea_strength'), alpha=0.9, stroke=0.7)
    + geom_text(
        aes(x='label_x', y='label_y', label='label'),
        inherit_aes=False,
        ha='left', va='top', size=8.2, lineheight=1.05, color='#202020',
    )
    + scale_color_manual(values=PROGRAM_COLORS, name='Program')
    + scale_size_continuous(
        range=(4.3, 9.0),
        breaks=[1.3, 2.0, 3.0],
        labels=['FDR 0.05', 'FDR 0.01', 'FDR 0.001'],
        name='GSEA significance',
    )
    + scale_x_continuous(breaks=np.arange(-1.5, 2.1, 0.5))
    + scale_y_continuous(breaks=np.arange(-0.5, 3.1, 0.5))
    + coord_cartesian(xlim=(-1.70, 2.45), ylim=(-0.65, 3.22))
    + labs(
        title='Independent human PAEC disease validation',
        subtitle='Concordant donor-score and GSEA directions indicate external human PAH support',
        x='Human PAH preranked GSEA NES',
        y='Human donor-score Hedges g (PAH - healthy)',
        caption='Point area reflects -log10(GSEA FDR). GSE185479 contains three PAH and three healthy PAEC donors.',
    )
    + guides(
        color=guide_legend(order=1, override_aes={'size': 4.2}),
        size=guide_legend(order=2),
    )
    + base_theme(8.3, 5.8, legend_position='right')
)
save_plot(p4, '04_human_PAH_program_validation', 8.3, 5.8)


# -----------------------------------------------------------------------------
# 5. Top PAH GWAS driver genes
# -----------------------------------------------------------------------------
p5_df = drivers.copy()
p5_df['program'] = 'P' + p5_df['program'].astype(str)
p5_df['human_symbol'] = p5_df['human_symbol'].astype(str).str.upper()
p5_df['program_gene'] = p5_df['program'] + '  ' + p5_df['human_symbol']
p5_df['gene_p_label'] = p5_df['gene_P'].apply(lambda x: f"p={p_label(float(x))}")

ordered_rows = []
positions = []
separators = []
current_y = 24.0
for program in PROGRAM_ORDER:
    subset = p5_df[p5_df['program'] == program].sort_values(
        'weighted_gene_Z_contribution', ascending=False
    )
    for _, row in subset.iterrows():
        out = row.to_dict()
        out['y'] = current_y
        ordered_rows.append(out)
        positions.append((current_y, out['program_gene']))
        current_y -= 1.0
    separators.append(current_y + 0.45)
    current_y -= 0.75
p5_plot = pd.DataFrame(ordered_rows)
separators = separators[:-1]
max_x = float(p5_plot['weighted_gene_Z_contribution'].max())
p5_plot['label_x'] = p5_plot['weighted_gene_Z_contribution'] + max_x * 0.018

p5 = (
    ggplot(p5_plot)
    + geom_segment(
        aes(
            x=0,
            xend='weighted_gene_Z_contribution',
            y='y',
            yend='y',
            color='program',
        ),
        size=6.2,
        lineend='butt',
    )
    + geom_text(
        aes(x='label_x', y='y', label='gene_p_label'),
        ha='left', va='center', size=7.2, color='#222222',
    )
    + geom_hline(yintercept=separators, color='#D9D9D9', size=0.5)
    + scale_color_manual(values=PROGRAM_COLORS, guide=None)
    + scale_y_continuous(
        breaks=[x[0] for x in positions],
        labels=[x[1] for x in positions],
        limits=(current_y + 0.35, 24.7),
    )
    + scale_x_continuous(
        breaks=np.arange(0, 0.051, 0.01),
        limits=(0, max_x * 1.28),
        expand=(0, 0),
    )
    + labs(
        title='Top PAH GWAS driver genes within selected programs',
        subtitle='Rank-weighted MAGMA gene-Z contributions after excluding experimentally targeted genes and extended-MHC genes',
        x='Rank-weighted MAGMA gene-Z contribution',
        y=None,
        caption='Bars show the five largest positive gene-level contributions within each selected program; labels report nominal MAGMA gene p-values.',
    )
    + base_theme(9.0, 7.6, legend_position='none')
    + theme(
        axis_text_y=element_text(size=8.4),
        panel_grid_major_x=element_line(color='#ECECEC', size=0.4),
        axis_ticks_major_y=element_blank(),
        axis_title_y=element_blank(),
    )
)
save_plot(p5, '05_top_GWAS_driver_genes', 9.0, 7.6)

# Combine all pages while preserving vector content.
try:
    from pypdf import PdfReader, PdfWriter
except Exception:
    from PyPDF2 import PdfReader, PdfWriter

individual_pdfs = [
    OUT / '01_integrated_program_evidence_editable.pdf',
    OUT / '02_nominal_in_vivo_perturbation_effects_editable.pdf',
    OUT / '03_common_variant_genetic_evidence_editable.pdf',
    OUT / '04_human_PAH_program_validation_editable.pdf',
    OUT / '05_top_GWAS_driver_genes_editable.pdf',
]
writer = PdfWriter()
for path in individual_pdfs:
    reader = PdfReader(str(path))
    for page in reader.pages:
        writer.add_page(page)
combined = OUT / 'selected_program_evidence_all_plots_editable.pdf'
with combined.open('wb') as handle:
    writer.write(handle)

readme = '''Selected PAH program plots - Illustrator-editable vector PDFs
================================================================

The five plots were recreated directly from the TSV files in the supplied
IP_GV_selected_program_plots.zip archive using plotnine/matplotlib.

Deliverables
------------
- 01_integrated_program_evidence_editable.pdf
- 02_nominal_in_vivo_perturbation_effects_editable.pdf
- 03_common_variant_genetic_evidence_editable.pdf
- 04_human_PAH_program_validation_editable.pdf
- 05_top_GWAS_driver_genes_editable.pdf
- selected_program_evidence_all_plots_editable.pdf (five-page combined PDF)

Editability
-----------
Lines, points, bars, labels, and axes are vector objects. Text uses embedded
TrueType fonts (PDF font type 42) rather than rasterized labels. Illustrator
may display subset font names; the source family is Arimo, an Arial-compatible
sans-serif font.

Notes
-----
The scientific content and input values are unchanged. Layout, typography,
colors, legends, and annotations were standardized for publication use.
'''
(OUT / 'README_vector_plots.txt').write_text(readme)
print(f'Created {len(individual_pdfs)} individual PDFs and {combined}')
