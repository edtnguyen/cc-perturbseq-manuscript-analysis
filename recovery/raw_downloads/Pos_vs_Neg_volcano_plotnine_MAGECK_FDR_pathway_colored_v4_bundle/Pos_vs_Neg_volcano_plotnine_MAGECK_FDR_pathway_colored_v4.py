from __future__ import annotations

from pathlib import Path
import zipfile

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from plotnine import (
    aes,
    arrow,
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
    labs,
    scale_color_identity,
    scale_x_continuous,
    scale_y_continuous,
    theme,
    theme_classic,
)

HERE = Path('/mnt/data')
INPUT_ZIP = HERE / 'Pos_vs_Neg.zip'
EXTRACT_DIR = HERE / '_Pos_vs_Neg_extracted_plotnine_pathway_v4'
OUTPUT_STEM = HERE / 'Pos_vs_Neg_volcano_plotnine_MAGECK_FDR_pathway_colored_v4'

PANEL_LABEL = 'A'
PLOT_TITLE = 'Positive versus negative fraction genome-wide screen'
HEADER_FILL = '#C8D7DF'
X_AXIS_LABEL = 'Log₂(positive / negative)'
Y_AXIS_LABEL = '−Log₁₀(MaGECK FDR)'
LEFT_ANNOTATION = 'Enriched in\nnegative fraction'
RIGHT_ANNOTATION = 'Enriched in\npositive fraction'

EFFECT_CUTOFF = 0.5
FDR_CUTOFF = 0.05
Y_CUTOFF = -np.log10(FDR_CUTOFF)
X_LIMITS = (-3.45, 2.55)
Y_LIMITS = (-0.05, 4.85)

PURPLE = '#7A43A2'
ORANGE = '#D97A34'
GREEN = '#059C48'
PINK = '#ED4A97'

PATHWAY_COLORS = {
    'SMAD1': PURPLE,
    'SMAD5': PURPLE,
    'BMPR1A': PURPLE,
    'MIXL1': ORANGE,
    'EOMES': ORANGE,
    'KDR': GREEN,
    'ETV2': GREEN,
    'PLCG1': GREEN,
    'KCTD10': GREEN,
    'PRKCD': GREEN,
    'NF2': PINK,
    'TAOK1': PINK,
    'LATS2': PINK,
    'AMOTL2': PINK,
}
HIGHLIGHT_GENES = list(PATHWAY_COLORS)

LABEL_POSITIONS = {
    'KDR': (-3.28, 4.15),
    'NF2': (-2.22, 4.56),
    'ETV2': (-2.22, 4.38),
    'PLCG1': (-2.22, 4.20),
    'SMAD1': (-2.82, 3.75),
    'SMAD5': (-2.76, 3.48),
    'BMPR1A': (-2.84, 3.21),
    'KCTD10': (-1.10, 4.60),
    'MIXL1': (-1.10, 4.38),
    'PRKCD': (-1.08, 4.18),
    'TAOK1': (0.10, 4.40),
    'LATS2': (0.10, 4.14),
    'AMOTL2': (0.10, 3.88),
    'EOMES': (-0.34, 2.18),
}

FONT_FAMILY = 'Liberation Sans'
mpl.rcParams.update({
    'font.family': FONT_FAMILY,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
    'axes.unicode_minus': True,
})

EXTRACT_DIR.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(INPUT_ZIP) as archive:
    archive.extractall(EXTRACT_DIR)

mageck_path = next(EXTRACT_DIR.rglob('Pos_vs_Neg.gene_summary.txt'))
utest_path = next(EXTRACT_DIR.rglob('combined_utest_P1P2_adj_test.txt'))

mageck = pd.read_csv(mageck_path, sep='\t')
utest = pd.read_csv(utest_path, sep='\t').rename(columns={
    'p_value': 'u_test_p_value',
    'log2FC': 'u_test_log2FC',
    'p_value_corrected': 'u_test_fdr',
})

plot_df = mageck.merge(utest, left_on='id', right_on='Target', how='outer')
plot_df['id'] = plot_df['id'].fillna(plot_df['Target'])
plot_df = plot_df.drop(columns=['Target'])
plot_df['gene'] = (
    plot_df['id'].astype(str)
    .str.replace(r'_-.*$', '', regex=True)
    .str.replace(r'-P1P2$', '', regex=True)
)
plot_df['x'] = pd.to_numeric(plot_df['u_test_log2FC'], errors='coerce')
neg_fdr = pd.to_numeric(plot_df['neg|fdr'], errors='coerce')
pos_fdr = pd.to_numeric(plot_df['pos|fdr'], errors='coerce')
plot_df['mageck_fdr'] = np.where(plot_df['x'] < 0, neg_fdr, pos_fdr)
raw_fdr = pd.to_numeric(plot_df['mageck_fdr'], errors='coerce').clip(lower=np.finfo(float).tiny)
plot_df['y'] = -np.log10(raw_fdr)
plot_df = plot_df.dropna(subset=['x', 'y']).copy()

control_names = {'non-targeting', 'safe-targeting'}
plot_df['category'] = 'Gene'
plot_df.loc[plot_df['gene'].isin(control_names), 'category'] = 'Control'
mageck_hit = (
    pd.to_numeric(plot_df['neg|fdr'], errors='coerce') < FDR_CUTOFF
) | (
    pd.to_numeric(plot_df['pos|fdr'], errors='coerce') < FDR_CUTOFF
)
plot_df.loc[mageck_hit & ~plot_df['gene'].isin(control_names), 'category'] = 'MAGeCK hit'
plot_df.loc[plot_df['gene'].isin(HIGHLIGHT_GENES), 'category'] = 'Gene highlight'

bg_df = plot_df[plot_df['category'] == 'Gene'].copy()
control_df = plot_df[plot_df['category'] == 'Control'].copy()
sig_df = plot_df[plot_df['category'] == 'MAGeCK hit'].copy()
highlight_df = plot_df[plot_df['category'] == 'Gene highlight'].copy()
for gene, color in PATHWAY_COLORS.items():
    highlight_df.loc[highlight_df['gene'] == gene, 'point_color'] = color

best_highlight = (
    highlight_df.sort_values('mageck_fdr')
    .drop_duplicates('gene')
    .set_index('gene')
)
labels_df = pd.DataFrame([
    {
        'gene': gene,
        'label_x': LABEL_POSITIONS[gene][0],
        'label_y': LABEL_POSITIONS[gene][1],
        'label_color': PATHWAY_COLORS[gene],
    }
    for gene in HIGHLIGHT_GENES if gene in best_highlight.index
])

x_min, x_max = X_LIMITS
arrow_df = pd.DataFrame([
    {'x': -1.35, 'y': 0.45, 'xend': x_min + 0.03, 'yend': 0.45},
    {'x': 0.75, 'y': 0.45, 'xend': x_max - 0.03, 'yend': 0.45},
])
arrow_text_df = pd.DataFrame([
    {'x': -2.40, 'y': 0.17, 'label': LEFT_ANNOTATION},
    {'x': 1.68, 'y': 0.17, 'label': RIGHT_ANNOTATION},
])

volcano = (
    ggplot()
    + geom_vline(xintercept=[-EFFECT_CUTOFF, EFFECT_CUTOFF], color='#E99A9A', linetype='dotted', size=1.0, alpha=0.95)
    + geom_hline(yintercept=Y_CUTOFF, color='#B3BAC3', linetype='dashed', size=1.1, alpha=1.0)
    + geom_point(bg_df, aes('x', 'y'), color='#D5D5D5', size=2.2, alpha=0.52, stroke=0, raster=True)
    + geom_point(control_df, aes('x', 'y'), color='#AEB3D2', size=2.6, alpha=0.72, stroke=0)
    + geom_point(sig_df, aes('x', 'y'), color='#111111', size=1.4, alpha=0.95, stroke=0)
    + geom_point(highlight_df, aes('x', 'y', color='point_color'), size=6.6, alpha=1.0, stroke=0)
    + geom_text(labels_df, aes(x='label_x', y='label_y', label='gene', color='label_color'), inherit_aes=False, family=FONT_FAMILY, fontstyle='italic', size=11.2)
    + geom_segment(arrow_df, aes(x='x', y='y', xend='xend', yend='yend'), inherit_aes=False, color='black', size=0.55, arrow=arrow(length=0.10, ends='last', type='closed'))
    + geom_text(arrow_text_df, aes(x='x', y='y', label='label'), inherit_aes=False, family=FONT_FAMILY, size=10.2, lineheight=0.95, color='#171717')
    + scale_color_identity()
    + scale_x_continuous(breaks=np.arange(-3, 3, 1), expand=(0, 0))
    + scale_y_continuous(breaks=np.arange(0, 5, 1), expand=(0, 0))
    + coord_cartesian(xlim=X_LIMITS, ylim=Y_LIMITS)
    + labs(x=X_AXIS_LABEL, y=Y_AXIS_LABEL)
    + theme_classic(base_family=FONT_FAMILY, base_size=14)
    + theme(
        figure_size=(12.0, 6.2),
        axis_title_x=element_text(size=17, margin={'t': 10}),
        axis_title_y=element_text(size=17, margin={'r': 10}),
        axis_text=element_text(size=13, color='#111111'),
        axis_ticks=element_line(color='#111111', size=0.8),
        axis_line=element_line(color='#111111', size=1.0),
        legend_position='none',
        panel_grid=element_blank(),
        plot_background=element_rect(fill='white', color='white'),
        panel_background=element_rect(fill='white', color='white'),
    )
)

figure = volcano.draw(show=False)
figure.subplots_adjust(top=0.82, left=0.12, right=0.985, bottom=0.17)
axis = figure.axes[0]
axis_box = axis.get_position()
header_y = axis_box.y1 + 0.02
header_height = 0.072

figure.add_artist(Rectangle((axis_box.x0, header_y), axis_box.width, header_height, transform=figure.transFigure, facecolor=HEADER_FILL, edgecolor='none', clip_on=False, zorder=2))
figure.text(axis_box.x0 + axis_box.width / 2, header_y + header_height / 2, PLOT_TITLE, ha='center', va='center', family=FONT_FAMILY, fontsize=20, color='#171717', zorder=3)
if PANEL_LABEL:
    figure.text(axis_box.x0 - 0.065, header_y + header_height / 2, PANEL_LABEL, ha='center', va='center', family=FONT_FAMILY, fontsize=22, fontweight='bold', color='#111111')

for suffix in ('png', 'pdf', 'svg'):
    output_path = OUTPUT_STEM.with_suffix(f'.{suffix}')
    figure.savefig(output_path, dpi=600, bbox_inches='tight', pad_inches=0.03, facecolor='white')
plt.close(figure)
