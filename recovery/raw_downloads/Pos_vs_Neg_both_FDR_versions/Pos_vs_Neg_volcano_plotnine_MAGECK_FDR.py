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
    geom_point,
    geom_segment,
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

HERE = Path('/mnt/data')
INPUT_ZIP = HERE / 'Pos_vs_Neg.zip'
EXTRACT_DIR = HERE / '_Pos_vs_Neg_extracted_MAGECK_FDR'
OUTPUT_STEM = HERE / 'Pos_vs_Neg_volcano_plotnine_MAGECK_FDR'

PANEL_LABEL = 'A'
PLOT_TITLE = 'Positive versus negative fraction genome-wide screen'
HEADER_FILL = '#C8D7DF'
X_AXIS_LABEL = 'Log₂(positive / negative)'
Y_AXIS_LABEL = '−Log₁₀(MAGeCK FDR)'
LEFT_ANNOTATION = 'Enriched in\nnegative fraction'
RIGHT_ANNOTATION = 'Enriched in\npositive fraction'

EFFECT_CUTOFF = 0.5
MAGECK_FDR_CUTOFF = 0.05
X_LIMITS = (-3.45, 2.55)
Y_LIMITS = (-0.05, 4.85)

HIGHLIGHT_GENES = [
    'SMAD1', 'SMAD5', 'ACVR2A', 'LEMD3', 'BMPR1A', 'PLCG1', 'CDH5',
    'NF2', 'ETV2', 'RUVBL2', 'CBLL1', 'WDR12', 'CENPC', 'PLRG1',
]

LABEL_POSITIONS = {
    'SMAD1': (-3.17, 4.15),
    'SMAD5': (-2.92, 3.86),
    'ACVR2A': (-2.72, 3.03),
    'LEMD3': (-2.60, 4.32),
    'BMPR1A': (-2.48, 3.52),
    'PLCG1': (-2.18, 4.03),
    'CDH5': (-1.94, 4.38),
    'ETV2': (-1.62, 3.84),
    'NF2': (-1.50, 3.20),
    'RUVBL2': (1.00, 2.92),
    'CBLL1': (1.12, 4.02),
    'WDR12': (1.55, 3.62),
    'CENPC': (2.02, 4.30),
    'PLRG1': (2.10, 2.75),
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
    pd.to_numeric(plot_df['neg|fdr'], errors='coerce') < MAGECK_FDR_CUTOFF
) | (
    pd.to_numeric(plot_df['pos|fdr'], errors='coerce') < MAGECK_FDR_CUTOFF
)
plot_df.loc[mageck_hit & ~plot_df['gene'].isin(control_names), 'category'] = 'MAGeCK hit'
plot_df.loc[plot_df['gene'].isin(HIGHLIGHT_GENES), 'category'] = 'Gene highlight'

draw_order = {'Gene': 0, 'Control': 1, 'MAGeCK hit': 2, 'Gene highlight': 3}
plot_df['_draw_order'] = plot_df['category'].map(draw_order)
plot_df = plot_df.sort_values(['_draw_order', 'y'])

best_highlight = (
    plot_df[plot_df['gene'].isin(HIGHLIGHT_GENES)]
    .sort_values('mageck_fdr')
    .drop_duplicates('gene')
    .set_index('gene')
)
label_rows = []
for gene in HIGHLIGHT_GENES:
    if gene not in best_highlight.index:
        continue
    row = best_highlight.loc[gene]
    if gene in LABEL_POSITIONS:
        label_x, label_y = LABEL_POSITIONS[gene]
    else:
        label_x = float(row['x'] + (-0.22 if row['x'] < 0 else 0.22))
        label_y = float(min(row['y'] + 0.20, Y_LIMITS[1] - 0.08))
    label_rows.append({'gene': gene, 'x': float(row['x']), 'y': float(row['y']), 'label_x': label_x, 'label_y': label_y})
labels_df = pd.DataFrame(label_rows)

x_min, x_max = X_LIMITS
arrow_df = pd.DataFrame([
    {'x': -1.35, 'y': 0.45, 'xend': x_min + 0.03, 'yend': 0.45},
    {'x': 0.75, 'y': 0.45, 'xend': x_max - 0.03, 'yend': 0.45},
])
arrow_text_df = pd.DataFrame([
    {'x': -2.40, 'y': 0.17, 'label': LEFT_ANNOTATION},
    {'x': 1.68, 'y': 0.17, 'label': RIGHT_ANNOTATION},
])

legend_order = ['Gene highlight', 'Gene', 'Control', 'MAGeCK hit']
colors = {
    'Gene highlight': '#8F0000',
    'Gene': '#D5D5D5',
    'Control': '#AEB3D2',
    'MAGeCK hit': '#111111',
}
sizes = {'Gene highlight': 4.0, 'Gene': 2.2, 'Control': 2.6, 'MAGeCK hit': 0.75}
alphas = {'Gene highlight': 1.0, 'Gene': 0.52, 'Control': 0.72, 'MAGeCK hit': 0.90}

volcano = (
    ggplot(plot_df, aes('x', 'y'))
    + geom_vline(xintercept=[-EFFECT_CUTOFF, EFFECT_CUTOFF], color='#E99A9A', linetype='dotted', size=1.0, alpha=0.95)
    + geom_point(aes(color='category', size='category', alpha='category'), stroke=0)
    + geom_segment(labels_df, aes(x='x', y='y', xend='label_x', yend='label_y'), inherit_aes=False, color='black', size=0.45)
    + geom_text(labels_df, aes(x='label_x', y='label_y', label='gene'), inherit_aes=False, family=FONT_FAMILY, fontstyle='italic', size=10.5, color='#151515')
    + geom_segment(arrow_df, aes(x='x', y='y', xend='xend', yend='yend'), inherit_aes=False, color='black', size=0.55, arrow=arrow(length=0.10, ends='last', type='closed'))
    + geom_text(arrow_text_df, aes(x='x', y='y', label='label'), inherit_aes=False, family=FONT_FAMILY, size=10.2, lineheight=0.95, color='#171717')
    + scale_color_manual(values=colors, breaks=legend_order)
    + scale_size_manual(values=sizes, breaks=legend_order, guide=None)
    + scale_alpha_manual(values=alphas, breaks=legend_order, guide=None)
    + scale_x_continuous(breaks=np.arange(-3, 3, 1), expand=(0, 0))
    + scale_y_continuous(breaks=np.arange(0, 5, 1), expand=(0, 0))
    + coord_cartesian(xlim=X_LIMITS, ylim=Y_LIMITS)
    + labs(x=X_AXIS_LABEL, y=Y_AXIS_LABEL, color='')
    + theme_classic(base_family=FONT_FAMILY, base_size=14)
    + theme(
        figure_size=(12.0, 6.2),
        axis_title_x=element_text(size=17, margin={'t': 10}),
        axis_title_y=element_text(size=17, margin={'r': 10}),
        axis_text=element_text(size=13, color='#111111'),
        axis_ticks=element_line(color='#111111', size=0.8),
        axis_line=element_line(color='#111111', size=1.0),
        legend_position=(0.83, 0.33),
        legend_background=element_rect(fill='white', color='#333333', size=0.5),
        legend_key=element_blank(),
        legend_text=element_text(size=10.5, color='#111111'),
        legend_title=element_blank(),
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

figure.add_artist(Rectangle((axis_box.x0, header_y), axis_box.width, header_height,
                             transform=figure.transFigure, facecolor=HEADER_FILL,
                             edgecolor='none', clip_on=False, zorder=2))
figure.text(axis_box.x0 + axis_box.width / 2, header_y + header_height / 2,
            PLOT_TITLE, ha='center', va='center', family=FONT_FAMILY,
            fontsize=20, color='#171717', zorder=3)
if PANEL_LABEL:
    figure.text(axis_box.x0 - 0.065, header_y + header_height / 2, PANEL_LABEL,
                ha='center', va='center', family=FONT_FAMILY, fontsize=22,
                fontweight='bold', color='#111111')

for suffix in ('png', 'pdf', 'svg'):
    output_path = OUTPUT_STEM.with_suffix(f'.{suffix}')
    figure.savefig(output_path, dpi=600, bbox_inches='tight', pad_inches=0.03, facecolor='white')
plt.close(figure)
