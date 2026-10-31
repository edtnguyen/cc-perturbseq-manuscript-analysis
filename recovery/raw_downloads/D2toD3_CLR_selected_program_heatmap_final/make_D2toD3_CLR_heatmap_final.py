from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

INPUT = Path('/mnt/data/D2toD3_CLR_selected_program_effects_source.tsv')
OUTDIR = Path('/mnt/data/D2toD3_CLR_selected_program_heatmap_final')
OUTDIR.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 9,
    'axes.titlesize': 14,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 10,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
})

PROGRAM_ORDER = ['P34', 'P27', 'P24', 'P2']
PROGRAM_LABELS = {
    'P34': 'P34\nYAP/TAZ\nmesenchyme',
    'P27': 'P27\nalt.',
    'P24': 'P24\nmesenchymal',
    'P2': 'P2\narterial EC',
}
ROW_ORDER = ['Control', 'TDI', 'IAG']


def star_from_q(q):
    if not np.isfinite(q):
        return ''
    if q < 0.001:
        return '***'
    if q < 0.01:
        return '**'
    if q < 0.05:
        return '*'
    return ''


df = pd.read_csv(INPUT, sep='\t')

# Build heatmap matrices. Control/control is a visual zero reference.
effect = pd.DataFrame(0.0, index=ROW_ORDER, columns=PROGRAM_ORDER)
qval = pd.DataFrame(np.nan, index=ROW_ORDER, columns=PROGRAM_ORDER)
pval = pd.DataFrame(np.nan, index=ROW_ORDER, columns=PROGRAM_ORDER)

for _, row in df.iterrows():
    if row['row'] in effect.index and row['program'] in effect.columns:
        effect.loc[row['row'], row['program']] = float(row['mean_diff_clr'])
        qval.loc[row['row'], row['program']] = float(row['qvalue_BH_8tests'])
        pval.loc[row['row'], row['program']] = float(row['pvalue'])

# Save source matrices.
effect.to_csv(OUTDIR / 'D2toD3_CLR_effect_matrix.tsv', sep='\t')
qval.to_csv(OUTDIR / 'D2toD3_CLR_qvalue_matrix.tsv', sep='\t')
pval.to_csv(OUTDIR / 'D2toD3_CLR_pvalue_matrix.tsv', sep='\t')

# Square-cell layout: figure width/height proportional to 4 columns / 3 rows.
fig, ax = plt.subplots(figsize=(6.6, 5.2))
mat = effect.to_numpy(dtype=float)
vmax = float(np.nanmax(np.abs(mat)))
if vmax == 0:
    vmax = 1.0

mesh = ax.pcolormesh(
    np.arange(mat.shape[1] + 1),
    np.arange(mat.shape[0] + 1),
    mat,
    cmap='RdBu_r',
    vmin=-vmax,
    vmax=vmax,
    shading='flat',
    edgecolors='white',
    linewidth=1.6,
)

ax.set_xlim(0, mat.shape[1])
ax.set_ylim(mat.shape[0], 0)
ax.set_aspect('equal')
ax.set_xticks(np.arange(mat.shape[1]) + 0.5)
ax.set_yticks(np.arange(mat.shape[0]) + 0.5)
ax.set_xticklabels([PROGRAM_LABELS[p] for p in PROGRAM_ORDER], rotation=42, ha='right')
ax.set_yticklabels(ROW_ORDER)
ax.set_xlabel('Projected perturb-seq cNMF program')
ax.set_ylabel('Treatment')
ax.set_title('D2→D3 treatment effect on selected D3 programs', pad=10, weight='bold')
ax.tick_params(length=0)
for spine in ax.spines.values():
    spine.set_visible(False)

# Annotate CLR effects and BH-adjusted significance stars.
for i, row_name in enumerate(ROW_ORDER):
    for j, program in enumerate(PROGRAM_ORDER):
        val = effect.loc[row_name, program]
        star = '' if row_name == 'Control' else star_from_q(qval.loc[row_name, program])
        label = f'{val:.2f}{star}'
        # White text on darker cells, black otherwise.
        text_color = 'white' if abs(val) > 0.60 * vmax else 'black'
        ax.text(j + 0.5, i + 0.5, label, ha='center', va='center', color=text_color, fontsize=10)

cbar = fig.colorbar(mesh, ax=ax, fraction=0.047, pad=0.04)
cbar.set_label('Δ CLR (treated − matched control)')

fig.text(
    0.99,
    0.01,
    'Paired t-test across 3 matched replicates; BH correction across 8 treatment–program tests.  *q<0.05',
    ha='right',
    va='bottom',
    fontsize=6.5,
)

fig.subplots_adjust(left=0.16, right=0.90, top=0.88, bottom=0.28)

stem = OUTDIR / 'D2toD3_TDI_IAG_selected_program_CLR_heatmap_control_square_BHstars'
for ext in ['svg', 'pdf', 'png']:
    fig.savefig(stem.with_suffix('.' + ext), dpi=300 if ext == 'png' else None, bbox_inches='tight')
plt.close(fig)

# Also create a compact statistics table for publication/source tracking.
stats_cols = [
    'row', 'treated_condition', 'control_condition', 'program', 'n_matched_pairs',
    'mean_control_clr', 'mean_treated_clr', 'mean_diff_clr', 'sd_diff_clr',
    'paired_t_statistic', 'pvalue', 'qvalue_BH_8tests',
    'n_positive_diffs', 'n_negative_diffs',
]
df[stats_cols].to_csv(OUTDIR / 'D2toD3_CLR_heatmap_statistics.tsv', sep='\t', index=False)

print('Wrote outputs to', OUTDIR)
print(effect)
print('\nMinimum q-value:', df['qvalue_BH_8tests'].min())
