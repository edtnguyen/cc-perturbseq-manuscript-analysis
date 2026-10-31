from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

IN = Path('/mnt/data/projected_cNMF_program_score_matched_contrasts.tsv')
OUTDIR = Path('/mnt/data/D2toD3_TDI_IAG_selected_program_heatmap')
OUTDIR.mkdir(exist_ok=True)

# Editable text in Illustrator
mpl.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 9,
    'axes.titlesize': 14,
    'axes.labelsize': 11,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'figure.dpi': 180,
    'savefig.dpi': 300,
})

df = pd.read_csv(IN, sep='\t')

contrasts = {
    'TDI': 'TDI_D2_to_D3_vs_D3_AE',
    'IAG': 'IAG_D2_to_D3_vs_D3_AE',
}
programs_requested = ['P2', 'P24', 'P27', 'P34']
programs_figure_like = ['P34', 'P27', 'P24', 'P2']
program_labels = {
    'P2': 'P2\narterial EC',
    'P24': 'P24\nmesenchymal',
    'P27': 'P27\nalt.',
    'P34': 'P34\nYAP/TAZ\nmesenchyme',
}

records = []
for treatment, contrast in contrasts.items():
    for program in programs_requested:
        row = df[(df['contrast_name'] == contrast) & (df['program'] == program)]
        if row.empty:
            raise ValueError(f'Missing {contrast} {program}')
        r = row.iloc[0].to_dict()
        records.append({
            'treatment': treatment,
            'contrast_name': contrast,
            'program': program,
            'program_label': program_labels[program].replace('\n', ' '),
            'mean_control': r['mean_control'],
            'mean_treated': r['mean_treated'],
            'delta_mean_usage': r['diff_treated_minus_control'],
            'delta_mean_usage_pct_points': 100 * r['diff_treated_minus_control'],
            'log2_ratio_treated_vs_control': r['log2_ratio_treated_vs_control'],
            'pvalue': r['pvalue'],
            'padj_bh_within_contrast': r['padj_bh_within_contrast'],
            'n_matched_pairs': r['n_matched_pairs'],
        })
source = pd.DataFrame(records)
source.to_csv(OUTDIR / 'D2toD3_TDI_IAG_selected_program_effects_source.tsv', sep='\t', index=False)

def star(padj, pval):
    # Use adjusted p-value for star; fall back to none if missing.
    if pd.notna(padj) and padj < 0.05:
        return '*'
    return ''

def make_matrix(value_col, program_order):
    mat = source.pivot(index='treatment', columns='program', values=value_col).loc[['TDI','IAG'], program_order]
    pmat = source.pivot(index='treatment', columns='program', values='pvalue').loc[['TDI','IAG'], program_order]
    qmat = source.pivot(index='treatment', columns='program', values='padj_bh_within_contrast').loc[['TDI','IAG'], program_order]
    return mat, pmat, qmat

def plot_heatmap(mat, pmat, qmat, program_order, title, cbar_label, stem, fmt='{:.2f}', symmetric=True):
    vals = mat.to_numpy(dtype=float)
    if symmetric:
        vmax = np.nanmax(np.abs(vals))
        vmax = vmax if np.isfinite(vmax) and vmax > 0 else 1
        vmin = -vmax
    else:
        vmin = np.nanmin(vals)
        vmax = np.nanmax(vals)
    fig_w = 1.25 * len(program_order) + 2.7
    fig_h = 2.7
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    mesh = ax.pcolormesh(
        np.arange(vals.shape[1] + 1),
        np.arange(vals.shape[0] + 1),
        vals,
        cmap='RdBu_r',
        vmin=vmin,
        vmax=vmax,
        shading='flat',
        linewidth=1,
        edgecolors='white',
    )
    ax.set_xlim(0, vals.shape[1])
    ax.set_ylim(vals.shape[0], 0)
    ax.set_xticks(np.arange(vals.shape[1]) + 0.5)
    ax.set_yticks(np.arange(vals.shape[0]) + 0.5)
    ax.set_xticklabels([program_labels[p] for p in program_order], rotation=45, ha='right')
    ax.set_yticklabels(mat.index.tolist())
    ax.set_title(title, fontweight='bold', pad=14)
    ax.set_xlabel('Projected perturb-seq cNMF program')
    ax.set_ylabel('Small molecule')
    for i, trt in enumerate(mat.index):
        for j, program in enumerate(program_order):
            val = mat.loc[trt, program]
            txt = fmt.format(val) + star(qmat.loc[trt, program], pmat.loc[trt, program])
            ax.text(j + 0.5, i + 0.5, txt, ha='center', va='center', fontsize=10, color='black')
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.03, fraction=0.05)
    cbar.set_label(cbar_label)
    fig.subplots_adjust(left=0.15, bottom=0.30, right=0.88, top=0.82)
    for ext in ['svg','pdf','png']:
        fig.savefig(OUTDIR / f'{stem}.{ext}', bbox_inches='tight')
    plt.close(fig)

# Recommended plot: delta projected usage in percentage points.
for order_name, order in [('requested_order', programs_requested), ('figure_like_order', programs_figure_like)]:
    mat, pmat, qmat = make_matrix('delta_mean_usage_pct_points', order)
    mat.to_csv(OUTDIR / f'D2toD3_TDI_IAG_delta_usage_pct_points_matrix_{order_name}.tsv', sep='\t')
    plot_heatmap(
        mat, pmat, qmat, order,
        title='D2→D3 treatment effect on selected D3 programs',
        cbar_label='Δ mean projected usage\npercentage points',
        stem=f'D2toD3_TDI_IAG_selected_program_heatmap_delta_usage_pct_points_{order_name}',
        fmt='{:.2f}',
        symmetric=True,
    )

# Secondary plot: log2 ratio; useful visually but can be unstable near zero.
for order_name, order in [('requested_order', programs_requested), ('figure_like_order', programs_figure_like)]:
    mat, pmat, qmat = make_matrix('log2_ratio_treated_vs_control', order)
    mat.to_csv(OUTDIR / f'D2toD3_TDI_IAG_log2_ratio_matrix_{order_name}.tsv', sep='\t')
    plot_heatmap(
        mat, pmat, qmat, order,
        title='D2→D3 treatment effect on selected D3 programs',
        cbar_label='log2 ratio\ntreated / control',
        stem=f'D2toD3_TDI_IAG_selected_program_heatmap_log2_ratio_{order_name}',
        fmt='{:.2f}',
        symmetric=True,
    )

# README
(OUTDIR / 'README.md').write_text('''# D2→D3 TDI/IAG selected-program heatmap\n\nInputs: projected_cNMF_program_score_matched_contrasts.tsv\n\nRows: TDI and IAG transition treatments compared to D3_AE matched control.\nColumns: P2, P24, P27, P34.\n\nRecommended main figure:\n- `D2toD3_TDI_IAG_selected_program_heatmap_delta_usage_pct_points_requested_order.svg`\n\nThis plots `100 * diff_treated_minus_control`, i.e. the change in mean projected normalized program usage in percentage points.\n\nAdditional log2-ratio versions are included for reference, but interpret them cautiously when control or treated usage is close to zero.\n\nAsterisks are drawn only for adjusted p-value < 0.05. No cells in this uploaded table pass that threshold for these programs.\n''')
print('Wrote outputs to', OUTDIR)
print(source)
