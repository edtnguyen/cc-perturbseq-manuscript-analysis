from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 8,
    'axes.titlesize': 11,
    'axes.labelsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'figure.dpi': 160,
    'savefig.dpi': 300,
})

IN = Path('/mnt/data/projected_cNMF_program_score_matched_contrasts.tsv')
OUT_DIR = Path('/mnt/data/selected_program_heatmap_D2toD3_TDI_IAG')
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Figure-like order from the example image. Numeric-order version is also exported below.
PROGRAMS_FIGURE_ORDER = ['P34', 'P27', 'P24', 'P2']
PROGRAMS_NUMERIC_ORDER = ['P2', 'P24', 'P27', 'P34']
PROGRAM_LABELS = {
    'P34': 'P34\nYAP/TAZ\nmesenchyme',
    'P27': 'P27\nalt.',
    'P24': 'P24\nmesenchymal',
    'P2':  'P2\narterial EC',
}
CONTRASTS = {
    'TDI': 'TDI_D2_to_D3_vs_D3_AE',
    'IAG': 'IAG_D2_to_D3_vs_D3_AE',
}

def star(p):
    if pd.isna(p): return ''
    if p < 0.001: return '***'
    if p < 0.01: return '**'
    if p < 0.05: return '*'
    return ''

def make_matrix(df, value_col, program_order):
    records=[]
    for row_label, contrast in CONTRASTS.items():
        for p in program_order:
            sub=df[(df['contrast_name']==contrast)&(df['program']==p)]
            if sub.empty:
                raise ValueError(f'Missing row for contrast={contrast}, program={p}')
            r=sub.iloc[0]
            records.append({
                'treatment': row_label,
                'contrast_name': contrast,
                'program': p,
                'program_label': PROGRAM_LABELS.get(p,p).replace('\n',' '),
                'mean_control': r['mean_control'],
                'mean_treated': r['mean_treated'],
                'diff_treated_minus_control': r['diff_treated_minus_control'],
                'log2_ratio_treated_vs_control': r['log2_ratio_treated_vs_control'],
                'pvalue': r['pvalue'],
                'padj_bh_within_contrast': r['padj_bh_within_contrast'],
                'star_pvalue': star(r['pvalue']),
                'star_padj': star(r['padj_bh_within_contrast']),
                'n_matched_pairs': r['n_matched_pairs'],
            })
    source=pd.DataFrame(records)
    mat=source.pivot(index='treatment', columns='program', values=value_col).loc[list(CONTRASTS.keys()), program_order]
    pmat=source.pivot(index='treatment', columns='program', values='pvalue').loc[list(CONTRASTS.keys()), program_order]
    return source, mat, pmat

def plot_heatmap(mat, pmat, program_order, title, cbar_label, out_prefix, fmt, symmetric=True, value_cap=None):
    vals=mat.to_numpy(dtype=float)
    plot_vals=vals.copy()
    if value_cap is not None:
        plot_vals=np.clip(plot_vals, -value_cap, value_cap)
        vmax=value_cap; vmin=-value_cap
    else:
        if symmetric:
            vmax=float(np.nanmax(np.abs(plot_vals)))
            if vmax == 0 or not np.isfinite(vmax): vmax=1e-9
            vmin=-vmax
        else:
            vmin=float(np.nanmin(plot_vals)); vmax=float(np.nanmax(plot_vals))
            if np.isclose(vmin, vmax): vmax = vmin + 1e-9
    fig, ax=plt.subplots(figsize=(7.1, 2.8))
    mesh=ax.pcolormesh(
        np.arange(mat.shape[1]+1), np.arange(mat.shape[0]+1), plot_vals,
        cmap='RdBu_r', vmin=vmin, vmax=vmax, shading='flat',
        edgecolors='white', linewidth=1.0
    )
    ax.set_xlim(0, mat.shape[1]); ax.set_ylim(mat.shape[0],0)
    ax.set_xticks(np.arange(mat.shape[1])+.5)
    ax.set_xticklabels([PROGRAM_LABELS.get(p,p) for p in program_order], rotation=40, ha='right')
    ax.set_yticks(np.arange(mat.shape[0])+.5)
    ax.set_yticklabels(mat.index.tolist(), fontsize=11)
    ax.set_xlabel('Program')
    ax.set_ylabel('Treatment')
    ax.set_title(title, pad=10, fontweight='bold')
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    # annotations show original unclipped values
    for i, treatment in enumerate(mat.index):
        for j, program in enumerate(program_order):
            v=mat.loc[treatment, program]
            pv=pmat.loc[treatment, program]
            txt=(fmt.format(v) if pd.notna(v) else 'NA') + star(pv)
            ax.text(j+.5, i+.5, txt, ha='center', va='center', fontsize=10, color='black')
    cbar=fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(cbar_label, fontsize=10)
    if value_cap is not None:
        cbar.ax.set_title(f'cap ±{value_cap:g}', fontsize=7, pad=6)
    fig.tight_layout()
    for ext in ['svg','pdf','png']:
        fig.savefig(OUT_DIR / f'{out_prefix}.{ext}', bbox_inches='tight')
    plt.close(fig)

# Load and create source tables
df=pd.read_csv(IN, sep='\t')

for program_order, order_name in [(PROGRAMS_FIGURE_ORDER,'figure_order'), (PROGRAMS_NUMERIC_ORDER,'numeric_order')]:
    source_diff, mat_diff, pmat = make_matrix(df, 'diff_treated_minus_control', program_order)
    source_log2, mat_log2, _ = make_matrix(df, 'log2_ratio_treated_vs_control', program_order)
    source_diff.to_csv(OUT_DIR / f'D2toD3_TDI_IAG_selected_program_effects_source_{order_name}.tsv', sep='\t', index=False)
    mat_diff.to_csv(OUT_DIR / f'D2toD3_TDI_IAG_selected_program_diff_matrix_{order_name}.tsv', sep='\t')
    mat_log2.to_csv(OUT_DIR / f'D2toD3_TDI_IAG_selected_program_log2ratio_matrix_{order_name}.tsv', sep='\t')

    plot_heatmap(
        mat_diff, pmat, program_order,
        title='D2→D3 treatment effects on selected D3 programs',
        cbar_label='Δ mean projected program usage\n(treated − D3 AE control)',
        out_prefix=f'D2toD3_TDI_IAG_selected_program_diff_heatmap_{order_name}',
        fmt='{:.3f}', symmetric=True, value_cap=None,
    )
    plot_heatmap(
        mat_log2, pmat, program_order,
        title='D2→D3 treatment effects on selected D3 programs',
        cbar_label='log2 ratio of projected usage\n(treated / D3 AE control)',
        out_prefix=f'D2toD3_TDI_IAG_selected_program_log2ratio_heatmap_{order_name}_capped',
        fmt='{:.2f}', symmetric=True, value_cap=8,
    )

print('Wrote outputs to', OUT_DIR)
print(source_diff.to_string(index=False))
