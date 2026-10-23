from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

mpl.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 7,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'figure.dpi': 150,
    'savefig.dpi': 300,
})

IN = Path('/mnt/data')
OUT = Path('/mnt/data/marker_gene_control_centered_plots')
OUT.mkdir(exist_ok=True)

rows = pd.read_csv(IN/'representative_gene_heatmap_rows.tsv', sep='\t')
rep_expr = pd.read_csv(IN/'representative_gene_log2CPM_by_sample.replicates.tsv', sep='\t', index_col=0)

row_order = rows['display_label'].tolist()
rep_expr = rep_expr.loc[row_order]
row_to_prog = rows.set_index('display_label')['program'].to_dict()
row_to_gene = rows.set_index('display_label')['gene'].to_dict()
program_labels = dict(rows[['program','program_annotation']].drop_duplicates().values)
program_order = ['P15','P23','P17','P2','P34','P24','P33']

# Contrasts are explicitly matched to same-stage control.
contrasts = [
    ('D2_LM+TDI','D2_LM','D2 TDI'),
    ('D2_LM+IAG','D2_LM','D2 IAG'),
    ('D3_AE+TDI','D3_AE','D3 TDI'),
    ('D3_AE+IAG','D3_AE','D3 IAG'),
    ('D2_LM+TDI->D3_AE+TDI','D3_AE','D2→D3 TDI'),
    ('D2_LM+IAG->D3_AE+IAG','D3_AE','D2→D3 IAG'),
]
contrast_order = [x[2] for x in contrasts]

# For condition-expression plots only, not z-scored.
focused_conditions = [
    'D2_LM','D2_LM+TDI','D2_LM+IAG',
    'D3_AE','D3_AE+TDI','D3_AE+IAG',
    'D2_LM+TDI->D3_AE+TDI','D2_LM+IAG->D3_AE+IAG'
]
focused_sample_cols = []
for cond in focused_conditions:
    for rep in [1,2,3]:
        col=f'{cond}_rep{rep}'
        if col in rep_expr.columns:
            focused_sample_cols.append(col)


def save(fig, name):
    for ext in ['svg','pdf','png']:
        fig.savefig(OUT/f'{name}.{ext}', bbox_inches='tight', dpi=300 if ext=='png' else None)
    plt.close(fig)


def symmetric_vlim(x):
    x=np.asarray(x,dtype=float)
    v=np.nanmax(np.abs(x)) if x.size else 1
    return v if np.isfinite(v) and v>0 else 1


def plot_heatmap(df, title, name, cmap='RdBu_r', center_zero=True, xlabel='', ylabel='', annotate=True, y_fontsize=7, x_rotation=45, figsize=None):
    mat=df.to_numpy(dtype=float)
    if figsize is None:
        figsize=(max(7, 0.6*df.shape[1]+3), max(4, 0.28*df.shape[0]+1.5))
    fig, ax=plt.subplots(figsize=figsize)
    if center_zero:
        v=symmetric_vlim(mat); vmin=-v; vmax=v
    else:
        vmin=np.nanmin(mat); vmax=np.nanmax(mat)
    mesh=ax.pcolormesh(np.arange(df.shape[1]+1), np.arange(df.shape[0]+1), mat, cmap=cmap, vmin=vmin, vmax=vmax, shading='flat', linewidth=0, antialiased=False)
    ax.set_xlim(0,df.shape[1]); ax.set_ylim(df.shape[0],0)
    ax.set_xticks(np.arange(df.shape[1])+0.5); ax.set_xticklabels(df.columns, rotation=x_rotation, ha='right')
    ax.set_yticks(np.arange(df.shape[0])+0.5); ax.set_yticklabels(df.index, fontsize=y_fontsize)
    ax.set_title(title); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    for sp in ax.spines.values(): sp.set_visible(False)
    if annotate:
        for i in range(df.shape[0]):
            for j in range(df.shape[1]):
                val=mat[i,j]
                if np.isfinite(val):
                    ax.text(j+0.5,i+0.5,f'{val:+.2f}',ha='center',va='center',fontsize=5)
    cbar=fig.colorbar(mesh, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label('Δ log2(CPM+1) vs matched control' if center_zero else 'log2(CPM+1)')
    fig.tight_layout()
    save(fig,name)

# ------------------------------------------------------------------
# Gene-level paired deltas: no z-score, directly treated - matching control.
# ------------------------------------------------------------------
gene_delta_records=[]
for treated, control, cname in contrasts:
    for label in row_order:
        prog=row_to_prog[label]
        gene=row_to_gene[label]
        for rep in [1,2,3]:
            ccol=f'{control}_rep{rep}'
            tcol=f'{treated}_rep{rep}'
            if ccol not in rep_expr.columns or tcol not in rep_expr.columns:
                continue
            cval=float(rep_expr.loc[label,ccol])
            tval=float(rep_expr.loc[label,tcol])
            gene_delta_records.append({
                'contrast': cname,
                'treated_condition': treated,
                'control_condition': control,
                'program': prog,
                'program_label': program_labels.get(prog,prog),
                'display_label': label,
                'gene': gene,
                'replicate': rep,
                'control_log2cpm': cval,
                'treated_log2cpm': tval,
                'delta_log2cpm': tval-cval,
            })
gene_delta_df=pd.DataFrame(gene_delta_records)
gene_delta_summary=gene_delta_df.groupby(['contrast','treated_condition','control_condition','program','program_label','display_label','gene'],as_index=False).agg(mean_delta_log2cpm=('delta_log2cpm','mean'),sd_delta_log2cpm=('delta_log2cpm','std'),n=('delta_log2cpm','size'))
gene_delta_df.to_csv(OUT/'marker_gene_delta_log2cpm_by_replicate_control_centered.tsv',sep='\t',index=False)
gene_delta_summary.to_csv(OUT/'marker_gene_delta_log2cpm_summary_control_centered.tsv',sep='\t',index=False)

# Gene-level summary heatmap: rows genes, columns contrasts, values mean paired Δlog2CPM
mat_gene=gene_delta_summary.pivot(index='display_label',columns='contrast',values='mean_delta_log2cpm').reindex(row_order)[contrast_order]
mat_gene.to_csv(OUT/'marker_gene_delta_log2cpm_summary_matrix_control_centered.tsv',sep='\t')
plot_heatmap(mat_gene, 'Representative marker genes: matched treatment effect\nvalues = mean paired Δ log2(CPM+1), no z-score', 'A_marker_gene_delta_heatmap_control_centered', ylabel='marker gene', annotate=True, figsize=(8.5,7.0))

# Gene-level replicate heatmap: columns contrast_repN, values paired deltas
cols=[]
for cname in contrast_order:
    for rep in [1,2,3]:
        cols.append(f'{cname} rep{rep}')
gene_rep_mat=pd.DataFrame(index=row_order, columns=cols, dtype=float)
for _,r in gene_delta_df.iterrows():
    col=f"{r['contrast']} rep{int(r['replicate'])}"
    gene_rep_mat.loc[r['display_label'], col]=r['delta_log2cpm']
gene_rep_mat.to_csv(OUT/'marker_gene_delta_log2cpm_by_replicate_matrix_control_centered.tsv',sep='\t')
plot_heatmap(gene_rep_mat, 'Representative marker genes by matched replicate\nvalues = treated_rep − control_rep log2(CPM+1), no z-score', 'B_marker_gene_delta_by_replicate_heatmap_control_centered', ylabel='marker gene', annotate=False, x_rotation=60, figsize=(11.5,7.0))

# ------------------------------------------------------------------
# Program marker modules: mean of gene-level paired deltas, not z-scored.
# ------------------------------------------------------------------
module_delta_df=gene_delta_df.groupby(['contrast','treated_condition','control_condition','program','program_label','replicate'],as_index=False).agg(mean_delta_log2cpm=('delta_log2cpm','mean'), n_genes=('gene','nunique'))
module_delta_summary=module_delta_df.groupby(['contrast','treated_condition','control_condition','program','program_label'],as_index=False).agg(mean_delta_log2cpm=('mean_delta_log2cpm','mean'),sd_delta_log2cpm=('mean_delta_log2cpm','std'),n_reps=('mean_delta_log2cpm','size'),n_genes=('n_genes','max'))
module_delta_df.to_csv(OUT/'marker_program_module_delta_log2cpm_by_replicate_control_centered.tsv',sep='\t',index=False)
module_delta_summary.to_csv(OUT/'marker_program_module_delta_log2cpm_summary_control_centered.tsv',sep='\t',index=False)

program_labels_pretty=[f"{p} — {program_labels.get(p,'')}" for p in program_order]
mat_mod=module_delta_summary.pivot(index='program',columns='contrast',values='mean_delta_log2cpm').reindex(program_order)[contrast_order]
mat_mod.index=program_labels_pretty
mat_mod.to_csv(OUT/'marker_program_module_delta_log2cpm_summary_matrix_control_centered.tsv',sep='\t')
plot_heatmap(mat_mod, 'Representative marker modules: matched treatment effect\nvalues = mean paired Δ log2(CPM+1), no z-score', 'C_marker_program_module_delta_heatmap_control_centered', ylabel='marker module', annotate=True, figsize=(8.2,4.7))

# Module replicate heatmap
mod_cols=cols
mod_rep_mat=pd.DataFrame(index=program_labels_pretty, columns=mod_cols, dtype=float)
prog_to_label=dict(zip(program_order,program_labels_pretty))
for _,r in module_delta_df.iterrows():
    col=f"{r['contrast']} rep{int(r['replicate'])}"
    mod_rep_mat.loc[prog_to_label[r['program']], col]=r['mean_delta_log2cpm']
mod_rep_mat.to_csv(OUT/'marker_program_module_delta_log2cpm_by_replicate_matrix_control_centered.tsv',sep='\t')
plot_heatmap(mod_rep_mat, 'Representative marker modules by matched replicate\nvalues = treated_rep − control_rep mean log2(CPM+1), no z-score', 'D_marker_program_module_delta_by_replicate_heatmap_control_centered', ylabel='marker module', annotate=True, x_rotation=60, figsize=(11.5,4.8))

# ------------------------------------------------------------------
# Paired slope plots with untransformed log2CPM module scores.
# Module expression score = mean log2CPM of representative genes per program/sample.
# ------------------------------------------------------------------
module_expr_records=[]
for prog in program_order:
    genes=[label for label in row_order if row_to_prog[label]==prog]
    if not genes: continue
    vals=rep_expr.loc[genes].mean(axis=0)
    for sample,v in vals.items():
        m=re.match(r'(.+)_rep(\d+)$', str(sample))
        if m:
            cond=m.group(1); rep=int(m.group(2))
        else:
            cond=str(sample); rep=np.nan
        module_expr_records.append({'program':prog,'program_label':program_labels.get(prog,prog),'sample_name':sample,'condition':cond,'replicate':rep,'mean_marker_log2cpm':float(v),'n_genes':len(genes)})
module_expr_df=pd.DataFrame(module_expr_records)
module_expr_df.to_csv(OUT/'marker_program_module_mean_log2cpm_by_sample.tsv',sep='\t',index=False)

key_panels=[
    ('P15','D2_LM','D2_LM+TDI','P15 D2 TDI'),
    ('P23','D2_LM','D2_LM+TDI','P23 D2 TDI'),
    ('P17','D2_LM','D2_LM+TDI','P17 D2 TDI'),
    ('P2','D3_AE','D3_AE+TDI','P2 D3 TDI'),
    ('P34','D3_AE','D3_AE+TDI','P34 D3 TDI'),
    ('P24','D3_AE','D2_LM+TDI->D3_AE+TDI','P24 D2→D3 TDI'),
]
fig, axes=plt.subplots(2,3,figsize=(8.4,5.2))
for ax,(prog,control,treated,title) in zip(axes.ravel(),key_panels):
    sub=module_expr_df[module_expr_df.program==prog]
    for rep in [1,2,3]:
        c=sub[(sub.condition==control)&(sub.replicate==rep)]
        t=sub[(sub.condition==treated)&(sub.replicate==rep)]
        if c.empty or t.empty: continue
        ax.plot([0,1],[c.mean_marker_log2cpm.iloc[0],t.mean_marker_log2cpm.iloc[0]], marker='o', lw=1.2, alpha=0.9)
    ax.set_xticks([0,1]); ax.set_xticklabels([control,treated], rotation=25, ha='right', fontsize=6)
    ax.set_ylabel('mean marker\nlog2(CPM+1)', fontsize=6)
    ax.set_title(f'{title}\n{program_labels.get(prog,prog)}', fontsize=7)
    ax.spines[['top','right']].set_visible(False)
fig.suptitle('Paired marker module expression, compared to matching control\nno across-sample z-score', y=1.02)
fig.tight_layout(); save(fig,'E_key_marker_module_paired_slopeplots_raw_log2cpm_control_centered')

# ------------------------------------------------------------------
# Collapsed category plot from control-centered deltas.
# ------------------------------------------------------------------
cat_map={'EC/progenitor markers':['P15','P23','P2'], 'Alternative/YAP markers':['P17','P34','P24','P33']}
cat_delta=[]
for cat,progs in cat_map.items():
    sub=module_delta_df[module_delta_df.program.isin(progs)]
    tmp=sub.groupby(['contrast','replicate'],as_index=False).mean_delta_log2cpm.mean()
    tmp['category']=cat
    cat_delta.append(tmp)
cat_delta=pd.concat(cat_delta,ignore_index=True)
cat_delta.to_csv(OUT/'marker_category_delta_log2cpm_by_replicate_control_centered.tsv',sep='\t',index=False)
fig, axes=plt.subplots(1,2,figsize=(8.6,3.4),sharey=True)
for ax,cat in zip(axes,cat_map):
    sub=cat_delta[cat_delta.category==cat]
    for i,cname in enumerate(contrast_order):
        vals=sub[sub.contrast==cname].mean_delta_log2cpm.values
        xs=np.full(len(vals),i)+(np.arange(len(vals))-1)*0.06
        ax.scatter(xs,vals,s=25,alpha=0.9)
        if len(vals): ax.plot([i-0.22,i+0.22],[vals.mean(),vals.mean()],color='black',lw=1.5)
    ax.axhline(0,color='0.6',lw=0.8)
    ax.set_title(cat)
    ax.set_xticks(range(len(contrast_order))); ax.set_xticklabels(contrast_order,rotation=45,ha='right')
    ax.spines[['top','right']].set_visible(False)
axes[0].set_ylabel('paired Δ mean marker log2(CPM+1)')
fig.suptitle('Collapsed marker classes: matched treatment effect\nno across-sample z-score')
fig.tight_layout(); save(fig,'F_collapsed_marker_category_delta_log2cpm_control_centered')

# README
(OUT/'README.md').write_text('''# Control-centered marker gene plots\n\nThese plots do not z-score across samples. They compare treated samples to the matching control at the same time point.\n\nGene-level effects are paired differences:\n  treated_rep - control_rep in log2(CPM+1).\n\nProgram marker module effects are the average of gene-level paired deltas for representative genes within each program.\n\nSVG files have editable text.\n''')
print(OUT)
