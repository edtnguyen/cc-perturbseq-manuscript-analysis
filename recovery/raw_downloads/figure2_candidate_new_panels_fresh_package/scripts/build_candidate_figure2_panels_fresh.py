import os, math, zipfile, shutil
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from matplotlib.lines import Line2D

BASE = '/mnt/data'
OUT = os.path.join(BASE, 'figure2_candidate_new_panels_fresh')
os.makedirs(OUT, exist_ok=True)
os.makedirs(os.path.join(OUT, 'source_data'), exist_ok=True)
os.makedirs(os.path.join(OUT, 'scripts'), exist_ok=True)
os.makedirs(os.path.join(OUT, 'individual_panels'), exist_ok=True)

XLSX = os.path.join(BASE, 'cNMF_50_0_2_032526_6co(1).xlsx')
ANN = os.path.join(BASE, 'figure1_annotation_literature_draft', 'screen_regulators_244_literature_annotation_draft.tsv')

mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'DejaVu Sans'

# Selected individual programs only; no composite/contrast scores.
programs = [
    ('D1', 30, 'P30'),
    ('D1', 13, 'P13'),
    ('D2', 15, 'P15'),
    ('D2', 17, 'P17'),
    ('D3', 3,  'P3'),
    ('D3', 34, 'P34'),
]

assoc = pd.read_excel(
    XLSX,
    sheet_name='Perturbation Association 1',
    usecols=['target_name','program_name','adj_pval','approx_log2FC','Sample']
)
ann = pd.read_csv(ANN, sep='\t')

# Map the literature annotation to concise mechanistic groups.
def class_group(row):
    primary = str(row['primary_category'])
    secondary = str(row['secondary_category'])
    if primary == 'Developmental signaling & lineage regulation':
        if 'Hippo' in secondary:
            return 'Hippo / YAP'
        if 'TGF' in secondary or 'BMP' in secondary:
            return 'BMP / TGFβ'
        return 'EC / mesoderm fate & signaling'
    if primary == 'Transcription & chromatin regulation':
        return 'Transcription / chromatin'
    if primary == 'RNA processing & surveillance':
        return 'RNA processing'
    if primary == 'Translation, ribosome & protein targeting':
        return 'Translation / protein targeting'
    if primary == 'Proteostasis & post-translational modification':
        return 'Proteostasis / PTM'
    if primary == 'Cell cycle, chromosome & DNA repair' or primary == 'Apoptosis & stress response':
        return 'Cell cycle / DNA / stress'
    if primary == 'Membrane trafficking, cytoskeleton & adhesion':
        return 'Trafficking / cytoskeleton / adhesion'
    if primary == 'Metabolism, mitochondria & glycosylation':
        return 'Metabolism / mitochondria'
    return 'Other'

ann['class_group'] = ann.apply(class_group, axis=1)
ann = ann[['gene','class_group']].drop_duplicates()

category_order = [
    'EC / mesoderm fate & signaling',
    'BMP / TGFβ',
    'Hippo / YAP',
    'Transcription / chromatin',
    'RNA processing',
    'Translation / protein targeting',
    'Proteostasis / PTM',
    'Cell cycle / DNA / stress',
    'Trafficking / cytoskeleton / adhesion',
    'Metabolism / mitochondria',
]

# Universe = all 244 screen-derived regulators with reviewed annotations.
universe = ann[ann['class_group'].isin(category_order)].copy()

rows = []
for day, prog, prog_label in programs:
    sub = assoc[(assoc['Sample'] == day) & (assoc['program_name'] == prog)].copy()
    sub = sub.merge(universe, left_on='target_name', right_on='gene', how='inner')
    sub['sig_positive'] = (sub['adj_pval'] < 0.05) & (sub['approx_log2FC'] > 0)
    sub['sig_negative'] = (sub['adj_pval'] < 0.05) & (sub['approx_log2FC'] < 0)
    for direction, sig_col in [('Positive effect','sig_positive'),('Negative effect','sig_negative')]:
        sig = sub[sig_col].astype(bool)
        for cat in category_order:
            in_cat = (sub['class_group'] == cat)
            a = int((in_cat & sig).sum())
            b = int((in_cat & ~sig).sum())
            c = int((~in_cat & sig).sum())
            d = int((~in_cat & ~sig).sum())
            _, p = fisher_exact([[a,b],[c,d]], alternative='greater')
            # Haldane–Anscombe correction for stable log2 OR
            or_corr = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
            rows.append({
                'day': day,
                'program': prog,
                'program_label': prog_label,
                'direction': direction,
                'category': cat,
                'category_n': int(in_cat.sum()),
                'significant_in_category': a,
                'significant_total': int(sig.sum()),
                'odds_ratio_corrected': or_corr,
                'log2_odds_ratio': np.log2(or_corr),
                'p_value': p,
                'n_universe': len(sub),
            })

enrich = pd.DataFrame(rows)
enrich['fdr'] = multipletests(enrich['p_value'], method='fdr_bh')[1]
enrich['minus_log10_fdr'] = -np.log10(enrich['fdr'].clip(lower=1e-300))
enrich['log2_or_plot'] = enrich['log2_odds_ratio'].clip(-2.5, 2.5)
enrich['size_plot'] = enrich['minus_log10_fdr'].clip(0, 6)
enrich['significant'] = enrich['fdr'] < 0.10

enrich.to_csv(os.path.join(OUT, 'source_data', 'Fig2_candidate_program_category_enrichment.tsv'), sep='\t', index=False)
universe.to_csv(os.path.join(OUT, 'source_data', 'screen_regulator_class_mapping.tsv'), sep='\t', index=False)

# ---- Plot ----
fig, axes = plt.subplots(1, 2, figsize=(11.4, 6.5), dpi=300, sharey=True)
fig.patch.set_facecolor('white')

prog_labels = [x[2] for x in programs]
xpos = {p:i for i,p in enumerate(prog_labels)}
ypos = {c:i for i,c in enumerate(category_order[::-1])}

cmap = mpl.colors.LinearSegmentedColormap.from_list('or_div', ['#3B78B5','#F7F7F7','#D95F4A'])
norm = mpl.colors.TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5)

for ax, direction in zip(axes, ['Positive effect','Negative effect']):
    ss = enrich[enrich['direction'] == direction]
    for _, r in ss.iterrows():
        size = 22 + 46 * r['size_plot']
        edge = '#1A1A1A' if r['significant'] else '#777777'
        lw = 1.1 if r['significant'] else 0.35
        alpha = 1.0 if r['significant'] else 0.68
        ax.scatter(xpos[r['program_label']], ypos[r['category']], s=size,
                   c=[cmap(norm(r['log2_or_plot']))], edgecolors=edge,
                   linewidths=lw, alpha=alpha, zorder=3)
    ax.set_xticks(range(len(prog_labels)))
    ax.set_xticklabels(prog_labels, fontsize=11, fontweight='bold')
    ax.set_xlim(-0.55, len(prog_labels)-0.45)
    ax.set_ylim(-0.55, len(category_order)-0.45)
    ax.set_title(direction, fontsize=13, fontweight='bold', pad=10)
    ax.grid(axis='both', color='#E7E7E7', linewidth=0.65)
    ax.set_axisbelow(True)
    for s in ['top','right','left','bottom']:
        ax.spines[s].set_visible(False)
    ax.tick_params(axis='x', length=0)
    ax.tick_params(axis='y', length=0)

axes[0].set_yticks(range(len(category_order)))
axes[0].set_yticklabels(category_order[::-1], fontsize=10)
axes[1].tick_params(labelleft=False)

fig.suptitle('Regulator classes are enriched among effects on individual programs',
             fontsize=17, fontweight='bold', y=0.985)
fig.text(0.5, 0.945,
         'Each column is one cNMF program; no program contrasts or composite scores are used',
         ha='center', va='top', fontsize=10.5, color='#555555')

# day brackets under x axis
for ax in axes:
    trans = ax.get_xaxis_transform()
    for start, end, label in [(-0.35,1.35,'D1'),(1.65,3.35,'D2'),(3.65,5.35,'D3')]:
        ax.plot([start,end],[-0.12,-0.12], transform=trans, color='#666666', lw=1.0, clip_on=False)
        ax.text((start+end)/2,-0.18,label,transform=trans,ha='center',va='top',fontsize=10,color='#444444')

# Colorbar
cax = fig.add_axes([0.37, 0.055, 0.26, 0.025])
cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')
cb.set_label('log2 odds ratio (class enrichment)', fontsize=9.5)
cb.ax.tick_params(labelsize=8.5, length=2)

# Size legend
legend_sizes = [0, 2, 4, 6]
handles=[]
labels=[]
for v in legend_sizes:
    handles.append(Line2D([0],[0], marker='o', linestyle='None', markersize=np.sqrt(22+46*v)/1.5,
                          markerfacecolor='white', markeredgecolor='#555555'))
    labels.append(f'{v:g}')
leg = fig.legend(handles, labels, title='−log10(FDR)', frameon=False,
                 loc='lower right', bbox_to_anchor=(0.95,0.035), ncol=4,
                 fontsize=8.5, title_fontsize=9.5, handletextpad=0.2, columnspacing=0.8)
fig.text(0.075,0.06,'Black outline: enrichment FDR < 0.10',fontsize=9.2,color='#444444')

fig.subplots_adjust(left=0.29, right=0.985, top=0.90, bottom=0.18, wspace=0.08)

panel_base = os.path.join(OUT, 'individual_panels', 'Fig2_candidate_panel_program_category_enrichment')
for ext in ['png','pdf','svg']:
    fig.savefig(panel_base + '.' + ext, dpi=350 if ext=='png' else None, bbox_inches='tight')
plt.close(fig)

# Missing-input note for exact fate-composition candidate
missing_note = '''Exact D3 fate-composition panel cannot be regenerated from the files currently present.

Please re-upload one of the following source tables:
1. target_D3_state_counts_and_fractions.csv (preferred), containing target, control_mode, group, state, n_cells, total_D3_cells, and fraction_among_all_D3; or
2. target_D3_endpoint_bias_FibroMesen_vs_artEC.csv, provided it contains target-level artEC and Fibro-Mesen fractions/counts.

The currently retained source table contains artEC fractions only, so Fibro-Mesen fractions cannot be reconstructed without inventing data.
'''
with open(os.path.join(OUT,'MISSING_INPUT_FOR_D3_FATE_COMPOSITION_PANEL.txt'),'w') as f:
    f.write(missing_note)

# README
readme = '''# Fresh Figure 2 candidate-panel rebuild

This directory was recreated in the active workspace.

Completed:
- `Fig2_candidate_panel_program_category_enrichment` in PNG, PDF, and SVG.
- Source tables and reproducible script.

Pending exact regeneration:
- D3 fate-composition map, because the target-level Fibro-Mesen fraction/count table is no longer present. See `MISSING_INPUT_FOR_D3_FATE_COMPOSITION_PANEL.txt`.

The enrichment panel tests each selected cNMF program independently (P30, P13, P15, P17, P3, P34). It uses no contrast or composite program score.
'''
with open(os.path.join(OUT,'README.md'),'w') as f:
    f.write(readme)

# copy script
shutil.copyfile(__file__, os.path.join(OUT,'scripts','build_candidate_figure2_panels_fresh.py'))

zip_path = os.path.join(BASE,'figure2_candidate_new_panels_fresh_package.zip')
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for root,_,files in os.walk(OUT):
        for fn in files:
            fp=os.path.join(root,fn)
            z.write(fp,arcname=os.path.relpath(fp,OUT))

print('Created:', panel_base+'.png', panel_base+'.pdf', panel_base+'.svg', zip_path, sep='\n')
