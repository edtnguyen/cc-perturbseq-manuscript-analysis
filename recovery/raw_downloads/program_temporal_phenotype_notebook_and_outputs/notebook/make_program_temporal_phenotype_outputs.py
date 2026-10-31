import os, zipfile, textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle

try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except Exception:
    HAS_ADJUST_TEXT = False

BASE = '/mnt/data'
ADV_OUT = os.path.join(BASE, 'outputs_program_first_advanced_unzipped')
OUT = os.path.join(BASE, 'program_temporal_phenotype_outputs')
os.makedirs(OUT, exist_ok=True)

EFFECT_CSV = os.path.join(ADV_OUT, 'target_by_curated_program_effect_matrix.csv')
FDR_CSV = os.path.join(ADV_OUT, 'target_by_curated_program_fdr_matrix.csv')
FIRST_DAY_CSV = os.path.join(ADV_OUT, 'target_first_affected_day_and_persistence.csv')
ANNOT_CSV = os.path.join('/mnt/data/outputs_program_level_temporal_map_v4_unzipped', 'curated_program_annotations_used.csv')

# -----------------------------
# Load data
# -----------------------------
eff = pd.read_csv(EFFECT_CSV)
fdr = pd.read_csv(FDR_CSV)
first = pd.read_csv(FIRST_DAY_CSV)
annot = pd.read_csv(ANNOT_CSV) if os.path.exists(ANNOT_CSV) else None

programs = [c for c in eff.columns if c != 'target']
eff = eff.set_index('target')
fdr = fdr.set_index('target')

# Program groups
D0_PLURI = ['P7', 'P22', 'P29']
D1_PROGRAMS = ['P28', 'P13', 'P30', 'P48', 'P33']
D2_PROGRAMS = ['P15', 'P14', 'P17']
D3_EC = ['P2', 'P3', 'P8', 'P11']
D3_ALT = ['P34', 'P27']

# Ensure present
for group in [D0_PLURI, D1_PROGRAMS, D2_PROGRAMS, D3_EC, D3_ALT]:
    for p in group:
        if p not in eff.columns:
            raise ValueError(f'Missing {p} from effect matrix')

# -----------------------------
# Phenotype feature table
# -----------------------------
features = pd.DataFrame(index=eff.index)
features['pluripotency_mean_effect'] = eff[D0_PLURI].mean(axis=1)
features['pluripotency_max_effect'] = eff[D0_PLURI].max(axis=1)
features['D1_primitive_streak_absmax_effect'] = eff[D1_PROGRAMS].abs().max(axis=1)
features['P48_neuroectoderm_effect'] = eff['P48']
features['P15_vascular_progenitor_effect'] = eff['P15']
features['P14_standard_LM_effect'] = eff['P14']
features['P17_YAP_LM_effect'] = eff['P17']
features['P34_YAP_mesenchyme_effect'] = eff['P34']
features['P27_ECM_mesenchyme_effect'] = eff['P27']
features['D3_EC_mean_effect'] = eff[D3_EC].mean(axis=1)
features['D3_EC_min_effect'] = eff[D3_EC].min(axis=1)
features['D3_EC_max_effect'] = eff[D3_EC].max(axis=1)
features['D3_ALT_mean_effect'] = eff[D3_ALT].mean(axis=1)
features['D3_ALT_max_effect'] = eff[D3_ALT].max(axis=1)

# significance helper
sig = fdr < 0.05
features['n_D0_pluri_up'] = ((eff[D0_PLURI] > 0) & sig[D0_PLURI]).sum(axis=1)
features['n_D1_programs_affected'] = sig[D1_PROGRAMS].sum(axis=1)
features['n_D2_programs_affected'] = sig[D2_PROGRAMS].sum(axis=1)
features['n_EC_programs_down'] = ((eff[D3_EC] < 0) & sig[D3_EC]).sum(axis=1)
features['n_EC_programs_up'] = ((eff[D3_EC] > 0) & sig[D3_EC]).sum(axis=1)
features['n_alt_programs_up'] = ((eff[D3_ALT] > 0) & sig[D3_ALT]).sum(axis=1)
features['n_alt_programs_down'] = ((eff[D3_ALT] < 0) & sig[D3_ALT]).sum(axis=1)
features['n_total_curated_sig_programs'] = sig[programs].sum(axis=1)
features['effect_spread_sd'] = eff[programs].std(axis=1)
features['strongest_abs_effect'] = eff[programs].abs().max(axis=1)

# Raw pairwise pattern evidence (not subtracting)
features['P15_down_sig'] = (eff['P15'] < 0) & sig['P15']
features['P15_up_sig'] = (eff['P15'] > 0) & sig['P15']
features['P17_up_sig'] = (eff['P17'] > 0) & sig['P17']
features['P17_down_sig'] = (eff['P17'] < 0) & sig['P17']
features['P14_up_sig'] = (eff['P14'] > 0) & sig['P14']
features['P14_down_sig'] = (eff['P14'] < 0) & sig['P14']
features['P34_up_sig'] = (eff['P34'] > 0) & sig['P34']
features['P27_up_sig'] = (eff['P27'] > 0) & sig['P27']
features['P48_up_sig'] = (eff['P48'] > 0) & sig['P48']
features['EC_endpoint_loss_sig'] = features['n_EC_programs_down'] >= 2
features['EC_endpoint_gain_sig'] = features['n_EC_programs_up'] >= 2
features['ALT_endpoint_gain_sig'] = features['n_alt_programs_up'] >= 1

# Composite tags: these are binary patterns requiring raw program evidence, not calibrated subtraction scores.
tags = pd.DataFrame(index=eff.index)
tags['pluripotency retention / delayed exit'] = features['n_D0_pluri_up'] >= 1
tags['D1 primitive-streak patterning'] = features['n_D1_programs_affected'] >= 2
tags['early neuroectoderm-like diversion (P48 up)'] = features['P48_up_sig']
tags['loss of P15 vascular progenitor'] = features['P15_down_sig']
tags['gain of P15 vascular progenitor'] = features['P15_up_sig']
tags['gain of standard LM P14'] = features['P14_up_sig']
tags['gain of YAP/TAZ LM P17'] = features['P17_up_sig']
tags['D2 split/diversion: P17 up + P15 down'] = features['P17_up_sig'] & features['P15_down_sig']
tags['D3 EC endpoint loss'] = features['EC_endpoint_loss_sig']
tags['D3 EC endpoint gain'] = features['EC_endpoint_gain_sig']
tags['gain of P34 YAP/TAZ mesenchyme'] = features['P34_up_sig']
tags['gain of P27 ECM mesenchyme'] = features['P27_up_sig']
tags['D3 diversion: P34 up + EC down'] = features['P34_up_sig'] & features['EC_endpoint_loss_sig']
tags['D3 ECM diversion: P27 up + EC down'] = features['P27_up_sig'] & features['EC_endpoint_loss_sig']
tags['broad / pleiotropic curated effect'] = features['n_total_curated_sig_programs'] >= 6
tags['persistent effect across >=2 days'] = first.set_index('target').reindex(eff.index)['n_days_affected'].fillna(0).astype(int) >= 2

# Dominant phenotype tag: prioritize interpretable composite tags
priority = [
    'D3 diversion: P34 up + EC down',
    'D3 ECM diversion: P27 up + EC down',
    'D2 split/diversion: P17 up + P15 down',
    'D3 EC endpoint loss',
    'gain of P34 YAP/TAZ mesenchyme',
    'gain of P27 ECM mesenchyme',
    'loss of P15 vascular progenitor',
    'gain of YAP/TAZ LM P17',
    'gain of standard LM P14',
    'early neuroectoderm-like diversion (P48 up)',
    'D1 primitive-streak patterning',
    'pluripotency retention / delayed exit',
    'broad / pleiotropic curated effect',
]

def assign_primary(row):
    for t in priority:
        if bool(row.get(t, False)):
            return t
    return 'no curated phenotype tag'

tags['primary_program_phenotype_tag'] = tags.apply(assign_primary, axis=1)

# Merge first affected day
first_idx = first.set_index('target').reindex(eff.index)
features['first_affected_day'] = first_idx['first_affected_day'].fillna('none')
features['n_days_affected'] = first_idx['n_days_affected'].fillna(0).astype(int)
features['primary_program_phenotype_tag'] = tags['primary_program_phenotype_tag']

# Strongest program labels
features['strongest_positive_program'] = eff[programs].idxmax(axis=1)
features['strongest_positive_effect'] = eff[programs].max(axis=1)
features['strongest_negative_program'] = eff[programs].idxmin(axis=1)
features['strongest_negative_effect'] = eff[programs].min(axis=1)

features.reset_index().rename(columns={'target':'target'}).to_csv(os.path.join(OUT, 'program_defined_temporal_phenotype_features.csv'), index=False)
tags.reset_index().rename(columns={'target':'target'}).to_csv(os.path.join(OUT, 'program_defined_temporal_phenotype_tags.csv'), index=False)

# -----------------------------
# Plot 1: tag counts
# -----------------------------
count_cols = [c for c in tags.columns if c != 'primary_program_phenotype_tag']
tag_counts = tags[count_cols].sum().sort_values(ascending=True)
fig, ax = plt.subplots(figsize=(7.2, 5.0), dpi=300)
y = np.arange(len(tag_counts))
ax.barh(y, tag_counts.values, color='#4E79A7')
ax.set_yticks(y)
ax.set_yticklabels(tag_counts.index, fontsize=8.5)
ax.set_xlabel('Number of target genes', fontsize=10.5, fontweight='bold')
ax.set_title('Program-defined temporal phenotype tags', fontsize=12.5, fontweight='bold')
for spine in ['top','right']:
    ax.spines[spine].set_visible(False)
ax.tick_params(axis='x', labelsize=9)
fig.tight_layout()
for ext in ['png','pdf','svg']:
    fig.savefig(os.path.join(OUT, f'program_defined_phenotype_tag_counts.{ext}'), bbox_inches='tight', dpi=300 if ext=='png' else None)
plt.close(fig)

# -----------------------------
# Plot 2: selected target raw program fingerprint features
# -----------------------------
selected = [
    'NF2','LATS2','SAV1','TAOK1','AMOTL2',
    'SMAD1','SMAD4','SMAD5','LEMD3','BMPR1A','ACVR2A','ACVR1B','ACVR2B',
    'PRKCD','KDR','ETV2','ERG','FLI1','TAL1','SOX17',
]
selected = [t for t in selected if t in features.index]
feature_cols = [
    'pluripotency_mean_effect','P48_neuroectoderm_effect','P15_vascular_progenitor_effect',
    'P14_standard_LM_effect','P17_YAP_LM_effect','P34_YAP_mesenchyme_effect',
    'P27_ECM_mesenchyme_effect','D3_EC_mean_effect'
]
pretty_cols = ['D0 pluri mean','P48 neuroectoderm','P15 vascular progenitor','P14 standard LM','P17 YAP/TAZ LM','P34 YAP mesenchyme','P27 ECM mesenchyme','D3 EC mean']
mat = features.loc[selected, feature_cols].copy()
# cap for plotting
plot_mat = mat.clip(-4,4)
fig, ax = plt.subplots(figsize=(9.2, max(3.6, 0.34*len(selected)+1.2)), dpi=300)
im = ax.imshow(plot_mat.values, aspect='auto', cmap='RdBu_r', vmin=-4, vmax=4)
ax.set_xticks(np.arange(len(feature_cols)))
ax.set_xticklabels(pretty_cols, rotation=45, ha='right', fontsize=8.5)
ax.set_yticks(np.arange(len(selected)))
ax.set_yticklabels(selected, fontsize=8.5)
ax.set_title('Raw program-effect fingerprints for selected regulators', fontsize=12.5, fontweight='bold')
for i in range(plot_mat.shape[0]):
    for j in range(plot_mat.shape[1]):
        val = mat.iloc[i,j]
        txt_color = 'white' if abs(plot_mat.iloc[i,j]) > 2.2 else 'black'
        ax.text(j, i, f'{val:.1f}', ha='center', va='center', fontsize=6.5, color=txt_color)
fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label='program effect')
fig.tight_layout()
for ext in ['png','pdf','svg']:
    fig.savefig(os.path.join(OUT, f'selected_regulator_program_temporal_fingerprints.{ext}'), bbox_inches='tight', dpi=300 if ext=='png' else None)
plt.close(fig)

# -----------------------------
# Plot 3: phenotype tag binary matrix for selected + top phenotype targets
# -----------------------------
extra_targets = []
for col in ['D2 split/diversion: P17 up + P15 down','D3 diversion: P34 up + EC down','early neuroectoderm-like diversion (P48 up)','loss of P15 vascular progenitor']:
    tmp = features.loc[tags[col].values].copy()
    if col == 'D2 split/diversion: P17 up + P15 down':
        rank = (tmp['P17_YAP_LM_effect'] - tmp['P15_vascular_progenitor_effect']).sort_values(ascending=False).head(8).index.tolist()
    elif col == 'D3 diversion: P34 up + EC down':
        rank = (tmp['P34_YAP_mesenchyme_effect'] - tmp['D3_EC_mean_effect']).sort_values(ascending=False).head(8).index.tolist()
    elif col == 'early neuroectoderm-like diversion (P48 up)':
        rank = tmp['P48_neuroectoderm_effect'].sort_values(ascending=False).head(6).index.tolist()
    else:
        rank = (-tmp['P15_vascular_progenitor_effect']).sort_values(ascending=False).head(6).index.tolist()
    extra_targets.extend(rank)
selected2 = list(dict.fromkeys(selected + extra_targets))[:55]
bin_cols = [
    'pluripotency retention / delayed exit',
    'D1 primitive-streak patterning',
    'early neuroectoderm-like diversion (P48 up)',
    'loss of P15 vascular progenitor',
    'gain of standard LM P14',
    'gain of YAP/TAZ LM P17',
    'D2 split/diversion: P17 up + P15 down',
    'D3 EC endpoint loss',
    'gain of P34 YAP/TAZ mesenchyme',
    'gain of P27 ECM mesenchyme',
    'D3 diversion: P34 up + EC down',
    'D3 ECM diversion: P27 up + EC down',
    'persistent effect across >=2 days',
]
binmat = tags.loc[selected2, bin_cols].astype(int)
fig, ax = plt.subplots(figsize=(9.5, max(4.0, 0.20*len(selected2)+1.6)), dpi=300)
im = ax.imshow(binmat.values, aspect='auto', cmap=ListedColormap(['#F2F2F2','#2F5597']), vmin=0, vmax=1)
ax.set_yticks(np.arange(len(selected2)))
ax.set_yticklabels(selected2, fontsize=6.6)
ax.set_xticks(np.arange(len(bin_cols)))
ax.set_xticklabels(bin_cols, rotation=45, ha='right', fontsize=7.0)
ax.set_title('Program-defined phenotype tags for selected targets', fontsize=12.5, fontweight='bold')
# grid lines
for i in range(len(selected2)+1):
    ax.axhline(i-0.5, color='white', lw=0.4)
for j in range(len(bin_cols)+1):
    ax.axvline(j-0.5, color='white', lw=0.4)
fig.tight_layout()
for ext in ['png','pdf','svg']:
    fig.savefig(os.path.join(OUT, f'program_defined_phenotype_tag_matrix_selected_targets.{ext}'), bbox_inches='tight', dpi=300 if ext=='png' else None)
plt.close(fig)

# -----------------------------
# Helper for quadrant scatter plots
# -----------------------------
def quadrant_scatter(xcol, ycol, xlabel, ylabel, title, fname, highlight_condition=None, label_targets=None):
    df = features.copy()
    x = df[xcol]
    y = df[ycol]
    fig, ax = plt.subplots(figsize=(5.8,5.2), dpi=300)
    if highlight_condition is None:
        highlight = pd.Series(False, index=df.index)
    else:
        highlight = highlight_condition(df)
    ax.scatter(x[~highlight], y[~highlight], s=18, color='#BDBDBD', alpha=0.45, edgecolors='none')
    ax.scatter(x[highlight], y[highlight], s=34, color='#D95F02', alpha=0.9, edgecolors='white', linewidth=0.4)
    ax.axhline(0, color='#777777', lw=0.85, ls='--')
    ax.axvline(0, color='#777777', lw=0.85, ls='--')
    ax.set_xlabel(xlabel, fontsize=10.5, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=10.5, fontweight='bold')
    ax.set_title(title, fontsize=12.0, fontweight='bold')
    for spine in ['top','right']:
        ax.spines[spine].set_visible(False)
    texts=[]
    if label_targets is None:
        label_targets = []
    for t in label_targets:
        if t in df.index:
            texts.append(ax.text(df.loc[t, xcol], df.loc[t, ycol], t, fontsize=8.0, color='black'))
    if HAS_ADJUST_TEXT and texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#777777', lw=0.45, alpha=0.7), time_lim=1.0)
    fig.tight_layout()
    for ext in ['png','pdf','svg']:
        fig.savefig(os.path.join(OUT, f'{fname}.{ext}'), bbox_inches='tight', dpi=300 if ext=='png' else None)
    plt.close(fig)

label_targets = ['NF2','LATS2','SAV1','TAOK1','AMOTL2','PRKCD','SMAD1','SMAD4','SMAD5','BMPR1A','ACVR2A','KDR','ETV2','ERG','FLI1','TAL1']
quadrant_scatter(
    'P15_vascular_progenitor_effect','P17_YAP_LM_effect',
    'raw effect on P15 vascular progenitor', 'raw effect on P17 YAP/TAZ LM',
    'D2 raw program quadrant: P15 vs P17',
    'D2_raw_P15_vs_P17_quadrant',
    highlight_condition=lambda df: (df['P17_up_sig'] & df['P15_down_sig']),
    label_targets=label_targets
)
quadrant_scatter(
    'D3_EC_mean_effect','P34_YAP_mesenchyme_effect',
    'mean raw effect on D3 EC programs', 'raw effect on P34 YAP/TAZ mesenchyme',
    'D3 raw endpoint quadrant: EC mean vs P34',
    'D3_raw_ECmean_vs_P34_quadrant',
    highlight_condition=lambda df: (df['P34_up_sig'] & df['EC_endpoint_loss_sig']),
    label_targets=label_targets
)
quadrant_scatter(
    'P17_YAP_LM_effect','P34_YAP_mesenchyme_effect',
    'raw effect on P17 YAP/TAZ LM', 'raw effect on P34 YAP/TAZ mesenchyme',
    'Raw program continuity: D2 P17 vs D3 P34',
    'raw_P17_vs_P34_continuity_scatter',
    highlight_condition=lambda df: (df['P17_up_sig'] & df['P34_up_sig']),
    label_targets=label_targets
)

# -----------------------------
# Plot 7: first affected day x primary phenotype tag counts heatmap
# -----------------------------
ct = features.reset_index().groupby(['first_affected_day','primary_program_phenotype_tag']).size().unstack(fill_value=0)
# order rows and columns
row_order = [x for x in ['D0','D1','D2','D3','none'] if x in ct.index] + [x for x in ct.index if x not in ['D0','D1','D2','D3','none']]
col_order = ct.sum(axis=0).sort_values(ascending=False).index.tolist()
ct = ct.loc[row_order, col_order]
fig, ax = plt.subplots(figsize=(max(8,0.45*len(col_order)+2), 3.4), dpi=300)
im = ax.imshow(ct.values, aspect='auto', cmap='Blues')
ax.set_yticks(np.arange(ct.shape[0]))
ax.set_yticklabels(ct.index, fontsize=9)
ax.set_xticks(np.arange(ct.shape[1]))
ax.set_xticklabels(ct.columns, rotation=45, ha='right', fontsize=7.6)
for i in range(ct.shape[0]):
    for j in range(ct.shape[1]):
        val = ct.iloc[i,j]
        if val>0:
            ax.text(j,i,str(val),ha='center',va='center',fontsize=7,color='black' if val<10 else 'white')
ax.set_title('First affected day by primary program phenotype tag', fontsize=12.5, fontweight='bold')
fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label='targets')
fig.tight_layout()
for ext in ['png','pdf','svg']:
    fig.savefig(os.path.join(OUT, f'first_affected_day_by_primary_phenotype_tag.{ext}'), bbox_inches='tight', dpi=300 if ext=='png' else None)
plt.close(fig)

# Package outputs
zip_path = '/mnt/data/program_temporal_phenotype_outputs_package.zip'
with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
    for fn in os.listdir(OUT):
        z.write(os.path.join(OUT, fn), arcname=fn)
print('Wrote outputs to', OUT)
print('ZIP:', zip_path)
