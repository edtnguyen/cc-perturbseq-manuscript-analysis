#!/usr/bin/env python3
"""
Start-from-scratch target-day module clustering for program-first Perturb-seq analysis.

Main idea:
1. Read curated target x program effects from outputs_program_first_advanced.zip.
2. Build target-day matrices from target x day x curated program effects.
3. Cluster target-day rows using significant program effects only (FDR < 0.05) for pattern modules.
4. Select/evaluate k, but generate k=12 modules because this was the data-driven sparse target-day recommendation and the M2/M3/M4 interpretation used this cut.
5. Annotate Hippo/pathway categories only after clustering.

Clustering unit:
    row = target_day perturbation profile, e.g. NF2_D3
    columns = expert-curated cNMF programs
    value = program effect if FDR < 0.05, else 0 for sparse pattern clustering

Distance/linkage:
    z-score columns -> cosine distance -> average-linkage hierarchical clustering.
"""

from __future__ import annotations
import os, re, zipfile, shutil, warnings, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, adjusted_rand_score, davies_bouldin_score
warnings.filterwarnings('ignore')

mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.family'] = 'DejaVu Sans'

BASE = Path('/mnt/data')
OUT = BASE / 'target_day_module_clustering_from_scratch'
OUT.mkdir(parents=True, exist_ok=True)
SRC = OUT / '_source'
SRC.mkdir(parents=True, exist_ok=True)

EFFECT_ZIP = BASE / 'outputs_program_first_advanced.zip'
ANNOTATION_CSV = BASE / 'program_first_advanced_analyses_notebook' / 'curated_program_annotations.csv'
K_MAIN = 12
FDR_CUTOFF = 0.05
RANDOM_STATE = 123

# ------------------------------
# Load source tables
# ------------------------------
needed = [
    'curated_program_effects_aggregated.csv',
    'target_by_curated_program_effect_matrix.csv',
    'target_by_curated_program_fdr_matrix.csv',
]
if EFFECT_ZIP.exists():
    with zipfile.ZipFile(EFFECT_ZIP) as z:
        for name in needed:
            if not (SRC / name).exists():
                z.extract(name, SRC)
else:
    raise FileNotFoundError(f"Missing {EFFECT_ZIP}. Please reupload outputs_program_first_advanced.zip or provide these CSVs: {needed}")

long = pd.read_csv(SRC / 'curated_program_effects_aggregated.csv')
target_mat = pd.read_csv(SRC / 'target_by_curated_program_effect_matrix.csv', index_col=0)
target_fdr = pd.read_csv(SRC / 'target_by_curated_program_fdr_matrix.csv', index_col=0)
ann = pd.read_csv(ANNOTATION_CSV)

# Normalize column names/fields
ann['program_label'] = ann['program_label'].astype(str)
ann['display_label'] = ann['program_label'] + '\n' + ann['expert_annotation'].astype(str)
program_order = [p for p in ann.sort_values('program_order')['program_label'].tolist() if p in target_mat.columns]
program_day = dict(zip(ann['program_label'], ann['day']))
program_group = dict(zip(ann['program_label'], ann['program_group']))
program_display = dict(zip(ann['program_label'], ann['display_label']))
DAYS = ['D0','D1','D2','D3']

# Save source metadata
ann.to_csv(OUT / 'curated_program_annotations_used.csv', index=False)

# ------------------------------
# Build target-day matrices
# ------------------------------
def build_target_day_matrix(long_df: pd.DataFrame, sig_only=True, fdr_cutoff=0.05):
    df = long_df.copy()
    if sig_only:
        df['value'] = np.where(df['fdr'] < fdr_cutoff, df['effect'], 0.0)
        value_label = f'sigFDR{fdr_cutoff}'
    else:
        df['value'] = df['effect']
        value_label = 'raw'
    rows = []
    for (target, day), sub in df.groupby(['target','day']):
        row = {'item': f'{target}_{day}', 'target': target, 'day': day}
        for p in program_order:
            row[p] = 0.0
        for _, r in sub.iterrows():
            p = str(r['program_label'])
            if p in program_order:
                row[p] = float(r['value'])
        rows.append(row)
    meta = pd.DataFrame(rows).set_index('item')
    X = meta[program_order].copy()
    # keep rows with at least one non-zero feature
    X = X.loc[X.abs().sum(axis=1) > 0]
    meta = meta.loc[X.index, ['target','day']]
    return X, meta, value_label

X_sig, meta_sig, _ = build_target_day_matrix(long, sig_only=True, fdr_cutoff=FDR_CUTOFF)
X_raw, meta_raw, _ = build_target_day_matrix(long, sig_only=False)
X_sig.to_csv(OUT / 'target_day_matrix_sigFDR0.05.csv')
meta_sig.to_csv(OUT / 'target_day_metadata_sigFDR0.05.csv')
X_raw.to_csv(OUT / 'target_day_matrix_raw.csv')
meta_raw.to_csv(OUT / 'target_day_metadata_raw.csv')

# ------------------------------
# Clustering utilities
# ------------------------------
def standardize_X(Xdf: pd.DataFrame):
    scaler = StandardScaler(with_mean=True, with_std=True)
    Xz = scaler.fit_transform(Xdf.values.astype(float))
    Xz = np.nan_to_num(Xz, nan=0.0, posinf=0.0, neginf=0.0)
    return Xz

def make_linkage(Xz):
    d = pdist(Xz, metric='cosine')
    if np.allclose(d, 0):
        rng = np.random.default_rng(RANDOM_STATE)
        d = pdist(Xz + rng.normal(0, 1e-8, Xz.shape), metric='cosine')
    return linkage(d, method='average')

def relabel_modules_by_dendrogram(labels, Z, index):
    order = dendrogram(Z, no_plot=True)['leaves']
    ordered_items = list(index[order])
    ser = pd.Series(labels, index=index)
    seen = []
    for item in ordered_items:
        m = int(ser.loc[item])
        if m not in seen:
            seen.append(m)
    remap = {m:i+1 for i,m in enumerate(seen)}
    return ser.map(remap).astype(int), ordered_items

def bootstrap_stability(Xz, Z, ks, n_boot=60):
    rng = np.random.default_rng(RANDOM_STATE)
    base = {k: fcluster(Z, t=k, criterion='maxclust') for k in ks}
    out = {k: [] for k in ks}
    n, p = Xz.shape
    for b in range(n_boot):
        cols = rng.choice(np.arange(p), size=p, replace=True)
        Xb = Xz[:, cols] + rng.normal(0, 0.015, size=(n, p))
        Zb = make_linkage(Xb)
        for k in ks:
            out[k].append(adjusted_rand_score(base[k], fcluster(Zb, t=k, criterion='maxclust')))
    return {k: float(np.mean(out[k])) for k in ks}

def evaluate_k(Xdf, label='target_day_sigFDR0.05', k_min=2, k_max=20):
    Xz = standardize_X(Xdf)
    Z = make_linkage(Xz)
    max_k = min(k_max, Xdf.shape[0]-1)
    ks = list(range(k_min, max_k + 1))
    stability = bootstrap_stability(Xz, Z, ks)
    D = squareform(pdist(Xz, metric='cosine'))
    rows = []
    for k in ks:
        lab = fcluster(Z, t=k, criterion='maxclust')
        sizes = pd.Series(lab).value_counts()
        try:
            sil = silhouette_score(D, lab, metric='precomputed')
        except Exception:
            sil = np.nan
        try:
            db = davies_bouldin_score(Xz, lab)
        except Exception:
            db = np.nan
        rows.append({
            'analysis': label,
            'k': k,
            'silhouette_cosine': sil,
            'bootstrap_ARI_mean': stability[k],
            'davies_bouldin': db,
            'min_module_size': int(sizes.min()),
            'n_small_modules_lt_5': int((sizes < 5).sum()),
            'frac_items_in_small_modules': float(sizes[sizes < 5].sum() / sizes.sum()) if (sizes < 5).any() else 0.0,
        })
    m = pd.DataFrame(rows)
    def scale(s, higher=True):
        v = s.astype(float).values
        if np.all(np.isnan(v)) or np.nanmax(v) == np.nanmin(v):
            x = np.ones_like(v) * 0.5
        else:
            x = (v - np.nanmin(v)) / (np.nanmax(v) - np.nanmin(v))
        return x if higher else 1-x
    m['score_silhouette'] = scale(m['silhouette_cosine'], True)
    m['score_stability'] = scale(m['bootstrap_ARI_mean'], True)
    m['score_db'] = scale(m['davies_bouldin'], False)
    m['composite_score'] = 0.45*m['score_silhouette'] + 0.40*m['score_stability'] + 0.15*m['score_db'] - 0.35*m['frac_items_in_small_modules']
    max_score = m['composite_score'].max()
    recommended_k = int(m.loc[m['composite_score'] >= 0.95*max_score, 'k'].min())
    return Xz, Z, m, recommended_k

# ------------------------------
# Evaluate k and cluster k=12
# ------------------------------
Xz_sig, Z_sig, k_metrics, recommended_k = evaluate_k(X_sig, label='target_day_sigFDR0.05')
k_metrics.to_csv(OUT / 'k_selection_metrics_target_day_sigFDR0.05.csv', index=False)

labels_raw = fcluster(Z_sig, t=K_MAIN, criterion='maxclust')
module_ser, ordered_items = relabel_modules_by_dendrogram(labels_raw, Z_sig, X_sig.index)
assign = meta_sig.copy()
assign['item'] = assign.index
assign['module'] = module_ser.loc[assign.index].values
assign['module_label'] = assign['module'].apply(lambda x: f'M{x}')
assign['is_hippo_robust'] = assign['target'].isin(['NF2','LATS2','TAOK1','AMOTL2'])
assign['is_hippo_extended'] = assign['target'].isin(['NF2','LATS2','TAOK1','AMOTL2','SAV1'])
assign.to_csv(OUT / f'target_day_module_assignments_k{K_MAIN}.csv', index=False)

# module centroids
centroids = []
for m in sorted(assign['module'].unique()):
    members = assign.loc[assign['module'] == m, 'item']
    row = X_sig.loc[members, program_order].mean(axis=0)
    row.name = f'M{m}\n(n={len(members)})'
    centroids.append(row)
centroid_df = pd.DataFrame(centroids)
centroid_df.to_csv(OUT / f'target_day_module_centroids_k{K_MAIN}.csv')

# temporal specificity
temporal = pd.crosstab(assign['module_label'], assign['day']).reindex(columns=DAYS, fill_value=0)
temporal_frac = temporal.div(temporal.sum(axis=1), axis=0)
temporal.to_csv(OUT / f'target_day_module_temporal_counts_k{K_MAIN}.csv')
temporal_frac.to_csv(OUT / f'target_day_module_temporal_fraction_k{K_MAIN}.csv')

# ------------------------------
# Post-hoc categories
# ------------------------------
CATEGORY_SETS = {
    'Hippo robust': ['NF2','LATS2','TAOK1','AMOTL2'],
    'SAV1': ['SAV1'],
    'BMP/TGF-SMAD': ['SMAD1','SMAD4','SMAD5','LEMD3','BMPR1A','ACVR2A','ACVR1B','ACVR2B'],
    'EC regulators': ['KDR','ETV2','ERG','FLI1','TAL1','SOX17','CDH5','PECAM1'],
    'Apoptosis': ['BAD','APAF1','BAX','CASP3','CYCS'],
    'Chromatin/cohesin/mediator': ['MED12','MED15','SETDB1','KDM1A','RAD21','STAG2','MAU2','NEDD8','CUL3','RNF20','TAF6L','TADA2B','SGF29','SMARCC1','SMARCE1','CCNC','CTR9'],
    'Signaling/candidate': ['PRKCD','PLCG1','PRKD1','PTEN','PPP2R1A','VHL','KDR'],
}
def target_category(t):
    cats = []
    for name, genes in CATEGORY_SETS.items():
        if t in genes:
            cats.append(name)
    if not cats:
        return 'Other'
    return '; '.join(cats)
assign['posthoc_category'] = assign['target'].map(target_category)
assign.to_csv(OUT / f'target_day_module_assignments_k{K_MAIN}_with_posthoc_categories.csv', index=False)
cat_counts = assign.groupby(['module_label','posthoc_category']).size().reset_index(name='n')
cat_counts.to_csv(OUT / f'posthoc_category_counts_by_module_k{K_MAIN}.csv', index=False)

# Summaries for modules
summary_rows=[]
for m in sorted(assign['module'].unique()):
    members = assign[assign['module']==m]
    centroid = centroid_df.loc[f'M{m}\n(n={len(members)})'] if f'M{m}\n(n={len(members)})' in centroid_df.index else X_sig.loc[members['item'], program_order].mean(axis=0)
    top_pos = centroid.sort_values(ascending=False).head(4)
    top_neg = centroid.sort_values(ascending=True).head(4)
    hippo = sorted(set(members.loc[members['is_hippo_extended'], 'item']))
    cats = members['posthoc_category'].value_counts().head(5).to_dict()
    summary_rows.append({
        'module': f'M{m}',
        'n_rows': len(members),
        'day_counts': json.dumps(members['day'].value_counts().to_dict()),
        'top_positive_programs': '; '.join([f'{p}:{v:.2f}' for p,v in top_pos.items()]),
        'top_negative_programs': '; '.join([f'{p}:{v:.2f}' for p,v in top_neg.items()]),
        'hippo_items': '; '.join(hippo),
        'top_posthoc_categories': json.dumps(cats),
    })
summary = pd.DataFrame(summary_rows)
summary.to_csv(OUT / f'module_summary_k{K_MAIN}.csv', index=False)

# ------------------------------
# Plot utilities
# ------------------------------
POS = '#E68613'
NEG = '#2F76B7'
PURPLE = '#7A4BA3'
GREY = '#BDBDBD'
RED = '#C43C39'

def savefig(fig, stem):
    for ext in ['png','pdf','svg']:
        fig.savefig(OUT/f'{stem}.{ext}', dpi=350 if ext == 'png' else None, bbox_inches='tight')
    plt.close(fig)

def program_xticklabels():
    labels=[]
    for p in program_order:
        nm = ann.set_index('program_label').loc[p, 'expert_annotation']
        labels.append(f'{p}\n{nm}')
    return labels

# ------------------------------
# Plot 1: k selection metrics
# ------------------------------
fig, axes = plt.subplots(2, 2, figsize=(9.4, 6.2), sharex=True)
axes = axes.flatten()
metrics_to_plot = [
    ('silhouette_cosine','Silhouette (cosine)','higher better'),
    ('bootstrap_ARI_mean','Bootstrap stability','higher better'),
    ('frac_items_in_small_modules','Fraction in small modules','lower better'),
    ('composite_score','Composite score','higher better'),
]
for ax,(col,title,sub) in zip(axes, metrics_to_plot):
    ax.plot(k_metrics['k'], k_metrics[col], marker='o', lw=1.6, color='#2F76B7')
    ax.axvline(K_MAIN, color=RED, ls='--', lw=1.2, label=f'k={K_MAIN}')
    ax.axvline(recommended_k, color='black', ls=':', lw=1.2, label=f'recommended={recommended_k}')
    ax.set_title(title, fontsize=10.5, weight='bold')
    ax.text(0.03,0.94,sub,transform=ax.transAxes,ha='left',va='top',fontsize=8,color='#555555')
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(labelsize=8.5)
axes[0].legend(frameon=False, fontsize=7.8, loc='best')
for ax in axes[2:]: ax.set_xlabel('number of modules k', fontsize=9)
fig.suptitle('Selecting module number for target-day sparse program-effect profiles', fontsize=12.5, weight='bold')
fig.tight_layout()
savefig(fig, '01_k_selection_metrics_target_day_sigFDR0.05')

# ------------------------------
# Plot 2: dendrogram with k=12
# ------------------------------
fig, ax = plt.subplots(figsize=(12, 4.8))
dendrogram(Z_sig, no_labels=True, color_threshold=None, above_threshold_color='#555555', ax=ax)
ax.set_title(f'Hierarchical clustering of target-day profiles (cosine + average linkage, k={K_MAIN})', fontsize=12, weight='bold')
ax.set_ylabel('cosine distance')
ax.set_xticks([])
ax.spines[['top','right','bottom']].set_visible(False)
savefig(fig, f'02_target_day_dendrogram_k{K_MAIN}')

# ------------------------------
# Plot 3: module centroid heatmap + temporal specificity side by side
# ------------------------------
fig = plt.figure(figsize=(12.6, 5.0))
gs = fig.add_gridspec(1, 2, width_ratios=[3.7, 1.4], wspace=0.35)
ax1 = fig.add_subplot(gs[0,0])
vals = centroid_df.values
lim = max(0.5, np.nanmax(np.abs(vals)))
im = ax1.imshow(vals, aspect='auto', cmap='RdBu_r', vmin=-lim, vmax=lim)
ax1.set_xticks(range(len(program_order)))
ax1.set_xticklabels(program_xticklabels(), rotation=45, ha='right', fontsize=7.5)
ax1.set_yticks(range(len(centroid_df)))
ax1.set_yticklabels(centroid_df.index, fontsize=8.5)
ax1.set_title('Module centroids: mean significant program effect', fontsize=11.5, weight='bold')
cb = fig.colorbar(im, ax=ax1, fraction=0.035, pad=0.02)
cb.set_label('mean effect (FDR < 0.05 only)', fontsize=8)
cb.ax.tick_params(labelsize=7)

ax2 = fig.add_subplot(gs[0,1])
tf = temporal_frac.reindex(index=[f'M{i}' for i in sorted(assign['module'].unique())])
im2 = ax2.imshow(tf.values, aspect='auto', cmap='Blues', vmin=0, vmax=1)
ax2.set_xticks(range(len(DAYS)))
ax2.set_xticklabels(DAYS, fontsize=8.5)
ax2.set_yticks(range(len(tf)))
ax2.set_yticklabels(tf.index, fontsize=8.5)
for i in range(tf.shape[0]):
    for j in range(tf.shape[1]):
        cnt = temporal.loc[tf.index[i], tf.columns[j]] if tf.index[i] in temporal.index else 0
        ax2.text(j,i,str(cnt),ha='center',va='center',fontsize=7.2,color='white' if tf.iloc[i,j]>0.5 else 'black')
ax2.set_title('Module day composition\n(counts inside cells)', fontsize=10.5, weight='bold')
cb2 = fig.colorbar(im2, ax=ax2, fraction=0.07, pad=0.02)
cb2.set_label('fraction', fontsize=8)
cb2.ax.tick_params(labelsize=7)
fig.suptitle(f'Target-day modules from sparse curated program-effect profiles (k={K_MAIN})', fontsize=13, weight='bold')
savefig(fig, f'03_module_centroids_and_temporal_specificity_k{K_MAIN}')

# ------------------------------
# Plot 4: post-hoc category counts for modules M2/M3/M4 and all modules
# ------------------------------
focus_mods = ['M2','M3','M4']
focus_counts = cat_counts[cat_counts['module_label'].isin(focus_mods) & (cat_counts['posthoc_category']!='Other')].copy()
# ensure categories sorted
cats = [c for c in ['Hippo robust','SAV1','BMP/TGF-SMAD','EC regulators','Apoptosis','Chromatin/cohesin/mediator','Signaling/candidate'] if c in focus_counts['posthoc_category'].unique()]
fig, ax = plt.subplots(figsize=(7.4, 3.8))
# stacked bars
bottom = np.zeros(len(focus_mods))
colors = {'Hippo robust':PURPLE,'SAV1':'#9C77B9','BMP/TGF-SMAD':POS,'EC regulators':'#2166AC','Apoptosis':'#C43C39','Chromatin/cohesin/mediator':'#7F7F7F','Signaling/candidate':'#4DAF4A'}
for cat in cats:
    vals=[]
    for m in focus_mods:
        val = focus_counts[(focus_counts['module_label']==m)&(focus_counts['posthoc_category']==cat)]['n']
        vals.append(int(val.iloc[0]) if len(val)>0 else 0)
    ax.bar(focus_mods, vals, bottom=bottom, color=colors.get(cat,'#AAAAAA'), label=cat, alpha=0.9)
    bottom += np.array(vals)
ax.set_ylabel('number of target-day rows')
ax.set_title('Post-hoc composition of selected modules (not used for clustering)', fontsize=12, weight='bold')
ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02,1), loc='upper left')
ax.spines[['top','right']].set_visible(False)
savefig(fig, f'04_M2_M3_M4_posthoc_category_composition_k{K_MAIN}')

# ------------------------------
# Plot 5: zoom heatmap for M2/M3/M4 rows with row labels
# ------------------------------
focus_assign = assign[assign['module_label'].isin(focus_mods)].copy()
# order by module then dendrogram order within global dendrogram
rank = {item:i for i,item in enumerate(ordered_items)}
focus_assign['dendro_rank'] = focus_assign['item'].map(rank)
focus_assign = focus_assign.sort_values(['module','dendro_rank'])
# cap number of rows maybe all rows M2-M4 manageable; if too many, show all because important
focus_X = X_sig.loc[focus_assign['item'], program_order]
focus_assign.to_csv(OUT / f'M2_M3_M4_target_day_members_k{K_MAIN}.csv', index=False)
row_labels = []
for _, r in focus_assign.iterrows():
    lab = f"{r['module_label']} | {r['target']}_{r['day']}"
    if r['target'] in ['NF2','LATS2','TAOK1','AMOTL2','SAV1']:
        lab = '★ ' + lab
    row_labels.append(lab)
fig, ax = plt.subplots(figsize=(9.8, max(5.2, 0.16*len(focus_X))))
vals = focus_X.values
lim = max(0.5, np.nanmax(np.abs(vals)))
im = ax.imshow(vals, aspect='auto', cmap='RdBu_r', vmin=-lim, vmax=lim)
ax.set_xticks(range(len(program_order)))
ax.set_xticklabels(program_xticklabels(), rotation=45, ha='right', fontsize=7.2)
ax.set_yticks(range(len(row_labels)))
ax.set_yticklabels(row_labels, fontsize=6.2)
# horizontal separators between modules
last = None
for i, m in enumerate(focus_assign['module_label'].tolist()):
    if last is not None and m != last:
        ax.axhline(i-0.5, color='black', lw=0.8)
    last = m
ax.set_title('Rows in modules M2/M3/M4 (★ Hippo/YAP genes; clustering was unsupervised)', fontsize=12, weight='bold')
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.015)
cb.set_label('significant effect (0 if FDR ≥ 0.05)', fontsize=8)
cb.ax.tick_params(labelsize=7)
fig.tight_layout()
savefig(fig, f'05_M2_M3_M4_target_day_effect_heatmap_k{K_MAIN}')

# ------------------------------
# Plot 6: concise module interpretation table as image
# ------------------------------
summary_focus = summary[summary['module'].isin(focus_mods)].copy()
fig, ax = plt.subplots(figsize=(11.5, 2.4))
ax.axis('off')
cols = ['module','n_rows','day_counts','top_positive_programs','top_negative_programs','hippo_items']
tab = ax.table(cellText=summary_focus[cols].values, colLabels=cols, cellLoc='left', loc='center')
tab.auto_set_font_size(False)
tab.set_fontsize(7.0)
tab.scale(1, 1.35)
for (r,c), cell in tab.get_celld().items():
    if r == 0:
        cell.set_text_props(weight='bold')
        cell.set_facecolor('#F0F0F0')
    else:
        cell.set_facecolor('white')
ax.set_title('Interpretation of selected modules M2/M3/M4', fontsize=12, weight='bold', pad=10)
savefig(fig, f'06_M2_M3_M4_interpretation_table_k{K_MAIN}')

# ------------------------------
# README text summary
# ------------------------------
readme = f"""
Target-day module clustering from scratch
=========================================

Input tables:
- {EFFECT_ZIP.name}: contains curated_program_effects_aggregated.csv, target_by_curated_program_effect_matrix.csv, target_by_curated_program_fdr_matrix.csv
- {ANNOTATION_CSV}

Clustering steps:
1. Build target-day profiles. Each row is one target at one day, e.g. NF2_D3.
2. Columns are expert-curated cNMF programs: {', '.join(program_order)}.
3. For the sparse pattern clustering used to interpret M2/M3/M4, values are program effect sizes only if FDR < {FDR_CUTOFF}; otherwise values are set to 0.
4. Columns are z-scored.
5. Rows are clustered by cosine distance with average-linkage hierarchical clustering.
6. A range of k values is evaluated by silhouette score, bootstrap stability, Davies-Bouldin score, and small-module penalty.
7. The main plots use k = {K_MAIN}, matching the sparse target-day clustering that generated M2/M3/M4.

Important caveat:
Module labels M1, M2, etc. are specific to this analysis and this k. They are not comparable to target-level modules from other scripts.

Core interpretation:
- M2: inspect module_summary_k{K_MAIN}.csv and M2_M3_M4_target_day_members_k{K_MAIN}.csv.
- M3: contains most D3 Hippo/YAP target-day entries in the k={K_MAIN} sparse target-day clustering.
- M4: broad EC endpoint-loss / P27-like module with BMP/SMAD, EC regulators, PRKCD/PLCG1 etc.

Post-hoc categories were added only after clustering and were not used to define modules.
"""
(OUT / 'README_clustering_from_beginning.txt').write_text(readme)

# Package
zip_path = BASE / 'target_day_module_clustering_from_scratch_package.zip'
with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
    for fp in OUT.glob('*'):
        if fp.is_file():
            z.write(fp, arcname=fp.name)

print('Done.')
print(f'Output directory: {OUT}')
print(f'Package: {zip_path}')
print(f'Recommended k by composite: {recommended_k}. Main k used: {K_MAIN}')
print(summary[summary['module'].isin(focus_mods)].to_string(index=False))
