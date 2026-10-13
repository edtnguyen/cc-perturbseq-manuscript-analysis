from pathlib import Path
import os, math, textwrap, zipfile
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'DejaVu Sans'

SOURCE = Path('/mnt/data/Figure1d_all244_redesign_v2/source_data/Fig1d_all244_dotstrip_source.tsv')
if not SOURCE.exists():
    raise FileNotFoundError(SOURCE)
OUT = Path('/mnt/data/Figure1d_compact_mainfigure_candidates_v4')
(OUT/'individual_panels').mkdir(parents=True, exist_ok=True)
(OUT/'source_data').mkdir(parents=True, exist_ok=True)
(OUT/'scripts').mkdir(parents=True, exist_ok=True)

df0 = pd.read_csv(SOURCE, sep='\t')
effect_col = 'u_test_effect_log2_CDH5pos_vs_CDH5neg' if 'u_test_effect_log2_CDH5pos_vs_CDH5neg' in df0.columns else 'plot_effect'
df = df0[['gene','category',effect_col]].copy()
df.columns = ['gene','category','effect']
df['effect'] = pd.to_numeric(df['effect'], errors='coerce')
df = df.dropna(subset=['effect']).drop_duplicates('gene')
df['direction'] = np.where(df['effect'] < 0, 'Required', 'Restricts')

# Cleaner names/order by total count descending for top-to-bottom
cat_counts = df.groupby('category', observed=True).agg(
    n_total=('gene','size'),
    n_required=('direction', lambda s: int((s=='Required').sum())),
    n_restricts=('direction', lambda s: int((s=='Restricts').sum()))
).reset_index().sort_values('n_total', ascending=False)
cat_counts['y'] = np.arange(len(cat_counts))[::-1]
cat_order = cat_counts['category'].tolist()
y_map = dict(zip(cat_counts['category'], cat_counts['y']))

nice_cat = {
    'Transcription / chromatin': 'Transcription / chromatin',
    'RNA processing': 'RNA processing',
    'Developmental signaling': 'Developmental signaling',
    'Protein homeostasis / trafficking': 'Proteostasis / PTM',
    'Protein homeostasis/trafficking': 'Proteostasis / PTM',
    'Trafficking / cytoskeleton / adhesion': 'Trafficking / cytoskeleton',
    'Metabolism / mitochondria': 'Metabolism / mitochondria',
    'Translation / ribosome': 'Translation / ribosome',
    'Cell cycle / DNA repair': 'Cell cycle / DNA repair',
    'Apoptosis / stress': 'Apoptosis / stress',
    'Other / unclear': 'Other / unclear',
    'Other / unclassified': 'Other / unclear',
}
def cat_label(c): return nice_cat.get(c, str(c))

def wrap(x, width=25): return textwrap.fill(str(x), width=width)

# Representative genes per category: prefer known anchors + strongest required/restrict hits
anchors = {'KDR','ETV2','CDH5','PLCG1','BMPR1A','SMAD1','SMAD4','SMAD5','ACVR2A','ACVR2B','ACVR1B','MIXL1','EOMES','TBXT','HAND1','VHL','PTEN','AMOTL2','LATS2','KCTD10','SMARCA4','SMARCC1','MED12'}
examples = []
for c in cat_order:
    sub = df[df['category']==c].copy()
    forced = sub[sub['gene'].isin(anchors)].copy()
    top_req = sub[sub['effect']<0].sort_values('effect').head(3)
    top_res = sub[sub['effect']>0].sort_values('effect', ascending=False).head(1)
    sel = pd.concat([forced, top_req, top_res]).drop_duplicates('gene')
    # Sort by effect direction/magnitude, cap to 5 for text column
    sel = sel.reindex(sel['effect'].abs().sort_values(ascending=False).index).head(5)
    examples.append({'category':c, 'examples': ', '.join(sel['gene'].tolist())})
examples = pd.DataFrame(examples)
cat_counts = cat_counts.merge(examples, on='category', how='left')

df.to_csv(OUT/'source_data'/'Fig1d_compact_all244_source.tsv', sep='\t', index=False)
cat_counts.to_csv(OUT/'source_data'/'Fig1d_compact_category_counts_examples.tsv', sep='\t', index=False)

BLUE = '#4C78A8'
ORANGE = '#E68613'
LIGHT_BLUE = '#AFC6DE'
LIGHT_ORANGE = '#F2C48A'
DARK = '#222222'

# --------------------
# C1: main-figure compact bars + readable examples
# --------------------
fig, ax = plt.subplots(figsize=(7.6, 4.15))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')
bar_h = 0.42
# alternating row backgrounds
for _, r in cat_counts.iterrows():
    if int(r['y']) % 2 == 0:
        ax.add_patch(Rectangle((-0.5, r['y']-0.44), 88, 0.88, facecolor='#FAFAFA', edgecolor='none', zorder=0))

ax.barh(cat_counts['y'], cat_counts['n_required'], height=bar_h, color=BLUE, alpha=0.94, edgecolor='none', zorder=2)
ax.barh(cat_counts['y'], cat_counts['n_restricts'], height=bar_h, left=cat_counts['n_required'], color=ORANGE, alpha=0.94, edgecolor='none', zorder=3)

# labels and counts
for _, r in cat_counts.iterrows():
    ax.text(-1.3, r['y'], f"{wrap(cat_label(r['category']), 23)}\n(n={int(r['n_total'])})", ha='right', va='center', fontsize=7.3, fontweight='bold', color=DARK)
    ax.text(r['n_total'] + 0.55, r['y'], f"{int(r['n_required'])}|{int(r['n_restricts'])}", ha='left', va='center', fontsize=7.0, color='#555555')
    ax.text(68.0, r['y'], str(r['examples']), ha='left', va='center', fontsize=6.9, color='#3B3B3B')

ax.text(68.0, cat_counts['y'].max()+0.52, 'Representative genes', ha='left', va='bottom', fontsize=7.8, fontweight='bold', color=DARK)
ax.text(0, cat_counts['y'].max()+0.52, 'Regulator count', ha='left', va='bottom', fontsize=7.8, fontweight='bold', color=DARK)
ax.text(51, -0.75, 'required | restricts', ha='left', va='center', fontsize=6.8, color='#666666')

ax.set_xlim(-22, 92)
ax.set_ylim(-0.75, cat_counts['y'].max()+0.9)
ax.set_yticks([])
ax.set_xlabel('Number of screen-derived regulators', fontsize=8.2)
ax.set_xticks([0,20,40,60])
ax.tick_params(axis='x', labelsize=7.2)
ax.grid(axis='x', color='#E0E0E0', linewidth=0.55)
ax.set_axisbelow(True)
for s in ['top','right','left']:
    ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color('#777777')
ax.spines['bottom'].set_linewidth(0.8)
ax.set_title('Functional classes of 244 screen-derived regulators', loc='left', fontsize=11.0, fontweight='bold', pad=6)
ax.text(0, 1.015, 'Compact category summary; gene names mark representative high-effect hits.', transform=ax.transAxes, ha='left', va='bottom', fontsize=7.4, color='#555555')
handles = [Line2D([0],[0], marker='s', color='none', markerfacecolor=BLUE, markersize=6.6, label='Required'),
           Line2D([0],[0], marker='s', color='none', markerfacecolor=ORANGE, markersize=6.6, label='Restricts')]
ax.legend(handles=handles, frameon=False, ncol=2, fontsize=7.2, loc='lower right', bbox_to_anchor=(0.99, 0.00), handletextpad=0.4, columnspacing=0.8)
fig.tight_layout()
base = OUT/'individual_panels'/'Fig1d_compact_category_summary_C1'
for ext in ['png','pdf','svg']:
    fig.savefig(base.with_suffix('.'+ext), dpi=350 if ext=='png' else None, bbox_inches='tight')
plt.close(fig)

# --------------------
# C2: compact all-244 barcode, labels reduced
# --------------------
plot_df = df.copy()
plot_df['y'] = plot_df['category'].map(y_map).astype(float)
rng = np.random.default_rng(123)
plot_df['yj'] = plot_df['y'] + rng.uniform(-0.19, 0.19, len(plot_df))
plot_df['color'] = np.where(plot_df['effect'] < 0, BLUE, ORANGE)
label_set = set()
for c in cat_order:
    sub = df[df['category']==c]
    # one strongest required and one strongest restrict per category if present
    if (sub['effect']<0).any(): label_set.add(sub[sub['effect']<0].sort_values('effect').iloc[0]['gene'])
    if (sub['effect']>0).any(): label_set.add(sub[sub['effect']>0].sort_values('effect', ascending=False).iloc[0]['gene'])
label_set |= anchors & set(df['gene'])
# cap labels to high effect anchors only
label_df = plot_df[plot_df['gene'].isin(label_set)].copy()
label_df = label_df.reindex(label_df['effect'].abs().sort_values(ascending=False).index).head(28)

fig2, ax2 = plt.subplots(figsize=(7.6, 4.15))
fig2.patch.set_facecolor('white')
ax2.set_facecolor('white')
for _, r in cat_counts.iterrows():
    if int(r['y']) % 2 == 0:
        ax2.add_patch(Rectangle((-4.2, r['y']-0.44), 8.4, 0.88, facecolor='#FAFAFA', edgecolor='none', zorder=0))
ax2.axvline(0, color='#555555', linewidth=0.8, zorder=1)
ax2.scatter(plot_df['effect'], plot_df['yj'], c=plot_df['color'], s=12, alpha=0.66, edgecolors='white', linewidths=0.18, zorder=2)

for _, r in label_df.iterrows():
    dx = -0.055 if r['effect'] < 0 else 0.055
    ha = 'right' if r['effect'] < 0 else 'left'
    ax2.text(r['effect'] + dx, r['yj'], r['gene'], ha=ha, va='center', fontsize=6.1, color=DARK,
             bbox=dict(boxstyle='round,pad=0.05', facecolor='white', edgecolor='none', alpha=0.70), zorder=3)

for _, r in cat_counts.iterrows():
    ax2.text(-4.35, r['y'], f"{wrap(cat_label(r['category']), 24)}\n(n={int(r['n_total'])})", ha='right', va='center', fontsize=7.3, fontweight='bold', color=DARK)

xlim = max(1.0, float(np.nanmax(np.abs(plot_df['effect'])))) * 1.22
ax2.set_xlim(-xlim, xlim)
ax2.set_ylim(-0.75, cat_counts['y'].max()+0.9)
ax2.set_yticks([])
ax2.set_xlabel('U-test log2 fold-change (CDH5+ vs CDH5-)', fontsize=8.2)
ax2.tick_params(axis='x', labelsize=7.2)
ax2.grid(axis='x', color='#E0E0E0', linewidth=0.55)
ax2.set_axisbelow(True)
for s in ['top','right','left']:
    ax2.spines[s].set_visible(False)
ax2.spines['bottom'].set_color('#777777')
ax2.spines['bottom'].set_linewidth(0.8)
ax2.set_title('All 244 regulators by functional class and screen effect', loc='left', fontsize=11.0, fontweight='bold', pad=6)
ax2.text(0, 1.015, 'Each point is one regulator; labels show selected high-effect examples.', transform=ax2.transAxes, ha='left', va='bottom', fontsize=7.4, color='#555555')
ax2.legend(handles=handles, frameon=False, ncol=2, fontsize=7.2, loc='lower right', bbox_to_anchor=(0.99, 0.00), handletextpad=0.4, columnspacing=0.8)
fig2.tight_layout()
base2 = OUT/'individual_panels'/'Fig1d_compact_all244_barcode_C2'
for ext in ['png','pdf','svg']:
    fig2.savefig(base2.with_suffix('.'+ext), dpi=350 if ext=='png' else None, bbox_inches='tight')
plt.close(fig2)

# Copy script and zip
try:
    import shutil
    shutil.copyfile(__file__, OUT/'scripts'/'build_compact_fig1d_candidates_v4.py')
except Exception:
    (OUT/'scripts'/'build_compact_fig1d_candidates_v4.py').write_text('# generated in runtime\n')
README = f"""Compact Figure 1d main-figure candidates v4\n\nInput: {SOURCE}\n\nC1 = compact category counts with representative gene names.\nC2 = compact barcode showing all 244 individual regulators as points.\n\nThe full 244-gene source and category counts are in source_data/.\n"""
(OUT/'README.md').write_text(README)
zip_path = Path('/mnt/data/Figure1d_compact_mainfigure_candidates_v4_package.zip')
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
    for p in OUT.rglob('*'):
        z.write(p, arcname=str(p.relative_to(OUT.parent)))
print('Created', zip_path)
