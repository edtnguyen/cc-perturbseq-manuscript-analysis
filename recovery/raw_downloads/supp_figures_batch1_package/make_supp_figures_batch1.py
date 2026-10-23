import os
import zipfile
import shutil
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import openpyxl

OUTDIR = '/mnt/data/supp_figures_batch1'
os.makedirs(OUTDIR, exist_ok=True)

CNMF_XLSX = '/mnt/data/cNMF_50_0_2_032526_6co.xlsx'
TARGET_KD_CSV = '/mnt/data/combined_TargetKD.csv'

mpl.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 7,
    'axes.labelsize': 7,
    'axes.titlesize': 8,
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5,
    'legend.fontsize': 6.2,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
    'axes.linewidth': 0.7,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
})


def read_xlsx_sheet(path, sheet_name):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet_name]
    rows = ws.iter_rows(values_only=True)
    header = next(rows)
    data = list(rows)
    df = pd.DataFrame(data, columns=header)
    if None in df.columns:
        df = df.drop(columns=[None])
    return df


def savefig(fig, stem):
    png = os.path.join(OUTDIR, f'{stem}.png')
    pdf = os.path.join(OUTDIR, f'{stem}.pdf')
    svg = os.path.join(OUTDIR, f'{stem}.svg')
    fig.savefig(png, dpi=450, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf, bbox_inches='tight', facecolor='white')
    fig.savefig(svg, bbox_inches='tight', facecolor='white')
    plt.close(fig)

# Read cNMF perturbation association table.
assoc = read_xlsx_sheet(CNMF_XLSX, 'Perturbation Association 1')
assoc['program_name'] = pd.to_numeric(assoc['program_name'], errors='coerce').astype('Int64')
assoc['approx_log2FC'] = pd.to_numeric(assoc['approx_log2FC'], errors='coerce')
assoc['adj_pval'] = pd.to_numeric(assoc['adj_pval'], errors='coerce')
assoc['Sample'] = assoc['Sample'].astype(str)
assoc['target_name'] = assoc['target_name'].astype(str)
assoc.to_csv(os.path.join(OUTDIR, 'source_cNMF_perturbation_association_clean.csv'), index=False)

# ------------------------------
# Panel 1: target knockdown effect by day
# ------------------------------
kd = pd.read_csv(TARGET_KD_CSV)
kd['day_label'] = kd['day'].map(lambda x: f'D{int(x)}' if pd.notnull(x) else str(x))
kd['log2FC'] = pd.to_numeric(kd['log_2_fold_change'], errors='coerce')
kd['significant'] = kd['significant'].astype(bool)
kd.to_csv(os.path.join(OUTDIR, 'source_knockdown_efficiency_by_day.csv'), index=False)

days = ['D0', 'D1', 'D2', 'D3']
positions = np.arange(len(days))
np.random.seed(1)
fig, ax = plt.subplots(figsize=(3.6, 2.75))
for i, d in enumerate(days):
    vals = kd.loc[kd['day_label'] == d, 'log2FC'].dropna().values
    if len(vals) == 0:
        continue
    q1, med, q3 = np.percentile(vals, [25, 50, 75])
    whisk_low, whisk_high = np.percentile(vals, [5, 95])
    ax.add_patch(Rectangle((i - 0.22, q1), 0.44, q3 - q1, facecolor='#E8E8E8', edgecolor='#333333', linewidth=0.7, zorder=2))
    ax.plot([i - 0.22, i + 0.22], [med, med], color='#333333', lw=1.0, zorder=3)
    ax.plot([i, i], [whisk_low, q1], color='#333333', lw=0.7)
    ax.plot([i, i], [q3, whisk_high], color='#333333', lw=0.7)
    ax.plot([i - 0.12, i + 0.12], [whisk_low, whisk_low], color='#333333', lw=0.7)
    ax.plot([i - 0.12, i + 0.12], [whisk_high, whisk_high], color='#333333', lw=0.7)
    xj = i + np.random.uniform(-0.16, 0.16, len(vals))
    ax.scatter(xj, vals, s=7, color='#4A6FA5', alpha=0.33, edgecolors='none', rasterized=True, zorder=1)
ax.axhline(0, color='#9A9A9A', ls='--', lw=0.7)
ax.set_xticks(positions)
ax.set_xticklabels(days)
ax.set_ylabel('Target expression log$_2$FC\n(target guide vs control)')
ax.set_title('Target knockdown across timepoints', weight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ymin, ymax = ax.get_ylim()
for i, d in enumerate(days):
    n = kd.loc[kd['day_label'] == d].shape[0]
    ax.text(i, ymin - 0.055 * (ymax - ymin), f'n={n}', ha='center', va='top', fontsize=6, color='#555555')
savefig(fig, 'supp_panel_knockdown_efficiency_by_day')

# ------------------------------
# Panel 2: regulatory burden per target by day
# ------------------------------
burden = (assoc.assign(significant=assoc['adj_pval'] < 0.05)
          .groupby(['target_name', 'Sample'], as_index=False)
          .agg(n_significant_programs=('significant', 'sum'),
               max_abs_effect=('approx_log2FC', lambda x: np.nanmax(np.abs(x)))))
burden.to_csv(os.path.join(OUTDIR, 'source_regulatory_burden_by_day.csv'), index=False)

fig, ax = plt.subplots(figsize=(3.6, 2.75))
np.random.seed(2)
for i, d in enumerate(days):
    vals = burden.loc[burden['Sample'] == d, 'n_significant_programs'].values
    if len(vals) == 0:
        continue
    q1, med, q3 = np.percentile(vals, [25, 50, 75])
    wh1, wh2 = np.percentile(vals, [5, 95])
    ax.add_patch(Rectangle((i - 0.22, q1), 0.44, q3 - q1, facecolor='#E8E8E8', edgecolor='#333333', linewidth=0.7, zorder=2))
    ax.plot([i - 0.22, i + 0.22], [med, med], color='#333333', lw=1)
    ax.plot([i, i], [wh1, q1], color='#333333', lw=0.7)
    ax.plot([i, i], [q3, wh2], color='#333333', lw=0.7)
    ax.plot([i - 0.12, i + 0.12], [wh1, wh1], color='#333333', lw=0.7)
    ax.plot([i - 0.12, i + 0.12], [wh2, wh2], color='#333333', lw=0.7)
    ax.scatter(i + np.random.uniform(-0.16, 0.16, len(vals)), vals, s=7, color='#8661A9', alpha=0.28, edgecolors='none', rasterized=True, zorder=1)
ax.set_xticks(positions)
ax.set_xticklabels(days)
ax.set_ylabel('Significant cNMF program effects\nper target (FDR < 0.05)')
ax.set_title('Regulatory burden across days', weight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
savefig(fig, 'supp_panel_regulatory_burden_by_day')

# ------------------------------
# Panel 3: active regulators by day
# ------------------------------
active = (burden.assign(active=burden['n_significant_programs'] > 0)
          .groupby('Sample', as_index=False)
          .agg(active_regulators=('active', 'sum'), total_regulators=('target_name', 'nunique')))
active['fraction_active'] = active['active_regulators'] / active['total_regulators']
active.to_csv(os.path.join(OUTDIR, 'source_active_regulators_by_day.csv'), index=False)
active_plot = active.set_index('Sample').reindex(days).reset_index()
fig, ax = plt.subplots(figsize=(3.2, 2.6))
ax.bar(np.arange(len(days)), active_plot['active_regulators'], color='#6F80B9', edgecolor='white', width=0.68)
for i, row in active_plot.iterrows():
    ax.text(i, row['active_regulators'] + 3, f"{int(row['active_regulators'])}", ha='center', va='bottom', fontsize=7)
ax.set_xticks(positions)
ax.set_xticklabels(days)
ax.set_ylabel('Targets with >=1 significant\nprogram effect')
ax.set_title('Active regulators by day', weight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
savefig(fig, 'supp_panel_active_regulators_by_day')

# ------------------------------
# Heatmap panels: raw selected program effects
# ------------------------------
program_selection = [
    ('D1 P13', 'D1', 13),
    ('D1 P30', 'D1', 30),
    ('D2 P15', 'D2', 15),
    ('D2 P17', 'D2', 17),
    ('D3 P3', 'D3', 3),
    ('D3 P34', 'D3', 34),
]

def build_selected_effect_matrix(targets):
    rows = []
    for target in targets:
        row = {'target': target}
        for label, sample, prog in program_selection:
            sub = assoc[(assoc['target_name'].str.upper() == target.upper()) & (assoc['Sample'] == sample) & (assoc['program_name'] == prog)]
            if len(sub) > 0:
                row[label] = float(sub['approx_log2FC'].iloc[0])
                row[label + '_fdr'] = float(sub['adj_pval'].iloc[0])
            else:
                row[label] = np.nan
                row[label + '_fdr'] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def plot_heatmap(mat, title, stem):
    labels = [x[0] for x in program_selection]
    mat.to_csv(os.path.join(OUTDIR, f'source_{stem}.csv'), index=False)
    values = mat[labels].astype(float).values
    fig_h = max(1.9, 0.32 * len(mat) + 0.75)
    fig, ax = plt.subplots(figsize=(4.25, fig_h))
    vmax = np.nanpercentile(np.abs(values), 95)
    vmax = max(0.4, min(2.2, vmax))
    cmap = LinearSegmentedColormap.from_list('effect_cmap', ['#2166AC', '#F7F7F7', '#D6604D'])
    im = ax.imshow(values, aspect='auto', cmap=cmap, vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha='right')
    ax.set_yticks(np.arange(len(mat)))
    ax.set_yticklabels(mat['target'])
    ax.set_title(title, weight='bold', pad=8)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for x in np.arange(-0.5, len(labels), 1):
        ax.axvline(x, color='white', lw=0.8)
    for y in np.arange(-0.5, len(mat), 1):
        ax.axhline(y, color='white', lw=0.8)
    for i in range(len(mat)):
        for j, label in enumerate(labels):
            f = mat.iloc[i][label + '_fdr']
            if pd.notnull(f) and f < 0.05:
                ax.text(j, i, '*', ha='center', va='center', fontsize=9, color='black', fontweight='bold')
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label('Effect size\n(approx. log$_2$FC)', rotation=270, labelpad=20)
    cbar.ax.tick_params(labelsize=6)
    fig.text(0.99, 0.02, '* FDR < 0.05', ha='right', va='bottom', fontsize=6)
    savefig(fig, stem)

hippo_targets = ['NF2', 'LATS2', 'SAV1', 'TAOK1', 'AMOTL2']
hippo = build_selected_effect_matrix(hippo_targets)
plot_heatmap(hippo, 'Hippo pathway effects on selected programs', 'supp_panel_Hippo_selected_program_effects_heatmap')

tgfb_targets = ['ACVR1B', 'ACVR2A', 'ACVR2B', 'BMPR1A', 'TGFBR2', 'SMAD1', 'SMAD4', 'SMAD5']
tgfb = build_selected_effect_matrix(tgfb_targets)
tgfb = tgfb[tgfb[[x[0] for x in program_selection]].notna().any(axis=1)].reset_index(drop=True)
plot_heatmap(tgfb, 'TGF-beta/BMP pathway effects on selected programs', 'supp_panel_TGFB_BMP_selected_program_effects_heatmap')

with open(os.path.join(OUTDIR, 'README.txt'), 'w') as f:
    f.write('Supplement figure batch 1 regenerated from current uploads.\n\n')
    f.write('Inputs:\n')
    f.write('- cNMF_50_0_2_032526_6co.xlsx\n')
    f.write('- combined_TargetKD.csv\n\n')
    f.write('Outputs are PDF/SVG/PNG panels plus source CSVs. SVG/PDF keep text editable in Illustrator.\n')
    f.write('Dense nonessential point layers are rasterized in some panels to keep vector files responsive.\n')

shutil.copyfile(__file__, os.path.join(OUTDIR, 'make_supp_figures_batch1.py'))
zip_path = '/mnt/data/supp_figures_batch1_package.zip'
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
    for fn in os.listdir(OUTDIR):
        fp = os.path.join(OUTDIR, fn)
        if os.path.isfile(fp):
            z.write(fp, arcname=fn)
print('Wrote', OUTDIR, zip_path)
