from pathlib import Path
import os, zipfile, shutil, re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy.stats import pearsonr

ROOT = Path('/mnt/data')
OUT = ROOT / 'Figure1_text_editable_1cde_fresh'
PAN = OUT / 'individual_panels'
SRC = OUT / 'source_data'
SCR = OUT / 'scripts'
for d in [OUT, PAN, SRC, SCR]:
    d.mkdir(parents=True, exist_ok=True)

# Illustrator-friendly vector text settings
mpl.rcParams.update({
    'font.family': 'Liberation Sans',
    'font.size': 9,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'axes.unicode_minus': False,
})

DARK = '#20262B'
MID = '#66707A'
BLUE = '#2E6F97'
LIGHT_BLUE = '#D9E6EE'
ORANGE = '#D87528'
PURPLE = '#74539A'
GREY = '#BFC4C7'

# Source files already present in the active runtime
C_SOURCE = ROOT / 'Figure1_restyled_1c1d_remade/source_data/Fig1c_nonHippo_known_regulators_dotmatrix.tsv'
D_SOURCE = ROOT / 'Figure1_restyled_1c1d_remade/source_data/Fig1d_functional_categories_restyled.tsv'
E_SOURCE = ROOT / 'Figure1_nature_style_remade_fresh/source_data/Fig1e_validation_source.tsv'
for p in [C_SOURCE, D_SOURCE, E_SOURCE]:
    if not p.exists():
        raise FileNotFoundError(p)

c = pd.read_csv(C_SOURCE, sep='\t')
d = pd.read_csv(D_SOURCE, sep='\t')
e = pd.read_csv(E_SOURCE, sep='\t')

c.to_csv(SRC / 'Fig1c_nonHippo_known_regulators_dotmatrix.tsv', sep='\t', index=False)
d.to_csv(SRC / 'Fig1d_functional_categories_restyled.tsv', sep='\t', index=False)
e.to_csv(SRC / 'Fig1e_validation_source.tsv', sep='\t', index=False)


def save_all(fig, stem, width, height):
    fig.set_size_inches(width, height)
    fig.savefig(PAN / f'{stem}.svg', bbox_inches='tight', facecolor='white')
    fig.savefig(PAN / f'{stem}.pdf', bbox_inches='tight', facecolor='white')
    fig.savefig(PAN / f'{stem}.png', dpi=420, bbox_inches='tight', facecolor='white')
    plt.close(fig)


# -----------------------------------------------------------------------------
# Figure 1c — dot matrix, no Hippo genes
# -----------------------------------------------------------------------------
assay_order = ['MaGECK', 'U-test', 'Validation']
gene_order = [
    'KDR','CDH5','ETV2','PLCG1','KCTD10','VHL','PTEN',
    'BMPR1A','SMAD1','SMAD5','SMAD4','ACVR2A',
    'EOMES','HAND1','MIXL1','TBXT'
]
c = c[c['gene'].isin(gene_order)].copy()
# ensure order top to bottom
c['gene'] = pd.Categorical(c['gene'], categories=gene_order[::-1], ordered=True)
c['assay'] = pd.Categorical(c['assay'], categories=assay_order, ordered=True)

fig, ax = plt.subplots(figsize=(5.35, 5.25))

# custom blue-to-white colormap used in prior panel
cmap = LinearSegmentedColormap.from_list('blue_white', ['#154C73', '#8FAFC1', '#EEF3F5'])
norm = Normalize(vmin=-4.5, vmax=0)

xmap = {a: i for i, a in enumerate(assay_order)}
ymap = {g: i for i, g in enumerate(gene_order[::-1])}

for _, r in c.iterrows():
    x = xmap[str(r['assay'])]
    y = ymap[str(r['gene'])]
    size = np.interp(float(r['minus_log10_fdr']), [0.4, 6.0], [45, 330])
    ax.scatter(x, y, s=size, c=[cmap(norm(float(r['effect'])))], edgecolor='#2F3438', linewidth=0.8, zorder=3)

# separators between the three gene groups
for sep in [len(gene_order)-7-0.5, len(gene_order)-12-0.5]:
    ax.axhline(sep, color='#D7DDE2', linewidth=1.1, zorder=1)

ax.set_xlim(-0.6, 2.6)
ax.set_ylim(-0.6, len(gene_order)-0.4)
ax.set_xticks(range(3), assay_order, fontsize=11, fontweight='bold', color='#273038')
ax.set_yticks(range(len(gene_order)), gene_order[::-1], fontsize=9.4, fontstyle='italic', color='#273038')
ax.tick_params(axis='both', length=0)
for s in ax.spines.values(): s.set_visible(False)
ax.set_title('Known EC regulators across assays', loc='left', fontsize=14.5, fontweight='bold', color='#171C21', pad=28)
ax.text(0, 1.01, 'Non-Hippo hits; color = effect, size = significance', transform=ax.transAxes,
        fontsize=10.2, color='#606B74', va='bottom')

# Effect colorbar
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, fraction=0.052, pad=0.06, aspect=11)
cbar.set_ticks([-4, -2, 0])
cbar.set_ticklabels(['-4','-2','0'])
cbar.ax.tick_params(labelsize=8, length=0)
cbar.outline.set_visible(False)
cbar.ax.set_title('Effect\nlog2FC', fontsize=9.2, fontweight='bold', pad=8)

# Size legend as editable text and vector circles
legend_vals = [1, 3, 5]
handles = []
for v in legend_vals:
    ms = np.sqrt(np.interp(v, [0.4, 6.0], [45, 330]))
    handles.append(Line2D([0],[0], marker='o', linestyle='None', color='#2F3438', markerfacecolor='#2F3438', markersize=ms/1.35, label=str(v)))
leg = ax.legend(handles=handles, title='-log10(FDR)', frameon=False, fontsize=8.3, title_fontsize=9.0,
                loc='lower left', bbox_to_anchor=(1.05, 0.02), borderaxespad=0, labelspacing=1.1)
leg.get_title().set_fontweight('bold')

fig.subplots_adjust(left=0.21, right=0.82, top=0.82, bottom=0.10)
save_all(fig, 'Fig1c_known_regulators_nonHippo_text_editable', 5.35, 5.25)


# -----------------------------------------------------------------------------
# Figure 1d — horizontal stacked bars
# -----------------------------------------------------------------------------
d = d.copy()
# order largest first in visual top-to-bottom
plot_d = d.sort_values('n_regulators', ascending=True).reset_index(drop=True)
y = np.arange(len(plot_d))

fig, ax = plt.subplots(figsize=(6.25, 5.25))
ax.barh(y, plot_d['n_validated'], color=BLUE, height=0.66, edgecolor='none', label='Validated')
ax.barh(y, plot_d['screen_only'], left=plot_d['n_validated'], color=LIGHT_BLUE, height=0.66, edgecolor='none', label='Not validated')

for i, r in plot_d.iterrows():
    ax.text(r['n_regulators'] + 0.8, i, r['validated_label'], va='center', ha='left', fontsize=7.8, color='#34404A')
    ax.text(plot_d['n_regulators'].max() + 8.2, i, str(int(r['n_restricting'])), va='center', ha='center', fontsize=8.2, fontweight='bold', color=ORANGE)

ax.set_yticks(y, plot_d['category_label'], fontsize=9.5, color='#273038')
ax.set_xlim(0, plot_d['n_regulators'].max() + 11.5)
ax.set_xticks([0,20,40,60])
ax.tick_params(axis='x', labelsize=8.8, colors='#5A646C')
ax.tick_params(axis='y', length=0)
ax.set_xlabel('Number of regulators', fontsize=10.5, fontweight='bold', color='#273038', labelpad=8)
ax.set_title('Functional categories and validation', loc='left', fontsize=14.5, fontweight='bold', color='#171C21', pad=28)
ax.text(0, 1.01, 'Dark = validated; label = validated/total; orange = restrictors', transform=ax.transAxes,
        fontsize=10.0, color='#606B74', va='bottom')
for s in ['top','right','left']:
    ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color('#4A535A')
ax.spines['bottom'].set_linewidth(0.8)
ax.grid(False)
ax.legend(frameon=False, ncol=2, loc='upper right', bbox_to_anchor=(0.88,1.13), fontsize=9.2)
# editable heading over restricting-count column
ax.text(plot_d['n_regulators'].max() + 8.2, len(plot_d)-0.05, 'Restricts', ha='center', va='bottom', fontsize=8.2, color=ORANGE, fontweight='bold')

fig.subplots_adjust(left=0.34, right=0.96, top=0.82, bottom=0.12)
save_all(fig, 'Fig1d_functional_categories_text_editable', 6.25, 5.25)


# -----------------------------------------------------------------------------
# Figure 1e — main validation scatter
# -----------------------------------------------------------------------------
class_map = {
    'genome_wide_hit': ('Screen-derived regulator', BLUE),
    'genome_wide_non_hit': ('Random non-hit', GREY),
    'tri_perturb_seq': ('Additional curated target', PURPLE),
}

fig, ax = plt.subplots(figsize=(5.2, 4.45))
for typ, (label, color) in class_map.items():
    sub = e[e['Type'] == typ]
    ax.scatter(sub['Original FC U Test'], sub['Average U Test FC'], s=23, c=color, alpha=0.78,
               edgecolors='white', linewidths=0.35, label=label)

x = e['Original FC U Test'].to_numpy(float)
yv = e['Average U Test FC'].to_numpy(float)
coef = np.polyfit(x, yv, 1)
xs = np.linspace(-4.2, 2.5, 120)
ax.plot(xs, coef[0]*xs + coef[1], color=DARK, linewidth=1.35)
ax.axhline(0, linestyle='--', linewidth=0.8, color='#9EA4A9')
ax.axvline(0, linestyle='--', linewidth=0.8, color='#9EA4A9')
r, p = pearsonr(x, yv)
ax.text(0.03, 0.96, f'n = {len(e)}\nPearson r = {r:.2f}', transform=ax.transAxes,
        ha='left', va='top', fontsize=8.5, color=DARK)

for g, dx, dy in [('KDR',-.18,.18),('CDH5',.10,.13),('ETV2',.10,-.15)]:
    rr = e[e['Gene'] == g]
    if len(rr):
        xx = float(rr['Original FC U Test'].iloc[0]); yy = float(rr['Average U Test FC'].iloc[0])
        ax.text(xx+dx, yy+dy, g, fontsize=8.0, fontstyle='italic', color=DARK)

hit = e[e['Type'] == 'genome_wide_hit']
sig = int((hit['FDR-corrected'] < .05).sum())
same = int(((hit['Original FC U Test'] * hit['Average U Test FC']) > 0).sum())
ax.set_title('Secondary screening validates effects', loc='left', fontsize=13.5, fontweight='bold', color='#171C21', pad=25)
ax.text(0, 1.01, f'{sig}/244 screen hits significant; {same}/{len(hit)} same direction', transform=ax.transAxes,
        fontsize=9.0, color='#606B74', va='bottom')
ax.set_xlabel('Genome-wide U-test log2FC', fontsize=10.0, fontweight='bold')
ax.set_ylabel('Validation U-test log2FC', fontsize=10.0, fontweight='bold')
ax.set_xlim(-4.1,2.5); ax.set_ylim(-4.2,4.3)
ax.tick_params(labelsize=8.5)
for s in ['top','right']:
    ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=7.3, loc='lower right')
fig.subplots_adjust(left=0.17, right=0.96, top=0.80, bottom=0.16)
save_all(fig, 'Fig1e_validation_main_text_editable', 5.2, 4.45)


# -----------------------------------------------------------------------------
# Figure 1e replicate concordance as a separate editable panel
# -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(4.2, 4.0))
xr = e['U Test FC BR1'].to_numpy(float)
yr = e['U Test FC BR2'].to_numpy(float)
ax.scatter(xr, yr, s=20, c='#9CA3A8', alpha=0.55, edgecolors='white', linewidths=0.25)
rr, rp = pearsonr(xr, yr)
cc = np.polyfit(xr, yr, 1)
xi = np.linspace(np.nanmin(xr), np.nanmax(xr), 100)
ax.plot(xi, cc[0]*xi + cc[1], color=DARK, linewidth=1.25)
ax.axhline(0, linestyle='--', linewidth=0.75, color='#AEB3B7')
ax.axvline(0, linestyle='--', linewidth=0.75, color='#AEB3B7')
ax.set_title('Validation replicate concordance', loc='left', fontsize=13.0, fontweight='bold', color='#171C21', pad=18)
ax.text(0.04, 0.95, f'Pearson r = {rr:.2f}', transform=ax.transAxes, fontsize=8.5, va='top')
ax.set_xlabel('Biological replicate 1 log2FC', fontsize=9.8, fontweight='bold')
ax.set_ylabel('Biological replicate 2 log2FC', fontsize=9.8, fontweight='bold')
ax.tick_params(labelsize=8.2)
for s in ['top','right']:
    ax.spines[s].set_visible(False)
fig.subplots_adjust(left=0.18, right=0.96, top=0.84, bottom=0.17)
save_all(fig, 'Fig1e_validation_replicates_text_editable', 4.2, 4.0)


# -----------------------------------------------------------------------------
# QC: confirm SVG contains text nodes
# -----------------------------------------------------------------------------
qc_rows = []
for svg in sorted(PAN.glob('*.svg')):
    txt = svg.read_text(errors='ignore')
    qc_rows.append({
        'file': svg.name,
        'text_element_count': len(re.findall(r'<text\b', txt)),
        'path_element_count': len(re.findall(r'<path\b', txt)),
        'contains_font_family': 'font-family' in txt or 'font:' in txt,
        'editable_text_qc_pass': len(re.findall(r'<text\b', txt)) > 0,
    })
qc = pd.DataFrame(qc_rows)
qc.to_csv(OUT / 'svg_text_qc.tsv', sep='\t', index=False)

# README
(OUT / 'README.md').write_text(
"""# Figure 1 panels 1c, 1d, and 1e — Illustrator-editable text exports\n\nThe SVG files preserve labels as SVG `<text>` elements (`matplotlib.rcParams['svg.fonttype'] = 'none'`). The PDF files use embedded TrueType fonts (`pdf.fonttype = 42`).\n\nOpen the SVG files first in Adobe Illustrator for the most reliable editable-text behavior. Font: Liberation Sans.\n\nIncluded panels:\n- 1c: non-Hippo known regulators across MaGECK, U-test, and validation.\n- 1d: functional categories and validation.\n- 1e main: genome-wide versus validation effects.\n- 1e replicates: biological replicate concordance.\n"""
)

shutil.copyfile(__file__, SCR / 'create_text_editable_panels_1cde_fresh.py')
zip_path = ROOT / 'Figure1_text_editable_1cde_fresh_package.zip'
with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
    for p in OUT.rglob('*'):
        if p.is_file():
            z.write(p, p.relative_to(OUT.parent))

print('Created:', zip_path)
print(qc.to_string(index=False))
