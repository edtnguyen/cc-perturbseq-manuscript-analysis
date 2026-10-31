import os, zipfile, math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.lines import Line2D

# Make SVG text editable in Illustrator.
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'DejaVu Sans'

BASE_V4 = '/mnt/data/outputs_program_level_temporal_map_v4_unzipped'
BASE_ADV = '/mnt/data/outputs_program_first_advanced_unzipped'
OUTDIR = '/mnt/data/program_first_hippo_discovery_figures'
os.makedirs(OUTDIR, exist_ok=True)

# Inputs.
program_corr = pd.read_csv(os.path.join(BASE_V4, 'curated_program_program_pearson_correlation.csv'), index_col=0)
program_corr.index = program_corr.index.astype(str)
program_corr.columns = program_corr.columns.astype(str)
reg = pd.read_csv(os.path.join(BASE_ADV, 'D3_programs_predicted_by_D2_programs_multivariable_regression.csv'))
conv = pd.read_csv(os.path.join(BASE_ADV, 'target_set_convergence_by_program.csv'))
effects_long = pd.read_csv(os.path.join(BASE_V4, 'curated_program_effects_long.csv'))
ann = pd.read_csv(os.path.join(BASE_V4, 'curated_program_annotations_used.csv'))
effect_mat = pd.read_csv(os.path.join(BASE_V4, 'target_by_curated_program_effect_matrix_labels.csv')).set_index('target')
fdr_mat = pd.read_csv(os.path.join(BASE_V4, 'target_by_curated_program_fdr_matrix_program_ids.csv')).set_index('target')
fdr_mat.columns = ['P'+str(c) for c in fdr_mat.columns]

# Colors.
pos_col = '#E68613'  # orange
neg_col = '#2F76B7'  # blue
purple = '#7A4BA3'
grey = '#BDBDBD'
red = '#C43C39'
black = '#222222'
light_grey = '#F2F2F2'

cmap_rb = LinearSegmentedColormap.from_list('orange_white_blue', [neg_col, '#F7F7F7', pos_col])
# note: negative blue, positive orange.
norm_effect = TwoSlopeNorm(vmin=-4, vcenter=0, vmax=4)
norm_corr = TwoSlopeNorm(vmin=-0.7, vcenter=0, vmax=0.7)

# ------------------------------
# Panel A: cleaned program coupling network.
# ------------------------------
def draw_network(ax):
    ax.set_xlim(-0.15, 4.95)
    ax.set_ylim(-1.55, 1.65)
    ax.axis('off')

    nodes = {
        'P30': {'xy': (0.45, 0.92), 'label': 'P30\nprimitive streak', 'edge': '#9B1B5A', 'fill': '#FFF2F6'},
        'P13': {'xy': (0.45, -0.92), 'label': 'P13\nposterior PS', 'edge': pos_col, 'fill': '#FFF7EB'},
        'P15': {'xy': (2.12, 0.92), 'label': 'P15\nvascular progenitor', 'edge': neg_col, 'fill': '#EFF6FF'},
        'P17': {'xy': (2.12, -0.92), 'label': 'P17\nYAP/TAZ LM', 'edge': pos_col, 'fill': '#FFF7EB'},
        'EC': {'xy': (4.00, 0.92), 'label': 'D3 EC\nP2/P3/P8/P11', 'edge': neg_col, 'fill': '#EFF6FF'},
        'P34': {'xy': (4.00, -0.40), 'label': 'P34\nYAP mesenchyme', 'edge': red, 'fill': '#FFF1F1'},
        'P27': {'xy': (4.00, -1.17), 'label': 'P27\nECM mesenchyme', 'edge': '#A55A2A', 'fill': '#FFF5E9'},
    }
    # Draw day labels.
    day_labels = [('D1',0.45),('D2',2.12),('D3',4.00)]
    for d,x in day_labels:
        ax.text(x, 1.45, d, ha='center', va='bottom', fontsize=12.5, fontweight='bold', color=black)
        ax.plot([x-0.55, x+0.55], [1.38, 1.38], color='#A0A0A0', lw=0.9)

    # Draw nodes.
    for n, info in nodes.items():
        x,y = info['xy']
        w = 0.9 if n not in ['P34','P27','EC'] else 1.15
        h = 0.43
        box = FancyBboxPatch((x-w/2,y-h/2), w, h, boxstyle='round,pad=0.02,rounding_size=0.04',
                             linewidth=1.2, edgecolor=info['edge'], facecolor=info['fill'])
        ax.add_patch(box)
        ax.text(x,y,info['label'], ha='center', va='center', fontsize=7.8, color=black)

    def r(a,b):
        return float(program_corr.loc[a,b])
    def edge(a,b,label=None, curvature=0.0, dashed=False, alpha=1.0):
        x1,y1 = nodes[a]['xy']; x2,y2 = nodes[b]['xy']
        val = r(a,b) if a != 'EC' and b != 'EC' else None
        # For EC node, use P15/P17 average correlations to D3 EC programs if requested separately.
        if label is None:
            val = r(a,b)
            label = f'{val:+.2f}'
        else:
            # parse numeric label
            try:
                val = float(label.replace('+',''))
            except Exception:
                val = 0
        color = pos_col if val >= 0 else neg_col
        lw = 0.8 + 4.0*abs(val)
        style = '--' if dashed else '-'
        arr = FancyArrowPatch((x1+0.48*np.sign(x2-x1), y1), (x2-0.54*np.sign(x2-x1), y2),
                              arrowstyle='-|>', mutation_scale=10, linewidth=lw,
                              color=color, alpha=alpha, linestyle=style,
                              connectionstyle=f'arc3,rad={curvature}')
        ax.add_patch(arr)
        xm=(x1+x2)/2; ym=(y1+y2)/2
        # perpendicular offset for label.
        ax.text(xm, ym + 0.12*np.sign(curvature if curvature!=0 else (y2-y1+0.01)), label,
                ha='center', va='center', fontsize=7.8, color=color, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.08', facecolor='white', edgecolor='none', alpha=0.8))

    # Core edges from correlation table.
    edge('P30','P15', curvature=0.10)
    edge('P30','P17', curvature=-0.08)
    edge('P13','P17', curvature=-0.10)
    edge('P13','P15', curvature=0.08, dashed=True, alpha=0.8)
    # D2 internal antagonism.
    x,y1 = nodes['P15']['xy']; _,y2 = nodes['P17']['xy']
    val = r('P15','P17')
    ax.plot([x,x], [y2+0.24,y1-0.24], color=neg_col, lw=0.8+4*abs(val), alpha=0.95)
    ax.text(x+0.34, (y1+y2)/2, f'{val:+.2f}', color=neg_col, fontsize=7.8, fontweight='bold',
            ha='left', va='center', bbox=dict(boxstyle='round,pad=0.08', facecolor='white', edgecolor='none', alpha=0.8))

    # D2 to D3 EC: use mean correlations across EC programs from raw correlation table.
    ec_prog = ['P2','P3','P8','P11']
    val_P15_EC = np.mean([r('P15',p) for p in ec_prog])
    val_P17_EC = np.mean([r('P17',p) for p in ec_prog])
    # Draw EC as mean edge label.
    edge('P15','EC', label=f'{val_P15_EC:+.2f}', curvature=0.04)
    edge('P17','EC', label=f'{val_P17_EC:+.2f}', curvature=-0.04, dashed=True, alpha=0.75)
    edge('P15','P34', curvature=-0.08)
    edge('P17','P34', curvature=0.08)
    edge('P15','P27', curvature=-0.12, dashed=True, alpha=0.8)
    edge('P17','P27', curvature=0.10, dashed=True, alpha=0.8)

    # Legend.
    ax.plot([0.05,0.35],[-1.42,-1.42], color=pos_col, lw=2.6)
    ax.text(0.40,-1.42,'positive r', fontsize=7.5, va='center', ha='left')
    ax.plot([1.02,1.32],[-1.42,-1.42], color=neg_col, lw=2.6)
    ax.text(1.37,-1.42,'negative r', fontsize=7.5, va='center', ha='left')
    ax.text(2.15,-1.42,'line width ~ |r|', fontsize=7.5, va='center', ha='left', color='#555555')

# ------------------------------
# Panel B: D2 -> D3 multivariable model coefficients.
# ------------------------------
def plot_multivariable(ax):
    coef_cols = ['beta_P15','beta_P14','beta_P17']
    row_order = ['P34','P27','P2','P3','P8','P11']
    df = reg.set_index('D3_program').loc[row_order]
    mat = df[coef_cols].values
    im = ax.imshow(mat, cmap=cmap_rb, norm=TwoSlopeNorm(vmin=-0.5, vcenter=0, vmax=0.5), aspect='auto')
    ax.set_xticks(range(3)); ax.set_xticklabels(['P15\nvascular', 'P14\nLM', 'P17\nYAP LM'], fontsize=7)
    ax.set_yticks(range(len(row_order))); ax.set_yticklabels(row_order, fontsize=8)
    ax.tick_params(length=0)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val=mat[i,j]
            ax.text(j,i,f'{val:+.2f}', ha='center', va='center', fontsize=7.5,
                    color='white' if abs(val)>0.32 else black, fontweight='bold')
    ax.set_title('D2 predictors of D3 programs', fontsize=10.5, fontweight='bold', pad=8)
    # Add adjusted R² as text in a mini column on right.
    for i, prog in enumerate(row_order):
        ax.text(3.05, i, f'adj. R²={df.loc[prog,"adj_r2"]:.2f}', fontsize=7.2, va='center', ha='left', color='#555555')
    ax.set_xlim(-0.5, 4.2)
    for spine in ax.spines.values(): spine.set_visible(False)
    return im

# ------------------------------
# Panel C: program convergence dot plot.
# ------------------------------
def plot_convergence(ax):
    target_sets = ['Hippo/YAP','BMP/TGF-SMAD','EC regulators','Notch/vascular','PKC/stress']
    programs = ['P48','P15','P14','P17','P34','P27','P2','P3','P8','P11']
    # Choose best direction per target_set/program by min fdr.
    rows=[]
    for ts in target_sets:
        for prog in programs:
            sub=conv[(conv.target_set==ts)&(conv.program_label==prog)].copy()
            if sub.empty: continue
            sub=sub.sort_values('fdr').iloc[0]
            rows.append(sub.to_dict())
    df=pd.DataFrame(rows)
    # only plot if fdr < 0.2 or key Hippo rows; faint dots otherwise with zero-ish size.
    xmap={p:i for i,p in enumerate(programs)}
    ymap={ts:i for i,ts in enumerate(target_sets[::-1])}
    for _,row in df.iterrows():
        x=xmap[row['program_label']]; y=ymap[row['target_set']]
        val=row['mean_effect_in_set']
        fdr=max(float(row['fdr']), 1e-300)
        size=20 + 45*min(-np.log10(fdr), 5)
        alpha=0.95 if fdr<0.1 else 0.15
        ax.scatter(x,y,s=size,c=[val],cmap=cmap_rb,norm=norm_effect,edgecolor='white',linewidth=0.45,alpha=alpha)
        if fdr<0.05:
            ax.text(x,y,str(int(row['overlap_n'])), ha='center', va='center', fontsize=6.8, color='white' if abs(val)>1.3 else black, fontweight='bold')
    ax.set_xticks(range(len(programs))); ax.set_xticklabels(programs, rotation=45, ha='right', fontsize=7.4)
    ax.set_yticks(range(len(target_sets))); ax.set_yticklabels(target_sets[::-1], fontsize=8)
    ax.set_title('Target-set convergence on curated programs', fontsize=10.5, fontweight='bold', pad=8)
    ax.set_xlim(-0.6,len(programs)-0.4); ax.set_ylim(-0.6,len(target_sets)-0.4)
    ax.grid(True, color='#E6E6E6', linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values(): spine.set_visible(False)

# ------------------------------
# Panel D: Hippo-focused program heatmap.
# ------------------------------
def plot_hippo_heatmap(ax):
    targets = ['NF2','LATS2','SAV1','TAOK1','AMOTL2']
    programs = ['P15','P14','P17','P34','P27','P2','P3','P8','P11']
    mat = effect_mat.loc[targets, programs].values
    im = ax.imshow(mat, cmap=cmap_rb, norm=norm_effect, aspect='auto')
    ax.set_xticks(range(len(programs))); ax.set_xticklabels(programs, fontsize=7.2, rotation=45, ha='right')
    ax.set_yticks(range(len(targets))); ax.set_yticklabels(targets, fontsize=8)
    ax.set_title('Hippo/YAP targets: raw program effects', fontsize=10.5, fontweight='bold', pad=8)
    for i,t in enumerate(targets):
        for j,p in enumerate(programs):
            val = mat[i,j]
            fdr = float(fdr_mat.loc[t,p]) if t in fdr_mat.index and p in fdr_mat.columns else 1.0
            label = f'{val:.1f}' if abs(val)>=1 else ''
            if fdr<0.05:
                ax.text(j,i,label if label else '•', ha='center', va='center', fontsize=7.5,
                        color='white' if abs(val)>2.0 else black, fontweight='bold')
    ax.tick_params(length=0)
    for spine in ax.spines.values(): spine.set_visible(False)
    return im

# Create individual panels and composite.
def save_individual():
    panels=[]
    fig, ax=plt.subplots(figsize=(7.2,3.6), dpi=350)
    draw_network(ax); fig.tight_layout();
    for ext in ['png','pdf','svg']:
        fig.savefig(os.path.join(OUTDIR, f'panelA_clean_lagged_program_coupling_network.{ext}'), bbox_inches='tight', dpi=350 if ext=='png' else None)
    plt.close(fig)
    panels.append('panelA_clean_lagged_program_coupling_network')

    fig, ax=plt.subplots(figsize=(5.2,3.7), dpi=350)
    im=plot_multivariable(ax); cbar=fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02); cbar.set_label('coefficient', fontsize=8); cbar.ax.tick_params(labelsize=7)
    fig.tight_layout();
    for ext in ['png','pdf','svg']:
        fig.savefig(os.path.join(OUTDIR, f'panelB_D2_to_D3_multivariable_coefficients.{ext}'), bbox_inches='tight', dpi=350 if ext=='png' else None)
    plt.close(fig)
    panels.append('panelB_D2_to_D3_multivariable_coefficients')

    fig, ax=plt.subplots(figsize=(6.3,3.6), dpi=350)
    plot_convergence(ax)
    # fake size legend and colorbar
    sm = mpl.cm.ScalarMappable(cmap=cmap_rb, norm=norm_effect); sm.set_array([])
    cbar=fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02); cbar.set_label('mean effect', fontsize=8); cbar.ax.tick_params(labelsize=7)
    fig.tight_layout();
    for ext in ['png','pdf','svg']:
        fig.savefig(os.path.join(OUTDIR, f'panelC_target_set_convergence_dotplot_simplified.{ext}'), bbox_inches='tight', dpi=350 if ext=='png' else None)
    plt.close(fig)
    panels.append('panelC_target_set_convergence_dotplot_simplified')

    fig, ax=plt.subplots(figsize=(5.2,2.9), dpi=350)
    im=plot_hippo_heatmap(ax); cbar=fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02); cbar.set_label('effect', fontsize=8); cbar.ax.tick_params(labelsize=7)
    fig.tight_layout();
    for ext in ['png','pdf','svg']:
        fig.savefig(os.path.join(OUTDIR, f'panelD_hippo_curated_program_heatmap_clean.{ext}'), bbox_inches='tight', dpi=350 if ext=='png' else None)
    plt.close(fig)
    panels.append('panelD_hippo_curated_program_heatmap_clean')
    return panels

# Composite figure.
fig = plt.figure(figsize=(13.5, 8.2), dpi=350)
gs = fig.add_gridspec(2, 2, width_ratios=[1.25,1.0], height_ratios=[1.0,1.0], wspace=0.25, hspace=0.35)
axA=fig.add_subplot(gs[0,0]); draw_network(axA)
axB=fig.add_subplot(gs[0,1]); imB=plot_multivariable(axB); cbarB=fig.colorbar(imB, ax=axB, fraction=0.035, pad=0.015); cbarB.ax.tick_params(labelsize=7); cbarB.set_label('coefficient', fontsize=8)
axC=fig.add_subplot(gs[1,0]); plot_convergence(axC); sm=mpl.cm.ScalarMappable(cmap=cmap_rb,norm=norm_effect); sm.set_array([]); cbarC=fig.colorbar(sm, ax=axC, fraction=0.03, pad=0.015); cbarC.ax.tick_params(labelsize=7); cbarC.set_label('mean effect', fontsize=8)
axD=fig.add_subplot(gs[1,1]); imD=plot_hippo_heatmap(axD); cbarD=fig.colorbar(imD, ax=axD, fraction=0.035, pad=0.015); cbarD.ax.tick_params(labelsize=7); cbarD.set_label('effect', fontsize=8)
# Panel letters
for letter, ax in zip(['A','B','C','D'], [axA, axB, axC, axD]):
    ax.text(-0.08, 1.08, letter, transform=ax.transAxes, fontsize=14, fontweight='bold', va='top', ha='left')
fig.suptitle('Program-first discovery of Hippo/YAP-associated diversion', fontsize=15.5, fontweight='bold', y=0.99)
fig.tight_layout(rect=[0,0,1,0.97])
for ext in ['png','pdf','svg']:
    fig.savefig(os.path.join(OUTDIR, f'program_first_hippo_discovery_composite.{ext}'), bbox_inches='tight', dpi=350 if ext=='png' else None)
plt.close(fig)

# Save individual panels.
save_individual()

# Create source summary CSVs used.
# Network edge table.
edge_rows=[]
def getr(a,b): return float(program_corr.loc[a,b])
ec_prog=['P2','P3','P8','P11']
for a,b in [('P30','P15'),('P30','P17'),('P13','P17'),('P13','P15'),('P15','P17'),('P15','P34'),('P17','P34'),('P15','P27'),('P17','P27')]:
    edge_rows.append({'from':a,'to':b,'pearson_r':getr(a,b)})
edge_rows.append({'from':'P15','to':'EC_mean','pearson_r':np.mean([getr('P15',p) for p in ec_prog])})
edge_rows.append({'from':'P17','to':'EC_mean','pearson_r':np.mean([getr('P17',p) for p in ec_prog])})
pd.DataFrame(edge_rows).to_csv(os.path.join(OUTDIR,'program_first_hippo_discovery_network_edges.csv'),index=False)
reg.to_csv(os.path.join(OUTDIR,'D2_to_D3_multivariable_regression_used.csv'),index=False)
conv.to_csv(os.path.join(OUTDIR,'target_set_convergence_used.csv'),index=False)
effect_mat.loc[['NF2','LATS2','SAV1','TAOK1','AMOTL2'], ['P15','P14','P17','P34','P27','P2','P3','P8','P11']].to_csv(os.path.join(OUTDIR,'hippo_heatmap_effects_used.csv'))
fdr_mat.loc[['NF2','LATS2','SAV1','TAOK1','AMOTL2'], ['P15','P14','P17','P34','P27','P2','P3','P8','P11']].to_csv(os.path.join(OUTDIR,'hippo_heatmap_fdr_used.csv'))

# Copy script to output.
import shutil
shutil.copyfile(__file__, os.path.join(OUTDIR,'make_program_first_hippo_discovery_figures.py'))

# Zip package.
zip_path='/mnt/data/program_first_hippo_discovery_figures_package.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for root,dirs,files in os.walk(OUTDIR):
        for fn in files:
            fp=os.path.join(root,fn)
            z.write(fp, arcname=os.path.basename(fp))
print('Created outputs in', OUTDIR)
print('ZIP', zip_path)
