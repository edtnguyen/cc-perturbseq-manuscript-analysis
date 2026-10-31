import os, zipfile, shutil, math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats import hypergeom

BASE = '/mnt/data/outputs_program_first_advanced_full'
OUT = '/mnt/data/hippo_identification_remade'
os.makedirs(OUT, exist_ok=True)

# ------------------------------
# Load source data
# ------------------------------
conv = pd.read_csv(os.path.join(BASE, 'target_set_convergence_by_program.csv'))
effects = pd.read_csv(os.path.join(BASE, 'target_by_curated_program_effect_matrix.csv')).set_index('target')
fdr_mat = pd.read_csv(os.path.join(BASE, 'target_by_curated_program_fdr_matrix.csv')).set_index('target')
program_meta = pd.read_csv(os.path.join(BASE, 'curated_program_effects_aggregated.csv'))[['program_label','program_id','day','program_group','expert_annotation']].drop_duplicates()
program_order = ['P7','P22','P29','P28','P13','P30','P48','P33','P15','P14','P17','P34','P27','P2','P3','P8','P11']
program_order = [p for p in program_order if p in effects.columns]

TARGET_SETS = {
    'Hippo/YAP': ['NF2','LATS2','SAV1','TAOK1','AMOTL2','WWC1','WWC2','MST1','STK4','STK3'],
    'BMP/TGF-SMAD': ['SMAD1','SMAD4','SMAD5','LEMD3','BMPR1A','ACVR2A','ACVR1B','ACVR2B','BMPR2','TGFBR1','TGFBR2'],
    'EC regulators': ['ETV2','KDR','ERG','FLI1','TAL1','SOX17','PECAM1','CDH5','VWF','ESAM','NOTCH1','DLL4'],
    'Notch/vascular': ['NOTCH1','NOTCH4','DLL4','JAG1','RBPJ','HEY1','HEY2','KDR','FLT1'],
    'Primitive streak/mesoderm': ['TBXT','MIXL1','EOMES','MESP1','HAND1','HAND2','CDX2','MSX1'],
    'PKC/stress': ['PRKCD','MAPK1','MAPK3','JUN','FOS','DUSP1','DUSP6'],
}
TARGET_SET_ORDER = ['Hippo/YAP','BMP/TGF-SMAD','EC regulators','Notch/vascular','Primitive streak/mesoderm','PKC/stress']

HIPPO = ['NF2','LATS2','SAV1','TAOK1','AMOTL2']

# --------------------------------
# Utility
# --------------------------------

def savefig(fig, stem):
    png = os.path.join(OUT, f'{stem}.png')
    pdf = os.path.join(OUT, f'{stem}.pdf')
    svg = os.path.join(OUT, f'{stem}.svg')
    fig.savefig(png, dpi=350, bbox_inches='tight')
    fig.savefig(pdf, bbox_inches='tight')
    fig.savefig(svg, bbox_inches='tight')
    return png, pdf, svg

# Colors
orange = '#E68613'
blue = '#2F76B7'
purple = '#7A1FA2'
red = '#D73027'
grey = '#BDBDBD'
dark = '#222222'

# ------------------------------
# Figure 1: convergence-first identification
# ------------------------------
selected_pairs = [
    ('P15', 'increased'), ('P17', 'increased'), ('P34', 'increased'), ('P27', 'increased'),
    ('P15', 'decreased'), ('P2', 'decreased'), ('P3', 'decreased'), ('P8', 'decreased'), ('P11', 'decreased'),
]
selected_labels = [f'{p}\n{d.replace("increased","↑").replace("decreased","↓")}' for p,d in selected_pairs]

dot_rows = []
for y, target_set in enumerate(TARGET_SET_ORDER):
    for x, (p, d) in enumerate(selected_pairs):
        row = conv[(conv.target_set==target_set)&(conv.program_label==p)&(conv.direction==d)]
        if len(row):
            r = row.iloc[0].to_dict()
            dot_rows.append({**r, 'x': x, 'y': y, 'column_label': selected_labels[x]})
dot_df = pd.DataFrame(dot_rows)
dot_df['fdr_plot'] = dot_df['fdr'].clip(lower=1e-300)
dot_df['sig'] = dot_df['fdr'] < 0.05
dot_df['size'] = 25 + 95 * (dot_df['overlap_n'] / dot_df['set_size'].replace(0, np.nan)).fillna(0)
dot_df['abs_mean'] = dot_df['mean_effect_in_set'].abs()

fig, ax = plt.subplots(figsize=(10.6, 4.8))
ax.set_facecolor('white')
norm = TwoSlopeNorm(vmin=-4.0, vcenter=0, vmax=4.0)
cmap = LinearSegmentedColormap.from_list('eff', ['#2C7BB6','#F7F7F7','#D7191C'])
# background grid
for x in range(len(selected_pairs)):
    ax.axvline(x, color='#F0F0F0', lw=0.8, zorder=0)
for y in range(len(TARGET_SET_ORDER)):
    ax.axhline(y, color='#F0F0F0', lw=0.8, zorder=0)
# plot non-sig then sig with outline
for sig, lw, ec in [(False,0.25,'#AAAAAA'),(True,1.1,'black')]:
    sub=dot_df[dot_df['sig']==sig]
    ax.scatter(sub['x'], sub['y'], s=sub['size'], c=sub['mean_effect_in_set'], cmap=cmap, norm=norm,
               edgecolors=ec, linewidths=lw, alpha=0.95, zorder=3 if sig else 2)
# mark hippo row
hippo_y = TARGET_SET_ORDER.index('Hippo/YAP')
ax.add_patch(Rectangle((-0.5, hippo_y-0.47), len(selected_pairs), 0.94, fill=False, lw=1.4, ec=purple, zorder=1))
# annotations for key Hippo convergences
for p,d,label in [('P34','increased','5/5'),('P2','decreased','5/5'),('P8','decreased','4/5'),('P17','increased','3/5')]:
    r=dot_df[(dot_df.target_set=='Hippo/YAP')&(dot_df.program_label==p)&(dot_df.direction==d)]
    if len(r):
        ax.text(r.iloc[0]['x'], r.iloc[0]['y']+0.28, label, ha='center', va='center', fontsize=8.5, color=purple, fontweight='bold')
ax.set_xticks(range(len(selected_pairs)))
ax.set_xticklabels(selected_labels, fontsize=10, fontweight='bold')
ax.set_yticks(range(len(TARGET_SET_ORDER)))
ax.set_yticklabels(TARGET_SET_ORDER, fontsize=10.5)
ax.set_xlim(-0.6, len(selected_pairs)-0.4)
ax.set_ylim(len(TARGET_SET_ORDER)-0.6, -0.6)
ax.set_title('Convergence analysis identifies Hippo/YAP perturbations\non P17/P34 increase and EC endpoint loss', fontsize=14, fontweight='bold', pad=12)
ax.set_xlabel('Curated program effect direction', fontsize=11.5, fontweight='bold')
for spine in ax.spines.values(): spine.set_visible(False)
# colorbar and legend
cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.035, pad=0.02)
cbar.set_label('Mean effect in regulator set', fontsize=10)
# size legend manually
legend_ax = fig.add_axes([0.78, 0.12, 0.18, 0.18]); legend_ax.axis('off')
legend_ax.text(0.0, 1.0, 'Dot size = fraction\nof set affected', fontsize=9.2, ha='left', va='top')
for i,(frac,lab) in enumerate([(0.2,'20%'),(0.6,'60%'),(1.0,'100%')]):
    legend_ax.scatter(0.10, 0.62-0.22*i, s=25+95*frac, c='white', edgecolors='black', lw=0.8)
    legend_ax.text(0.25, 0.62-0.22*i, lab, va='center', fontsize=8.5)
legend_ax.text(0.0, -0.12, 'Black outline = FDR < 0.05', fontsize=8.5, ha='left')
fig.tight_layout()
savefig(fig, 'hippo_identification_convergence_first')
plt.close(fig)

dot_df.to_csv(os.path.join(OUT, 'source_convergence_first_dotplot.csv'), index=False)

# ------------------------------
# Module analysis: target modules from curated program effects
# ------------------------------
# Cluster targets by curated program effects.
mat = effects[program_order].fillna(0)
X = mat.values
Xz = (X - X.mean(0)) / X.std(0)
Z = linkage(Xz, method='ward')
modules = fcluster(Z, 4, criterion='maxclust')
module_df = pd.DataFrame({'target': mat.index, 'module': modules})
# Name modules by centroid signatures.
centroids = mat.groupby(modules).mean()
# stable order by manually sorting modules: BMP module, Hippo-diversion, mild EC, broad/background etc.
# We'll keep numeric but append simple labels based on centroids.
module_labels = {}
for m in sorted(set(modules)):
    c = centroids.loc[m]
    if c.get('P48',0)>1 and c.get('P14',0)<-1 and c.get('P2',0)<-1:
        name = 'BMP/SMAD\nearly-patterning'
    elif c.get('P34',0)>1.0 and c.get('P17',0)>0.8 and c.get('P2',0)<-1:
        name = 'P17/P34\nEC-loss module'
    elif c.get('P27',0)>0.3 and c.get('P2',0)<-0.2:
        name = 'D3 alt/EC-loss\nmodule'
    else:
        name = 'weak/mixed\nprogram effects'
    module_labels[m] = f'M{m}\n{ name }\n(n={int((modules==m).sum())})'

# Module enrichment hypergeometric
N = len(mat.index)
module_enrich_rows=[]
for m in sorted(set(modules)):
    genes=set(module_df[module_df['module']==m]['target'])
    M=len(genes)
    for set_name,set_genes in TARGET_SETS.items():
        bg=set([g for g in set_genes if g in mat.index])
        K=len(bg)
        x=len(genes & bg)
        p=hypergeom.sf(x-1,N,K,M) if K>0 else 1
        module_enrich_rows.append({'module':m,'module_label':module_labels[m], 'target_set':set_name, 'overlap_n':x, 'set_size':K, 'module_size':M, 'p_value':p, 'genes':';'.join(sorted(genes&bg))})
module_enrich = pd.DataFrame(module_enrich_rows)
# BH across all module x set tests
order=np.argsort(module_enrich['p_value'].values)
pvals=module_enrich['p_value'].values
q=np.empty_like(pvals, dtype=float)
# standard BH monotonic
rank=np.arange(1,len(pvals)+1)
q_temp=pvals[order]*len(pvals)/rank
q_temp=np.minimum.accumulate(q_temp[::-1])[::-1]
q[order]=np.minimum(q_temp,1)
module_enrich['fdr']=q
module_enrich['neglog10_fdr']=-np.log10(module_enrich['fdr'].clip(lower=1e-300))
module_enrich.to_csv(os.path.join(OUT, 'target_module_enrichment_source.csv'), index=False)
module_df.to_csv(os.path.join(OUT, 'target_module_assignments_source.csv'), index=False)
centroids.rename(index=module_labels).to_csv(os.path.join(OUT, 'target_module_centroids_source.csv'))

# Figure 2: module first, two panels
fig, axes = plt.subplots(1,2,figsize=(13,5.4), gridspec_kw={'width_ratios':[1.45,1.0]})
ax=axes[0]
cent = centroids.loc[sorted(centroids.index), program_order]
y_labels=[module_labels[m] for m in cent.index]
im=ax.imshow(cent.values, aspect='auto', cmap=cmap, norm=TwoSlopeNorm(vmin=-4, vcenter=0, vmax=4))
ax.set_xticks(range(len(program_order))); ax.set_xticklabels(program_order, rotation=45, ha='right', fontsize=9, fontweight='bold')
ax.set_yticks(range(len(y_labels))); ax.set_yticklabels(y_labels, fontsize=9.5)
ax.set_title('Target modules from curated program effects', fontsize=13, fontweight='bold')
# group day lines
for x in [2.5, 7.5, 10.5, 12.5]: ax.axvline(x, color='white', lw=2)
for y in np.arange(-0.5, len(cent.index),1): ax.axhline(y, color='white', lw=0.6)
# mark hippo module row
for idx,m in enumerate(cent.index):
    if 'P17/P34' in module_labels[m]:
        ax.add_patch(Rectangle((-0.5,idx-0.5),len(program_order),1,fill=False,edgecolor=purple,lw=2.2))
# text values sparse
for i in range(cent.shape[0]):
    for j in range(cent.shape[1]):
        v=cent.values[i,j]
        if abs(v)>=1.0:
            ax.text(j,i,f'{v:.1f}', ha='center', va='center', fontsize=7.2, color='white' if abs(v)>2.0 else 'black')
cbar=fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
cbar.set_label('Mean target effect in module', fontsize=9)

ax=axes[1]
# dotplot module enrichments: rows target sets, cols modules
mod_order=sorted(set(modules))
mod_x={m:i for i,m in enumerate(mod_order)}
set_y={s:i for i,s in enumerate(TARGET_SET_ORDER)}
plot_en=module_enrich.copy()
plot_en['x']=plot_en['module'].map(mod_x); plot_en['y']=plot_en['target_set'].map(set_y)
plot_en['size']=25+110*(plot_en['overlap_n']/plot_en['set_size'].replace(0,np.nan)).fillna(0)
plot_en['sig']=plot_en['fdr']<0.05
for sig, lw, ec in [(False,0.2,'#AAAAAA'),(True,1.1,'black')]:
    sub=plot_en[plot_en['sig']==sig]
    ax.scatter(sub['x'], sub['y'], s=sub['size'], c=sub['neglog10_fdr'], cmap='Purples', vmin=0, vmax=max(3,plot_en['neglog10_fdr'].max()), edgecolors=ec, linewidths=lw)
# label overlaps if significant
for _,r in plot_en[plot_en['sig']].iterrows():
    if r['overlap_n']>0:
        ax.text(r['x'], r['y']+0.28, f"{int(r['overlap_n'])}/{int(r['set_size'])}", ha='center', fontsize=8, color=dark, fontweight='bold')
ax.set_xticks(range(len(mod_order))); ax.set_xticklabels([f'M{m}' for m in mod_order], fontsize=10, fontweight='bold')
ax.set_yticks(range(len(TARGET_SET_ORDER))); ax.set_yticklabels(TARGET_SET_ORDER, fontsize=9.5)
ax.set_xlim(-0.5,len(mod_order)-0.5); ax.set_ylim(len(TARGET_SET_ORDER)-0.6,-0.6)
for x in range(len(mod_order)): ax.axvline(x,color='#F0F0F0',lw=0.8,zorder=0)
for y in range(len(TARGET_SET_ORDER)): ax.axhline(y,color='#F0F0F0',lw=0.8,zorder=0)
ax.set_title('Pathway enrichment in target modules', fontsize=13, fontweight='bold')
for spine in ax.spines.values(): spine.set_visible(False)
cbar=fig.colorbar(plt.cm.ScalarMappable(cmap='Purples', norm=plt.Normalize(0,max(3,plot_en['neglog10_fdr'].max()))), ax=ax, fraction=0.055, pad=0.03)
cbar.set_label('-log10(FDR)', fontsize=9)
fig.suptitle('Target-module analysis nominates a Hippo/YAP-enriched P17/P34 diversion module', fontsize=15, fontweight='bold', y=1.03)
fig.tight_layout()
savefig(fig, 'hippo_identification_module_first')
plt.close(fig)

# ------------------------------
# Figure 3: compact discovery story
# ------------------------------
fig = plt.figure(figsize=(12.5,7.5))
gs = fig.add_gridspec(2, 3, height_ratios=[0.95, 1.05], width_ratios=[1.1,1.25,1.25], hspace=0.35, wspace=0.35)
ax_net = fig.add_subplot(gs[0,0])
ax_dot = fig.add_subplot(gs[0,1:])
ax_heat = fig.add_subplot(gs[1,:])

# network simplified
ax_net.axis('off')
positions = {
    'P15': (0.1,0.70), 'P17': (0.1,0.25), 'P34': (0.72,0.25), 'EC': (0.72,0.70)
}
node_labels = {'P15':'P15\nvascular\nprogenitor','P17':'P17\nYAP/TAZ-LM','P34':'P34\nYAP/TAZ\nmesenchyme','EC':'D3 EC\nP2/P3/P8/P11'}
node_colors = {'P15':'#EAF3FC','P17':'#FDEAEA','P34':'#FDEAEA','EC':'#EAF3FC'}
edge_defs=[('P15','EC',0.40,orange,'+ EC programs'),('P15','P34',-0.42,blue,'- P34'),('P17','P34',0.45,orange,'+ P34'),('P17','EC',-0.18,blue,'- EC weak')]
for n,(x,y) in positions.items():
    ax_net.add_patch(FancyBboxPatch((x-0.09,y-0.08),0.18,0.16,boxstyle='round,pad=0.02',fc=node_colors[n],ec=('#2F76B7' if n in ['P15','EC'] else '#D95F4A'),lw=1.4))
    ax_net.text(x,y,node_labels[n],ha='center',va='center',fontsize=8.7,fontweight='bold')
for a,b,r,col,lab in edge_defs:
    xa,ya=positions[a]; xb,yb=positions[b]
    lw=1.4+3*abs(r)
    style='solid' if abs(r)>=0.3 else (0,(3,3))
    arr=FancyArrowPatch((xa+0.11,ya),(xb-0.11,yb),arrowstyle='-|>',mutation_scale=12,lw=lw,color=col,linestyle=style,connectionstyle='arc3,rad=0.05')
    ax_net.add_patch(arr)
    ax_net.text((xa+xb)/2, (ya+yb)/2+0.06, f'r={r:.2f}', color=col, fontsize=9, fontweight='bold', ha='center')
ax_net.set_title('Program-effect\ncoupling model', fontsize=12, fontweight='bold')
ax_net.set_xlim(0,1); ax_net.set_ylim(0,1)

# convergence dotplot zoom for Hippo row + selected target sets
sel_sets = ['Hippo/YAP','BMP/TGF-SMAD','EC regulators','PKC/stress']
zoom = dot_df[dot_df['target_set'].isin(sel_sets)].copy()
set_order=sel_sets
set_y={s:i for i,s in enumerate(set_order)}
zoom['y']=zoom['target_set'].map(set_y)
for sig,lw,ec in [(False,0.2,'#AAAAAA'),(True,1.1,'black')]:
    sub=zoom[zoom['sig']==sig]
    ax_dot.scatter(sub['x'], sub['y'], s=sub['size'], c=sub['mean_effect_in_set'], cmap=cmap, norm=norm, edgecolors=ec, linewidths=lw, alpha=0.95)
ax_dot.add_patch(Rectangle((-0.5,set_y['Hippo/YAP']-0.45),len(selected_pairs),0.9,fill=False,ec=purple,lw=1.4))
ax_dot.set_xticks(range(len(selected_pairs))); ax_dot.set_xticklabels(selected_labels, rotation=0, fontsize=8.5, fontweight='bold')
ax_dot.set_yticks(range(len(set_order))); ax_dot.set_yticklabels(set_order, fontsize=9.5)
ax_dot.set_xlim(-0.6,len(selected_pairs)-0.4); ax_dot.set_ylim(len(set_order)-0.6,-0.6)
for spine in ax_dot.spines.values(): spine.set_visible(False)
for x in range(len(selected_pairs)): ax_dot.axvline(x,color='#F0F0F0',lw=0.8,zorder=0)
for y in range(len(set_order)): ax_dot.axhline(y,color='#F0F0F0',lw=0.8,zorder=0)
ax_dot.set_title('Regulator-set convergence highlights Hippo/YAP', fontsize=12, fontweight='bold')

# heatmap Hippo effect
hippo_effects = effects.loc[[g for g in HIPPO if g in effects.index], [p for p in ['P15','P14','P17','P34','P27','P2','P3','P8','P11'] if p in effects.columns]]
hippo_fdr = fdr_mat.loc[hippo_effects.index, hippo_effects.columns]
im=ax_heat.imshow(hippo_effects.values, aspect='auto', cmap=cmap, norm=TwoSlopeNorm(vmin=-5, vcenter=0, vmax=5))
ax_heat.set_xticks(range(hippo_effects.shape[1])); ax_heat.set_xticklabels(hippo_effects.columns, fontsize=10, fontweight='bold')
ax_heat.set_yticks(range(hippo_effects.shape[0])); ax_heat.set_yticklabels(hippo_effects.index, fontsize=10, fontweight='bold')
for i in range(hippo_effects.shape[0]):
    for j in range(hippo_effects.shape[1]):
        v=hippo_effects.values[i,j]
        sig=hippo_fdr.values[i,j]<0.05
        txt=f'{v:.1f}' + ('*' if sig else '')
        ax_heat.text(j,i,txt,ha='center',va='center',fontsize=9,color='white' if abs(v)>2.7 else 'black',fontweight='bold' if sig else 'normal')
for x in np.arange(-0.5,hippo_effects.shape[1],1): ax_heat.axvline(x,color='white',lw=1)
for y in np.arange(-0.5,hippo_effects.shape[0],1): ax_heat.axhline(y,color='white',lw=1)
ax_heat.set_title('Hippo/YAP targets share a P17/P34↑ and EC-program↓ signature', fontsize=12.5, fontweight='bold')
cbar=fig.colorbar(im, ax=ax_heat, fraction=0.018, pad=0.02)
cbar.set_label('Target effect', fontsize=9)
fig.suptitle('Program-first identification of Hippo/YAP as a convergent alternative-fate module', fontsize=16, fontweight='bold', y=0.995)
fig.tight_layout(rect=[0,0,1,0.97])
savefig(fig, 'hippo_identification_compact_story')
plt.close(fig)

# README
with open(os.path.join(OUT, 'README.txt'), 'w') as f:
    f.write('Hippo identification figures remade from program-first output tables.\n')
    f.write('Figures:\n')
    f.write('1. hippo_identification_convergence_first: regulator-set convergence analysis.\n')
    f.write('2. hippo_identification_module_first: unsupervised target-module analysis + module enrichment.\n')
    f.write('3. hippo_identification_compact_story: compact discovery narrative.\n')
    f.write('SVG/PDF/PNG are included. SVG text is preserved where possible (matplotlib svg.fonttype=none).\n')

# copy script
shutil.copy(__file__, os.path.join(OUT, 'remake_hippo_identification_figures.py'))

# Package
zip_path='/mnt/data/hippo_identification_figures_options_REM ADE.zip'.replace(' ','_')
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
    for root, _, files in os.walk(OUT):
        for fn in files:
            fp=os.path.join(root, fn)
            z.write(fp, arcname=os.path.relpath(fp, OUT))
print(zip_path)
