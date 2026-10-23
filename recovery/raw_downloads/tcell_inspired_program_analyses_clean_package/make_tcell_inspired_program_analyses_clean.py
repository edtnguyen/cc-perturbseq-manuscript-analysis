#!/usr/bin/env python3
from __future__ import annotations
import os, zipfile, shutil, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr, mannwhitneyu
warnings.filterwarnings('ignore')

mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.family'] = 'DejaVu Sans'

BASE = Path('/mnt/data')
SRC = BASE / '_analysis_sources'
ADV = SRC / 'adv'
V4 = SRC / 'v4'
PHENO = BASE / 'program_temporal_phenotype_outputs'
OUT = BASE / 'tcell_inspired_program_analyses_clean'
OUT.mkdir(exist_ok=True, parents=True)

# If source CSVs not available, extract zips.
def unzip_if_needed(zip_path, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    if not any(outdir.glob('*.csv')) and Path(zip_path).exists():
        with zipfile.ZipFile(zip_path) as z:
            z.extractall(outdir)

unzip_if_needed(BASE/'outputs_program_first_advanced.zip', ADV)
unzip_if_needed(BASE/'outputs_program_level_temporal_map_v4.zip', V4)
if not (PHENO/'program_defined_temporal_phenotype_features.csv').exists() and (BASE/'program_temporal_phenotype_outputs_package.zip').exists():
    with zipfile.ZipFile(BASE/'program_temporal_phenotype_outputs_package.zip') as z:
        z.extractall(PHENO)

long = pd.read_csv(ADV/'curated_program_effects_aggregated.csv')
mat = pd.read_csv(ADV/'target_by_curated_program_effect_matrix.csv', index_col=0)
fdr = pd.read_csv(ADV/'target_by_curated_program_fdr_matrix.csv', index_col=0)
ann = pd.read_csv(V4/'curated_program_annotations_used.csv')
features = pd.read_csv(PHENO/'program_defined_temporal_phenotype_features.csv')

ann['program_label'] = ann['program_label'].astype(str)
program_order = [p for p in ann['program_label'].tolist() if p in mat.columns]
mat = mat.loc[:, program_order]
fdr = fdr.loc[:, program_order]
ann_map = ann.set_index('program_label')
program_display = {p: f"{p}\n{ann_map.loc[p, 'Name']}" for p in program_order}
program_day = {p: str(ann_map.loc[p, 'curated_day']) for p in program_order}

DAYS = ['D0','D1','D2','D3']
POS = '#E68613'
NEG = '#2F76B7'
PURPLE = '#7A4BA3'
GREY = '#BDBDBD'
RED = '#C43C39'
ECBLUE = '#2166AC'

CORE_SETS = {
    'Hippo/YAP': ['NF2','LATS2','SAV1','TAOK1','AMOTL2'],
    'BMP/TGF-SMAD': ['SMAD1','SMAD4','SMAD5','LEMD3','BMPR1A','ACVR2A','ACVR1B','ACVR2B'],
    'EC regulators': ['KDR','ETV2','ERG','FLI1','TAL1','SOX17','CDH5','PECAM1'],
    'Example': ['PRKCD'],
}
CORE_COLORS = {'Hippo/YAP': PURPLE, 'BMP/TGF-SMAD': POS, 'EC regulators': ECBLUE, 'Example': '#555555'}

# Utility
def savefig(fig, stem):
    for ext in ['png','pdf','svg']:
        fig.savefig(OUT/f'{stem}.{ext}', dpi=350 if ext=='png' else None, bbox_inches='tight')
    plt.close(fig)

def plot_heatmap(ax, arr, cmap='RdBu_r', vmin=None, vmax=None, cbar_label='', xticks=None, yticks=None, xticklabels=None, yticklabels=None, fontsize=7):
    im=ax.imshow(arr, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    if xticks is not None:
        ax.set_xticks(xticks); ax.set_xticklabels(xticklabels, rotation=45, ha='right', fontsize=fontsize)
    if yticks is not None:
        ax.set_yticks(yticks); ax.set_yticklabels(yticklabels, fontsize=fontsize)
    cb=plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    if cbar_label:
        cb.set_label(cbar_label, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    return im

# =====================================================================
# A. Top positive/negative regulators per program + volcano grid
# =====================================================================
selected_programs = [p for p in ['P15','P17','P34','P27','P3','P8','P48'] if p in mat.columns]
top_rows=[]
for p in selected_programs:
    tmp = pd.DataFrame({'target': mat.index, 'effect': mat[p].values, 'fdr': fdr[p].values})
    pos = tmp.sort_values('effect', ascending=False).head(10).assign(program=p, direction='positive')
    neg = tmp.sort_values('effect', ascending=True).head(10).assign(program=p, direction='negative')
    top_rows.append(pd.concat([pos, neg]))
top_reg = pd.concat(top_rows, ignore_index=True)
top_reg.to_csv(OUT/'top_positive_negative_regulators_per_selected_program.csv', index=False)

fig, axes = plt.subplots(len(selected_programs), 1, figsize=(8.2, 1.25*len(selected_programs)), sharex=False)
if len(selected_programs)==1: axes=[axes]
for ax,p in zip(axes, selected_programs):
    t=top_reg[top_reg['program']==p].copy().sort_values('effect')
    colors=[NEG if v<0 else POS for v in t['effect']]
    ax.barh(t['target'], t['effect'], color=colors, alpha=0.9)
    ax.axvline(0, color='black', lw=0.8)
    ax.set_title(f'{p}: {ann_map.loc[p,"Name"]}', loc='left', fontsize=9, weight='bold')
    ax.tick_params(axis='y', labelsize=6)
    ax.tick_params(axis='x', labelsize=7)
    ax.spines[['top','right']].set_visible(False)
axes[-1].set_xlabel('program effect size')
fig.suptitle('Top positive and negative regulators of selected programs', fontsize=12, weight='bold', y=1.01)
fig.tight_layout()
savefig(fig,'fig2_style_top_regulators_per_program')

# Volcano grid for same selected programs
fig, axes = plt.subplots(2, int(np.ceil(len(selected_programs)/2)), figsize=(12, 6.2))
axes = axes.flatten()
for ax,p in zip(axes, selected_programs):
    df = pd.DataFrame({'target':mat.index, 'effect':mat[p].values, 'fdr':fdr[p].values})
    df['neglog10fdr'] = -np.log10(df['fdr'].clip(lower=1e-300))
    sig = df['fdr'] < 0.05
    ax.scatter(df.loc[~sig,'effect'], df.loc[~sig,'neglog10fdr'], s=8, color=GREY, alpha=0.4, edgecolors='none')
    ax.scatter(df.loc[sig,'effect'], df.loc[sig,'neglog10fdr'], s=14, color=POS, alpha=0.8, edgecolors='none')
    lab = pd.concat([df.sort_values('effect').head(4), df.sort_values('effect',ascending=False).head(4)]).drop_duplicates('target')
    for _,r in lab.iterrows():
        ax.text(r['effect'], r['neglog10fdr'], r['target'], fontsize=6, ha='center')
    ax.axvline(0, color='#888888', ls='--', lw=0.7)
    ax.axhline(-np.log10(0.05), color='#888888', ls='--', lw=0.7)
    ax.set_title(f'{p} {ann_map.loc[p,"Name"]}', fontsize=8, weight='bold')
    ax.tick_params(labelsize=7)
    ax.spines[['top','right']].set_visible(False)
for ax in axes[len(selected_programs):]: ax.axis('off')
fig.supxlabel('target effect on program', fontsize=10)
fig.supylabel('-log10(FDR)', fontsize=10)
fig.suptitle('Program-specific regulator volcano plots', fontsize=12, weight='bold')
fig.tight_layout()
savefig(fig,'fig2_style_selected_program_volcano_grid')

# =====================================================================
# B. Target-day modules / temporal specificity (Fig 3 style)
# =====================================================================
long_sig = long.copy()
long_sig['sig_effect'] = np.where(long_sig['fdr'] < 0.05, long_sig['effect'], 0.0)
rows=[]
for (target, day), sub in long_sig.groupby(['target','day']):
    row={'target':target,'day':day,'target_day':f'{target}_{day}'}
    row.update({p:0.0 for p in program_order})
    for _,r in sub.iterrows():
        p=r['program_label']
        if p in program_order: row[p]=r['sig_effect']
    rows.append(row)
td = pd.DataFrame(rows).set_index('target_day')
td_features=td[program_order]
td_nonzero=td_features.loc[td_features.abs().sum(axis=1)>0]
# standardize columns; cosine cluster
X=td_nonzero.values
sd=X.std(axis=0); sd[sd==0]=1
Xz=(X-X.mean(axis=0))/sd
Z=linkage(pdist(Xz, metric='cosine'), method='average')
order=dendrogram(Z, no_plot=True)['leaves']
ordered=td_nonzero.index[order]
mods=fcluster(Z, t=10, criterion='maxclust')
mod_ser=pd.Series(mods,index=td_nonzero.index)
# reorder module labels by dendrogram appearance
seen=[]
for idx in ordered:
    m=mod_ser.loc[idx]
    if m not in seen: seen.append(m)
mod_map={m:i+1 for i,m in enumerate(seen)}
mod_ser=mod_ser.map(mod_map)
assign=td.loc[td_nonzero.index, ['target','day']].copy()
assign['target_day']=assign.index; assign['module']=mod_ser
assign.to_csv(OUT/'target_day_unsupervised_module_assignments.csv', index=False)

# module centroids
centroids=[]
for m in sorted(mod_ser.unique()):
    members=mod_ser[mod_ser==m].index
    centroids.append(td_features.loc[members].mean(axis=0).rename(f'M{m}\n(n={len(members)})'))
centroid_df=pd.DataFrame(centroids)
centroid_df.to_csv(OUT/'target_day_module_centroid_program_effects.csv')
fig,ax=plt.subplots(figsize=(8.8,4.2))
lim=np.nanmax(np.abs(centroid_df.values))
plot_heatmap(ax, centroid_df.values, vmin=-lim, vmax=lim, cbar_label='mean significant effect', xticks=range(len(program_order)), xticklabels=[program_display[p] for p in program_order], yticks=range(len(centroid_df)), yticklabels=centroid_df.index, fontsize=7)
ax.set_title('Unsupervised target-day modules: program-effect centroids', fontsize=12, weight='bold')
savefig(fig,'target_day_module_centroid_heatmap')

spec=pd.crosstab(assign['module'], assign['day']).reindex(columns=DAYS, fill_value=0)
spec_frac=spec.div(spec.sum(axis=1), axis=0)
spec_frac.to_csv(OUT/'target_day_module_temporal_specificity_fraction.csv')
fig,ax=plt.subplots(figsize=(4.8,3.6))
im=ax.imshow(spec_frac.values, aspect='auto', cmap='Blues', vmin=0, vmax=1)
for i in range(spec.shape[0]):
    for j in range(spec.shape[1]):
        ax.text(j,i,str(spec.iloc[i,j]),ha='center',va='center',fontsize=8,color='white' if spec_frac.iloc[i,j]>0.55 else 'black')
ax.set_xticks(range(len(DAYS))); ax.set_xticklabels(DAYS, fontsize=9)
ax.set_yticks(range(spec.shape[0])); ax.set_yticklabels([f'M{i}' for i in spec.index], fontsize=9)
ax.set_title('Temporal specificity of target-day modules', fontsize=12, weight='bold')
cb=plt.colorbar(im, ax=ax, fraction=0.05, pad=0.02); cb.set_label('fraction of module', fontsize=8); cb.ax.tick_params(labelsize=7)
savefig(fig,'target_day_module_temporal_specificity')

# =====================================================================
# C. Program-space EC-vs-alternative signature regulator model (Fig 4 style)
# =====================================================================
weights=pd.Series(0.0,index=program_order)
for p in ['P2','P3','P8','P11']:
    if p in weights.index: weights[p]=1.0
if 'P15' in weights.index: weights['P15']=0.7
if 'P14' in weights.index: weights['P14']=0.2
if 'P17' in weights.index: weights['P17']=-0.7
for p in ['P34','P27']:
    if p in weights.index: weights[p]=-1.0
weights.to_csv(OUT/'program_space_EC_vs_alternative_signature_weights.csv')
Mz=(mat-mat.mean(axis=0))/mat.std(axis=0).replace(0,1)
w=weights.values; w=w/np.sqrt(np.sum(w*w))
sig_score=pd.Series(Mz.values@w,index=Mz.index,name='EC_vs_alt_signature_similarity')
rank_df=pd.DataFrame({'target':sig_score.index,'EC_vs_alt_signature_similarity':sig_score.values})
for col in ['P15','P17','P34','P27']:
    if col in mat.columns: rank_df[col]=mat[col].values
EC_cols=[p for p in ['P2','P3','P8','P11'] if p in mat.columns]
rank_df['EC_endpoint_mean']=mat[EC_cols].mean(axis=1).values
rank_df.to_csv(OUT/'EC_vs_alternative_program_signature_target_scores.csv',index=False)

fig,ax=plt.subplots(figsize=(8.4,4.2))
plot=rank_df.sort_values('EC_vs_alt_signature_similarity').reset_index(drop=True)
plot['rank']=np.arange(1,len(plot)+1)
ax.plot(plot['rank'],plot['EC_vs_alt_signature_similarity'],color='black',lw=1.1)
ax.axhline(0,color='#999999',ls='--',lw=0.8)
for setname,genes in CORE_SETS.items():
    sub=plot[plot['target'].isin(genes)]
    ax.scatter(sub['rank'],sub['EC_vs_alt_signature_similarity'],s=38,color=CORE_COLORS[setname],label=setname,zorder=4)
    for _,r in sub.iterrows(): ax.text(r['rank'],r['EC_vs_alt_signature_similarity'],r['target'],fontsize=7,ha='center',va='bottom')
ax.text(0.02,0.05,'alternative-like\n(P34/P27 up, EC down)',transform=ax.transAxes,fontsize=8,color=PURPLE,ha='left')
ax.text(0.98,0.95,'EC-like\n(EC/P15 up)',transform=ax.transAxes,fontsize=8,color=ECBLUE,ha='right',va='top')
ax.set_xlabel('targets ranked by program-space signature similarity')
ax.set_ylabel('EC-like  ←  similarity  →  alternative-like', fontsize=9)
ax.set_title('EC-vs-alternative endpoint signature model in curated-program space',fontsize=12,weight='bold')
ax.legend(frameon=False,fontsize=8,ncols=2,loc='lower right')
ax.spines[['top','right']].set_visible(False)
savefig(fig,'EC_vs_alternative_signature_ranked_targets')

# Ridge reconstruction, with caution because program-space has few features
try:
    from sklearn.linear_model import RidgeCV
    X=Mz.T.values # programs x targets
    y=weights.values
    model=RidgeCV(alphas=np.logspace(-2,3,60)).fit(X,y)
    coef=pd.DataFrame({'target':Mz.index,'ridge_coef':model.coef_}).sort_values('ridge_coef')
    coef.to_csv(OUT/'EC_vs_alternative_signature_ridge_reconstruction_coefficients.csv',index=False)
    top=pd.concat([coef.head(18),coef.tail(18)]).sort_values('ridge_coef')
    fig,ax=plt.subplots(figsize=(7.4,6.2))
    ax.barh(top['target'],top['ridge_coef'],color=[PURPLE if x<0 else ECBLUE for x in top['ridge_coef']])
    ax.axvline(0,color='black',lw=0.8)
    ax.set_xlabel('ridge coefficient')
    ax.set_title('Program-space reconstruction of EC-vs-alternative signature',fontsize=12,weight='bold')
    ax.text(0.02,0.02,f'Exploratory ridge model; alpha={model.alpha_:.2g}\npositive = reconstructs EC-like signature',transform=ax.transAxes,fontsize=8,ha='left',va='bottom')
    ax.spines[['top','right']].set_visible(False)
    savefig(fig,'EC_vs_alternative_signature_ridge_reconstruction_top_coefficients')
except Exception as e:
    print('Ridge skipped',e)

# =====================================================================
# D. Divergent-regulator quadrants
# =====================================================================
def quadrant(xcol,ycol,xlab,ylab,title,stem):
    df=features.dropna(subset=[xcol,ycol]).copy()
    r,p=pearsonr(df[xcol],df[ycol])
    fig,ax=plt.subplots(figsize=(5.2,4.8))
    ax.scatter(df[xcol],df[ycol],s=22,color=GREY,alpha=0.55,edgecolors='none')
    for setname,genes in CORE_SETS.items():
        sub=df[df['target'].isin(genes)]
        ax.scatter(sub[xcol],sub[ycol],s=50,color=CORE_COLORS[setname],edgecolors='white',lw=0.4,label=setname)
        for _,row in sub.iterrows(): ax.text(row[xcol],row[ycol],row['target'],fontsize=7,ha='left',va='bottom')
    ax.axhline(0,color='#888888',ls='--',lw=0.8); ax.axvline(0,color='#888888',ls='--',lw=0.8)
    ax.set_xlabel(xlab); ax.set_ylabel(ylab)
    ax.set_title(title,fontsize=11.5,weight='bold')
    ax.text(0.03,0.97,f'Pearson r={r:.2f}\nP={p:.1e}',transform=ax.transAxes,ha='left',va='top',fontsize=8,bbox=dict(boxstyle='round,pad=0.25',fc='white',ec='#ccc',lw=0.5))
    ax.legend(frameon=False,fontsize=7,loc='best')
    ax.spines[['top','right']].set_visible(False)
    savefig(fig,stem)

quadrant('D3_EC_mean_effect','P34_YAP_mesenchyme_effect','D3 EC endpoint mean effect','P34 YAP/mesenchyme effect','Divergent regulators: EC endpoint vs P34','quadrant_D3_ECmean_vs_P34')
quadrant('P15_vascular_progenitor_effect','P17_YAP_LM_effect','P15 vascular progenitor effect','P17 YAP/TAZ-LM effect','Divergent regulators: P15 vs P17','quadrant_D2_P15_vs_P17')
quadrant('P34_YAP_mesenchyme_effect','P27_ECM_mesenchyme_effect','P34 YAP/mesenchyme effect','P27 ECM mesenchyme effect','D3 alternative endpoints: P34 vs P27','quadrant_D3_P34_vs_P27')

# =====================================================================
# E. Core regulator similarity analysis
# =====================================================================
profile_z=Mz.copy()
core_sim=pd.DataFrame(index=profile_z.index)
for setname,genes in CORE_SETS.items():
    present=[g for g in genes if g in profile_z.index]
    vals=[]
    for target in profile_z.index:
        refs=[g for g in present if g!=target] or present
        corrs=[np.corrcoef(profile_z.loc[target],profile_z.loc[g])[0,1] for g in refs]
        vals.append(np.nanmean(corrs) if len(corrs) else np.nan)
    core_sim[setname]=vals
core_sim=core_sim.reset_index().rename(columns={'index':'target'})
core_sim.to_csv(OUT/'target_similarity_to_core_regulator_sets.csv',index=False)
core_cols=list(CORE_SETS.keys())
core_sim['max_abs_core_similarity']=core_sim[core_cols].abs().max(axis=1)
top=core_sim.sort_values('max_abs_core_similarity',ascending=False).head(60)
fig,ax=plt.subplots(figsize=(5.4,8.0))
arr=top[core_cols].values
im=ax.imshow(arr,aspect='auto',cmap='RdBu_r',vmin=-1,vmax=1)
ax.set_xticks(range(len(core_cols))); ax.set_xticklabels(core_cols,rotation=45,ha='right',fontsize=8)
ax.set_yticks(range(len(top))); ax.set_yticklabels(top['target'],fontsize=6)
ax.set_title('Targets most similar to core regulator profiles',fontsize=12,weight='bold')
cb=plt.colorbar(im,ax=ax,fraction=0.04,pad=0.02); cb.set_label('mean profile correlation',fontsize=8); cb.ax.tick_params(labelsize=7)
savefig(fig,'core_regulator_similarity_heatmap_top_targets')

if 'Hippo/YAP' in core_sim.columns and 'EC regulators' in core_sim.columns:
    fig,ax=plt.subplots(figsize=(5.3,4.8))
    ax.scatter(core_sim['EC regulators'],core_sim['Hippo/YAP'],s=22,color=GREY,alpha=0.55)
    for setname,genes in CORE_SETS.items():
        sub=core_sim[core_sim['target'].isin(genes)]
        ax.scatter(sub['EC regulators'],sub['Hippo/YAP'],s=55,color=CORE_COLORS[setname],edgecolors='white',lw=0.4,label=setname)
        for _,row in sub.iterrows(): ax.text(row['EC regulators'],row['Hippo/YAP'],row['target'],fontsize=7,ha='left',va='bottom')
    ax.axhline(0,color='#888888',ls='--',lw=0.8); ax.axvline(0,color='#888888',ls='--',lw=0.8)
    ax.set_xlabel('similarity to EC regulator profiles')
    ax.set_ylabel('similarity to Hippo/YAP profiles')
    ax.set_title('Post-hoc similarity to core regulator profiles',fontsize=12,weight='bold')
    ax.legend(frameon=False,fontsize=7)
    ax.spines[['top','right']].set_visible(False)
    savefig(fig,'core_similarity_EC_vs_Hippo_scatter')

# =====================================================================
# F. Core-set rank/rug enrichment on program metrics
# =====================================================================
metric_df=rank_df.merge(features,on='target',how='left')
metrics=[('EC_vs_alt_signature_similarity','EC-vs-alt signature similarity'),
         ('P34_YAP_mesenchyme_effect','P34 effect'),
         ('P17_YAP_LM_effect','P17 effect'),
         ('P15_vascular_progenitor_effect','P15 effect'),
         ('D3_EC_mean_effect','D3 EC mean effect'),
         ('P27_ECM_mesenchyme_effect','P27 effect')]
metrics=[m for m in metrics if m[0] in metric_df.columns]
rank_stats=[]
fig,axes=plt.subplots(len(metrics),1,figsize=(8.4,1.45*len(metrics)),sharex=False)
if len(metrics)==1: axes=[axes]
for ax,(col,label) in zip(axes,metrics):
    df=metric_df.dropna(subset=[col]).sort_values(col).reset_index(drop=True)
    df['rank']=np.arange(1,len(df)+1)
    ax.plot(df['rank'],df[col],color='black',lw=0.8)
    ax.axhline(0,color='#aaa',ls='--',lw=0.7)
    ymin,ymax=df[col].min(),df[col].max(); yr=ymin-0.08*(ymax-ymin if ymax>ymin else 1)
    for setname,genes in CORE_SETS.items():
        sub=df[df['target'].isin(genes)]
        ax.scatter(sub['rank'],[yr]*len(sub),marker='|',s=75,color=CORE_COLORS[setname],label=setname if ax is axes[0] else None)
        if len(sub)>=2:
            other=df.loc[~df['target'].isin(genes),col]
            p=mannwhitneyu(sub[col],other,alternative='two-sided').pvalue
            rank_stats.append({'metric':col,'metric_label':label,'target_set':setname,'n':len(sub),'mean_metric':sub[col].mean(),'median_metric':sub[col].median(),'mannwhitney_two_sided_p':p})
    ax.set_ylabel(label,fontsize=8)
    ax.tick_params(labelsize=7)
    ax.spines[['top','right']].set_visible(False)
axes[0].legend(frameon=False,fontsize=7,ncols=4,loc='upper center',bbox_to_anchor=(0.5,1.45))
axes[-1].set_xlabel('targets ranked by metric',fontsize=9)
fig.suptitle('Post-hoc core-set positions in ranked program-effect metrics',fontsize=12,weight='bold',y=1.01)
fig.tight_layout()
savefig(fig,'core_sets_rank_rug_across_program_metrics')
pd.DataFrame(rank_stats).to_csv(OUT/'core_set_rank_enrichment_statistics.csv',index=False)

# Optional single-cell template for adata_int, not run.
template = r'''# Optional single-cell signature distributions for selected perturbations
# Requires adata_int loaded with per-cell program usages and target carrier labels.
import scanpy as sc
SELECTED_TARGETS = ['NF2','LATS2','SAV1','TAOK1','AMOTL2','PRKCD','SMAD1','SMAD4','SMAD5','KDR','ETV2']
CONTROL_LABEL = 'NT_control'
carrier_col = 'target_carrier'  # replace with your column
EC_PROGRAMS = ['Program 2','Program 3','Program 8','Program 11']
ALT_PROGRAMS = ['Program 34','Program 27']
adata_int.obs['EC_endpoint_score'] = adata_int.obs[EC_PROGRAMS].mean(axis=1)
adata_int.obs['ALT_endpoint_score'] = adata_int.obs[ALT_PROGRAMS].mean(axis=1)
adata_int.obs['P34_score'] = adata_int.obs['Program 34']
adata_int.obs['P17_score'] = adata_int.obs['Program 17']
adata_int.obs['P15_score'] = adata_int.obs['Program 15']
for score in ['EC_endpoint_score','ALT_endpoint_score','P34_score','P17_score','P15_score']:
    sc.pl.violin(adata_int[adata_int.obs[carrier_col].isin(SELECTED_TARGETS+[CONTROL_LABEL])],
                 keys=score, groupby=carrier_col, rotation=90, stripplot=False)
'''
(OUT/'optional_single_cell_signature_distribution_template.py').write_text(template)

# README
readme = '''T-cell-paper-inspired program-first analyses\n\nGenerated analyses:\n1. Fig2-style top regulators and volcano plots for selected curated programs.\n2. Fig3-style target-day perturbation modules and temporal specificity.\n3. Fig4-style EC-vs-alternative signature scoring in curated program space.\n4. Divergent regulator quadrant plots.\n5. Core-regulator similarity and ranked/rug enrichment analyses.\n6. Optional single-cell signature distribution template requiring adata_int.\n\nImportant limitation: gene-level target perturbation DE tables were not available here, so analyses that parse which genes within the EC-vs-P34 signature are controlled by each regulator are included only conceptually/template.\n'''
(OUT/'README_tcell_inspired_program_analyses.txt').write_text(readme)

# copy this script
shutil.copyfile(__file__, OUT/'make_tcell_inspired_program_analyses_clean.py')

# package
zip_path=BASE/'tcell_inspired_program_analyses_clean_package.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for fp in OUT.iterdir():
        if fp.is_file(): z.write(fp, arcname=fp.name)
print('Wrote', OUT)
print('Package', zip_path)
