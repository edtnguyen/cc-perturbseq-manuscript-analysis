#!/usr/bin/env python
import os, re, zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from scipy.spatial.distance import pdist
from scipy.stats import zscore

IN_MATRIX='/mnt/data/target_day_clustering_full_heatmap_modules_v2/target_day_input_matrix_effect_if_FDR_lt_0.05_else_0.csv'
OUTDIR='/mnt/data/target_level_zoom_M4_M5_k15_downloadable'
os.makedirs(OUTDIR, exist_ok=True)

PROGRAMS=["P7","P22","P29","P28","P13","P30","P48","P33","P15","P14","P17","P34","P27","P2","P3","P8","P11"]
PROGRAM_LABELS={
    "P7":"P7\nPluri", "P22":"P22\nPluri", "P29":"P29\nPluri",
    "P28":"P28\nPrimitive\nstreak", "P13":"P13\nPost-\nPS", "P30":"P30\nAnt-\nPS", "P48":"P48\nNeuro-\nectoderm", "P33":"P33\nAnti-\nmesoderm?",
    "P15":"P15\nVascular\nprogenitor", "P14":"P14\nLateral\nmesoderm", "P17":"P17\nYAP/TAZ\nLM",
    "P34":"P34\nYAP/TAZ\nmesenchyme", "P27":"P27\nECM\nmesenchyme",
    "P2":"P2\nArtery\nEC 1", "P3":"P3\nEC\nfunction 1", "P8":"P8\nArtery\nEC 2", "P11":"P11\nEC\nfunction 2",
}
DAY_GROUPS={
    'D0':['P7','P22','P29'],
    'D1':['P28','P13','P30','P48','P33'],
    'D2':['P15','P14','P17'],
    'D3':['P34','P27','P2','P3','P8','P11'],
}
DAY_COLORS={'D0':'#D9D9D9','D1':'#F4B183','D2':'#9DC3E6','D3':'#C6E0B4'}
HIPPO=['NF2','LATS2','SAV1','TAOK1','AMOTL2']

# Reconstruct target-level FDR-masked matrix from available target-day matrix.
td=pd.read_csv(IN_MATRIX,index_col=0)
def target_from_item(x):
    return re.sub(r'_D[0-3]$','',str(x))
td['target']=[target_from_item(i) for i in td.index]
mat=td.groupby('target')[PROGRAMS].sum()
mat=mat.loc[(mat.abs().sum(axis=1)>0)].copy()

# Cluster target-level profiles at k=15.
X=mat.values.astype(float)
Xz=zscore(X,axis=0,nan_policy='omit')
Xz=np.nan_to_num(Xz,nan=0.0,posinf=0.0,neginf=0.0)
Z=linkage(pdist(Xz,metric='cosine'),method='average')
leaf_order=leaves_list(Z)
ordered_targets=mat.index[leaf_order].tolist()
raw_modules=pd.Series(fcluster(Z,t=15,criterion='maxclust'),index=mat.index,name='raw_module')
module_order=[]
for t in ordered_targets:
    m=raw_modules.loc[t]
    if m not in module_order:
        module_order.append(m)
module_map={old:i+1 for i,old in enumerate(module_order)}
modules=raw_modules.map(module_map).rename('module')

assign=pd.DataFrame({'target':mat.index,'module':modules.values,'raw_module':raw_modules.values})
assign['is_hippo']=assign['target'].isin(HIPPO)
assign.to_csv(os.path.join(OUTDIR,'target_level_module_assignments_k15.csv'),index=False)

selected_targets=[t for t in ordered_targets if int(modules.loc[t]) in [4,5]]
plot_mat=mat.loc[selected_targets,PROGRAMS]
row_modules=modules.loc[selected_targets]

zoom_df=pd.DataFrame({'target':selected_targets,'module':[int(modules.loc[t]) for t in selected_targets],'is_hippo':[t in HIPPO for t in selected_targets]})
zoom_df=pd.concat([zoom_df.reset_index(drop=True), plot_mat.reset_index(drop=True)], axis=1)
zoom_df.to_csv(os.path.join(OUTDIR,'zoom_target_level_M4_M5_members_k15.csv'),index=False)

centroids=mat.join(modules).groupby('module')[PROGRAMS].mean()
centroids.to_csv(os.path.join(OUTDIR,'target_level_module_centroids_k15.csv'))
summary=[]
for m in sorted(modules.unique()):
    mem=modules[modules==m].index.tolist()
    summary.append({'module':int(m),'n_targets':len(mem),'hippo_targets':';'.join([g for g in mem if g in HIPPO])})
pd.DataFrame(summary).to_csv(os.path.join(OUTDIR,'target_level_module_summary_k15.csv'),index=False)

# Plot using separate axes for heatmap/text/colorbar.
n_rows=len(selected_targets)
n_cols=len(PROGRAMS)
fig_h=max(9.0,0.205*n_rows+2.4)
fig=plt.figure(figsize=(13.4,fig_h))
gs=fig.add_gridspec(nrows=1,ncols=3,width_ratios=[1.0,0.24,0.035],wspace=0.08)
ax=fig.add_subplot(gs[0,0])
ax_text=fig.add_subplot(gs[0,1])
cax=fig.add_subplot(gs[0,2])

clip_val=4.2
norm=TwoSlopeNorm(vmin=-clip_val,vcenter=0,vmax=clip_val)
im=ax.imshow(plot_mat.clip(-clip_val,clip_val).values,aspect='auto',cmap='RdBu_r',norm=norm,interpolation='nearest')

ax.set_xticks(np.arange(n_cols))
ax.set_xticklabels([PROGRAM_LABELS.get(p,p) for p in PROGRAMS],rotation=55,ha='right',fontsize=8)
row_labels=[]
for t in selected_targets:
    star='★ ' if t in HIPPO else '  '
    row_labels.append(f"{star}M{int(modules.loc[t])} | {t}")
ax.set_yticks(np.arange(n_rows))
ax.set_yticklabels(row_labels,fontsize=7)

# module separators and outlines
for idx in range(1,n_rows):
    if row_modules.iloc[idx] != row_modules.iloc[idx-1]:
        ax.axhline(idx-0.5,color='black',lw=1.0)
ax.axhline(-0.5,color='black',lw=0.8)
ax.axhline(n_rows-0.5,color='black',lw=0.8)

# day top bars and vertical boundaries
for day,plist in DAY_GROUPS.items():
    start=PROGRAMS.index(plist[0])
    end=PROGRAMS.index(plist[-1])
    ax.add_patch(Rectangle((start-0.5,-2.40),end-start+1,0.60,facecolor=DAY_COLORS[day],edgecolor='none',clip_on=False))
    ax.text((start+end)/2,-2.88,day,ha='center',va='center',fontsize=10.5,weight='bold',clip_on=False)
    ax.axvline(start-0.5,color='black',lw=0.8)
ax.axvline(n_cols-0.5,color='black',lw=0.8)

# subtle gridlines
for x in np.arange(-0.5,n_cols,1):
    ax.axvline(x,color='white',lw=0.35,alpha=0.55)
for y in np.arange(-0.5,n_rows,1):
    ax.axhline(y,color='white',lw=0.25,alpha=0.40)

ax.set_ylim(n_rows-0.5,-3.25)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.8)

# Title across figure
fig.suptitle('Zoom on target-level modules M4/M5 from the k=15 clustered heatmap',fontsize=16,weight='bold',y=0.988)
fig.text(0.40,0.952,'Rows are target-level profiles; ★ marks Hippo/YAP genes. Modules were defined without pathway labels.',ha='center',va='top',fontsize=10,color='#555555')

# Right interpretation panel
ax_text.axis('off')
ax_text.text(0.0,0.98,'Interpretation',ha='left',va='top',fontsize=12,weight='bold')
interp=(
    'M4: P27-up / modest P34-up\n'
    '    includes weaker SAV1 profile\n\n'
    'M5: EC-down + P27/P34-up\n'
    '    contains NF2, LATS2,\n'
    '    TAOK1, AMOTL2\n\n'
    'Rows are whole target profiles,\n'
    'not target-day profiles.\n\n'
    'Displayed value: effect retained\n'
    'only when FDR < 0.05.'
)
ax_text.text(0.0,0.92,interp,ha='left',va='top',fontsize=9.5,linespacing=1.35)

cbar=fig.colorbar(im,cax=cax)
cbar.set_label('effect\nFDR < 0.05',fontsize=9)
cbar.ax.tick_params(labelsize=8)

fig.tight_layout(rect=[0.02,0.02,0.98,0.94])
for ext in ['png','pdf','svg']:
    path=os.path.join(OUTDIR,f'zoom_M4_M5_target_level_profiles_k15.{ext}')
    fig.savefig(path,dpi=350 if ext=='png' else None,bbox_inches='tight')
plt.close(fig)

# README
with open(os.path.join(OUTDIR,'README.txt'),'w') as f:
    f.write('Target-level M4/M5 zoom heatmap at k=15.\n')
    f.write('Rows are target genes, columns are expert-curated programs.\n')
    f.write('The matrix was reconstructed from the target-day FDR-masked effect matrix by summing target-day rows to target-level profiles.\n')
    f.write('Displayed value = effect size if FDR < 0.05, otherwise 0.\n')
    f.write('Clustering: column z-score -> cosine distance -> average linkage -> k=15, modules renumbered by dendrogram order.\n')

zip_path='/mnt/data/target_level_zoom_M4_M5_k15_downloadable_package.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    z.write(__file__,arcname='make_zoom_M4_M5_target_level_profiles_k15.py')
    for fn in os.listdir(OUTDIR):
        fp=os.path.join(OUTDIR,fn)
        if os.path.isfile(fp):
            z.write(fp,arcname=fn)
print('Wrote',OUTDIR)
print('Package',zip_path)
