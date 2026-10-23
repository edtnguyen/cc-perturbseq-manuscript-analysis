import os, zipfile, shutil
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['axes.unicode_minus'] = False

BASE='/mnt/data'
OUT=os.path.join(BASE,'Figure1c_redesign_v8_fresh')
PANEL=os.path.join(OUT,'individual_panels')
SRC=os.path.join(OUT,'source_data')
SCRIPTS=os.path.join(OUT,'scripts')
for d in [OUT,PANEL,SRC,SCRIPTS]: os.makedirs(d, exist_ok=True)

source_path='/mnt/data/Figure1_alternative_v5/source_data/Fig1c_displayed_gene_selection_v5.tsv'
if not os.path.exists(source_path):
    source_path='/mnt/data/Figure1_restyled_1c1d_v2/source_data/Fig1c_nonHippo_known_regulators_source.tsv'

df=pd.read_csv(source_path, sep='\t')
# standardize
if 'display_gene_symbol' not in df.columns and 'gene' in df.columns: df['display_gene_symbol']=df['gene']
if 'effect' not in df.columns:
    for col in ['u_test_effect_log2_CDH5pos_vs_CDH5neg','u_test_log2FC','screen_effect']:
        if col in df.columns: df['effect']=df[col]; break
if 'direction' not in df.columns and 'screen_direction' in df.columns: df['direction']=df['screen_direction']
if 'method_support' not in df.columns:
    if 'mageck_hit' in df.columns and 'u_test_hit' in df.columns:
        def ms(r):
            if bool(r['mageck_hit']) and bool(r['u_test_hit']): return 'MaGECK + U-test'
            if bool(r['mageck_hit']): return 'MaGECK only'
            if bool(r['u_test_hit']): return 'U-test only'
            return 'not called'
        df['method_support']=df.apply(ms,axis=1)
    else: df['method_support']='screen hit'
hippo={'NF2','LATS1','LATS2','SAV1','TAOK1','TAOK2','TAOK3','AMOT','AMOTL1','AMOTL2','WWC1','WWC2','WWC3'}
df=df[~df['display_gene_symbol'].astype(str).str.replace(r'[^A-Za-z0-9].*','',regex=True).isin(hippo)].copy()
required=df[df['direction'].eq('Required for EC differentiation')].sort_values('effect', ascending=True).head(15).copy()
restrict=df[df['direction'].eq('Restricts EC differentiation')].sort_values('effect', ascending=False).copy()
plot_df=pd.concat([required, restrict], ignore_index=True)
plot_df['abs_effect']=plot_df['effect'].abs()
plot_df['direction_short']=np.where(plot_df['direction'].eq('Required for EC differentiation'), 'Required', 'Restricts')
if 'known_EC_regulator' in plot_df.columns:
    plot_df['is_ec_lit']=plot_df['known_EC_regulator'].astype(str).str.lower().isin(['true','1','yes'])
else:
    plot_df['is_ec_lit']=False
plot_df.to_csv(os.path.join(SRC,'Fig1c_redesigned_nonHippo_source_top15plus8.tsv'), sep='\t', index=False)

blue='#2B6F9E'; blue_fill='#D8E7F0'; orange='#C96533'; orange_fill='#F1D8C8'; black='#222222'; gray='#666666'; grid='#E2E2E2'
method_markers={'MaGECK + U-test':'o','MaGECK only':'^','U-test only':'s','screen hit':'o','not called':'x'}
method_labels=['MaGECK + U-test','MaGECK only','U-test only']

# Helper for final axis styling
def style_axis(ax):
    ax.grid(axis='x', color=grid, lw=0.65, alpha=0.9)
    for sp in ['top','right','left']:
        ax.spines[sp].set_visible(False)
    ax.spines['bottom'].set_color('#4A4A4A')
    ax.spines['bottom'].set_linewidth(0.8)
    ax.tick_params(axis='y', length=0, pad=3)
    ax.tick_params(axis='x', length=2.2, width=0.8, labelsize=7.8)

# -------------------------
# Recommended split-bar panel v8
# -------------------------
req=plot_df[plot_df['direction_short']=='Required'].sort_values('abs_effect', ascending=True).copy()
res=plot_df[plot_df['direction_short']=='Restricts'].sort_values('abs_effect', ascending=True).copy()
fig=plt.figure(figsize=(6.65,4.15), dpi=350)
gs=fig.add_gridspec(1,2, width_ratios=[1.60,1.06], wspace=0.42, left=0.115, right=0.975, top=0.81, bottom=0.245)
ax1=fig.add_subplot(gs[0,0]); ax2=fig.add_subplot(gs[0,1])
for ax, sub, col, fill, title, subtitle in [
    (ax1, req, blue, blue_fill, 'Required for EC differentiation', '15 strongest non-Hippo hits'),
    (ax2, res, orange, orange_fill, 'Restricts EC differentiation', 'all 8 hits')]:
    y=np.arange(len(sub))
    maxx=np.nanmax(sub['abs_effect'])*1.15
    ax.barh(y, sub['abs_effect'], color=fill, height=0.58, edgecolor='none', zorder=1)
    # endpoints + method overlay
    ax.scatter(sub['abs_effect'], y, s=40, color=col, edgecolor='white', lw=0.85, zorder=3)
    for method, marker in method_markers.items():
        ss=sub[sub['method_support'].eq(method)]
        if not ss.empty:
            yy=[list(sub.index).index(i) for i in ss.index]
            ax.scatter(ss['abs_effect'], yy, s=48, facecolor='none', edgecolor='#222222', lw=0.75, marker=marker, zorder=4)
    labels=sub['display_gene_symbol'].tolist()
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.8, color=black)
    for tick, (_,r) in zip(ax.get_yticklabels(), sub.iterrows()):
        if r['is_ec_lit']:
            tick.set_color('#1E4E79')
            tick.set_fontweight('bold')
    ax.set_xlim(0, maxx)
    ax.set_ylim(-0.75, len(sub)-0.25)
    ax.set_xlabel('|U-test log2FC|', fontsize=8.5, labelpad=5)
    ax.set_title(title+'\n'+subtitle, fontsize=8.9, fontweight='bold', color=col, pad=5)
    style_axis(ax)

# Clean header, no overlap
fig.text(0.115,0.965,'Known non-Hippo regulators recovered in the screen', fontsize=10.7, fontweight='bold', ha='left', va='top', color=black)
fig.text(0.115,0.922,'Ranked U-test effects; endpoint symbols indicate hit-calling support.', fontsize=7.5, ha='left', va='top', color=gray)
# Legends lower, compact
handles=[Line2D([0],[0], marker=method_markers[m], color='none', markeredgecolor='#222222', markerfacecolor='none', lw=0, markersize=5.4, label=m) for m in method_labels]
handles += [Line2D([0],[0], marker='o', color='none', markeredgecolor='none', markerfacecolor='#1E4E79', lw=0, markersize=5.2, label='bold label = literature-supported EC/vascular regulator')]
fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.54,0.060), ncol=2, frameon=False, fontsize=7.1, columnspacing=1.4, handletextpad=0.45)
for ext in ['png','pdf','svg']:
    fig.savefig(os.path.join(PANEL,f'Fig1c_splitbars_publication_v8.{ext}'), dpi=350 if ext=='png' else None, bbox_inches='tight')
plt.close(fig)

# -------------------------
# Alternate compact centered bar panel v8
# -------------------------
water=plot_df.copy()
water['plot_effect']=np.where(water['direction_short'].eq('Required'), -water['abs_effect'], water['abs_effect'])
req2=water[water['direction_short']=='Required'].sort_values('plot_effect', ascending=True)
res2=water[water['direction_short']=='Restricts'].sort_values('plot_effect', ascending=True)
water2=pd.concat([req2,res2], ignore_index=True)
water2['y']=np.arange(len(water2))[::-1]
fig,ax=plt.subplots(figsize=(5.95,4.85), dpi=350)
for _,r in water2.iterrows():
    col=blue if r['direction_short']=='Required' else orange
    fill=blue_fill if r['direction_short']=='Required' else orange_fill
    ax.plot([0,r['plot_effect']],[r['y'],r['y']], color=fill, lw=4.2, solid_capstyle='round', zorder=1)
    ax.scatter([r['plot_effect']],[r['y']], s=39, color=col, edgecolor='white', lw=0.8, zorder=3)
    marker=method_markers.get(r['method_support'],'o')
    ax.scatter([r['plot_effect']],[r['y']], s=47, facecolor='none', edgecolor='#222222', lw=0.75, marker=marker, zorder=4)
    ha='right' if r['plot_effect']<0 else 'left'
    tx=r['plot_effect'] - 0.07 if r['plot_effect']<0 else r['plot_effect']+0.07
    ax.text(tx, r['y'], r['display_gene_symbol'], ha=ha, va='center', fontsize=7.3, color='#1E4E79' if r['is_ec_lit'] else black, fontweight='bold' if r['is_ec_lit'] else 'normal')
ax.axvline(0, color='#555555', lw=0.75)
lim=max(abs(water2['plot_effect']).max()*1.33,3.1)
ax.set_xlim(-lim,lim*0.82); ax.set_ylim(-1, water2['y'].max()+1.6)
ax.set_yticks([]); ax.grid(axis='x', color=grid, lw=0.6)
ax.set_xlabel('U-test log2FC in CDH5+ vs CDH5- cells', fontsize=8.6, labelpad=5)
ax.set_title('Bidirectional non-Hippo screen regulators', fontsize=10.8, fontweight='bold', pad=9)
ax.text(-lim*0.50, water2['y'].max()+0.83, 'Required for EC differentiation', color=blue, ha='center', va='center', fontsize=8.1, fontweight='bold')
ax.text(lim*0.37, water2['y'].max()+0.83, 'Restricts EC differentiation', color=orange, ha='center', va='center', fontsize=8.1, fontweight='bold')
for sp in ['top','right','left']:
    ax.spines[sp].set_visible(False)
ax.spines['bottom'].set_color('#4A4A4A'); ax.spines['bottom'].set_linewidth(0.8)
ax.tick_params(axis='x', labelsize=7.8)
handles=[Line2D([0],[0], marker=method_markers[m], color='none', markeredgecolor='#222222', markerfacecolor='none', lw=0, markersize=5.2, label=m) for m in method_labels]
handles += [Line2D([0],[0], marker='o', color='#1E4E79', lw=0, markersize=4.8, label='EC/vascular regulator label')]
ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.52,-0.235), ncol=2, frameon=False, fontsize=6.7, handletextpad=0.4, columnspacing=1.1)
for ext in ['png','pdf','svg']:
    fig.savefig(os.path.join(PANEL,f'Fig1c_centered_publication_v8.{ext}'), dpi=350 if ext=='png' else None, bbox_inches='tight')
plt.close(fig)

# QC / README / package
shutil.copy(__file__, os.path.join(SCRIPTS,'build_figure1c_redesign_v8.py'))
qc=[]
for fn in os.listdir(PANEL):
    if fn.endswith('.svg'):
        txt=open(os.path.join(PANEL,fn), errors='ignore').read()
        qc.append({'file':fn,'text_tag_count':txt.count('<text'),'path_tag_count':txt.count('<path')})
pd.DataFrame(qc).to_csv(os.path.join(OUT,'svg_text_qc.tsv'), sep='\t', index=False)
with open(os.path.join(OUT,'README.md'),'w') as f:
    f.write('''Figure 1c redesign v8 fresh\n\nRecommended file: individual_panels/Fig1c_splitbars_publication_v8.svg/pdf/png.\nAlternate file: individual_panels/Fig1c_centered_publication_v8.svg/pdf/png.\n\nDisplay: top 15 non-Hippo hits required for EC differentiation plus all 8 hits that restrict EC differentiation. Hippo genes are filtered out. Text is kept editable in SVG using Matplotlib svg.fonttype=none.\n''')
zip_path=os.path.join(BASE,'Figure1c_redesign_v8_fresh_package.zip')
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for root, dirs, files in os.walk(OUT):
        for fn in files:
            fp=os.path.join(root,fn)
            z.write(fp, arcname=os.path.relpath(fp, OUT))
print(zip_path)
