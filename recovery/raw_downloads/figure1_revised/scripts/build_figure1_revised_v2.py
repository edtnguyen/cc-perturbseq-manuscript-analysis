import os, re, zipfile, shutil, textwrap
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
from scipy.stats import pearsonr, spearmanr
try:
    from adjustText import adjust_text
except Exception:
    adjust_text=None

ROOT = Path('/mnt/data/figure1_revised_v2')
for sub in ['individual_panels','assembled_figures','source_data','scripts','data_audit','legends','candidate_panels','cleaned_data','supplementary_panels']:
    (ROOT/sub).mkdir(parents=True, exist_ok=True)
BASE = Path('/mnt/data/figure1_revised')
# copy source/audit data from prior build when available
for sub in ['cleaned_data','source_data','data_audit']:
    src=BASE/sub
    dst=ROOT/sub
    if src.exists():
        for f in src.glob('*'):
            if f.is_file(): shutil.copy(f,dst/f.name)

mpl.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':0.8})
COL_REQ='#3B73B9'; COL_RES='#D96B27'; COL_GRAY='#C9C9C9'; COL_DARK='#222222'; COL_ORANGE=COL_RES; COL_BLUE=COL_REQ

def save(fig, path, dpi=450):
    path=Path(path)
    fig.savefig(path.with_suffix('.svg'), bbox_inches='tight')
    fig.savefig(path.with_suffix('.pdf'), bbox_inches='tight')
    fig.savefig(path.with_suffix('.png'), dpi=dpi, bbox_inches='tight')

def neglog10(x): return -np.log10(np.clip(np.asarray(x,dtype=float),1e-300,None))
# Load cleaned data
gene = pd.read_csv(BASE/'cleaned_data/genomewide_results_unique_gene_level.tsv', sep='\t')
hits = pd.read_csv(BASE/'cleaned_data/screen_hits_244_unique_gene_level.tsv', sep='\t')
val = pd.read_csv(BASE/'cleaned_data/validation_results_harmonized.tsv', sep='\t')
# normalize validation column names across previous builds
if 'normalized_gene_symbol' not in val.columns:
    val['normalized_gene_symbol'] = val['gene'] if 'gene' in val.columns else val.get('Gene')
if 'screen_effect_mageck_lfc' not in val.columns and 'genome_wide_screen_log2fc' in val.columns:
    val['screen_effect_mageck_lfc'] = val['genome_wide_screen_log2fc']
if 'validation_utest_lfc' not in val.columns and 'validation_screen_log2fc' in val.columns:
    val['validation_utest_lfc'] = val['validation_screen_log2fc']
if 'validation_rep1_lfc' not in val.columns and 'U Test FC BR1' in val.columns:
    val['validation_rep1_lfc'] = val['U Test FC BR1']
if 'validation_rep2_lfc' not in val.columns and 'U Test FC BR2' in val.columns:
    val['validation_rep2_lfc'] = val['U Test FC BR2']
if 'validation_utest_fdr' not in val.columns and 'FDR-corrected' in val.columns:
    val['validation_utest_fdr'] = val['FDR-corrected']
if 'target_group' not in val.columns and 'Type' in val.columns:
    val['target_group'] = val['Type']
comp = pd.read_csv(BASE/'source_data/Fig1d_component_effects.tsv', sep='\t')
# Build if missing columns
if 'screen_effect_mageck_lfc' not in gene.columns:
    gene['screen_effect_mageck_lfc'] = gene['screen_effect'] if 'screen_effect' in gene.columns else gene['mageck_lfc']
if 'mageck_direction_fdr' not in gene.columns:
    gene['mageck_direction_fdr'] = gene['best_mageck_fdr'] if 'best_mageck_fdr' in gene.columns else gene[['neg|fdr','pos|fdr']].min(axis=1)
if 'method_support' not in gene.columns:
    gene['method_support'] = gene.get('hit_support', 'Non-hit')
if 'minus_log10_mageck_direction_fdr' not in gene.columns:
    gene['minus_log10_mageck_direction_fdr']=neglog10(gene['mageck_direction_fdr'])
# final quantities
N_full=gene['normalized_gene_symbol'].nunique(); N_hit=int(gene['union_hit'].sum()); N_req=int((hits['screen_direction']=='Required for EC differentiation').sum()); N_res=int((hits['screen_direction']=='Restricts EC differentiation').sum())
# label selections
labels_1b=['CDH5','KDR','ETV2','PLCG1','SMAD1','SMAD5','ACVR2A','ACVR2B','BMPR1A','CBLL1','SMC2']
labels_1d=['CDH5','KDR','ETV2','PLCG1','SMAD1','SMAD4','SMAD5','BMPR1A','ACVR2A','ACVR2B','NF2','LATS2','AMOTL2','SMARCA4','SMARCC1','MED12']
labels_1e=['CDH5','KDR','ETV2','PLCG1','ACVR2A','ACVR2B','SMAD1','SMAD5','BMPR1A']

# ---------------- Panel A ----------------
def fig1a():
    fig, ax = plt.subplots(figsize=(5.4,2.55))
    ax.set_xlim(0,10); ax.set_ylim(0,5.2); ax.axis('off')
    ax.text(0,5.05,'Genome-wide CRISPRi screen design',fontsize=12.5,fontweight='bold',ha='left',va='top')
    boxes=[(0.9,3.15,'CRISPRi\nhiPSCs'),(2.35,3.15,'Genome-wide\ngRNA library'),(3.8,3.15,'Dox +\nselection'),(5.25,3.15,'Directed\ndifferentiation'),(7.0,3.85,'CDH5+\nEC fraction'),(7.0,2.45,'CDH5−\nnon-EC'),(8.75,3.15,'gRNA\nsequencing')]
    fills=['#F7F7F7','#F7F7F7','#F7F7F7','#EFF6FF','#E6F1FF','#FFF1E5','#F7F7F7']
    for (x,y,t),fc in zip(boxes,fills):
        ax.add_patch(FancyBboxPatch((x-0.55,y-0.38),1.1,0.76,boxstyle='round,pad=0.06,rounding_size=0.06',fc=fc,ec='#4A4A4A',lw=0.8))
        ax.text(x,y,t,fontsize=7.6,ha='center',va='center')
    def arrow(x1,y1,x2,y2,c='#555',lw=0.9): ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=8,lw=lw,color=c))
    arrow(1.48,3.15,1.78,3.15); arrow(2.93,3.15,3.23,3.15); arrow(4.38,3.15,4.68,3.15)
    arrow(5.85,3.34,6.42,3.72,COL_REQ,1.0); arrow(5.85,2.96,6.42,2.58,COL_RES,1.0); arrow(7.58,3.85,8.22,3.35); arrow(7.58,2.45,8.22,2.95)
    ax.text(5.25,1.65,'D0 hiPSC → D1 primitive streak → D2 lateral mesoderm → D3 endpoint',fontsize=7.7,ha='center')
    ax.text(0.9,1.05,f'{N_full:,} genes tested\n{N_hit} screen-derived regulators',ha='center',fontsize=7.2)
    ax.text(7.0,1.05,'Phenotype: guide enrichment\nin CDH5+ vs CDH5− cells',ha='center',fontsize=7.2,fontweight='bold')
    save(fig, ROOT/'individual_panels/Fig1a_screen_workflow_v2')
    plt.close(fig)

# ---------------- Panel B ----------------
def fig1b():
    df=gene.copy()
    fig, ax = plt.subplots(figsize=(6.6,4.2))
    # nonhits grey
    bg=df[~df['union_hit']]
    ax.scatter(bg['screen_effect_mageck_lfc'], bg['minus_log10_mageck_direction_fdr'], s=5, color='#D4D4D4', alpha=0.32, edgecolors='none', rasterized=True)
    hit=df[df['union_hit']]
    colors=np.where(hit['screen_direction']=='Required for EC differentiation',COL_REQ,COL_RES)
    ax.scatter(hit['screen_effect_mageck_lfc'], hit['minus_log10_mageck_direction_fdr'], s=20, color=colors, alpha=0.9, edgecolors='none')
    both=hit[hit['method_support'].isin(['Both methods','Both'])]
    ax.scatter(both['screen_effect_mageck_lfc'], both['minus_log10_mageck_direction_fdr'], s=36, facecolors='none', edgecolors='black', lw=0.45)
    # labels with cutoff, no top pile
    lab=df[df['normalized_gene_symbol'].isin(labels_1b)].copy()
    lab=lab[lab['union_hit']]
    texts=[]
    for _,r in lab.iterrows():
        texts.append(ax.text(r['screen_effect_mageck_lfc'], r['minus_log10_mageck_direction_fdr'], r['normalized_gene_symbol'], fontsize=7.0, color='#111'))
        ax.scatter([r['screen_effect_mageck_lfc']],[r['minus_log10_mageck_direction_fdr']],s=48,facecolors='none',edgecolors='black',lw=0.6)
    if adjust_text and texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-',lw=0.35,color='#777'), expand=(1.1,1.2), time_lim=1.0, max_move=25)
    ax.axvline(0,color='#9A9A9A',ls='--',lw=0.8)
    ax.axhline(-np.log10(0.05),color='#9A9A9A',ls=':',lw=0.8)
    ax.set_xlabel('MaGECK log2 fold change (CDH5+ vs CDH5−)',fontsize=8.5)
    ax.set_ylabel('−log10 direction-specific MaGECK FDR',fontsize=8.5)
    ax.set_title('Genome-wide screen identifies 244 regulators',fontsize=11,fontweight='bold')
    ax.text(0.02,0.98,f'{N_req} required for EC differentiation\n{N_res} restrict EC differentiation',transform=ax.transAxes,ha='left',va='top',fontsize=7.3)
    handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=COL_REQ,label='Required for EC differentiation',markersize=6),Line2D([0],[0],marker='o',color='none',markerfacecolor=COL_RES,label='Restricts EC differentiation',markersize=6),Line2D([0],[0],marker='o',color='black',markerfacecolor='none',label='Both methods',markersize=6)]
    ax.legend(handles=handles,frameon=False,fontsize=6.8,loc='lower left')
    ax.spines[['top','right']].set_visible(False)
    save(fig, ROOT/'individual_panels/Fig1b_screen_result_v2')
    # source data
    df.to_csv(ROOT/'source_data/Fig1b_source_data_v2.tsv',sep='\t',index=False)
    plt.close(fig)

# ---------------- Panel D ----------------
def fig1d():
    df=comp[comp['gene_in_screen']].copy()
    # friendly row names
    rename={'VEGF / endothelial signaling':'EC / VEGF specification','TGF-beta / BMP signaling':'TGFβ / BMP signaling','SWI/SNF and chromatin':'SWI/SNF and chromatin'}
    df['pathway_display_name']=df['pathway_display_name'].replace(rename)
    order=['EC / VEGF specification','TGFβ / BMP signaling','Hippo / YAP signaling','SWI/SNF and chromatin','Mediator / transcription']
    order=[o for o in order if o in set(df['pathway_display_name'])]
    if len(order)==0:
        order=list(df['pathway_display_name'].dropna().unique())[:5]
    df=df[df['pathway_display_name'].isin(order)].copy()
    df['row']=df['pathway_display_name'].map({p:i for i,p in enumerate(order[::-1])})
    fig, ax = plt.subplots(figsize=(6.6,3.35))
    for p,i in {p:i for i,p in enumerate(order[::-1])}.items():
        sub=df[df['pathway_display_name']==p]
        ax.axhline(i,color='#E0E0E0',lw=0.7,zorder=0)
        # background strip extent
        ax.text(-4.05,i+0.23,p,fontsize=7.0,fontweight='bold',ha='left',va='bottom')
        n_hit=int(sub['gene_is_final_hit'].sum()); n_all=len(sub)
        ax.text(2.05,i+0.23,f'{n_hit}/{n_all} hits',fontsize=6.2,ha='right',va='bottom',color='#555')
    # points
    colors=np.where(df['screen_effect']<0,COL_REQ,COL_RES)
    ax.scatter(df['screen_effect'], df['row'], s=np.where(df['gene_is_final_hit'],48,22), c=colors, alpha=np.where(df['gene_is_final_hit'],0.95,0.28), edgecolors=np.where(df['gene_is_final_hit'],'black','#999'), lw=np.where(df['gene_is_final_hit'],0.45,0.2), zorder=2)
    lab=df[(df['member_gene'].isin(labels_1d)) & (df['gene_is_final_hit'])].copy()
    texts=[]
    for _,r in lab.iterrows():
        texts.append(ax.text(r['screen_effect'], r['row']-0.22, r['member_gene'], fontsize=6.0, ha='center'))
    if adjust_text and texts:
        adjust_text(texts,ax=ax,only_move={'text':'x'},arrowprops=dict(arrowstyle='-',lw=0.25,color='#777'),time_lim=0.8)
    ax.axvline(0,color='#9A9A9A',ls='--',lw=0.8)
    ax.set_xlim(-4.15,2.2)
    ax.set_yticks([])
    ax.set_xlabel('MaGECK log2 fold change (CDH5+ vs CDH5−)',fontsize=8.2)
    ax.set_title('Screen hits converge on pathway components',fontsize=11,fontweight='bold')
    ax.text(0.01,0.02,'Faint points: tested component members; outlined points: final screen regulators',transform=ax.transAxes,fontsize=6.0,ha='left',va='bottom')
    ax.spines[['top','right','left']].set_visible(False)
    df.to_csv(ROOT/'source_data/Fig1d_component_effects_v2.tsv',sep='\t',index=False)
    save(fig,ROOT/'individual_panels/Fig1d_pathway_components_v2')
    plt.close(fig)

# ---------------- Panel E ----------------
def fig1e():
    df=val.dropna(subset=['screen_effect_mageck_lfc','validation_utest_lfc']).copy()
    fig=plt.figure(figsize=(6.6,4.2))
    gs=fig.add_gridspec(2,2,width_ratios=[2.25,1.28],height_ratios=[1,1],wspace=0.42,hspace=0.45)
    ax=fig.add_subplot(gs[:,0]); ax2=fig.add_subplot(gs[0,1]); ax3=fig.add_subplot(gs[1,1])
    group_colors={'genome_wide_hit':COL_REQ,'non_hit':'#9E9E9E','candidate':COL_RES,'curated_candidate':COL_RES,'unknown':'#CFCFCF'}
    labels={'genome_wide_hit':'Screen hit','non_hit':'Random non-hit','candidate':'Curated candidate','curated_candidate':'Curated candidate','unknown':'unknown'}
    for grp,sub in df.groupby('target_group'):
        ax.scatter(sub['screen_effect_mageck_lfc'],sub['validation_utest_lfc'],s=28,c=group_colors.get(grp,'#CCC'),alpha=0.75,edgecolors='white',lw=0.3,label=labels.get(grp,grp))
    ax.axhline(0,color='#999',ls='--',lw=0.8); ax.axvline(0,color='#999',ls='--',lw=0.8)
    x=df['screen_effect_mageck_lfc']; y=df['validation_utest_lfc']; rr,pp=pearsonr(x,y); sp,_=spearmanr(x,y)
    xx=np.linspace(x.min(),x.max(),120); cf=np.polyfit(x,y,1); ax.plot(xx,cf[0]*xx+cf[1],color='black',lw=1.1)
    lab=df[df['normalized_gene_symbol'].isin(labels_1e)]
    texts=[]
    for _,r in lab.iterrows(): texts.append(ax.text(r['screen_effect_mageck_lfc'],r['validation_utest_lfc'],r['normalized_gene_symbol'],fontsize=6.3))
    if adjust_text and texts: adjust_text(texts,ax=ax,arrowprops=dict(arrowstyle='-',lw=0.25,color='#777'),time_lim=0.8)
    ax.set_xlabel('Genome-wide screen effect',fontsize=8.5); ax.set_ylabel('Validation screen effect',fontsize=8.5)
    ax.set_title('Secondary screen validation',fontsize=11,fontweight='bold')
    ax.text(0.02,0.98,f'n={len(df)} targets\nPearson r={rr:.2f}\nSpearman ρ={sp:.2f}',transform=ax.transAxes,va='top',fontsize=7.0)
    ax.legend(frameon=False,fontsize=6.3,loc='lower right')
    ax.spines[['top','right']].set_visible(False)
    reps=val.dropna(subset=['validation_rep1_lfc','validation_rep2_lfc']).copy()
    ax2.scatter(reps['validation_rep1_lfc'],reps['validation_rep2_lfc'],s=10,color='#737373',alpha=0.6)
    ax2.axhline(0,color='#AAA',ls='--',lw=0.6); ax2.axvline(0,color='#AAA',ls='--',lw=0.6)
    r2,_=pearsonr(reps['validation_rep1_lfc'],reps['validation_rep2_lfc']); xx=np.linspace(reps['validation_rep1_lfc'].min(),reps['validation_rep1_lfc'].max(),100); c2=np.polyfit(reps['validation_rep1_lfc'],reps['validation_rep2_lfc'],1); ax2.plot(xx,c2[0]*xx+c2[1],color='black',lw=0.8)
    ax2.set_title('Replicates',fontsize=8,fontweight='bold'); ax2.set_xlabel('rep 1',fontsize=6.5); ax2.set_ylabel('rep 2',fontsize=6.5); ax2.tick_params(labelsize=6); ax2.text(0.05,0.95,f'r={r2:.2f}',transform=ax2.transAxes,va='top',fontsize=6.2); ax2.spines[['top','right']].set_visible(False)
    sm=[]
    for grp,sub in val.groupby('target_group'):
        sm.append({'group':grp,'label':labels.get(grp,grp),'n':len(sub),'validated':int((sub['validation_utest_fdr']<0.05).sum()),'frac':(sub['validation_utest_fdr']<0.05).mean()})
    sm=pd.DataFrame(sm); sm=sm[sm['group'].isin(['genome_wide_hit','non_hit','candidate','curated_candidate'])]
    ax3.barh(np.arange(len(sm)), sm['frac'], color=[group_colors.get(g,'#CCC') for g in sm['group']], alpha=0.85)
    ax3.set_yticks(np.arange(len(sm))); ax3.set_yticklabels(sm['label'],fontsize=6.2); ax3.set_xlim(0,1); ax3.set_xlabel('fraction',fontsize=6.5); ax3.set_title('Validation rate',fontsize=8,fontweight='bold'); ax3.tick_params(labelsize=6)
    for i,r in sm.iterrows(): ax3.text(r['frac']+0.02, list(sm.index).index(i), f"{int(r['validated'])}/{int(r['n'])}",va='center',fontsize=6)
    ax3.spines[['top','right']].set_visible(False)
    val.to_csv(ROOT/'source_data/Fig1e_validation_source_v2.tsv',sep='\t',index=False)
    save(fig,ROOT/'individual_panels/Fig1e_validation_v2')
    plt.close(fig)

fig1a(); fig1b(); fig1d(); fig1e()
# assemble four-panel fallback
import matplotlib.image as mpimg
def assemble():
    fig=plt.figure(figsize=(12.6,8.0))
    gs=fig.add_gridspec(2,2,width_ratios=[0.38,0.62],height_ratios=[0.44,0.56],wspace=0.06,hspace=0.13)
    imgs=[('a',ROOT/'individual_panels/Fig1a_screen_workflow_v2.png',gs[0,0]),('b',ROOT/'individual_panels/Fig1b_screen_result_v2.png',gs[0,1]),('c',ROOT/'individual_panels/Fig1d_pathway_components_v2.png',gs[1,0]),('d',ROOT/'individual_panels/Fig1e_validation_v2.png',gs[1,1])]
    for letter,path,g in imgs:
        ax=fig.add_subplot(g); ax.imshow(mpimg.imread(path)); ax.axis('off'); ax.text(-0.02,1.03,letter,transform=ax.transAxes,fontsize=16,fontweight='bold',va='top')
    fig.suptitle('Genome-wide CRISPRi screen identifies validated regulators of endothelial differentiation',fontsize=15,fontweight='bold',y=0.99)
    fig.savefig(ROOT/'assembled_figures/Figure1_revised_fallback_v2.png',dpi=600,bbox_inches='tight')
    fig.savefig(ROOT/'assembled_figures/Figure1_revised_fallback_v2.pdf',bbox_inches='tight')
    fig.savefig(ROOT/'assembled_figures/Figure1_revised_fallback_v2.svg',bbox_inches='tight')
    fig.savefig(ROOT/'assembled_figures/Figure1_revised_fallback_v2_600dpi.tiff',dpi=600,bbox_inches='tight')
    plt.close(fig)
assemble()
# placeholder layout optional
# copy artifacts from previous build needed
for f in ['file_inventory.tsv','analysis_version_audit.tsv','file_to_panel_map.tsv','consolidated_missing_input_request.md']:
    src=Path('/mnt/data/figure1_revised/data_audit')/f
    if src.exists(): shutil.copy(src, ROOT/'data_audit'/f)
for f in ['unresolved_issues.md','final_QC_checklist.md']:
    src=Path('/mnt/data/figure1_revised')/f
    if src.exists(): shutil.copy(src, ROOT/f)
# legend/README
Path(ROOT/'README_figure1_revised_v2.md').write_text('''# Figure 1 revised v2\n\nThis build uses the revised main-text concept from Pasted text(4). No standalone ranked-hit panel is included. Panel 1c is absent in the selected fallback because no D0-versus-D3/baseline-to-endpoint abundance table was found after auditing available files.\n\nSelected panels:\n- 1a: genome-wide CRISPRi workflow and phenotype definition\n- 1b: full gene-level MaGECK screen result with U-test support encoded by outline\n- 1c in assembled fallback: pathway component effects (corresponds to requested 1d)\n- 1d in assembled fallback: validation composite (corresponds to requested 1e)\n\nThe preferred five-panel version can be made after uploading baseline-to-endpoint abundance effects.\n''')
Path(ROOT/'legends/Figure1_revised_legend_v2.md').write_text('''Figure 1. Genome-wide CRISPRi screening identifies regulators of endothelial differentiation.\n\n(a) Schematic of the genome-wide CRISPRi screen. Inducible CRISPRi hiPSCs were transduced with a genome-wide guide library, selected, induced, differentiated to the D3 endothelial endpoint, separated into CDH5-positive and CDH5-negative fractions, and analyzed by guide sequencing.\n\n(b) Full genome-wide gene-level screen result. Each point represents one harmonized biological gene. The x-axis shows MaGECK log2 fold change for CDH5-positive versus CDH5-negative fractions, and the y-axis shows the direction-specific MaGECK FDR. Blue points are required-for-EC-differentiation regulators; orange points are regulators whose knockdown restricts EC differentiation; gray points are non-hits. Open outlines mark genes supported by both MaGECK and U-test.\n\n(c) Selected pathway component effects. Each point is a tested component of a manuscript-relevant pathway group. Faint points are tested members; outlined points are final screen regulators. This panel is a selected component display, not a formal enrichment test.\n\n(d) Validation composite. The main scatter compares genome-wide and focused validation-screen effects. Inset shows validation replicate concordance. Bar plot shows validation rate by target group.\n''')
# packages
zip_path='/mnt/data/figure1_revised_v2_package.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for p in ROOT.rglob('*'):
        if p.is_file(): z.write(p, arcname=str(p.relative_to(ROOT)))
print('done', ROOT, zip_path)
