from pathlib import Path
import os, shutil, zipfile, textwrap
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import pearsonr

matplotlib.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 6.0,
    'axes.linewidth': 0.65,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'savefig.facecolor': 'white',
})

BASE = Path('/mnt/data')
SRC = BASE / 'Figure1_alternative_refined'
OUT = BASE / 'Figure1_alternative_v5'
for d in ['individual_panels','candidate_panels','assembled_figures','source_data','cleaned_data','data_audit','legends','scripts']:
    (OUT/d).mkdir(parents=True, exist_ok=True)

# Palette: restrained, color-blind-aware
BLUE = '#3B7EA1'
BLUE_DARK = '#27627E'
ORANGE = '#D98232'
GREY = '#D8DADD'
GREY_LIGHT = '#ECEDEE'
MID = '#6C747B'
DARK = '#202326'
PURPLE = '#76549A'
GREEN = '#688E55'

# Load existing audited/cleaned source data
u = pd.read_csv(SRC/'source_data/Fig1b_candidate_A_u_test_source.tsv', sep='\t')
m = pd.read_csv(SRC/'source_data/Fig1b_candidate_B_mageck_source.tsv', sep='\t')
r = pd.read_csv(SRC/'source_data/Fig1b_candidate_C_ranked_effect_source.tsv', sep='\t')
hits = pd.read_csv(SRC/'cleaned_data/screen_hits_unique_gene_level.tsv', sep='\t')
cat = pd.read_csv(SRC/'source_data/Fig1d_source_data.tsv', sep='\t')
val = pd.read_csv(SRC/'source_data/Fig1f_source_data.tsv', sep='\t')

for df in [u,m,r]:
    if 'display_gene_symbol' not in df.columns:
        df['display_gene_symbol'] = df['normalized_gene_symbol']

# Preserve source tables in new build
for p in (SRC/'source_data').glob('*'):
    if p.is_file(): shutil.copy2(p, OUT/'source_data'/p.name)
for p in (SRC/'cleaned_data').glob('*'):
    if p.is_file(): shutil.copy2(p, OUT/'cleaned_data'/p.name)
for p in (SRC/'data_audit').glob('*'):
    if p.is_file(): shutil.copy2(p, OUT/'data_audit'/p.name)

# ---------- general drawing helpers ----------
def clean_axis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.65)
    ax.spines['bottom'].set_linewidth(0.65)
    ax.tick_params(labelsize=5.5, width=0.6, length=2.3, pad=2)
    ax.grid(False)


def panel_spec(fig, spec, letter, title, subtitle=None, title_ratio=0.16):
    sub = spec.subgridspec(2, 1, height_ratios=[title_ratio, 1-title_ratio], hspace=0.02)
    tax = fig.add_subplot(sub[0,0])
    tax.axis('off')
    tax.text(0.0, 0.88, letter, fontsize=9.5, fontweight='bold', ha='left', va='top')
    tax.text(0.075, 0.88, title, fontsize=8.2, fontweight='bold', ha='left', va='top', color=DARK)
    if subtitle:
        tax.text(0.075, 0.30, subtitle, fontsize=5.5, ha='left', va='top', color=MID)
    return tax, sub[1,0]


def save_fig(fig, stem, width, height):
    fig.set_size_inches(width, height)
    for ext in ['png','pdf','svg']:
        fig.savefig(stem.with_suffix('.'+ext), dpi=500 if ext=='png' else None,
                    bbox_inches='tight', facecolor='white')
    plt.close(fig)

# ---------- Panel a: workflow ----------
def draw_workflow(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
    # starting CRISPRi state
    ax.add_patch(FancyBboxPatch((0.025,0.57),0.205,0.19,boxstyle='round,pad=0.012,rounding_size=0.018',
                                facecolor='#F4F0F8',edgecolor=PURPLE,linewidth=0.8))
    ax.text(0.1275,0.685,'CRISPRi hiPSCs',ha='center',va='center',fontsize=5.2,fontweight='bold',color=PURPLE)
    ax.text(0.1275,0.625,'genome-wide guide library',ha='center',va='center',fontsize=4.25,color=MID)

    # timeline
    y=0.51
    xs=[0.27,0.38,0.50,0.60,0.70,0.80]
    ax.add_patch(FancyArrowPatch((0.21,y),(0.82,y),arrowstyle='-|>',mutation_scale=8,linewidth=0.75,color='#7B7F83'))
    labels=[('Day -4','infect'),('Day -2','induce'),('D0','hiPSC'),('D1','primitive\nstreak'),('D2','lateral\nmesoderm'),('D3','EC')]
    cols=[MID,MID,'#777777','#9A834B','#6B8796','#4D8A7A']
    for x,(top,bottom),c in zip(xs,labels,cols):
        ax.add_patch(Circle((x,y),0.011,facecolor='white',edgecolor=c,linewidth=0.8,zorder=3))
        ax.text(x,y+0.055,top,ha='center',va='bottom',fontsize=4.7,fontweight='bold',color=DARK)
        ax.text(x,y-0.053,bottom,ha='center',va='top',fontsize=4.15,color=MID,linespacing=1.0)
    ax.text(0.325,0.405,'puromycin selection + doxycycline CRISPRi',ha='center',va='top',fontsize=4.0,color=MID)

    # D3 split
    ax.text(0.81,0.57,'MACS',fontsize=4.2,color=MID,ha='left')
    ec_box=(0.85,0.66,0.12,0.13); non_box=(0.85,0.35,0.12,0.13)
    for (x0,y0,w,h),lab,edge,fc in [(ec_box,'CDH5+\nEC',BLUE,'#EAF3F7'),(non_box,'CDH5-\nnon-EC',ORANGE,'#FBF0E6')]:
        ax.add_patch(FancyBboxPatch((x0,y0),w,h,boxstyle='round,pad=.008,rounding_size=.012',facecolor=fc,edgecolor=edge,linewidth=.8))
        ax.text(x0+w/2,y0+h/2,lab,ha='center',va='center',fontsize=5.0,fontweight='bold',color=edge,linespacing=.95)
    ax.add_patch(FancyArrowPatch((0.81,0.51),(0.85,0.715),arrowstyle='-|>',mutation_scale=8,linewidth=.8,color=BLUE,connectionstyle='arc3,rad=-.12'))
    ax.add_patch(FancyArrowPatch((0.81,0.51),(0.85,0.415),arrowstyle='-|>',mutation_scale=8,linewidth=.8,color=ORANGE,connectionstyle='arc3,rad=.12'))

    # sequencing and regulators
    ax.add_patch(FancyBboxPatch((0.70,0.10),0.27,0.14,boxstyle='round,pad=.01,rounding_size=.015',facecolor='#F7F7F7',edgecolor='#777777',linewidth=.7))
    ax.text(0.835,0.178,'guide sequencing',ha='center',va='center',fontsize=4.8,fontweight='bold',color=DARK)
    ax.text(0.835,0.132,'CDH5+ / CDH5- abundance',ha='center',va='center',fontsize=3.55,color=MID)
    for y0 in [0.66,0.35]:
        ax.add_patch(FancyArrowPatch((0.91,y0),(0.84,0.24),arrowstyle='-|>',mutation_scale=7,linewidth=.55,color='#777777',connectionstyle='arc3,rad=.12'))
    ax.add_patch(FancyBboxPatch((0.30,0.10),0.33,0.14,boxstyle='round,pad=.01,rounding_size=.015',facecolor='#EFF6F9',edgecolor=BLUE,linewidth=.8))
    ax.text(0.465,0.18,'244 unique regulators',ha='center',va='center',fontsize=5.0,fontweight='bold',color=BLUE_DARK)
    ax.text(0.465,0.132,'236 required  |  8 restricting',ha='center',va='center',fontsize=3.7,color=MID)
    ax.add_patch(FancyArrowPatch((0.70,0.17),(0.63,0.17),arrowstyle='-|>',mutation_scale=8,linewidth=.75,color=BLUE))
    ax.text(0.03,0.03,'129,638 guides; 20,528 genes; 6 sgRNAs per promoter',fontsize=4.0,color=MID,ha='left')

# ---------- Panel b: main screen ----------
def draw_screen(ax, source='mageck'):
    # use coherent MaGECK quantities
    df=m
    non=df[~df.final_union_hit]
    ax.scatter(non.mageck_lfc, non.mageck_minus_log10_fdr_plot, s=2.0, c=GREY, alpha=0.30, edgecolors='none', zorder=1)
    # hits by method and direction
    shapes={'MaGECK + U-test':'o','MaGECK only':'s','U-test only':'D'}
    for support,marker in shapes.items():
        for direction,color in [('Required for EC differentiation',BLUE),('Restricts EC differentiation',ORANGE)]:
            d=df[(df.final_union_hit)&(df.method_support==support)&(df.screen_direction==direction)]
            if d.empty: continue
            if support=='U-test only':
                ax.scatter(d.mageck_lfc,d.mageck_minus_log10_fdr_plot,s=18,marker=marker,facecolors='none',edgecolors=DARK,linewidths=.65,zorder=4)
                ax.scatter(d.mageck_lfc,d.mageck_minus_log10_fdr_plot,s=7,marker='o',color=color,edgecolors='none',zorder=5)
            else:
                ax.scatter(d.mageck_lfc,d.mageck_minus_log10_fdr_plot,s=9 if support=='MaGECK + U-test' else 10,
                           marker=marker,color=color,alpha=.90,edgecolors='white',linewidths=.2,zorder=3)
    ax.axhline(-np.log10(.05),color='#999999',lw=.55,ls='--',zorder=0)
    ax.axvline(0,color='#999999',lw=.55,ls='--',zorder=0)
    ax.set_xlabel('MaGECK log2 fold change',fontsize=6.0)
    ax.set_ylabel('-log10(minimum directional FDR)',fontsize=6.0)
    ax.set_ylim(-.15,4.85)
    ax.set_xlim(df.mageck_lfc.min()-.18, df.mageck_lfc.max()+.18)
    ax.text(.98,.96,'20,655 promoter targets\n244 unique regulator genes',transform=ax.transAxes,ha='right',va='top',fontsize=5.2,color=DARK)
    # fixed, well-separated top labels; no Hippo labels
    label_order=['KDR','SMAD1','BMPR1A','CDH5','PLCG1','ETV2']
    label_x={'KDR':-3.52,'SMAD1':-3.05,'BMPR1A':-2.70,'CDH5':-2.34,'PLCG1':-2.02,'ETV2':-1.70}
    for g in label_order:
        q=df[df.normalized_gene_symbol==g].iloc[0]
        lx=label_x[g]
        label_y={'KDR':4.56,'SMAD1':4.70,'BMPR1A':4.48,'CDH5':4.62,'PLCG1':4.43,'ETV2':4.70}
        ly=label_y[g]
        ax.annotate(g,(q.mageck_lfc,q.mageck_minus_log10_fdr_plot),xytext=(lx,ly),fontsize=5.05,fontstyle='italic',ha='center',va='center',
                    arrowprops=dict(arrowstyle='-',lw=.35,color='#777777',shrinkA=0,shrinkB=2),zorder=6)
    # concise shape legend inside lower-left
    handles=[Line2D([0],[0],marker='o',linestyle='none',markerfacecolor='#777777',markeredgecolor='white',markersize=3.7,label='both methods'),
             Line2D([0],[0],marker='s',linestyle='none',markerfacecolor='#777777',markeredgecolor='white',markersize=3.7,label='MaGECK only'),
             Line2D([0],[0],marker='D',linestyle='none',markerfacecolor='none',markeredgecolor=DARK,markersize=4,label='U-test only')]
    ax.legend(handles=handles,frameon=False,fontsize=4.4,loc='lower left',ncol=1,handletextpad=.2,labelspacing=.2,borderpad=.1)
    clean_axis(ax)

# Candidate A and C for comparison
def draw_utest(ax):
    df=u; non=df[~df.final_union_hit]
    ax.scatter(non.u_test_log2FC,non.u_test_minus_log10_fdr_plot,s=2,c=GREY,alpha=.3,edgecolors='none')
    for direction,color in [('Required for EC differentiation',BLUE),('Restricts EC differentiation',ORANGE)]:
        d=df[(df.final_union_hit)&(df.screen_direction==direction)]
        ax.scatter(d.u_test_log2FC,d.u_test_minus_log10_fdr_plot,s=8,color=color,alpha=.9,edgecolors='white',linewidths=.2)
    ut=df[df.method_support=='U-test only'];ax.scatter(ut.u_test_log2FC,ut.u_test_minus_log10_fdr_plot,s=18,facecolors='none',edgecolors=DARK,marker='D',linewidths=.6)
    ax.axhline(-np.log10(.05),color='#999',lw=.55,ls='--');ax.axvline(0,color='#999',lw=.55,ls='--')
    ax.set_xlabel('U-test log2 fold change',fontsize=6);ax.set_ylabel('-log10(adjusted P)',fontsize=6);ax.set_ylim(-.05,1.95);clean_axis(ax)


def draw_ranked(ax):
    df=r
    non=df[~df.final_union_hit]; req=df[(df.final_union_hit)&(df.screen_direction=='Required for EC differentiation')]; res=df[(df.final_union_hit)&(df.screen_direction=='Restricts EC differentiation')]
    ax.scatter(non['rank'],non.u_test_log2FC,s=1.3,color=GREY,alpha=.35,edgecolors='none')
    ax.scatter(req['rank'],req.u_test_log2FC,s=6,color=BLUE,alpha=.9,edgecolors='white',linewidths=.15)
    ax.scatter(res['rank'],res.u_test_log2FC,s=9,color=ORANGE,alpha=.95,edgecolors='white',linewidths=.2)
    ax.axhline(0,color='#999',lw=.55,ls='--');ax.set_xlabel('Ranked promoter targets',fontsize=6);ax.set_ylabel('U-test log2 fold change',fontsize=6);clean_axis(ax)

# ---------- Panel c: bidirectional ranked hits ----------
HIPPO={'NF2','LATS1','LATS2','SAV1','AMOTL1','AMOTL2','TAOK1','TAOK2','TAOK3','WWC1','WWC2','WWC3','PTPN14'}
req_all=hits[hits.direction=='Required for EC differentiation'].sort_values('effect')
req=req_all[~req_all.display_gene_symbol.isin(HIPPO)].head(20).copy()
res=hits[hits.direction=='Restricts EC differentiation'].sort_values('effect',ascending=False).copy()
req['selection_reason']='top 20 non-Hippo required hits by U-test effect'
res['selection_reason']='all genes in smaller hit direction'
sel_c=pd.concat([req,res],ignore_index=True)
sel_c.to_csv(OUT/'source_data'/'Fig1c_displayed_gene_selection_v5.tsv',sep='\t',index=False)
hits.to_csv(OUT/'source_data'/'Fig1c_source_data_v5.tsv',sep='\t',index=False)

MARK={'MaGECK + U-test':'o','MaGECK only':'s','U-test only':'D'}
def draw_lollipop(ax,df,color,left=True):
    # put strongest at bottom to read from weak-to-strong vertically
    dd=df.sort_values('effect',ascending=False if left else True).reset_index(drop=True)
    yy=np.arange(len(dd))
    for i,q in dd.iterrows():
        ax.hlines(i,0,q.effect,color=color,lw=.72,zorder=1)
        ax.scatter(q.effect,i,s=16,color=color,edgecolors='white',linewidths=.25,marker=MARK.get(q.method_support,'o'),zorder=2)
    ax.set_yticks(yy); ax.set_yticklabels(dd.display_gene_symbol,fontsize=5.0,fontstyle='italic')
    ax.axvline(0,color='#777777',lw=.55);ax.tick_params(axis='y',length=0,pad=2)
    ax.set_xlabel('U-test effect',fontsize=5.7)
    clean_axis(ax)
    if left: ax.set_xlim(min(-3.45,dd.effect.min()-.12),.10)
    else: ax.set_xlim(-.05,max(2.55,dd.effect.max()+.15))
    return dd


def draw_panel_c(fig,spec,letter='c'):
    tax,content=panel_spec(fig,spec,letter,'Bidirectional ranked-hit map','Top non-Hippo required hits by U-test effect; all 8 restricting hits',title_ratio=.22)
    sub=content.subgridspec(1,2,wspace=.58,width_ratios=[1.22,1.0])
    a1=fig.add_subplot(sub[0,0]);a2=fig.add_subplot(sub[0,1])
    draw_lollipop(a1,req,BLUE,True);draw_lollipop(a2,res,ORANGE,False)
    a1.set_title('Required for EC differentiation',fontsize=5.7,fontweight='bold',color=BLUE_DARK,loc='left',pad=2)
    a2.set_title('Restricts EC differentiation',fontsize=5.7,fontweight='bold',color=ORANGE,loc='left',pad=2)
    # method support legend in title band
    handles=[Line2D([0],[0],marker=MARK[k],linestyle='none',markerfacecolor='#777777',markeredgecolor='white',markersize=3.5,label=k) for k in MARK]
    tax.legend(handles=handles,loc='upper right',ncol=3,frameon=False,fontsize=4.2,handletextpad=.15,columnspacing=.55,borderaxespad=0.0,bbox_to_anchor=(1,.96))
    return tax,a1,a2

# ---------- Panel d ----------
SHORT={
'Transcription & chromatin regulation':'Transcription / chromatin',
'RNA processing & surveillance':'RNA processing / surveillance',
'Developmental signaling & lineage regulation':'Developmental signaling / lineage',
'Proteostasis & post-translational modification':'Proteostasis / PTM',
'Membrane trafficking, cytoskeleton & adhesion':'Trafficking / cytoskeleton / adhesion',
'Metabolism, mitochondria & glycosylation':'Metabolism / mitochondria / glycosylation',
'Translation, ribosome & protein targeting':'Translation / protein targeting',
'Cell cycle, chromosome & DNA repair':'Cell cycle / DNA repair',
'Apoptosis & stress response':'Apoptosis / stress',
'Other / insufficiently characterized':'Other / insufficiently characterized'
}
cat_plot=cat.copy();cat_plot['short_category']=cat_plot.primary_category.map(lambda x:SHORT.get(x,x));cat_plot=cat_plot.sort_values('total',ascending=True)
cat_plot.to_csv(OUT/'source_data'/'Fig1d_source_data_v5.tsv',sep='\t',index=False)

def draw_panel_d(ax):
    yy=np.arange(len(cat_plot))
    left=-cat_plot['Required for EC differentiation'].astype(float)
    right=cat_plot['Restricts EC differentiation'].astype(float)
    ax.barh(yy,left,color=BLUE,height=.64,label='Required')
    ax.barh(yy,right,color=ORANGE,height=.64,label='Restricts')
    for i,l,rr in zip(yy,-left,right):
        if l>0: ax.text(-l-.8,i,str(int(l)),ha='right',va='center',fontsize=4.8,color=BLUE_DARK,fontweight='bold')
        if rr>0: ax.text(rr+.22,i,str(int(rr)),ha='left',va='center',fontsize=4.8,color=ORANGE,fontweight='bold')
    ax.axvline(0,color='#777777',lw=.55)
    ax.set_yticks(yy);ax.set_yticklabels(cat_plot.short_category,fontsize=4.75)
    ax.tick_params(axis='y',length=0,pad=2)
    ax.set_xlim(-68,10);ax.set_xlabel('Number of regulators',fontsize=5.7)
    ax.legend(frameon=False,fontsize=4.5,loc='lower right',handlelength=1.2,labelspacing=.2)
    clean_axis(ax)

# ---------- Panel e placeholder ----------
def draw_panel_e(ax):
    ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
    items=[(.04,.55,.92,.30,'A  Stage-biased expression','Unperturbed/control RNA time course','gene x day x biological replicate'),
           (.04,.13,.92,.30,'B  Differentiation specificity versus broad abundance','D0-to-D3 guide or target abundance','gene-level abundance effect and significance')]
    for x,y,w,h,tt,sub,fields in items:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.012,rounding_size=.018',fc='#F7F7F7',ec='#B8B8B8',lw=.65))
        ax.text(x+.03,y+h-.07,tt,fontsize=5.9,fontweight='bold',ha='left',va='top',color=DARK)
        ax.text(x+.03,y+h-.16,sub,fontsize=5.2,ha='left',va='top',color=MID)
        ax.text(x+.03,y+.06,fields,fontsize=4.35,ha='left',va='bottom',color=MID)
    ax.text(.5,.02,'No expression-timing or broad-fitness conclusion was inferred from unrelated data.',ha='center',va='bottom',fontsize=4.6,color=MID,style='italic')

# ---------- Panel f validation ----------
TYPES={'genome_wide_hit':(BLUE,'screen-derived hit'),'genome_wide_non_hit':('#B9B9B9','random non-hit'),'tri_perturb_seq':(ORANGE,'additional curated')}
VAL_LABS=['KDR','CDH5','ETV2','PLCG1','BMPR1A','SMAD1']
VAL_OFF={'KDR':(-.32,.05),'CDH5':(-.12,.15),'ETV2':(.22,.02),'PLCG1':(.12,-.12),'BMPR1A':(.14,-.09),'SMAD1':(-.18,-.10)}

def draw_panel_f(ax):
    for typ,(col,lab) in TYPES.items():
        d=val[val.Type==typ]
        sig=d.validation_screen_hit.fillna(False).astype(bool)
        ax.scatter(d.loc[~sig,'Original FC Mageck'],d.loc[~sig,'Average U Test FC'],s=7,facecolors='white',edgecolors=col,linewidths=.45,alpha=.78,zorder=1)
        ax.scatter(d.loc[sig,'Original FC Mageck'],d.loc[sig,'Average U Test FC'],s=9,color=col,edgecolors='white',linewidths=.18,alpha=.88,label=lab,zorder=2)
    mm=val.dropna(subset=['Original FC Mageck','Average U Test FC'])
    coef=np.polyfit(mm['Original FC Mageck'],mm['Average U Test FC'],1)
    xx=np.linspace(mm['Original FC Mageck'].min(),mm['Original FC Mageck'].max(),100)
    ax.plot(xx,coef[0]*xx+coef[1],color='black',lw=.75,zorder=0)
    ax.axhline(0,color='#999999',lw=.5,ls='--');ax.axvline(0,color='#999999',lw=.5,ls='--')
    rr,pp=pearsonr(mm['Original FC Mageck'],mm['Average U Test FC'])
    for _,q in val[val.Gene.isin(VAL_LABS)].iterrows():
        if pd.isna(q['Original FC Mageck']): continue
        dx,dy=VAL_OFF[q.Gene]
        ax.annotate(q.Gene,(q['Original FC Mageck'],q['Average U Test FC']),xytext=(q['Original FC Mageck']+dx,q['Average U Test FC']+dy),fontsize=4.7,fontstyle='italic',ha='center',
                    arrowprops=dict(arrowstyle='-',lw=.28,color='#777777',shrinkA=0,shrinkB=2))
    ax.text(.03,.96,f'n = {len(mm)} targets\nPearson r = {rr:.2f}\n190/244 screen hits significant\nin the same direction',transform=ax.transAxes,ha='left',va='top',fontsize=4.8)
    ax.set_xlabel('Genome-wide MaGECK effect',fontsize=5.7);ax.set_ylabel('Validation U-test effect',fontsize=5.7)
    clean_axis(ax)
    # replicate inset
    ins=inset_axes(ax,width='36%',height='36%',loc='lower right',borderpad=.65)
    rep=val.dropna(subset=['U Test FC BR1','U Test FC BR2'])
    ins.scatter(rep['U Test FC BR1'],rep['U Test FC BR2'],s=3,color='#777777',alpha=.45,edgecolors='none')
    rrep,_=pearsonr(rep['U Test FC BR1'],rep['U Test FC BR2'])
    cf=np.polyfit(rep['U Test FC BR1'],rep['U Test FC BR2'],1);xs=np.linspace(rep['U Test FC BR1'].min(),rep['U Test FC BR1'].max(),100)
    ins.plot(xs,cf[0]*xs+cf[1],color='black',lw=.5);ins.axhline(0,color='#aaa',lw=.3,ls='--');ins.axvline(0,color='#aaa',lw=.3,ls='--')
    ins.set_title(f'Replicates r = {rrep:.2f}',fontsize=4.3,pad=2);ins.set_xlabel('BR1',fontsize=3.8,labelpad=1);ins.set_ylabel('BR2',fontsize=3.8,labelpad=1)
    ins.tick_params(labelsize=3.5,length=1,pad=1);ins.spines['top'].set_visible(False);ins.spines['right'].set_visible(False)
    ax.legend(frameon=False,fontsize=4.2,loc='upper right',bbox_to_anchor=(1,.81),handletextpad=.18,labelspacing=.2)

# ---------- Individual panels ----------
def individual_panel(letter,title,subtitle,draw_func,stem,width,height, title_ratio=.16):
    fig=plt.figure(figsize=(width,height),facecolor='white')
    gs=fig.add_gridspec(1,1,left=.03,right=.98,bottom=.06,top=.98)
    tax,content=panel_spec(fig,gs[0,0],letter,title,subtitle,title_ratio=title_ratio)
    ax=fig.add_subplot(content)
    draw_func(ax)
    save_fig(fig,OUT/'individual_panels'/stem,width,height)

individual_panel('a','Genome-wide CRISPRi screen','Endothelial differentiation and CDH5-based endpoint sorting',draw_workflow,'Fig1a_workflow',3.0,2.35,.18)
individual_panel('b','Genome-wide screen result','Final union hits highlighted; point shape indicates hit-calling support',draw_screen,'Fig1b_genomewide_screen',4.6,3.25,.15)
# panel c custom
fig=plt.figure(figsize=(5.0,3.55),facecolor='white');gs=fig.add_gridspec(1,1,left=.03,right=.98,bottom=.07,top=.98);draw_panel_c(fig,gs[0,0],'c');save_fig(fig,OUT/'individual_panels'/'Fig1c_bidirectional_ranked_hits',5.0,3.55)
individual_panel('d','Functional organization','Literature-informed primary categories for 244 genes',draw_panel_d,'Fig1d_functional_organization',3.6,3.25,.16)
individual_panel('e','Biological context - input required','One of two prespecified source datasets is needed',draw_panel_e,'Fig1e_biological_context_input_required',4.7,2.45,.18)
individual_panel('f','Secondary validation','Genome-wide versus focused-screen effects',draw_panel_f,'Fig1f_validation_reproducibility',3.5,2.75,.17)

# candidate panels b
individual_panel('','U-test screen result','Coherent effect and adjusted P value',draw_utest,'Fig1b_candidate_A_u_test',4.6,3.25,.15)
individual_panel('','MaGECK screen result','Coherent effect and directional FDR',draw_screen,'Fig1b_candidate_B_mageck',4.6,3.25,.15)
individual_panel('','Ranked genome-wide effects','All promoter targets; final union hits highlighted',draw_ranked,'Fig1b_candidate_C_ranked',4.6,3.25,.15)
for fn in (OUT/'individual_panels').glob('Fig1b_candidate_*'):
    shutil.copy2(fn,OUT/'candidate_panels'/fn.name)

# ---------- Assembled figures ----------
def assemble_main(stem,compact=False):
    if not compact:
        fig=plt.figure(figsize=(7.2,8.05),facecolor='white')
        gs=fig.add_gridspec(3,10,left=.035,right=.985,bottom=.035,top=.985,hspace=.30,wspace=.34,height_ratios=[1.00,1.27,.95])
        # row1
        tax_a,cont_a=panel_spec(fig,gs[0,:4],'a','Genome-wide CRISPRi screen','Endothelial differentiation and CDH5-based endpoint sorting',.16);ax_a=fig.add_subplot(cont_a);draw_workflow(ax_a)
        tax_b,cont_b=panel_spec(fig,gs[0,4:],'b','Genome-wide screen result','Final union hits highlighted; shape indicates method support',.14);ax_b=fig.add_subplot(cont_b);draw_screen(ax_b)
        # row2
        draw_panel_c(fig,gs[1,:6],'c')
        tax_d,cont_d=panel_spec(fig,gs[1,6:],'d','Functional organization','Literature-informed categories for 244 genes',.15);ax_d=fig.add_subplot(cont_d);draw_panel_d(ax_d)
        # row3
        tax_e,cont_e=panel_spec(fig,gs[2,:6],'e','Biological context - input required','One of two prespecified source datasets is needed',.17);ax_e=fig.add_subplot(cont_e);draw_panel_e(ax_e)
        tax_f,cont_f=panel_spec(fig,gs[2,6:],'f','Secondary validation','Genome-wide versus focused-screen effects',.16);ax_f=fig.add_subplot(cont_f);draw_panel_f(ax_f)
    else:
        # compact alternative: c and d aligned, biological context remains full-width bottom-left
        fig=plt.figure(figsize=(7.2,7.5),facecolor='white')
        gs=fig.add_gridspec(3,10,left=.035,right=.985,bottom=.035,top=.985,hspace=.30,wspace=.34,height_ratios=[1.0,1.08,.88])
        tax_a,cont_a=panel_spec(fig,gs[0,:4],'a','Genome-wide CRISPRi screen','Endothelial differentiation and CDH5-based endpoint sorting',.16);ax_a=fig.add_subplot(cont_a);draw_workflow(ax_a)
        tax_b,cont_b=panel_spec(fig,gs[0,4:],'b','Genome-wide screen result','Final union hits highlighted; shape indicates method support',.14);ax_b=fig.add_subplot(cont_b);draw_screen(ax_b)
        # compressed c: rank strips without labels + full d
        tax_c,cont_c=panel_spec(fig,gs[1,:5],'c','Ranked screen hits','Top required hits and all restricting hits',.16)
        sub=cont_c.subgridspec(1,2,wspace=.34);a1=fig.add_subplot(sub[0,0]);a2=fig.add_subplot(sub[0,1])
        dd1=req.sort_values('effect');a1.scatter(dd1.effect,np.arange(len(dd1)),s=9,color=BLUE);a1.axvline(0,color='#777',lw=.5);a1.set_yticks([]);a1.set_xlabel('Required effect',fontsize=5.5);clean_axis(a1)
        dd2=res.sort_values('effect');a2.scatter(dd2.effect,np.arange(len(dd2)),s=11,color=ORANGE);a2.axvline(0,color='#777',lw=.5);a2.set_yticks([]);a2.set_xlabel('Restricting effect',fontsize=5.5);clean_axis(a2)
        tax_d,cont_d=panel_spec(fig,gs[1,5:],'d','Functional organization','Literature-informed categories for 244 genes',.15);ax_d=fig.add_subplot(cont_d);draw_panel_d(ax_d)
        tax_e,cont_e=panel_spec(fig,gs[2,:6],'e','Biological context - input required','One of two prespecified source datasets is needed',.17);ax_e=fig.add_subplot(cont_e);draw_panel_e(ax_e)
        tax_f,cont_f=panel_spec(fig,gs[2,6:],'f','Secondary validation','Genome-wide versus focused-screen effects',.16);ax_f=fig.add_subplot(cont_f);draw_panel_f(ax_f)
    for ext in ['png','pdf','svg']:
        fig.savefig(OUT/'assembled_figures'/f'{stem}.{ext}',dpi=500 if ext=='png' else None,bbox_inches='tight',facecolor='white')
    plt.close(fig)

assemble_main('Figure1_alternative_provisional_v5',False)
assemble_main('Figure1_alternative_compact_v5',True)

# ---------- README, legend, missing-input request ----------
readme=f"""# Alternative Figure 1 v5

This build implements the data-rich alternative Figure 1 specified in `Pasted text(3).txt`.

## Completed panels

- 1a genome-wide CRISPRi workflow
- 1b coherent MaGECK effect-versus-directional-FDR result plot
- 1c bidirectional ranked-hit map, excluding Hippo-pathway genes from the displayed required-hit subset
- 1d literature-informed functional organization of 244 unique screen-derived regulators
- 1f genome-wide versus focused-validation effects with replicate inset

## Panel 1e status

Panel 1e is blocked because neither required biological-context dataset was present. The provisional assembly contains a labeled input-request panel rather than an inferred biological conclusion.

## Main 1b selection

The MaGECK panel was selected because both axes come from the same MaGECK analysis, the full tested set is displayed, U-test-only union hits can be shown with outlined diamonds, and the panel complements the U-test-ranked panel 1c. U-test volcano and ranked-effect alternatives are exported in `candidate_panels`.

## Counts recalculated from current files

- 20,655 MaGECK promoter-target records
- 20,480 U-test promoter-target records
- 245 final hit target records representing 244 unique biological genes
- 88 hits supported by both methods
- 151 MaGECK-only hits
- 5 U-test-only hits
- 236 required for EC differentiation
- 8 restricting EC differentiation
- validation library: 297 biological genes
- 193/244 screen-derived regulators significant in validation
- 190/244 significant in the same direction

## Vector text

SVG text is preserved as text (`svg.fonttype = none`). PDF fonts are embedded as TrueType (`pdf.fonttype = 42`).
"""
(OUT/'README_figure1_alternative_v5.md').write_text(readme)

missing="""# Missing input for panel 1e

Upload either one of the following source datasets.

## Option A: stage-biased expression

Preferred file: unperturbed/control RNA time course, ideally the uniformly processed multiome RNA modality.

Required columns:
- `gene`
- `day` (D0, D1, D2, D3; D4 optional)
- `biological_replicate`
- `expression_value`
- `expression_unit` (for example pseudobulk normalized expression, mean log-normalized expression, or fraction detected)

This blocks panel 1e-A. A provisional placeholder is included; no timing conclusion can be made without this file.

## Option B: differentiation specificity versus broad abundance

Preferred file: D0-to-D3 guide- or target-abundance analysis.

Required columns:
- `gene`
- `D0_to_D3_abundance_effect`
- `abundance_p_value` and/or `abundance_fdr`
- biological replicate fields, when available

The differentiation-screen effect can be joined from the existing MaGECK/U-test tables. This blocks panel 1e-B. A provisional placeholder is included; no broad-fitness conclusion can be made without this file.

Only one of Option A or Option B is needed to complete the main alternative assembly.
"""
(OUT/'data_audit'/'consolidated_missing_input_request_v5.md').write_text(missing)

legend="""# Draft legend: alternative Figure 1

**Figure 1. Genome-wide CRISPRi screening identifies validated regulators of endothelial differentiation and organizes them by effect and functional context.**

**a,** Schematic of the genome-wide CRISPRi differentiation screen. Inducible CRISPRi hiPSCs were infected with the genome-wide guide library, selected and induced before differentiation through D0 hiPSC, D1 primitive streak, D2 lateral mesoderm and D3 endothelial stages. D3 cells were separated into CDH5-positive EC and CDH5-negative non-EC populations, and guide abundance was quantified to identify 244 unique regulators.

**b,** MaGECK log2 fold change versus -log10 minimum directional FDR for 20,655 promoter targets. Final union hits are colored by biological direction. Point shape indicates support by MaGECK and U-test, MaGECK only or U-test only. Selected established endothelial and developmental regulators are labeled.

**c,** Bidirectional ranked-hit map showing the 20 strongest non-Hippo genes required for EC differentiation by U-test effect and all eight genes whose knockdown restricts EC differentiation. Point shape indicates hit-calling support. The displayed subset does not alter the full counts of 236 required and 8 restricting regulators.

**d,** Functional organization of the 244 unique screen-derived regulators using literature-informed primary categories. Bars show the number of regulators required for EC differentiation and the number that restrict EC differentiation.

**e,** Placeholder identifying the source data required for a biological-context panel. No expression-timing or broad-abundance conclusion was inferred in the absence of the prespecified source data.

**f,** Comparison of genome-wide MaGECK effects and focused-validation U-test effects. Filled points denote significant validation effects; colors indicate target-design group. The inset shows validation replicate concordance. Statistics were calculated from harmonized biological targets.
"""
(OUT/'legends'/'Figure1_alternative_legend_v5.md').write_text(legend)

# candidate comparison
pd.DataFrame([
 {'candidate':'A U-test volcano','selected':False,'reason':'Statistically coherent but adjusted-P dynamic range is compressed and the panel duplicates U-test ranking in 1c.'},
 {'candidate':'B MaGECK result','selected':True,'reason':'Statistically coherent, full tested set shown, union hits and U-test-only support are legible, and complements 1c.'},
 {'candidate':'C ranked effect','selected':False,'reason':'Very readable but too visually similar to the ranked-hit logic of panel 1c.'}
]).to_csv(OUT/'data_audit'/'Fig1b_candidate_comparison_v5.tsv',sep='\t',index=False)

# copy script
shutil.copy2(Path(__file__),OUT/'scripts'/'build_figure1_alternative_v5.py')

# package
zip_path=BASE/'Figure1_alternative_v5_package.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for fp in OUT.rglob('*'):
        if fp.is_file(): z.write(fp,fp.relative_to(OUT))
print('Built:',OUT)
print('Package:',zip_path)
