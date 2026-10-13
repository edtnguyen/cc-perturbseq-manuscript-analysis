import os, re, math, json, shutil, zipfile, base64
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D
from scipy.stats import pearsonr
from plotnine import (
    ggplot, aes, geom_point, geom_segment, geom_hline, geom_vline,
    geom_text, geom_col, coord_flip, geom_abline, geom_smooth,
    scale_color_manual, scale_fill_manual, scale_shape_manual,
    scale_size_manual, labs, theme_classic, theme, element_text,
    element_blank, element_line, element_rect, guides, guide_legend,
    facet_wrap, coord_cartesian, annotate
)
import cairosvg
import xml.etree.ElementTree as ET

BASE = Path('/mnt/data')
OUT = BASE/'Figure1_nature_style_remade'
PANELS = OUT/'individual_panels'
ASSEMBLED = OUT/'assembled_figures'
SOURCE = OUT/'source_data'
SCRIPTS = OUT/'scripts'
for d in [PANELS, ASSEMBLED, SOURCE, SCRIPTS]: d.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Palette / style
# -----------------------------
BLUE = '#2F6F9F'       # required / EC-related
ORANGE = '#D9772A'     # restrict / alternative
PURPLE = '#74539A'
TEAL = '#3B8D88'
GREEN = '#5F8B4C'
GRAY = '#C9CDD0'
DARK = '#22272B'
MID = '#66707A'
LIGHT = '#F5F6F7'

mpl.rcParams.update({
    'font.family':'DejaVu Sans',
    'font.size':8,
    'pdf.fonttype':42,
    'ps.fonttype':42,
    'svg.fonttype':'none',
})

# -----------------------------
# Helpers
# -----------------------------
def parse_gene(x):
    s = str(x)
    if '_-' in s:
        return s.split('_-')[0]
    if re.search(r'-P\d', s):
        return re.split(r'-P\d', s)[0]
    return s.split('_')[0]

def save_mpl(fig, stem, width=None, height=None):
    if width and height:
        fig.set_size_inches(width, height)
    for ext in ['png','pdf','svg']:
        fig.savefig(PANELS/f'{stem}.{ext}', dpi=500 if ext=='png' else None,
                    bbox_inches='tight', facecolor='white')
    plt.close(fig)

def save_plotnine(p, stem, width, height):
    svg = PANELS/f'{stem}.svg'
    pdf = PANELS/f'{stem}.pdf'
    png = PANELS/f'{stem}.png'
    p.save(svg, width=width, height=height, verbose=False)
    p.save(png, width=width, height=height, dpi=400, verbose=False)
    cairosvg.svg2pdf(url=str(svg), write_to=str(pdf))

# -----------------------------
# Load and clean data
# -----------------------------
full = pd.read_excel(BASE/'posneg_unpacked/Pos_vs_Neg/combined_mageck_u_test_pos_v_neg.xlsx')
ann = pd.read_csv(BASE/'figure_build/functional_annotation/screen_regulators_244_literature_annotation_FINAL_DRAFT.tsv', sep='\t')
known = pd.read_csv(BASE/'figure_build/functional_annotation/known_EC_regulators_literature_review_FINAL_DRAFT.tsv', sep='\t')
val = pd.read_csv(BASE/'validation_unpacked/validation screen results/validation_screen_d3pos_d3neg_mageckLFC.csv')
val_type = pd.read_excel(BASE/'validation_unpacked/validation screen results/validation_genes_type.xlsx')
val = val.merge(val_type, left_on='Gene', right_on='Target', how='left')

alias_path = BASE/'figure_build/cleaned_data/gene_alias_mapping.tsv'
alias = pd.read_csv(alias_path, sep='\t') if alias_path.exists() else pd.DataFrame(columns=['original_gene','gene_mapped'])
alias_map = dict(zip(alias.original_gene, alias.gene_mapped))

full['original_gene'] = full['id'].map(parse_gene)
full['gene'] = full['original_gene'].map(lambda g: alias_map.get(g,g))
full['mageck_lfc'] = full['neg|lfc'].combine_first(full['pos|lfc'])
full['mageck_min_fdr'] = full[['neg|fdr','pos|fdr']].min(axis=1, skipna=True)
full['minus_log10_mageck_fdr'] = -np.log10(full['mageck_min_fdr'].clip(lower=1e-300))

ann_join = ann[['gene','final_union_hit','screen_direction','hit_calling_method','recommended_must_label','recommended_label_priority','known_EC_regulator_tier','figure1d_category']].copy()
full = full.merge(ann_join, on='gene', how='left')
full['final_union_hit'] = full['final_union_hit'].fillna(False).astype(bool)
full['screen_direction'] = full['screen_direction'].fillna('Non-hit')
full['hit_calling_method'] = full['hit_calling_method'].fillna('Non-hit')
full['point_group'] = np.where(~full.final_union_hit, 'Non-hit', full.screen_direction)
full['method_group'] = np.where(~full.final_union_hit, 'Non-hit', full.hit_calling_method)

full.to_csv(SOURCE/'Fig1b_full_screen_source.tsv', sep='\t', index=False)
ann.to_csv(SOURCE/'Fig1d_annotation_source.tsv', sep='\t', index=False)
val.to_csv(SOURCE/'Fig1e_validation_source.tsv', sep='\t', index=False)
known.to_csv(SOURCE/'Fig1c_known_regulators_source.tsv', sep='\t', index=False)

# -----------------------------
# Panel A — workflow schematic (vector)
# -----------------------------
def build_panel_a():
    fig, ax = plt.subplots(figsize=(4.05,3.05))
    ax.set_xlim(0, 12); ax.set_ylim(0, 8); ax.axis('off')
    ax.text(0.05,7.72,'Genome-wide CRISPRi screen design', fontsize=11.3, fontweight='bold', color=DARK, va='top')
    ax.text(0.05,7.20,'129,638 guides targeting 20,528 genes', fontsize=7.7, color=MID, va='top')

    # input cells/library
    cell_box = FancyBboxPatch((0.35,4.72),2.15,1.35,boxstyle='round,pad=0.12,rounding_size=0.14',
                              facecolor='#F7F7F7',edgecolor='#7B848C',linewidth=0.9)
    ax.add_patch(cell_box)
    ax.text(1.42,5.66,'Doxycycline-inducible',ha='center',fontsize=7.4,fontweight='bold',color=DARK)
    ax.text(1.42,5.34,'CRISPRi hiPSCs',ha='center',fontsize=8.3,fontweight='bold',color=BLUE)
    ax.text(1.42,4.96,'genome-wide guide library',ha='center',fontsize=6.8,color=MID)
    for j in range(5):
        x=0.75+0.32*j
        ax.plot([x,x+0.12],[4.58,4.42],color=PURPLE,lw=1.1)
        ax.add_patch(Circle((x-0.02,4.60),0.045,color=PURPLE))

    # timeline
    y=3.48
    xs=[3.05,4.20,5.35,6.50,7.65,8.80]
    labs=[('Day −4','infect'),('Day −2','induce'),('D0','hiPSC'),('D1','primitive\nstreak'),('D2','lateral\nmesoderm'),('D3','EC endpoint')]
    ax.plot([xs[0],xs[-1]],[y,y],color='#626A72',lw=1.3)
    for i,(x,(day,lab)) in enumerate(zip(xs,labs)):
        fc=['#ECECEC','#ECECEC','#E9DFF0','#F4E4C5','#DDEBF4','#DCEFE6'][i]
        ax.add_patch(Circle((x,y),0.16,facecolor=fc,edgecolor='#626A72',lw=0.9,zorder=3))
        ax.text(x,y+0.55,day,ha='center',fontsize=7.2,fontweight='bold',color=DARK)
        ax.text(x,y-0.55,lab,ha='center',va='top',fontsize=6.3,color=MID)
    ax.annotate('',xy=(9.20,y),xytext=(8.80,y),arrowprops=dict(arrowstyle='-|>',lw=1.1,color='#626A72'))

    # MACS split and sequencing
    ax.text(9.35,4.24,'MACS',ha='center',fontsize=6.8,fontweight='bold',color=MID)
    ax.annotate('',xy=(10.08,5.05),xytext=(8.95,3.53),arrowprops=dict(arrowstyle='-|>',lw=1.3,color=BLUE,connectionstyle='arc3,rad=-0.15'))
    ax.annotate('',xy=(10.08,2.22),xytext=(8.95,3.43),arrowprops=dict(arrowstyle='-|>',lw=1.3,color=ORANGE,connectionstyle='arc3,rad=0.15'))
    for yy,lab,fc,ec in [(5.05,'CDH5+\nEC','#E7F2F8',BLUE),(2.22,'CDH5−\nnon-EC','#FBEDE5',ORANGE)]:
        b=FancyBboxPatch((10.08,yy-0.48),1.28,0.96,boxstyle='round,pad=0.08,rounding_size=0.12',facecolor=fc,edgecolor=ec,lw=1.0)
        ax.add_patch(b); ax.text(10.72,yy,lab,ha='center',va='center',fontsize=7.4,fontweight='bold',color=ec)
    seq=FancyBboxPatch((9.58,0.36),2.0,0.95,boxstyle='round,pad=0.08,rounding_size=0.12',facecolor='white',edgecolor='#7B848C',lw=0.9)
    ax.add_patch(seq); ax.text(10.58,0.94,'guide sequencing',ha='center',fontsize=7.4,fontweight='bold',color=DARK)
    ax.text(10.58,0.59,'244 unique regulators',ha='center',fontsize=7.4,fontweight='bold',color=PURPLE)
    ax.annotate('',xy=(10.58,1.31),xytext=(10.70,1.75),arrowprops=dict(arrowstyle='-|>',lw=0.9,color='#6F777E'))
    ax.annotate('',xy=(10.58,1.31),xytext=(10.70,4.55),arrowprops=dict(arrowstyle='-|>',lw=0.9,color='#6F777E',connectionstyle='arc3,rad=.18'))
    ax.text(0.35,0.55,'6 sgRNAs/promoter  •  2,680 safe-targeting  •  3,790 non-targeting',fontsize=6.1,color=MID)
    save_mpl(fig,'Fig1a_workflow_nature_style',4.05,3.05)

# -----------------------------
# Panel B — full MaGECK volcano
# -----------------------------
def build_panel_b():
    nonhits=full[~full.final_union_hit].copy()
    hits=full[full.final_union_hit].copy()
    label_genes=['KDR','CDH5','ETV2','PLCG1','BMPR1A','SMAD1','SMAD5','NF2']
    label_df=(hits[hits.gene.isin(label_genes)]
              .sort_values('mageck_min_fdr')
              .drop_duplicates('gene'))
    offsets={
        'KDR':(-0.22,0.18),'CDH5':(0.12,0.14),'ETV2':(0.10,-0.18),'PLCG1':(-0.10,-0.16),
        'BMPR1A':(-0.10,0.18),'SMAD1':(0.12,0.12),'SMAD5':(-0.18,-0.20),'NF2':(0.16,0.16)
    }
    label_df=label_df.copy()
    label_df['label_x']=[r.mageck_lfc+offsets.get(r.gene,(0.1,0.1))[0] for r in label_df.itertuples()]
    label_df['label_y']=[r.minus_log10_mageck_fdr+offsets.get(r.gene,(0.1,0.1))[1] for r in label_df.itertuples()]
    label_df['label']=label_df['gene']

    p=(ggplot()
       + geom_point(nonhits,aes('mageck_lfc','minus_log10_mageck_fdr'),color='#D7D9DB',size=0.8,alpha=0.35,raster=True)
       + geom_point(hits,aes('mageck_lfc','minus_log10_mageck_fdr',fill='point_group',shape='method_group'),
                    color='#3F464C',stroke=0.25,size=2.05,alpha=0.92)
       + geom_hline(yintercept=-np.log10(0.05),linetype='dashed',color='#9AA0A5',size=0.45)
       + geom_vline(xintercept=0,linetype='dashed',color='#9AA0A5',size=0.45)
       + geom_segment(label_df,aes(x='mageck_lfc',y='minus_log10_mageck_fdr',xend='label_x',yend='label_y'),color='#737A80',size=0.35)
       + geom_text(label_df,aes('label_x','label_y',label='label'),size=7.2,fontstyle='italic',color=DARK)
       + scale_fill_manual(values={'Required for EC differentiation':BLUE,'Restricts EC differentiation':ORANGE,'Non-hit':'#D7D9DB'},
                           name='Screen direction')
       + scale_shape_manual(values={'MaGECK + U-test':'o','MaGECK only':'^','U-test only':'s','Non-hit':'o'},name='Hit-calling method')
       + labs(title='Genome-wide screen identifies 244 regulators',
              subtitle=f'{full["id"].notna().sum():,} target records; 244 unique screen-derived regulators',
              x='MaGECK log2 fold change\n← required for EC differentiation     restricts EC differentiation →',
              y='−log10(minimum MaGECK FDR)')
       + coord_cartesian(xlim=(-4.25,2.75),ylim=(-0.15,4.65))
       + theme_classic(base_size=9)
       + theme(figure_size=(6.65,3.05),
               plot_title=element_text(size=12,weight='bold',color=DARK),
               plot_subtitle=element_text(size=8.3,color=MID),
               axis_title=element_text(size=8.4,weight='bold'),
               axis_text=element_text(size=7.5),
               legend_title=element_text(size=7.6,weight='bold'),
               legend_text=element_text(size=7.1),
               legend_position='right',
               panel_grid_major=element_blank(),panel_grid_minor=element_blank(),
               axis_line=element_line(color='#4D545A',size=0.55)))
    save_plotnine(p,'Fig1b_genomewide_screen_nature_style',6.65,3.05)

# -----------------------------
# Panel C — literature-supported known regulators
# -----------------------------
def build_panel_c():
    selected=['KDR','SMAD4','ACVR2A','SMAD1','SMAD5','BMPR1A','CDH5','NF2','PLCG1','ETV2','VHL','KCTD10','LATS2','AMOTL2']
    d=known[known.gene.isin(selected)].copy()
    d['tier_short']=np.where(d.known_EC_regulator_tier.str.startswith('Tier 1'),'Tier 1: direct / established','Tier 2: pathway-supported')
    d['gene_order']=pd.Categorical(d.gene,categories=d.sort_values('mageck_lfc').gene,ordered=True)
    d['validated']=np.where(d.validation_screen_hit,'validated','not significant in validation')
    d['zero']=0
    p=(ggplot(d,aes(y='gene_order'))
       + geom_segment(aes(x='zero',xend='mageck_lfc',yend='gene_order'),color='#AEB4B8',size=0.75)
       + geom_point(aes(x='mageck_lfc',fill='tier_short',shape='validated'),color='#3E454A',stroke=0.35,size=3.1)
       + geom_vline(xintercept=0,linetype='dashed',color='#9AA0A5',size=0.45)
       + scale_fill_manual(values={'Tier 1: direct / established':BLUE,'Tier 2: pathway-supported':PURPLE},name='Literature evidence')
       + scale_shape_manual(values={'validated':'o','not significant in validation':'^'},name='Validation')
       + labs(title='Known endothelial and vascular regulators are recovered',
              subtitle='Literature-supported regulators selected before plotting',
              x='MaGECK log2 fold change',y='')
       + coord_cartesian(xlim=(-3.85,0.2))
       + theme_classic(base_size=8)
       + theme(figure_size=(3.25,3.15),
               plot_title=element_text(size=10.1,weight='bold',color=DARK),
               plot_subtitle=element_text(size=7.2,color=MID),
               axis_title_x=element_text(size=7.6,weight='bold'),
               axis_text_y=element_text(size=7.5,style='italic'),
               axis_text_x=element_text(size=6.8),
               legend_title=element_text(size=6.7,weight='bold'),
               legend_text=element_text(size=6.2),
               legend_position='bottom',
               panel_grid_major=element_blank(),panel_grid_minor=element_blank(),
               axis_line_y=element_blank(),axis_ticks_major_y=element_blank()))
    save_plotnine(p,'Fig1c_known_regulators_nature_style',3.25,3.15)

# -----------------------------
# Panel D — functional categories diverging bars
# -----------------------------
def build_panel_d():
    cat=(ann.groupby(['figure1d_category','screen_direction']).size().reset_index(name='count'))
    req_name='Required for EC differentiation'
    res_name='Restricts EC differentiation'
    cats=ann.groupby('figure1d_category').size().sort_values(ascending=True).index.tolist()
    cat['category']=pd.Categorical(cat.figure1d_category,categories=cats,ordered=True)
    cat['signed_count']=np.where(cat.screen_direction==req_name,-cat['count'],cat['count'])
    cat['label']=cat['count'].astype(str)
    cat['label_x']=np.where(cat.signed_count<0,cat.signed_count-2,cat.signed_count+1.0)
    cat['direction_short']=np.where(cat.screen_direction==req_name,'Required for EC differentiation','Restricts EC differentiation')
    # add zeros for missing direction/category for stable bars
    complete=pd.MultiIndex.from_product([cats,[req_name,res_name]],names=['figure1d_category','screen_direction']).to_frame(index=False)
    cat=complete.merge(cat[['figure1d_category','screen_direction','count']],how='left').fillna({'count':0})
    cat['category']=pd.Categorical(cat.figure1d_category,categories=cats,ordered=True)
    cat['signed_count']=np.where(cat.screen_direction==req_name,-cat['count'],cat['count'])
    cat['label']=np.where(cat['count']>0,cat['count'].astype(int).astype(str),'')
    cat['label_x']=np.where(cat.signed_count<0,cat.signed_count-2,cat.signed_count+1.2)
    cat['direction_short']=np.where(cat.screen_direction==req_name,'Required for EC differentiation','Restricts EC differentiation')
    cat.to_csv(SOURCE/'Fig1d_functional_categories_source.tsv',sep='\t',index=False)
    p=(ggplot(cat,aes('category','signed_count',fill='direction_short'))
       + geom_col(width=0.72)
       + geom_text(aes(y='label_x',label='label'),size=7.0,color=DARK)
       + geom_hline(yintercept=0,color='#5A6065',size=0.55)
       + coord_flip()
       + scale_fill_manual(values={'Required for EC differentiation':BLUE,'Restricts EC differentiation':ORANGE},name='')
       + labs(title='Screen-derived regulators span lineage and core cellular processes',
              subtitle='244 unique regulators; bars show gene counts by screen direction',
              x='',y='Number of regulators  ← required     restricts →')
       + coord_cartesian(ylim=(-68,8))
       + theme_classic(base_size=8)
       + theme(figure_size=(4.45,3.15),
               plot_title=element_text(size=10.2,weight='bold',color=DARK),
               plot_subtitle=element_text(size=7.2,color=MID),
               axis_title_x=element_text(size=7.4,weight='bold'),
               axis_text_y=element_text(size=6.7),
               axis_text_x=element_text(size=6.6),
               legend_position='bottom',legend_text=element_text(size=6.4),
               panel_grid_major_y=element_blank(),panel_grid_minor=element_blank(),
               axis_line_y=element_blank(),axis_ticks_major_y=element_blank()))
    save_plotnine(p,'Fig1d_functional_categories_nature_style',4.45,3.15)

# -----------------------------
# Panel E — validation and replicate concordance
# -----------------------------
def build_panel_e():
    # main comparison
    main=val[['Gene','Original FC U Test','Average U Test FC','FDR-corrected','Type']].dropna().copy()
    main.columns=['gene','x','y','validation_fdr','target_class']
    class_map={'genome_wide_hit':'Screen-derived regulator','genome_wide_non_hit':'Random screen non-hit','tri_perturb_seq':'Additional curated target'}
    main['target_class']=main.target_class.map(class_map)
    label_genes=['KDR','CDH5','ETV2','ACVR2A','SMAD1','SMAD5','BMPR1A','PLCG1']
    lab=main[main.gene.isin(label_genes)].copy()
    off={'KDR':(-0.05,0.18),'CDH5':(0.12,0.10),'ETV2':(-0.12,-0.14),'ACVR2A':(0.12,-0.15),'SMAD1':(-0.10,0.12),'SMAD5':(0.14,0.10),'BMPR1A':(-0.12,-0.12),'PLCG1':(0.12,0.15)}
    lab['lx']=[r.x+off.get(r.gene,(.1,.1))[0] for r in lab.itertuples()]
    lab['ly']=[r.y+off.get(r.gene,(.1,.1))[1] for r in lab.itertuples()]
    r,pv=pearsonr(main.x,main.y)
    same=((main.x*main.y)>0).sum(); sig=(main.validation_fdr<.05).sum()
    hit=main[main.target_class=='Screen-derived regulator']
    hit_sig=(hit.validation_fdr<.05).sum(); hit_same=((hit.x*hit.y)>0).sum()
    main_stats=f'n = {len(main)}; Pearson r = {r:.2f}; R² = {r*r:.2f}\n{hit_sig}/{len(hit)} screen hits significant; {hit_same}/{len(hit)} same direction'
    p1=(ggplot(main,aes('x','y',color='target_class'))
        + geom_point(size=1.8,alpha=0.82,stroke=0.15)
        + geom_abline(slope=1,intercept=0,linetype='dashed',color='#9AA0A5',size=0.45)
        + geom_smooth(method='lm',se=False,color=DARK,size=0.75)
        + geom_hline(yintercept=0,linetype='dotted',color='#AEB3B7',size=0.4)
        + geom_vline(xintercept=0,linetype='dotted',color='#AEB3B7',size=0.4)
        + geom_segment(lab,aes(x='x',y='y',xend='lx',yend='ly'),inherit_aes=False,color='#777E83',size=0.3)
        + geom_text(lab,aes('lx','ly',label='gene'),inherit_aes=False,size=6.2,fontstyle='italic',color=DARK)
        + annotate('text',x=-3.7,y=4.0,label=main_stats,ha='left',va='top',size=6.6,color=DARK)
        + scale_color_manual(values={'Screen-derived regulator':BLUE,'Random screen non-hit':'#BFC4C7','Additional curated target':ORANGE},name='Target class')
        + labs(title='Secondary screening validates genome-wide effects',
               x='Genome-wide U-test log2 fold change',y='Validation U-test log2 fold change')
        + coord_cartesian(xlim=(-4.0,2.4),ylim=(-4.2,4.35))
        + theme_classic(base_size=8)
        + theme(figure_size=(3.75,3.15),plot_title=element_text(size=10.2,weight='bold',color=DARK),
                axis_title=element_text(size=7.4,weight='bold'),axis_text=element_text(size=6.7),
                legend_title=element_text(size=6.6,weight='bold'),legend_text=element_text(size=6.2),
                legend_position='bottom',panel_grid_major=element_blank(),panel_grid_minor=element_blank()))
    save_plotnine(p1,'Fig1e_validation_main_nature_style',3.75,3.15)

    rep=val[['Gene','U Test FC BR1','U Test FC BR2','Type']].dropna().copy()
    rep.columns=['gene','x','y','target_class']
    rep['target_class']=rep.target_class.map(class_map)
    rr,rp=pearsonr(rep.x,rep.y)
    p2=(ggplot(rep,aes('x','y',color='target_class'))
        + geom_point(size=1.55,alpha=0.78,stroke=0.1)
        + geom_abline(slope=1,intercept=0,linetype='dashed',color='#9AA0A5',size=0.4)
        + geom_smooth(method='lm',se=False,color=DARK,size=0.65)
        + annotate('text',x=-3.8,y=4.0,label=f'n = {len(rep)}\nr = {rr:.2f}; R² = {rr*rr:.2f}',ha='left',va='top',size=6.5,color=DARK)
        + scale_color_manual(values={'Screen-derived regulator':BLUE,'Random screen non-hit':'#BFC4C7','Additional curated target':ORANGE},name='')
        + labs(title='Validation replicate concordance',x='Replicate 1 log2 fold change',y='Replicate 2 log2 fold change')
        + coord_cartesian(xlim=(-4.1,2.5),ylim=(-4.1,4.3))
        + theme_classic(base_size=8)
        + theme(figure_size=(2.55,3.15),plot_title=element_text(size=9.1,weight='bold',color=DARK),
                axis_title=element_text(size=7.0,weight='bold'),axis_text=element_text(size=6.4),
                legend_position='none',panel_grid_major=element_blank(),panel_grid_minor=element_blank()))
    save_plotnine(p2,'Fig1e_validation_replicates_nature_style',2.55,3.15)

# -----------------------------
# SVG assembly helper
# -----------------------------
def parse_viewbox(root):
    vb=root.attrib.get('viewBox')
    if vb:
        vals=list(map(float,vb.replace(',',' ').split()))
        return vals
    w=float(re.sub('[^0-9.]','',root.attrib.get('width','100')))
    h=float(re.sub('[^0-9.]','',root.attrib.get('height','100')))
    return [0,0,w,h]

def nest_svg(parent, path, x, y, w, h):
    tree=ET.parse(path); root=tree.getroot(); vb=parse_viewbox(root)
    ns='http://www.w3.org/2000/svg'
    nested=ET.SubElement(parent,f'{{{ns}}}svg',{'x':str(x),'y':str(y),'width':str(w),'height':str(h),'viewBox':' '.join(map(str,vb)),'preserveAspectRatio':'xMidYMid meet'})
    for child in list(root):
        nested.append(child)

# -----------------------------
# Build + assemble
# -----------------------------
build_panel_a(); build_panel_b(); build_panel_c(); build_panel_d(); build_panel_e()

ns='http://www.w3.org/2000/svg'; ET.register_namespace('',ns)
W,H=1800,1420
root=ET.Element(f'{{{ns}}}svg',{'width':str(W),'height':str(H),'viewBox':f'0 0 {W} {H}'})
ET.SubElement(root,f'{{{ns}}}rect',{'x':'0','y':'0','width':str(W),'height':str(H),'fill':'white'})
# title
ET.SubElement(root,f'{{{ns}}}text',{'x':'55','y':'58','font-family':'DejaVu Sans','font-size':'34','font-weight':'700','fill':DARK}).text='Genome-wide CRISPRi screen identifies regulators of endothelial differentiation'
ET.SubElement(root,f'{{{ns}}}text',{'x':'55','y':'91','font-family':'DejaVu Sans','font-size':'16','fill':MID}).text='244 unique screen-derived regulators; functional annotation and secondary validation'
# panels
nest_svg(root,PANELS/'Fig1a_workflow_nature_style.svg',55,125,650,500)
nest_svg(root,PANELS/'Fig1b_genomewide_screen_nature_style.svg',735,125,1010,500)
nest_svg(root,PANELS/'Fig1c_known_regulators_nature_style.svg',55,700,430,600)
nest_svg(root,PANELS/'Fig1d_functional_categories_nature_style.svg',500,700,670,600)
nest_svg(root,PANELS/'Fig1e_validation_main_nature_style.svg',1185,700,400,600)
nest_svg(root,PANELS/'Fig1e_validation_replicates_nature_style.svg',1540,700,240,600)
# panel letters
for letter,x,y in [('a',42,134),('b',720,134),('c',42,709),('d',488,709),('e',1173,709)]:
    ET.SubElement(root,f'{{{ns}}}text',{'x':str(x),'y':str(y),'font-family':'DejaVu Sans','font-size':'31','font-weight':'700','fill':DARK}).text=letter

svg_path=ASSEMBLED/'Figure1_nature_style_remade.svg'
ET.ElementTree(root).write(svg_path,encoding='utf-8',xml_declaration=True)
cairosvg.svg2pdf(url=str(svg_path),write_to=str(ASSEMBLED/'Figure1_nature_style_remade.pdf'),output_width=W,output_height=H)
cairosvg.svg2png(url=str(svg_path),write_to=str(ASSEMBLED/'Figure1_nature_style_remade.png'),output_width=W*2,output_height=H*2)

# -----------------------------
# README / legend / QC
# -----------------------------
val_hit=val[val.Type=='genome_wide_hit']
r_main=pearsonr(val['Original FC U Test'],val['Average U Test FC'])[0]
r_rep=pearsonr(val['U Test FC BR1'],val['U Test FC BR2'])[0]
readme=f'''# Figure 1 nature-style remake\n\nAll quantitative panels were regenerated from the uploaded Batch A source files.\n\n- Fig. 1a: conceptual vector redraw of the screen workflow.\n- Fig. 1b: full MaGECK target-level result plot using 20,655 MaGECK target records; 244 unique final regulators are colored.\n- Fig. 1c: literature-supported EC/vascular regulators selected from the reviewed annotation table.\n- Fig. 1d: diverging category counts from the 244-gene literature-informed annotation draft.\n- Fig. 1e: genome-wide U-test effect versus validation U-test effect plus validation replicate concordance.\n\nValidation summary: {int((val_hit['FDR-corrected']<.05).sum())}/{len(val_hit)} screen hits significant in validation; {int(((val_hit['Original FC U Test']*val_hit['Average U Test FC'])>0).sum())}/{len(val_hit)} same direction. Genome-wide versus validation Pearson r={r_main:.3f}; replicate Pearson r={r_rep:.3f}.\n\nIndividual SVG/PDF panels are editable in Adobe Illustrator. The assembled SVG preserves vector panel contents.\n'''
(OUT/'README.md').write_text(readme)
legend='''**Figure 1 | A genome-wide CRISPRi screen identifies regulators of endothelial differentiation.** **a,** Doxycycline-inducible CRISPRi hiPSCs were infected with a genome-wide guide library, differentiated through primitive streak and lateral mesoderm to the day-3 endothelial endpoint, separated into CDH5-positive and CDH5-negative fractions, and analyzed by guide sequencing. **b,** Target-level MaGECK results for the full tested set. Final screen-derived regulators are colored according to whether knockdown was required for or restricted EC differentiation; point shape denotes hit-calling method. **c,** Literature-supported endothelial and vascular regulators recovered among screen hits. **d,** Functional organization of the 244 unique regulators. Counts are split by screen direction. Categories are literature-informed draft annotations and should be author-reviewed. **e,** Comparison of genome-wide and focused validation effects, with validation replicate concordance shown at right.\n'''
(OUT/'Figure1_legend_draft.md').write_text(legend)

# Copy this script
shutil.copy(__file__,SCRIPTS/'build_figure1_nature_style_remade.py')

# Package
zip_path=BASE/'Figure1_nature_style_remade_package.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for p in OUT.rglob('*'):
        if p.is_file(): z.write(p,p.relative_to(OUT.parent))
print('Created',zip_path)
