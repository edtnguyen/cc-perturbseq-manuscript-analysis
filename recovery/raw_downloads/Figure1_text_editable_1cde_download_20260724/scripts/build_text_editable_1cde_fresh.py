from pathlib import Path
import os, zipfile, shutil, re
import numpy as np
import pandas as pd
import matplotlib as mpl

# Preserve live text in Illustrator-friendly vector outputs.
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'DejaVu Sans'

from scipy.stats import pearsonr
from plotnine import (
    ggplot, aes, geom_point, geom_hline, geom_vline, geom_abline, geom_smooth,
    geom_segment, geom_text, geom_col, scale_fill_gradient, scale_size_continuous,
    scale_fill_manual, scale_color_manual, scale_y_continuous, scale_x_discrete,
    scale_y_discrete, labs, theme, theme_classic, element_text, element_blank,
    element_line, element_rect, coord_flip, coord_cartesian
)

ROOT = Path('/mnt/data')
OUT = ROOT / 'Figure1_text_editable_1cde_fresh'
PAN = OUT / 'individual_panels'
SRC = OUT / 'source_data'
SCR = OUT / 'scripts'
for d in (OUT, PAN, SRC, SCR): d.mkdir(parents=True, exist_ok=True)

# Input tables already available in this active workspace.
C_SRC = ROOT / 'Figure1_restyled_1c1d_remade/source_data/Fig1c_nonHippo_known_regulators_dotmatrix.tsv'
D_STACK = ROOT / 'Figure1_restyled_1c1d_remade/source_data/Fig1d_functional_categories_stacked.tsv'
D_SUM = ROOT / 'Figure1_restyled_1c1d_remade/source_data/Fig1d_functional_categories_restyled.tsv'
E_SRC = ROOT / 'Figure1_nature_style_remade/source_data/Fig1e_validation_source.tsv'
for p in (C_SRC, D_STACK, D_SUM, E_SRC):
    if not p.exists(): raise FileNotFoundError(p)

c_long = pd.read_csv(C_SRC, sep='\t')
bars = pd.read_csv(D_STACK, sep='\t')
d = pd.read_csv(D_SUM, sep='\t')
val = pd.read_csv(E_SRC, sep='\t')

# Copy plotting source tables into the new package.
shutil.copy2(C_SRC, SRC / C_SRC.name)
shutil.copy2(D_STACK, SRC / D_STACK.name)
shutil.copy2(D_SUM, SRC / D_SUM.name)
shutil.copy2(E_SRC, SRC / E_SRC.name)

DARK='#171C21'; MID='#606B74'; BLUE='#2E6F97'; ORANGE='#D87528'; PURPLE='#7C4D8E'

# -----------------------------
# Panel 1c — non-Hippo known regulators across assays
# -----------------------------
ASSAY_ORDER=['MaGECK','U-test','Validation']
c_long['assay']=pd.Categorical(c_long['assay'], categories=ASSAY_ORDER, ordered=True)
# Derive stable gene order from y position already saved.
gene_order=(c_long[['gene','y_pos']].drop_duplicates().sort_values('y_pos',ascending=False)['gene'].tolist())
n_genes=len(gene_order)
bottom_order=list(reversed(gene_order))
# group boundaries used by the original design
block_boundaries=[]
for _,grp in c_long[['gene','group','y_pos']].drop_duplicates().sort_values('y_pos',ascending=False).groupby('group',sort=False):
    pass
# hard-coded boundaries from original 3-group design, safe for 16 genes
if n_genes >= 12:
    block_boundaries=[n_genes-7+0.5, n_genes-12+0.5]

p_c=(
    ggplot(c_long, aes('assay','y_pos'))
    + geom_hline(yintercept=block_boundaries, color='#D7DDE2', size=0.55)
    + geom_point(aes(fill='effect', size='minus_log10_fdr'), shape='o', color='#2F3438', stroke=0.35, alpha=0.98)
    + scale_fill_gradient(low='#154C73', high='#EEF3F5', limits=(-4.5,0), breaks=[-4,-2,0], name='Effect\nlog2FC')
    + scale_size_continuous(range=(2.5,8.0), breaks=[1,3,5], limits=(0.4,6.0), name='−log10(FDR)')
    + scale_x_discrete(limits=ASSAY_ORDER)
    + scale_y_continuous(breaks=list(range(1,n_genes+1)), labels=bottom_order, limits=(0.5,n_genes+0.5))
    + labs(title='Known EC regulators across assays', subtitle='Non-Hippo regulators; color = effect, size = significance', x='', y='')
    + theme_classic(base_size=11)
    + theme(
        figure_size=(5.35,5.25), axis_line=element_blank(), axis_ticks=element_blank(),
        axis_text_x=element_text(size=10.2,weight='bold',color='#273038'),
        axis_text_y=element_text(size=9.0,style='italic',color='#273038'),
        axis_title=element_blank(), plot_title=element_text(size=12.6,weight='bold',color=DARK,margin={'b':4}),
        plot_subtitle=element_text(size=8.7,color=MID,margin={'b':8}),
        legend_title=element_text(size=8.7,weight='bold'), legend_text=element_text(size=7.8),
        legend_position='right', panel_grid=element_blank()
    )
)

# -----------------------------
# Panel 1d — functional category validation composition
# -----------------------------
d = d.sort_values('n_regulators', ascending=True).reset_index(drop=True)
d['category_label']=pd.Categorical(d['category_label'],categories=d['category_label'].tolist(),ordered=True)
bars['category_label']=pd.Categorical(bars['category_label'],categories=d['category_label'].tolist(),ordered=True)
bars['segment']=pd.Categorical(bars['segment'],categories=['Validated','Not validated'],ordered=True)
max_total=int(d['n_regulators'].max()); restrict_x=max_total+8
d['restrict_x']=restrict_x

p_d=(
    ggplot()
    + geom_col(data=bars,mapping=aes('category_label','n',fill='segment'),position='stack',width=0.66)
    + geom_text(data=d,mapping=aes('category_label','n_regulators + 1.0',label='validated_label'),ha='left',size=8.4/2.834,color='#34404A')
    + geom_text(data=d,mapping=aes('category_label','restrict_x',label='n_restricting'),color=ORANGE,size=9.0/2.834,fontweight='bold')
    + scale_fill_manual(values={'Validated':BLUE,'Not validated':'#D9E6EE'},breaks=['Validated','Not validated'],name='')
    + scale_y_continuous(limits=(0,restrict_x+3),breaks=[0,20,40,60],expand=(0,0))
    + coord_flip()
    + labs(title='Functional categories and validation',subtitle='Dark = validated; label = validated/total; orange = restrictors',x='',y='Number of regulators')
    + theme_classic(base_size=11)
    + theme(
        figure_size=(6.25,5.25), axis_line_y=element_blank(), axis_ticks_y=element_blank(),
        axis_line_x=element_line(color='#4A535A',size=0.6), axis_ticks_x=element_line(color='#4A535A',size=0.5),
        axis_text_y=element_text(size=8.8,color='#273038'), axis_text_x=element_text(size=8.5,color='#5A646C'),
        axis_title_x=element_text(size=9.5,weight='bold',color='#273038',margin={'t':6}), axis_title_y=element_blank(),
        plot_title=element_text(size=12.8,weight='bold',color=DARK,margin={'b':4}),
        plot_subtitle=element_text(size=8.8,color=MID,margin={'b':8}), legend_position='top',
        legend_title=element_blank(), legend_text=element_text(size=8.4),
        legend_key=element_rect(fill='white',color='white'), panel_grid=element_blank()
    )
)

# -----------------------------
# Panel 1e — validation effects and replicate concordance
# -----------------------------
class_map={'genome_wide_hit':'Screen-derived regulator','genome_wide_non_hit':'Random screen non-hit','tri_perturb_seq':'Additional curated target'}
main=val[['Gene','Original FC U Test','Average U Test FC','FDR-corrected','Type']].dropna().copy()
main.columns=['gene','x','y','validation_fdr','target_class']
main['target_class']=main.target_class.map(class_map)
label_genes=['KDR','CDH5','ETV2','ACVR2A','SMAD1','SMAD5','BMPR1A','PLCG1']
lab=main[main.gene.isin(label_genes)].copy()
off={'KDR':(-0.05,0.18),'CDH5':(0.12,0.10),'ETV2':(-0.12,-0.14),'ACVR2A':(0.12,-0.15),'SMAD1':(-0.10,0.12),'SMAD5':(0.14,0.10),'BMPR1A':(-0.12,-0.12),'PLCG1':(0.12,0.15)}
lab['lx']=[r.x+off.get(r.gene,(.1,.1))[0] for r in lab.itertuples()]
lab['ly']=[r.y+off.get(r.gene,(.1,.1))[1] for r in lab.itertuples()]
r,pv=pearsonr(main.x,main.y)
hit=main[main.target_class=='Screen-derived regulator']
hit_sig=(hit.validation_fdr<.05).sum(); hit_same=((hit.x*hit.y)>0).sum()
main_stats=f'n = {len(main)}; Pearson r = {r:.2f}; R² = {r*r:.2f}\n{hit_sig}/{len(hit)} screen hits significant; {hit_same}/{len(hit)} same direction'

p_e_main=(
    ggplot(main,aes('x','y',color='target_class'))
    + geom_point(size=1.8,alpha=0.82,stroke=0.15)
    + geom_abline(slope=1,intercept=0,linetype='dashed',color='#9AA0A5',size=0.45)
    + geom_smooth(method='lm',se=False,color=DARK,size=0.75)
    + geom_hline(yintercept=0,linetype='dotted',color='#AEB3B7',size=0.4)
    + geom_vline(xintercept=0,linetype='dotted',color='#AEB3B7',size=0.4)
    + geom_segment(lab,aes(x='x',y='y',xend='lx',yend='ly'),inherit_aes=False,color='#777E83',size=0.3)
    + geom_text(lab,aes('lx','ly',label='gene'),inherit_aes=False,size=6.2,fontstyle='italic',color=DARK)
    + geom_text(pd.DataFrame({'x':[-3.7],'y':[4.0],'label':[main_stats]}),aes('x','y',label='label'),inherit_aes=False,ha='left',va='top',size=6.6,color=DARK)
    + scale_color_manual(values={'Screen-derived regulator':BLUE,'Random screen non-hit':'#BFC4C7','Additional curated target':ORANGE},name='Target class')
    + labs(title='Secondary screening validates genome-wide effects',x='Genome-wide U-test log2 fold change',y='Validation U-test log2 fold change')
    + coord_cartesian(xlim=(-4.0,2.4),ylim=(-4.2,4.35))
    + theme_classic(base_size=8)
    + theme(figure_size=(3.75,3.15),plot_title=element_text(size=10.2,weight='bold',color=DARK),
            axis_title=element_text(size=7.4,weight='bold'),axis_text=element_text(size=6.7),
            legend_title=element_text(size=6.6,weight='bold'),legend_text=element_text(size=6.2),
            legend_position='bottom',panel_grid_major=element_blank(),panel_grid_minor=element_blank())
)

rep=val[['Gene','U Test FC BR1','U Test FC BR2','Type']].dropna().copy()
rep.columns=['gene','x','y','target_class']
rep['target_class']=rep.target_class.map(class_map)
rr,rp=pearsonr(rep.x,rep.y)
rep_stats=f'n = {len(rep)}\nr = {rr:.2f}; R² = {rr*rr:.2f}'
p_e_rep=(
    ggplot(rep,aes('x','y',color='target_class'))
    + geom_point(size=1.55,alpha=0.78,stroke=0.1)
    + geom_abline(slope=1,intercept=0,linetype='dashed',color='#9AA0A5',size=0.4)
    + geom_smooth(method='lm',se=False,color=DARK,size=0.65)
    + geom_text(pd.DataFrame({'x':[-3.8],'y':[4.0],'label':[rep_stats]}),aes('x','y',label='label'),inherit_aes=False,ha='left',va='top',size=6.5,color=DARK)
    + scale_color_manual(values={'Screen-derived regulator':BLUE,'Random screen non-hit':'#BFC4C7','Additional curated target':ORANGE},name='')
    + labs(title='Validation replicate concordance',x='Replicate 1 log2 fold change',y='Replicate 2 log2 fold change')
    + coord_cartesian(xlim=(-4.1,2.5),ylim=(-4.1,4.3))
    + theme_classic(base_size=8)
    + theme(figure_size=(2.55,3.15),plot_title=element_text(size=9.1,weight='bold',color=DARK),
            axis_title=element_text(size=7.0,weight='bold'),axis_text=element_text(size=6.4),
            legend_position='none',panel_grid_major=element_blank(),panel_grid_minor=element_blank())
)

plots=[
    ('Fig1c_known_regulators_nonHippo_text_editable',p_c,5.35,5.25),
    ('Fig1d_functional_categories_text_editable',p_d,6.25,5.25),
    ('Fig1e_validation_main_text_editable',p_e_main,3.75,3.15),
    ('Fig1e_validation_replicates_text_editable',p_e_rep,2.55,3.15),
]

for name,p,w,h in plots:
    p.save(PAN/f'{name}.png',width=w,height=h,dpi=420,verbose=False)
    p.save(PAN/f'{name}.pdf',width=w,height=h,verbose=False)
    p.save(PAN/f'{name}.svg',width=w,height=h,verbose=False)

# Verify SVGs contain live <text> elements.
qc=[]
for name,_,_,_ in plots:
    svg=PAN/f'{name}.svg'
    txt=svg.read_text(errors='ignore')
    qc.append({'panel':name,'svg_file':svg.name,'n_text_elements':len(re.findall(r'<text\b',txt)),'contains_font_family':('font-family' in txt)})
pd.DataFrame(qc).to_csv(OUT/'svg_text_qc.tsv',sep='\t',index=False)

readme='''# Figure 1 panels 1c, 1d and 1e — Illustrator text-editable exports\n\nThese files were freshly regenerated in the active workspace.\n\n- SVG export uses `matplotlib.rcParams["svg.fonttype"] = "none"`, preserving supported text as live SVG text.\n- PDF export uses TrueType font embedding (`pdf.fonttype = 42`).\n- Open the SVG versions first in Adobe Illustrator for the most reliable editable text.\n- Panel 1c excludes Hippo genes.\n- Panel 1e is supplied as the main validation panel and the replicate-concordance panel.\n'''
(OUT/'README.md').write_text(readme)
shutil.copy2(__file__,SCR/'build_text_editable_1cde_fresh.py')

zip_path=ROOT/'Figure1_text_editable_1cde_fresh_package.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for p in OUT.rglob('*'):
        if p.is_file(): z.write(p,arcname=str(p.relative_to(OUT)))
print('Created',OUT)
print('Package',zip_path)
print(pd.DataFrame(qc).to_string(index=False))
