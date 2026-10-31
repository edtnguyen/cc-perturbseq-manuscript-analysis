import os
import io
import base64
import textwrap
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageChops
import matplotlib.pyplot as plt

from plotnine import (
    ggplot, aes, geom_point, geom_col, geom_text, geom_hline, coord_flip,
    facet_grid, scale_fill_gradient2, scale_size_continuous,
    scale_fill_manual, scale_x_continuous, scale_y_continuous, labs, coord_cartesian,
    theme_classic, theme, element_text, element_blank, element_line,
    guides, guide_colorbar, guide_legend
)

ROOT = Path('/mnt/data/Figure1_restyled_1c1d_v2')
PANELS = ROOT / 'individual_panels'
ASSEMBLED = ROOT / 'assembled_figures'
SOURCE = ROOT / 'source_data'
SCRIPTS = ROOT / 'scripts'
for d in [PANELS, ASSEMBLED, SOURCE, SCRIPTS]:
    d.mkdir(parents=True, exist_ok=True)

# --------------------------
# Shared visual settings
# --------------------------
BLUE = '#2B6E9E'
BLUE_DARK = '#174A73'
ORANGE = '#D87625'
PURPLE = '#74539A'
GREY = '#AEB5BA'
TEXT = '#20262B'
MUTED = '#68727A'

# --------------------------
# Panel 1c: non-Hippo known regulators
# --------------------------
c_src = Path('/mnt/data/Figure1_restyled_1c1d_remade/source_data/Fig1c_nonHippo_known_regulators_dotmatrix.tsv')
if c_src.exists():
    c = pd.read_csv(c_src, sep='\t')
else:
    # Rebuild from the complete known-regulator table if the convenience file is missing.
    full = pd.read_csv('/mnt/data/Figure1_nature_style_remade_fresh/source_data/Fig1c_known_regulators_source.tsv', sep='\t')
    hippo = {'NF2','SAV1','LATS2','AMOTL2','TAOK1','WWC3','YAP1','WWTR1'}
    keep = ['CDH5','ETV2','KDR','KCTD10','PLCG1','PTEN','VHL',
            'BMPR1A','SMAD1','SMAD5','ACVR2A','SMAD4',
            'EOMES','HAND1','MIXL1','TBXT']
    full = full[(full['gene'].isin(keep)) & (~full['gene'].isin(hippo))].copy()
    group_map = {
        'CDH5':'Endothelial identity','ETV2':'Endothelial identity','KDR':'Endothelial identity',
        'KCTD10':'Endothelial identity','PLCG1':'Endothelial identity','PTEN':'Endothelial identity','VHL':'Endothelial identity',
        'BMPR1A':'BMP/TGF signaling','SMAD1':'BMP/TGF signaling','SMAD5':'BMP/TGF signaling','ACVR2A':'BMP/TGF signaling','SMAD4':'BMP/TGF signaling',
        'EOMES':'Primitive streak / mesoderm','HAND1':'Primitive streak / mesoderm','MIXL1':'Primitive streak / mesoderm','TBXT':'Primitive streak / mesoderm'
    }
    rows=[]
    for _,r in full.iterrows():
        vals=[('MaGECK',r['mageck_lfc'],min(r['mageck_neg_fdr'],r['mageck_pos_fdr'])),
              ('U-test',r['u_test_effect_log2_CDH5pos_vs_CDH5neg'],r['u_test_fdr']),
              ('Validation',r['validation_screen_log2fc'], np.nan if not r['validation_screen_hit'] else 1e-6)]
        for assay,effect,fdr in vals:
            rows.append({'gene':r['gene'],'group':group_map[r['gene']], 'assay':assay,
                         'effect':effect,'fdr':fdr,
                         'minus_log10_fdr':np.nan if pd.isna(fdr) else min(6,-np.log10(max(fdr,1e-6)))})
    c=pd.DataFrame(rows)

# Hard guard: no Hippo genes in panel 1c.
HIPPO = {'NF2','SAV1','LATS2','AMOTL2','TAOK1','WWC3','YAP1','WWTR1'}
c = c[~c['gene'].isin(HIPPO)].copy()

# order groups and genes by mean effect (stronger depletion first)
group_order = ['Endothelial identity','BMP/TGF signaling','Primitive streak / mesoderm']
assay_order = ['MaGECK','U-test','Validation']
mean_effect = c.groupby('gene')['effect'].mean().sort_values()
ordered_genes=[]
for g in group_order:
    genes = c[c['group']==g]['gene'].unique().tolist()
    genes = sorted(genes, key=lambda x: mean_effect.get(x,0))
    ordered_genes.extend(genes)
# Reverse levels so first gene appears at top in plotnine
c['gene'] = pd.Categorical(c['gene'], categories=list(reversed(ordered_genes)), ordered=True)
c['group'] = pd.Categorical(c['group'], categories=group_order, ordered=True)
c['assay_display'] = c['assay'].replace({'Validation':'Validation\nscreen'})
c['assay_display'] = pd.Categorical(c['assay_display'], categories=['MaGECK','U-test','Validation\nscreen'], ordered=True)
c['assay'] = pd.Categorical(c['assay'], categories=assay_order, ordered=True)
c['minus_log10_fdr_plot'] = c['minus_log10_fdr'].clip(lower=0.4, upper=6).fillna(0.4)

c.to_csv(SOURCE/'Fig1c_nonHippo_known_regulators_source.tsv', sep='\t', index=False)

p_c = (
    ggplot(c, aes('assay_display','gene'))
    + geom_point(
        aes(size='minus_log10_fdr_plot', fill='effect'),
        shape='o', color='#3F464B', stroke=0.38, alpha=0.98
    )
    + facet_grid('group ~ .', scales='free_y', space='free_y')
    + scale_fill_gradient2(
        low='#1F5A8A', mid='#F7F7F7', high='#C8582B', midpoint=0,
        limits=(-4.3,0.5), breaks=[-4,-2,0],
        name='Effect\nlog2FC'
    )
    + scale_size_continuous(
        range=(2.3,8.0), limits=(0,6), breaks=[1,3,5],
        name='-log10(FDR)'
    )
    + labs(
        tag='c',
        title='Known non-Hippo regulators recovered',
        subtitle='Consistent effects across genome-wide and validation screens',
        x='', y=''
    )
    + guides(fill=guide_colorbar(order=1), size=guide_legend(order=2))
    + theme_classic(base_size=11)
    + theme(
        figure_size=(6.8,5.2),
        axis_text_x=element_text(size=9.4, weight='bold', linespacing=0.9),
        axis_text_y=element_text(size=10, style='italic'),
        axis_ticks_major_y=element_blank(),
        axis_ticks_major_x=element_blank(),
        axis_line=element_blank(),
        panel_grid_major_x=element_line(color='#E6E9EB', size=0.55),
        panel_grid_minor=element_blank(),
        strip_background=element_blank(),
        strip_text_y=element_text(size=9.5, weight='bold', color=MUTED, angle=0),
        plot_title=element_text(size=14.6, weight='bold', color=TEXT),
        plot_subtitle=element_text(size=9.6, color=MUTED),
        plot_tag=element_text(size=23, weight='bold', color=TEXT),
        legend_title=element_text(size=9.5, weight='bold'),
        legend_text=element_text(size=8.8),
        legend_position='right',
    )
)

for ext in ['png','pdf','svg']:
    p_c.save(PANELS/f'Fig1c_known_regulators_nonHippo_restyled.{ext}', width=6.8, height=5.2, dpi=350 if ext=='png' else 300, verbose=False)

# --------------------------
# Panel 1d: functional organization
# --------------------------
d_src = pd.read_csv('/mnt/data/Figure1_nature_style_remade_fresh/source_data/Fig1d_annotation_source.tsv', sep='\t')
summary = (d_src.groupby(['figure1d_category','screen_direction']).size()
           .unstack(fill_value=0).reset_index())
summary.columns.name=None
summary = summary.rename(columns={
    'Required for EC differentiation':'Required',
    'Restricts EC differentiation':'Restricts'
})
summary['Validated'] = summary['figure1d_category'].map(
    d_src.groupby('figure1d_category')['validation_screen_hit'].sum().to_dict()
).fillna(0).astype(int)
summary['Total'] = summary[['Required','Restricts']].sum(axis=1)
summary = summary.sort_values('Total', ascending=False).reset_index(drop=True)
summary['category_label'] = summary['figure1d_category'].map(lambda s: textwrap.fill(s.replace(' and ',' & '), width=29))
summary['category_label'] = pd.Categorical(summary['category_label'], categories=list(reversed(summary['category_label'].tolist())), ordered=True)
summary.to_csv(SOURCE/'Fig1d_functional_categories_restyled_source.tsv', sep='\t', index=False)

long = summary.melt(
    id_vars=['figure1d_category','category_label','Validated','Total'],
    value_vars=['Required','Restricts'],
    var_name='Direction', value_name='Count'
)
long['Direction'] = pd.Categorical(long['Direction'], categories=['Required','Restricts'], ordered=True)

# position totals slightly after bar; validated circle uses same count axis
p_d = (
    ggplot(long, aes('category_label','Count', fill='Direction'))
    + geom_col(position='stack', width=0.62, color='white', size=0.25)
    + geom_point(
        data=summary,
        mapping=aes(x='category_label', y='Validated'),
        inherit_aes=False,
        shape='o', size=3.7, stroke=0.9, fill='white', color='#222222'
    )
    + geom_text(
        data=summary,
        mapping=aes(x='category_label', y='Total+1.0', label='Total'),
        inherit_aes=False,
        size=3.2, va='center', color=TEXT
    )
    + scale_fill_manual(values={'Required':BLUE,'Restricts':ORANGE}, name='Screen effect')
    + scale_y_continuous(limits=(0,66), breaks=[0,20,40,60], expand=(0,0))
    + coord_flip()
    + labs(
        tag='d',
        title='Functional classes of 244 regulators',
        subtitle='Bars = screen hits; circles = validation hits',
        x='', y='Number of regulators'
    )
    + theme_classic(base_size=11)
    + theme(
        figure_size=(6.4,5.2),
        axis_text_y=element_text(size=9.2, color=TEXT),
        axis_text_x=element_text(size=9.5),
        axis_title_x=element_text(size=11, weight='bold'),
        axis_title_y=element_blank(),
        axis_ticks_major_y=element_blank(),
        axis_line_y=element_blank(),
        panel_grid_major_x=element_line(color='#E6E9EB', size=0.55),
        panel_grid_minor=element_blank(),
        plot_title=element_text(size=14.6, weight='bold', color=TEXT),
        plot_subtitle=element_text(size=9.6, color=MUTED),
        plot_tag=element_text(size=23, weight='bold', color=TEXT),
        legend_title=element_text(size=9.5, weight='bold'),
        legend_text=element_text(size=9),
        legend_position='top',
        legend_direction='horizontal',
    )
)

for ext in ['png','pdf','svg']:
    p_d.save(PANELS/f'Fig1d_functional_categories_restyled.{ext}', width=6.4, height=5.2, dpi=350 if ext=='png' else 300, verbose=False)

# --------------------------
# Assemble updated Figure 1 (raster-backed PDF/PNG, individual panels remain vector SVG/PDF)
# --------------------------

def crop_white(path, pad=10):
    im=Image.open(path).convert('RGB')
    bg=Image.new('RGB', im.size, 'white')
    diff=ImageChops.difference(im,bg)
    bbox=diff.getbbox()
    if bbox is None:
        return im
    x0,y0,x1,y1=bbox
    x0=max(0,x0-pad); y0=max(0,y0-pad); x1=min(im.width,x1+pad); y1=min(im.height,y1+pad)
    return im.crop((x0,y0,x1,y1))

panel_paths={
    'a':'/mnt/data/Figure1_nature_style_remade_fresh/individual_panels/Fig1a_workflow.png',
    'b':'/mnt/data/Figure1_nature_style_remade_fresh/individual_panels/Fig1b_genomewide_screen.png',
    'c':str(PANELS/'Fig1c_known_regulators_nonHippo_restyled.png'),
    'd':str(PANELS/'Fig1d_functional_categories_restyled.png'),
    'e':'/mnt/data/Figure1_nature_style_remade_fresh/individual_panels/Fig1e_validation.png',
}
images={k:crop_white(v, pad=12) for k,v in panel_paths.items()}

fig=plt.figure(figsize=(18,12), facecolor='white')
fig.text(0.035,0.975,'Genome-wide CRISPRi screen identifies regulators of endothelial differentiation',
         ha='left',va='top',fontsize=25,fontweight='bold',color=TEXT)
fig.text(0.035,0.945,'Screen design, genome-wide discovery, literature-supported recovery, functional organization and secondary validation',
         ha='left',va='top',fontsize=12.5,color=MUTED)

positions={
    'a':[0.03,0.51,0.35,0.40],
    'b':[0.39,0.51,0.58,0.40],
    'c':[0.03,0.055,0.31,0.40],
    'd':[0.35,0.055,0.30,0.40],
    'e':[0.66,0.055,0.31,0.40],
}
for k,pos in positions.items():
    ax=fig.add_axes(pos)
    ax.imshow(images[k])
    ax.axis('off')

png_out=ASSEMBLED/'Figure1_restyled_1c1d_v2.png'
pdf_out=ASSEMBLED/'Figure1_restyled_1c1d_v2.pdf'
fig.savefig(png_out,dpi=350,bbox_inches='tight',facecolor='white')
fig.savefig(pdf_out,bbox_inches='tight',facecolor='white')
plt.close(fig)

# Copy script into package
import shutil
shutil.copy(__file__, SCRIPTS/'build_restyled_1c1d_v2.py')

readme = f'''# Figure 1 restyled panels 1c and 1d\n\n- Panel 1c excludes all Hippo-pathway genes: {', '.join(sorted(HIPPO))}.\n- Panel 1c is a three-assay dot matrix: fill = effect size; size = -log10(FDR).\n- Panel 1d is a stacked functional-category bar chart: blue = required, orange = restricts, open circle = validation hits.\n- Individual panels are available as editable SVG and vector PDF.\n- The assembled PDF/PNG uses the existing A, B and E panels with the new C and D panels.\n'''
(ROOT/'README.md').write_text(readme)

# zip final package
zip_path=Path('/mnt/data/Figure1_restyled_1c1d_v2_package.zip')
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for fp in ROOT.rglob('*'):
        if fp.is_file():
            z.write(fp,arcname=str(fp.relative_to(ROOT)))
print('Created', ROOT)
print('ZIP', zip_path)
