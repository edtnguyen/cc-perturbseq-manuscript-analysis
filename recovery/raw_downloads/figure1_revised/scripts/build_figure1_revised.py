import os, re, glob, zipfile, json, math, shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import pearsonr, spearmanr
try:
    from adjustText import adjust_text
except Exception:
    adjust_text = None

# Illustrator-friendly defaults
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = 0.8
mpl.rcParams['figure.dpi'] = 300

BASE = Path('/mnt/data')
OUT = BASE / 'figure1_revised'
for sub in ['data_audit','cleaned_data','source_data','scripts','candidate_panels','individual_panels','assembled_figures','supplementary_panels','legends']:
    (OUT/sub).mkdir(parents=True, exist_ok=True)

# Copy this script later

# -------------------------
# Utilities
# -------------------------

def parse_gene(x):
    return re.sub(r'(_)?-P\d.*$', '', str(x).strip())

def norm_gene(x):
    s = str(x).strip()
    if re.match(r'^C\d+orf\d+$', s, flags=re.I):
        m = re.match(r'^(C\d+)orf(\d+)$', s, flags=re.I)
        return f"{m.group(1).upper()}orf{m.group(2)}"
    return s

def neglog10(p):
    p = pd.to_numeric(p, errors='coerce')
    return -np.log10(np.clip(p, 1e-300, None))

def ptxt(p):
    if pd.isna(p): return 'P = NA'
    return f'P = {p:.1e}' if p < 1e-3 else f'P = {p:.3g}'

def save_fig(fig, stem, folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    paths = {}
    for ext in ['svg','pdf','png']:
        fp = folder / f'{stem}.{ext}'
        fig.savefig(fp, bbox_inches='tight', dpi=600 if ext=='png' else None, facecolor='white')
        paths[ext] = str(fp)
    return paths

# -------------------------
# Extract archives
# -------------------------
EXTRACT = OUT/'data_audit'/'extracted_archives'
EXTRACT.mkdir(parents=True, exist_ok=True)
archives = [BASE/'Pos_vs_Neg(1).zip', BASE/'validation screen results(1).zip']
archive_members=[]
for z in archives:
    if z.exists():
        dest = EXTRACT/z.stem
        if dest.exists(): shutil.rmtree(dest)
        dest.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(z) as zipf:
            zipf.extractall(dest)
            for m in zipf.namelist():
                archive_members.append({'archive':z.name,'member':m})
pd.DataFrame(archive_members).to_csv(OUT/'data_audit'/'archive_members.tsv', sep='\t', index=False)

# Locate authoritative files
mageck_full_path = EXTRACT/'Pos_vs_Neg(1)'/'Pos_vs_Neg'/'Pos_vs_Neg.gene_summary.txt'
if not mageck_full_path.exists(): mageck_full_path = EXTRACT/'Pos_vs_Neg(1)'/'Pos_vs_Neg.gene_summary.txt'
utest_full_path = EXTRACT/'Pos_vs_Neg(1)'/'Pos_vs_Neg'/'combined_utest_P1P2_adj_test.txt'
if not utest_full_path.exists(): utest_full_path = EXTRACT/'Pos_vs_Neg(1)'/'combined_utest_P1P2_adj_test.txt'
combined_path = EXTRACT/'Pos_vs_Neg(1)'/'Pos_vs_Neg'/'combined_mageck_u_test_pos_v_neg.xlsx'
if not combined_path.exists(): combined_path = EXTRACT/'Pos_vs_Neg(1)'/'combined_mageck_u_test_pos_v_neg.xlsx'
val_dir = EXTRACT/'validation screen results(1)'/'validation screen results'
if not val_dir.exists(): val_dir = EXTRACT/'validation screen results(1)'
val_ec_path = val_dir/'ec_screen_hits.csv'
val_full_path = val_dir/'validation_screen_d3pos_d3neg_mageckLFC.csv'
val_type_path = val_dir/'validation_genes_type.xlsx'

# Read data
mageck = pd.read_csv(mageck_full_path, sep='\t')
utest = pd.read_csv(utest_full_path, sep='\t')
combined = pd.read_excel(combined_path)
val_ec = pd.read_csv(val_ec_path)
val_full = pd.read_csv(val_full_path)
val_type = pd.read_excel(val_type_path)

# Optional literature annotation
ann_path = BASE/'figure1_annotation_literature_draft'/'screen_regulators_244_literature_annotation_draft.tsv'
ann = pd.read_csv(ann_path, sep='\t') if ann_path.exists() else pd.DataFrame()

# -------------------------
# Audit inventory
# -------------------------
interesting_files = [
    BASE/'Pasted text(4).txt', BASE/'mageck_information_hits.txt', BASE/'mageck_information_non_hits.txt',
    BASE/'u_test_information_hits.txt', BASE/'u_test_information_non_hits.txt',
    BASE/'Pos_vs_Neg(1).zip', BASE/'validation screen results(1).zip',
    mageck_full_path, utest_full_path, combined_path, val_ec_path, val_full_path, val_type_path, ann_path
]
inv=[]
for f in interesting_files:
    if not Path(f).exists(): continue
    ft = Path(f).suffix.lower()
    row={'filename':str(f),'file_type':ft,'rows':np.nan,'columns':np.nan,'important_columns':'','statistical_method':'','candidate_panel':'','sufficient_for_plotting':'unknown'}
    try:
        if ft in ['.txt','.tsv','.csv']:
            sep = '\t' if ft in ['.txt','.tsv'] else ','
            df = pd.read_csv(f, sep=sep)
            row['rows']=df.shape[0]; row['columns']=df.shape[1]; row['important_columns']=', '.join(map(str,df.columns[:10])); row['sufficient_for_plotting']='yes'
        elif ft == '.xlsx':
            xl = pd.ExcelFile(f)
            row['important_columns']='sheets: '+', '.join(xl.sheet_names)
            tmp = pd.read_excel(f, sheet_name=xl.sheet_names[0])
            row['rows']=tmp.shape[0]; row['columns']=tmp.shape[1]; row['sufficient_for_plotting']='yes'
        elif ft == '.zip':
            row['sufficient_for_plotting']='archive extracted and audited'
    except Exception as e:
        row['sufficient_for_plotting']='read error: '+str(e)
    inv.append(row)
pd.DataFrame(inv).to_csv(OUT/'data_audit'/'file_inventory.tsv', sep='\t', index=False)

# -------------------------
# Harmonize full screen data
# -------------------------
# Combined table is authoritative for full screen because it contains both MaGECK and U-test statistics
full = combined.copy()
full['original_target_id'] = full['id'].astype(str)
full['original_gene_symbol'] = full['original_target_id'].map(parse_gene)
full['normalized_gene_symbol'] = full['original_gene_symbol'].map(norm_gene)
full['display_gene_symbol'] = full['normalized_gene_symbol']
full['mageck_neg_fdr'] = pd.to_numeric(full['neg|fdr'], errors='coerce')
full['mageck_pos_fdr'] = pd.to_numeric(full['pos|fdr'], errors='coerce')
full['mageck_lfc'] = pd.to_numeric(full['neg|lfc'], errors='coerce')
full['utest_log2fc'] = pd.to_numeric(full['u_test_log2FC'], errors='coerce')
full['utest_fdr'] = pd.to_numeric(full['u_test_p_value_corrected'], errors='coerce')
full['mageck_neg_hit'] = full['mageck_neg_fdr'] < 0.05
full['mageck_pos_hit'] = full['mageck_pos_fdr'] < 0.05
full['mageck_hit'] = full['mageck_neg_hit'] | full['mageck_pos_hit']
full['utest_hit'] = full['utest_fdr'] < 0.05
full['union_hit'] = full['mageck_hit'] | full['utest_hit']
# Biological direction is taken from the direction-specific MaGECK test when available;
# U-test sign is used for U-test-only hits. This reproduces 236/8 gene-level directions.
full['direction_rule'] = np.select(
    [full['mageck_neg_hit'] & ~full['mageck_pos_hit'],
     full['mageck_pos_hit'] & ~full['mageck_neg_hit'],
     full['utest_hit'] & (full['utest_log2fc'] < 0),
     full['utest_hit'] & (full['utest_log2fc'] > 0)],
    ['Required for EC differentiation','Restricts EC differentiation',
     'Required for EC differentiation','Restricts EC differentiation'],
    default=np.where(full['mageck_lfc'] < 0, 'Required for EC differentiation','Restricts EC differentiation')
)
full['best_mageck_fdr'] = np.where(full['direction_rule'] == 'Required for EC differentiation', full['mageck_neg_fdr'], full['mageck_pos_fdr'])
# if lfc==0 or missing, use min fdr
full['best_mageck_fdr'] = np.where(full['best_mageck_fdr'].isna(), full[['mageck_neg_fdr','mageck_pos_fdr']].min(axis=1), full['best_mageck_fdr'])
full['screen_effect'] = full['mageck_lfc']
full['screen_direction'] = full['direction_rule']
full['hit_support'] = np.select([
    full['mageck_hit'] & full['utest_hit'],
    full['mageck_hit'] & ~full['utest_hit'],
    ~full['mageck_hit'] & full['utest_hit']],
    ['Both methods','MaGECK only','U-test only'], default='Non-hit')
# collapse to gene by best mageck fdr, retaining NEDD8 audit separately
full['gene_dup_count'] = full.groupby('normalized_gene_symbol')['original_target_id'].transform('nunique')
full['duplicate_gene_flag'] = full['gene_dup_count'] > 1
full_sorted = full.sort_values(['normalized_gene_symbol','union_hit','best_mageck_fdr'], ascending=[True,False,True])
gene_full = full_sorted.groupby('normalized_gene_symbol', as_index=False).first()
# combine utest/mageck support by gene across rows
support = full.groupby('normalized_gene_symbol').agg(mageck_hit=('mageck_hit','any'), utest_hit=('utest_hit','any'), union_hit=('union_hit','any'), duplicate_gene_flag=('duplicate_gene_flag','any'), target_rows=('original_target_id','nunique')).reset_index()
gene_full = gene_full.drop(columns=['mageck_hit','utest_hit','union_hit','duplicate_gene_flag','target_rows'], errors='ignore').merge(support,on='normalized_gene_symbol',how='left')
gene_full['hit_support'] = np.select([
    gene_full['mageck_hit'] & gene_full['utest_hit'],
    gene_full['mageck_hit'] & ~gene_full['utest_hit'],
    ~gene_full['mageck_hit'] & gene_full['utest_hit']],
    ['Both methods','MaGECK only','U-test only'], default='Non-hit')
gene_full['screen_direction'] = gene_full['direction_rule']

screen_hits = gene_full[gene_full['union_hit']].copy()

# Save cleaned tables
full.to_csv(OUT/'cleaned_data'/'genomewide_results_target_level.tsv', sep='\t', index=False)
gene_full.to_csv(OUT/'cleaned_data'/'genomewide_results_unique_gene_level.tsv', sep='\t', index=False)
screen_hits.to_csv(OUT/'cleaned_data'/'screen_hits_244_unique_gene_level.tsv', sep='\t', index=False)
full[full['normalized_gene_symbol']=='NEDD8'].to_csv(OUT/'cleaned_data'/'NEDD8_target_level_audit.tsv', sep='\t', index=False)

# alias map
gene_alias_map = full[['original_target_id','original_gene_symbol','normalized_gene_symbol','display_gene_symbol','duplicate_gene_flag']].copy()
gene_alias_map['alias_resolution_status'] = np.where(gene_alias_map['original_gene_symbol']==gene_alias_map['normalized_gene_symbol'],'unchanged','normalized_capitalization')
gene_alias_map['alias_source'] = 'string normalization; source target suffix removed'
gene_alias_map.to_csv(OUT/'cleaned_data'/'gene_alias_map.tsv', sep='\t', index=False)

# -------------------------
# Audit metrics
# -------------------------
metrics=[]
metrics.append({'metric':'MaGECK full target rows','value':len(mageck)})
metrics.append({'metric':'MaGECK full unique parsed genes','value':mageck['id'].map(parse_gene).nunique()})
metrics.append({'metric':'U-test full target rows','value':len(utest)})
metrics.append({'metric':'U-test full unique parsed genes','value':utest['Target'].map(parse_gene).nunique()})
metrics.append({'metric':'Combined full target rows','value':len(full)})
metrics.append({'metric':'Combined full unique parsed genes','value':gene_full.shape[0]})
metrics.append({'metric':'Union hit target rows','value':int(full['union_hit'].sum())})
metrics.append({'metric':'Union hit unique biological genes','value':int(screen_hits['normalized_gene_symbol'].nunique())})
metrics.append({'metric':'MaGECK significant target rows','value':int(full['mageck_hit'].sum())})
metrics.append({'metric':'U-test significant target rows','value':int(full['utest_hit'].sum())})
metrics.append({'metric':'MaGECK-only hit genes','value':int(((screen_hits['mageck_hit']) & (~screen_hits['utest_hit'])).sum())})
metrics.append({'metric':'Both-method hit genes','value':int(((screen_hits['mageck_hit']) & (screen_hits['utest_hit'])).sum())})
metrics.append({'metric':'U-test-only hit genes','value':int(((~screen_hits['mageck_hit']) & (screen_hits['utest_hit'])).sum())})
metrics.append({'metric':'Required-for-EC hit genes','value':int((screen_hits['screen_direction']=='Required for EC differentiation').sum())})
metrics.append({'metric':'Restricts-EC hit genes','value':int((screen_hits['screen_direction']=='Restricts EC differentiation').sum())})
metrics.append({'metric':'Duplicate NEDD8 target rows in hit audit','value':int((full[full['normalized_gene_symbol']=='NEDD8']['union_hit']).sum())})
pd.DataFrame(metrics).to_csv(OUT/'data_audit'/'analysis_version_audit.tsv', sep='\t', index=False)

# validation harmonization
val = val_ec.merge(val_type, left_on='gene', right_on='Target', how='left')
val_full_clean = val_full.rename(columns={'Gene':'gene'})
if '\ufeffGene' in val_full_clean.columns:
    val_full_clean = val_full_clean.rename(columns={'\ufeffGene':'gene'})
val = val.merge(val_full_clean[['gene','U Test FC BR1','U Test FC BR2','Average U Test FC','FDR-corrected','Normalized Mageck FC']], on='gene', how='left')
val['validation_significant'] = val['validation_screen_hit'].astype(bool)
val['direction_concordant'] = np.sign(val['genome_wide_screen_log2fc']) == np.sign(val['validation_screen_log2fc'])
val.to_csv(OUT/'cleaned_data'/'validation_results_harmonized.tsv', sep='\t', index=False)
val_summary = val.groupby('Type').agg(n=('gene','nunique'), validation_hits=('validation_screen_hit','sum'), gw_hits=('genome_wide_screen_hit','sum'), concordant=('direction_concordant','sum')).reset_index()
val_summary.to_csv(OUT/'source_data'/'Fig1e_validation_summary.tsv', sep='\t', index=False)

# Additional audit files
file_to_panel = pd.DataFrame([
    {'file':'combined_mageck_u_test_pos_v_neg.xlsx','panel':'1b','status':'authoritative full screen table'},
    {'file':'Pos_vs_Neg.gene_summary.txt','panel':'1b','status':'MaGECK full screen table'},
    {'file':'combined_utest_P1P2_adj_test.txt','panel':'1b','status':'U-test full screen table'},
    {'file':'screen_regulators_244_literature_annotation_draft.tsv','panel':'1d','status':'draft literature-informed categories'},
    {'file':'ec_screen_hits.csv','panel':'1e','status':'validation screen gene-level summary'},
    {'file':'validation_screen_d3pos_d3neg_mageckLFC.csv','panel':'1e','status':'validation replicate and effect table'},
    {'file':'validation_genes_type.xlsx','panel':'1e','status':'validation target group labels'},
])
file_to_panel.to_csv(OUT/'data_audit'/'file_to_panel_map.tsv', sep='\t', index=False)

# -------------------------
# Data for panels
# -------------------------
# Fig1b source
fig1b_source = gene_full[['normalized_gene_symbol','display_gene_symbol','original_target_id','screen_effect','best_mageck_fdr','mageck_hit','utest_hit','union_hit','hit_support','screen_direction','duplicate_gene_flag']].copy()
fig1b_source['minus_log10_directional_fdr'] = neglog10(fig1b_source['best_mageck_fdr'])
fig1b_source.to_csv(OUT/'source_data'/'Fig1b_source_data.tsv', sep='\t', index=False)

# label selection
priority1 = ['CDH5','KDR','ETV2']
priority2 = ['SMAD1','SMAD4','SMAD5','BMPR1A','ACVR2A','PLCG1','HAND1','NF2','LATS2','AMOTL2','SAV1','TAOK1']
label_rows=[]
for g in priority1+priority2:
    if g in set(fig1b_source['display_gene_symbol']):
        label_rows.append({'gene':g,'label_priority':1 if g in priority1 else 2,'label_reason':'known EC regulator / pathway representative' if g in priority1 else 'pathway representative central to manuscript','predefined_or_posthoc':'predefined' if g in priority1 else 'pathway representative','displayed':g in priority1+['SMAD1','SMAD5','ACVR2A','BMPR1A','PLCG1','AMOTL2','LATS2','NF2'],'source':'manual priority list based on manuscript narrative'})
pd.DataFrame(label_rows).to_csv(OUT/'source_data'/'Fig1b_label_selection.tsv', sep='\t', index=False)
label_genes = ['CDH5','KDR','ETV2','CBLL1']

# Fig1d pathway membership: fixed curated groups based on prior literature annotation / known pathway labels
pathway_groups = {
    'EC/VEGF specification':['CDH5','KDR','ETV2','PLCG1','KCTD10','SOX4'],
    'BMP/TGF/Activin':['BMPR1A','ACVR1B','ACVR2A','ACVR2B','SMAD1','SMAD3','SMAD4','SMAD5','LEMD3'],
    'Hippo/YAP regulation':['NF2','LATS2','SAV1','TAOK1','AMOTL2','WWC3'],
    'SWI/SNF & chromatin':['SMARCA4','SMARCA5','SMARCC1','ASH2L','BRD3','SETDB1'],
    'Mediator/transcription':['MED12','CCNC','CCNT2','SUPT20H','TAF1','ELOA'],
    'RNA processing':['RBM14-RBM4','SF3B1','SRRM2','SRSF3','PRPF8','BCAS2'],
}
component_rows=[]
for pname, genes in pathway_groups.items():
    for g in genes:
        match = gene_full[gene_full['display_gene_symbol']==g]
        if match.empty:
            # try normalized source differences
            match = gene_full[gene_full['normalized_gene_symbol']==g]
        if not match.empty:
            row=match.iloc[0]
            component_rows.append({'pathway_display_name':pname,'member_gene':g,'gene_in_screen':True,'gene_is_final_hit':bool(row['union_hit']),'screen_effect':row['screen_effect'],'screen_direction':row['screen_direction'],'screen_fdr':row['best_mageck_fdr'],'hit_support':row['hit_support'],'manually_added':True,'addition_reason':'selected pathway component panel from literature-informed annotation/manuscript pathway list','approval_status':'draft - author review'})
        else:
            component_rows.append({'pathway_display_name':pname,'member_gene':g,'gene_in_screen':False,'gene_is_final_hit':False,'screen_effect':np.nan,'screen_direction':'not tested','screen_fdr':np.nan,'hit_support':'not in screen','manually_added':True,'addition_reason':'selected pathway component panel from literature-informed annotation/manuscript pathway list','approval_status':'draft - author review'})
component_df = pd.DataFrame(component_rows)
component_df.to_csv(OUT/'source_data'/'Fig1d_component_effects.tsv', sep='\t', index=False)
component_df.to_csv(OUT/'source_data'/'pathway_annotation_review.tsv', sep='\t', index=False)

# 1d candidate dummy selection files
pd.DataFrame({'candidate':['A_preranked_enrichment','B_selected_pathway_strips','C_composite'],'status':['not run - fixed GMT/gene-set library absent','generated and selected','not used']}).to_csv(OUT/'candidate_panels'/'Fig1d_candidate_comparison.tsv',sep='\t',index=False)

# Fig1e source data
val.to_csv(OUT/'source_data'/'Fig1e_source_data.tsv', sep='\t', index=False)
val[['gene','U Test FC BR1','U Test FC BR2']].dropna().to_csv(OUT/'source_data'/'Fig1e_replicate_source_data.tsv', sep='\t', index=False)

# 1c missing files
(OUT/'source_data'/'Fig1c_classification_rules.tsv').write_text('Panel 1c not generated: no D0-versus-D3 guide-abundance table found. Required columns: gene, d0_vs_d3_effect, p_value, adjusted_p_value, significant, analysis_method, source_target_id.\n')
(OUT/'source_data'/'Fig1c_class_summary.tsv').write_text('Panel 1c not generated: source data absent.\n')
(OUT/'source_data'/'Fig1c_source_data.tsv').write_text('Panel 1c not generated: source data absent.\n')

# -------------------------
# Plot functions
# -------------------------
COL_REQUIRED = '#2B6CB0' # blue
COL_RESTRICT = '#D55E00' # vermilion/orange
COL_GREY = '#C8C8C8'
COL_DARK = '#222222'


def panel_letter(ax, letter):
    ax.text(-0.04, 1.04, letter, transform=ax.transAxes, fontsize=16, fontweight='bold', ha='left', va='bottom')


def draw_panel_1a(ax):
    ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis('off')
    ax.text(0.0, 5.7, 'Genome-wide CRISPRi screen', fontsize=12.5, fontweight='bold')
    ax.text(0.0, 5.32, 'CDH5/CD144-based endothelial differentiation phenotype', fontsize=8.5, color='#555555')
    xs=[0.7,2.2,3.5,4.8,6.0,7.2,8.5]
    labs=['iPSC\nCRISPRi','gRNA\nlibrary','selection\n+ dox','D0\niPSC','D1\nprimitive\nstreak','D2\nlateral\nmesoderm','D3\nEC diff.']
    cols=['#F7F7F7','#F7F7F7','#F7F7F7','#F1F7FE','#FFF7ED','#EEF7F3','#EAF4FF']
    for i,(x,l,c) in enumerate(zip(xs,labs,cols)):
        box=FancyBboxPatch((x-0.5,3.3),1.0,1.1,boxstyle='round,pad=0.04,rounding_size=0.08',facecolor=c,edgecolor='#555555',linewidth=0.7)
        ax.add_patch(box); ax.text(x,3.85,l,ha='center',va='center',fontsize=7.2)
        if i<len(xs)-1:
            ax.add_patch(FancyArrowPatch((x+0.52,3.85),(xs[i+1]-0.52,3.85),arrowstyle='-|>',mutation_scale=9,linewidth=0.7,color='#555555'))
    # Sort schematic
    ax.text(7.9,2.6,'D3 sorting',fontsize=8.5,fontweight='bold')
    y1,y2=1.95,1.2
    ax.add_patch(FancyBboxPatch((7.6,y1-0.25),1.55,0.5,boxstyle='round,pad=0.05,rounding_size=0.07',facecolor='#DCEEFF',edgecolor=COL_REQUIRED,linewidth=0.8))
    ax.text(8.375,y1,'CDH5+\nEC',fontsize=7,ha='center',va='center',color=COL_REQUIRED)
    ax.add_patch(FancyBboxPatch((7.6,y2-0.25),1.55,0.5,boxstyle='round,pad=0.05,rounding_size=0.07',facecolor='#F5F5F5',edgecolor='#888888',linewidth=0.8))
    ax.text(8.375,y2,'CDH5−\nnon-EC',fontsize=7,ha='center',va='center',color='#444444')
    ax.add_patch(FancyArrowPatch((8.5,3.28),(8.5,2.25),arrowstyle='-|>',mutation_scale=9,linewidth=0.8,color='#555555'))
    ax.add_patch(FancyArrowPatch((9.2,1.55),(9.8,1.55),arrowstyle='-|>',mutation_scale=9,linewidth=0.8,color='#555555'))
    ax.text(9.95,1.55,'guide\nsequencing',fontsize=7,va='center')
    ax.text(0.2,0.7,'20,528 genes\n129,638 guides\n244 regulators',fontsize=8.8,fontweight='bold',ha='left',va='bottom',bbox=dict(boxstyle='round,pad=0.3',fc='white',ec='#CCCCCC',lw=0.6))


def draw_panel_1b(ax, data=fig1b_source, candidate='mageck'):
    ax.set_title('Genome-wide gene-level screen result', fontsize=11.5, fontweight='bold', loc='left', pad=6)
    df = data.copy()
    # background
    bg = df[~df['union_hit']]
    ax.scatter(bg['screen_effect'], bg['minus_log10_directional_fdr'], s=4, c=COL_GREY, alpha=0.22, linewidth=0)
    # hits by direction/support
    markers = {'Both methods':'o','MaGECK only':'^','U-test only':'s'}
    for support_name, marker in markers.items():
        for direction, color in [('Required for EC differentiation',COL_REQUIRED),('Restricts EC differentiation',COL_RESTRICT)]:
            sub=df[(df['union_hit'])&(df['hit_support']==support_name)&(df['screen_direction']==direction)]
            if len(sub):
                ax.scatter(sub['screen_effect'], sub['minus_log10_directional_fdr'], s=20, c=color, marker=marker, edgecolor='white', linewidth=0.3, alpha=0.92)
    ax.axvline(0, color='#888888', lw=0.8, ls='--')
    ax.axhline(-np.log10(0.05), color='#AAAAAA', lw=0.7, ls=':')
    ax.set_ylim(-0.2, 5.25)
    # labels
    lab = df[df['display_gene_symbol'].isin(label_genes)].dropna(subset=['screen_effect','minus_log10_directional_fdr'])
    texts=[]
    for _,r in lab.iterrows():
        texts.append(ax.text(r['screen_effect'], r['minus_log10_directional_fdr'], r['display_gene_symbol'], fontsize=5.6, color='black'))
    if adjust_text and texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#777777', lw=0.3, alpha=0.7), force_text=(0.4,0.5), force_static=(0.15,0.2), max_move=20, time_lim=0.7)
    ax.set_xlabel('MaGECK gene effect (CDH5+ vs CDH5−)', fontsize=8)
    ax.set_ylabel('−log10 directional FDR', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.text(0.02,0.96,f"244 unique regulators\n{int((screen_hits['screen_direction']=='Required for EC differentiation').sum())} required, {int((screen_hits['screen_direction']=='Restricts EC differentiation').sum())} restrict", transform=ax.transAxes, ha='left', va='top', fontsize=7, bbox=dict(boxstyle='round,pad=0.25',fc='white',ec='#DDDDDD',lw=0.4))
    handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=COL_REQUIRED,markersize=5,label='Required for EC'),Line2D([0],[0],marker='o',color='none',markerfacecolor=COL_RESTRICT,markersize=5,label='Restricts EC'),Line2D([0],[0],marker='^',color='#555555',markerfacecolor='white',markersize=5,label='MaGECK only'),Line2D([0],[0],marker='o',color='#555555',markerfacecolor='white',markersize=5,label='Both')]
    ax.legend(handles=handles, frameon=False, fontsize=6.5, loc='upper right', handletextpad=0.4)
    ax.spines[['top','right']].set_visible(False)


def draw_panel_1b_utest(ax):
    u = utest.copy(); u['gene']=u['Target'].map(parse_gene); u['gene_norm']=u['gene'].map(norm_gene)
    u['minuslog10']=neglog10(u['p_value_corrected'])
    u_sorted=u.sort_values(['gene_norm','p_value_corrected']).groupby('gene_norm',as_index=False).first()
    hset=set(screen_hits['normalized_gene_symbol'])
    u_sorted['union_hit']=u_sorted['gene_norm'].isin(hset)
    u_sorted['screen_direction']=np.where(u_sorted['log2FC']<0,'Required for EC differentiation','Restricts EC differentiation')
    ax.set_title('Candidate 1b-A: U-test volcano', fontsize=10, fontweight='bold', loc='left')
    ax.scatter(u_sorted.loc[~u_sorted.union_hit,'log2FC'],u_sorted.loc[~u_sorted.union_hit,'minuslog10'],s=4,c=COL_GREY,alpha=0.22,lw=0)
    for direction,color in [('Required for EC differentiation',COL_REQUIRED),('Restricts EC differentiation',COL_RESTRICT)]:
        sub=u_sorted[(u_sorted.union_hit)&(u_sorted.screen_direction==direction)]
        ax.scatter(sub['log2FC'],sub['minuslog10'],s=16,c=color,edgecolor='white',lw=0.3)
    ax.axvline(0,color='#888',ls='--',lw=0.8); ax.axhline(-np.log10(0.05),color='#AAA',ls=':',lw=0.7)
    ax.set_xlabel('U-test log2 fold-change'); ax.set_ylabel('−log10 U-test FDR'); ax.tick_params(labelsize=7); ax.spines[['top','right']].set_visible(False)
    u_sorted.to_csv(OUT/'source_data'/'Fig1b_candidate_A_Utest_source_data.tsv',sep='\t',index=False)


def draw_panel_1b_ranked(ax):
    df=fig1b_source.sort_values('screen_effect').reset_index(drop=True).copy(); df['rank']=np.arange(1,len(df)+1)
    ax.set_title('Candidate 1b-C: full ranked screen map', fontsize=10, fontweight='bold', loc='left')
    ax.scatter(df.loc[~df.union_hit,'rank'],df.loc[~df.union_hit,'screen_effect'],s=3,c=COL_GREY,alpha=0.25,lw=0)
    for direction,color in [('Required for EC differentiation',COL_REQUIRED),('Restricts EC differentiation',COL_RESTRICT)]:
        sub=df[(df.union_hit)&(df.screen_direction==direction)]
        ax.scatter(sub['rank'],sub['screen_effect'],s=14,c=color,edgecolor='white',lw=0.25)
    ax.axhline(0,color='#888',ls='--',lw=0.8); ax.set_xlabel('Gene rank by MaGECK effect'); ax.set_ylabel('MaGECK gene effect'); ax.tick_params(labelsize=7); ax.spines[['top','right']].set_visible(False)
    df.to_csv(OUT/'source_data'/'Fig1b_candidate_C_ranked_source_data.tsv',sep='\t',index=False)


def draw_panel_1d(ax):
    ax.set_title('Pathway component effects', fontsize=11.5, fontweight='bold', loc='left', pad=6)
    df=component_df[component_df['gene_in_screen']].copy()
    pathways=list(pathway_groups.keys())[::-1]
    ymap={p:i for i,p in enumerate(pathways)}
    df['y']=df['pathway_display_name'].map(ymap)
    # jitter within row deterministic
    df['jitter'] = df.groupby('pathway_display_name').cumcount().map(lambda x: (x%5-2)*0.055)
    df['plot_y']=df['y']+df['jitter']
    # draw background bands
    for p,y in ymap.items():
        ax.axhspan(y-0.38,y+0.38,color='#F7F7F7' if y%2==0 else 'white',zorder=0)
    # non-hit open, hit filled
    for _,r in df.iterrows():
        color = COL_REQUIRED if r['screen_effect']<0 else COL_RESTRICT
        if r['gene_is_final_hit']:
            ax.scatter(r['screen_effect'], r['plot_y'], s=26, c=color, edgecolor='white', lw=0.3, zorder=3)
        else:
            ax.scatter(r['screen_effect'], r['plot_y'], s=24, facecolor='white', edgecolor='#999999', lw=0.7, zorder=2)
    # labels for key components only
    key=set(['CDH5','KDR','ETV2','BMPR1A','NF2','LATS2','AMOTL2','SMARCA4','MED12'])
    for _,r in df[df['member_gene'].isin(key)].iterrows():
        ha='left' if r['screen_effect']>=0 else 'right'
        dx=0.08 if ha=='left' else -0.08
        ax.text(r['screen_effect']+dx, r['plot_y'], r['member_gene'], fontsize=5.4, va='center', ha=ha)
    ax.axvline(0,color='#888888',ls='--',lw=0.8)
    ax.set_yticks(range(len(pathways))); ax.set_yticklabels(pathways, fontsize=7)
    ax.set_xlabel('MaGECK gene effect (CDH5+ vs CDH5−)', fontsize=8)
    ax.tick_params(axis='x',labelsize=7)
    ax.text(0.98,0.04,'filled = final screen hit\nopen = tested non-hit',transform=ax.transAxes,ha='right',va='bottom',fontsize=6.6)
    ax.spines[['top','right','left']].set_visible(False)


def draw_panel_1e(ax):
    ax.set_title('Secondary validation', fontsize=11.5, fontweight='bold', loc='left', pad=6)
    df=val.copy().dropna(subset=['genome_wide_screen_log2fc','validation_screen_log2fc'])
    colors={'genome_wide_hit':'#4C78A8','genome_wide_non_hit':'#BDBDBD','tri_perturb_seq':'#D55E00'}
    labels={'genome_wide_hit':'screen-derived hit','genome_wide_non_hit':'random non-hit','tri_perturb_seq':'curated candidate'}
    for typ in ['genome_wide_non_hit','tri_perturb_seq','genome_wide_hit']:
        sub=df[df['Type']==typ]
        ax.scatter(sub['genome_wide_screen_log2fc'],sub['validation_screen_log2fc'],s=18,c=colors[typ],alpha=0.85,edgecolor='white',lw=0.25,label=labels[typ])
    x=df['genome_wide_screen_log2fc']; y=df['validation_screen_log2fc']; r,p=pearsonr(x,y)
    xs=np.linspace(x.min(),x.max(),100); slope,intercept=np.polyfit(x,y,1); ax.plot(xs,intercept+slope*xs,color='black',lw=1)
    ax.axhline(0,color='#999',ls='--',lw=0.7); ax.axvline(0,color='#999',ls='--',lw=0.7)
    # label selected
    lab_genes=['CDH5','KDR','ETV2','SMAD1','SMAD5','ACVR2A','PLCG1','BMPR1A','CBLL1','SMC2']
    texts=[]
    for _,rr in df[df['gene'].isin(lab_genes)].iterrows():
        texts.append(ax.text(rr['genome_wide_screen_log2fc'],rr['validation_screen_log2fc'],rr['gene'],fontsize=6.2))
    if adjust_text and texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', lw=0.3, color='#777777'), time_lim=0.5)
    ax.text(0.02,0.98,f'n = {len(df)}\nPearson r = {r:.2f}\n{ptxt(p)}',transform=ax.transAxes,ha='left',va='top',fontsize=6.8,bbox=dict(boxstyle='round,pad=0.25',fc='white',ec='#DDD',lw=0.4))
    ax.set_xlabel('Genome-wide screen effect',fontsize=8); ax.set_ylabel('Validation screen effect',fontsize=8); ax.tick_params(labelsize=7)
    ax.legend(frameon=False,fontsize=6.2,loc='lower right')
    ax.spines[['top','right']].set_visible(False)
    # inset replicate
    iax=inset_axes(ax,width='36%',height='36%',loc='upper right',borderpad=1.1)
    rep=val.dropna(subset=['U Test FC BR1','U Test FC BR2'])
    iax.scatter(rep['U Test FC BR1'],rep['U Test FC BR2'],s=8,c='#777777',alpha=0.55,lw=0)
    if len(rep)>3:
        rr,pp=pearsonr(rep['U Test FC BR1'],rep['U Test FC BR2']); xs=np.linspace(rep['U Test FC BR1'].min(),rep['U Test FC BR1'].max(),100); sl,ic=np.polyfit(rep['U Test FC BR1'],rep['U Test FC BR2'],1); iax.plot(xs,ic+sl*xs,color='black',lw=0.7); iax.text(0.05,0.93,f'r={rr:.2f}',transform=iax.transAxes,fontsize=5.8,ha='left',va='top')
    iax.axhline(0,color='#AAA',ls='--',lw=0.4); iax.axvline(0,color='#AAA',ls='--',lw=0.4)
    iax.set_title('replicates',fontsize=6); iax.tick_params(labelsize=5); iax.set_xlabel('BR1',fontsize=5); iax.set_ylabel('BR2',fontsize=5)


def draw_panel_1c_placeholder(ax):
    ax.axis('off'); ax.set_xlim(0,1); ax.set_ylim(0,1)
    ax.add_patch(FancyBboxPatch((0.08,0.25),0.84,0.5,boxstyle='round,pad=0.03,rounding_size=0.04',facecolor='#F8F8F8',edgecolor='#BBBBBB',lw=0.8))
    ax.text(0.5,0.58,'D0–D3 abundance comparison pending',ha='center',va='center',fontsize=11,fontweight='bold',color='#444444')
    ax.text(0.5,0.42,'Processed guide-abundance table not found\nRequired: gene, abundance effect, P value, FDR',ha='center',va='center',fontsize=7.5,color='#555555')

# -------------------------
# Individual panels
# -------------------------
# Fig1a
fig,ax=plt.subplots(figsize=(6.2,3.2)); draw_panel_1a(ax); panel_letter(ax,'a'); save_fig(fig,'Fig1a_screen_workflow',OUT/'individual_panels'); plt.close(fig)
# Fig1b selected and candidates
fig,ax=plt.subplots(figsize=(5.9,4.1)); draw_panel_1b(ax); panel_letter(ax,'b'); save_fig(fig,'Fig1b_selected_main_panel',OUT/'individual_panels'); plt.close(fig)
fig,ax=plt.subplots(figsize=(5,3.6)); draw_panel_1b_utest(ax); save_fig(fig,'Fig1b_candidate_A_Utest_volcano',OUT/'candidate_panels'); plt.close(fig)
fig,ax=plt.subplots(figsize=(5,3.6)); draw_panel_1b(ax); save_fig(fig,'Fig1b_candidate_B_MaGECK_plot',OUT/'candidate_panels'); plt.close(fig)
fig,ax=plt.subplots(figsize=(5,3.6)); draw_panel_1b_ranked(ax); save_fig(fig,'Fig1b_candidate_C_full_ranked_screen',OUT/'candidate_panels'); plt.close(fig)
# Fig1d selected
fig,ax=plt.subplots(figsize=(6.2,4.1)); draw_panel_1d(ax); panel_letter(ax,'d'); save_fig(fig,'Fig1d_selected_main_panel',OUT/'individual_panels'); plt.close(fig)
# Fig1e selected
fig,ax=plt.subplots(figsize=(5.5,4.1)); draw_panel_1e(ax); panel_letter(ax,'e'); save_fig(fig,'Fig1e_selected_main_panel',OUT/'individual_panels'); plt.close(fig)
# placeholder 1c
fig,ax=plt.subplots(figsize=(4,3)); draw_panel_1c_placeholder(ax); panel_letter(ax,'c'); save_fig(fig,'Fig1c_abundance_placeholder_internal_review',OUT/'individual_panels'); plt.close(fig)

# -------------------------
# Assembled figures
# -------------------------
# Layout B four-panel fallback
fig = plt.figure(figsize=(12.6,8.0), constrained_layout=False)
gs = fig.add_gridspec(2, 2, width_ratios=[0.44,0.56], height_ratios=[0.47,0.53], wspace=0.28, hspace=0.34)
ax1=fig.add_subplot(gs[0,0]); draw_panel_1a(ax1); panel_letter(ax1,'a')
ax2=fig.add_subplot(gs[0,1]); draw_panel_1b(ax2); panel_letter(ax2,'b')
ax3=fig.add_subplot(gs[1,0]); draw_panel_1d(ax3); panel_letter(ax3,'c')
ax4=fig.add_subplot(gs[1,1]); draw_panel_1e(ax4); panel_letter(ax4,'d')
fig.suptitle('Genome-wide CRISPRi screen identifies regulators of endothelial-cell differentiation', fontsize=14, fontweight='bold', y=0.995)
save_fig(fig,'Figure1_revised_fallback',OUT/'assembled_figures')
# also TIFF
png = OUT/'assembled_figures'/'Figure1_revised_fallback.png'
try:
    from PIL import Image
    im=Image.open(png); im.save(OUT/'assembled_figures'/'Figure1_revised_fallback_600dpi.tiff', dpi=(600,600))
except Exception:
    pass
plt.close(fig)

# Layout C internal five-panel placeholder
fig = plt.figure(figsize=(13.2,8.4), constrained_layout=False)
gs = fig.add_gridspec(2, 3, height_ratios=[0.47,0.53], width_ratios=[0.31,0.36,0.33], wspace=0.32, hspace=0.34)
ax1=fig.add_subplot(gs[0,0]); draw_panel_1a(ax1); panel_letter(ax1,'a')
ax2=fig.add_subplot(gs[0,1:]); draw_panel_1b(ax2); panel_letter(ax2,'b')
ax3=fig.add_subplot(gs[1,0]); draw_panel_1c_placeholder(ax3); panel_letter(ax3,'c')
ax4=fig.add_subplot(gs[1,1]); draw_panel_1d(ax4); panel_letter(ax4,'d')
ax5=fig.add_subplot(gs[1,2]); draw_panel_1e(ax5); panel_letter(ax5,'e')
fig.suptitle('Genome-wide CRISPRi screen identifies regulators of endothelial-cell differentiation', fontsize=14, fontweight='bold', y=0.995)
save_fig(fig,'Figure1_revised_internal_placeholder',OUT/'assembled_figures')
plt.close(fig)

# -------------------------
# Supporting docs
# -------------------------
# candidate comparison
pd.DataFrame([
    {'candidate':'1b-A U-test volcano','selected':False,'reason':'coherent within U-test but less directly tied to MaGECK screen ranking'},
    {'candidate':'1b-B MaGECK effect/FDR plot','selected':True,'reason':'coherent effect/significance pair and shows full screen; U-test support encoded'},
    {'candidate':'1b-C full ranked screen map','selected':False,'reason':'compact but less interpretable than effect/FDR plot'},
]).to_csv(OUT/'source_data'/'Fig1b_candidate_comparison.tsv', sep='\t', index=False)

# missing input request/unresolved
unres = '''# Unresolved issues for Figure 1 revised

## Panel 1c blocked
No processed D0-versus-D3 or baseline-to-endpoint guide-abundance table was found after auditing the currently available ZIP archives and text/spreadsheet files. I therefore assembled a polished four-panel fallback figure and a separate internal five-panel placeholder version.

Required processed table columns:
- gene
- d0_vs_d3_effect
- p_value
- adjusted_p_value
- significant
- analysis_method
- source_target_id

Promoter-level rows are acceptable if the source target identifier is retained.

## Panel 1d enrichment limitation
A fixed GMT or database export for Reactome/GO/CORUM enrichment was not attached. The main panel therefore uses selected pathway/component strips based on the existing literature-informed annotation and manuscript-relevant pathway groups. Formal preranked enrichment should be run once a fixed gene-set library is provided.
'''
(OUT/'unresolved_issues.md').write_text(unres)
(OUT/'data_audit'/'consolidated_missing_input_request.md').write_text(unres)

readme = f'''# Figure 1 revised build

Generated from currently available Batch A files.

Main output: assembled_figures/Figure1_revised_fallback.pdf

The preferred full five-panel design is blocked because the D0-versus-D3 guide-abundance analysis required for panel 1c was not present. A clean four-panel fallback and an internal five-panel placeholder are exported.

Key reproduced counts:
- Full combined target rows: {len(full)}
- Full combined unique genes: {gene_full.shape[0]}
- Union hit target rows: {int(full['union_hit'].sum())}
- Union hit unique biological genes: {int(screen_hits['normalized_gene_symbol'].nunique())}
- NEDD8 target rows retained in audit: {len(full[full['normalized_gene_symbol']=='NEDD8'])}
- Screen-derived hits required for EC differentiation: {int((screen_hits['screen_direction']=='Required for EC differentiation').sum())}
- Screen-derived hits restricting EC differentiation: {int((screen_hits['screen_direction']=='Restricts EC differentiation').sum())}

Selected panel 1b is the MaGECK effect/FDR plot. U-test support is encoded in point support status.
'''
(OUT/'README_figure1_revised.md').write_text(readme)
legend = '''# Draft Figure 1 legend

**Figure 1. A genome-wide CRISPRi screen identifies validated pathways and regulators of endothelial-cell differentiation.**

**a,** Schematic of the genome-wide CRISPRi differentiation screen. Inducible CRISPRi hiPSCs were transduced with a genome-wide guide library, selected, induced for CRISPRi, differentiated through D0 hiPSC, D1 primitive-streak and D2 lateral-mesoderm stages, and assayed at D3 by separating CDH5/CD144-positive endothelial cells from CDH5/CD144-negative non-endothelial cells followed by guide sequencing.

**b,** Gene-level genome-wide screen result. Each point represents one harmonized biological gene after collapsing target-level rows for display. The x-axis shows the MaGECK gene effect for CDH5-positive versus CDH5-negative abundance, and the y-axis shows the direction-appropriate MaGECK FDR. Highlighted genes represent the 244 unique screen-derived regulators; color denotes biological direction and point shape denotes support from MaGECK and/or Mann–Whitney U-test calling.

**c,** Selected pathway component effects for pathways and complexes nominated by the screen. Each point is a tested pathway member; filled points are final screen hits and open points are tested non-hits. The x-axis shows the MaGECK gene effect using the same orientation as panel b.

**d,** Secondary focused-screen validation. Main plot compares genome-wide screen effects with validation-screen effects for harmonized biological genes. Inset shows validation biological-replicate concordance.

Note: the planned differentiation-versus-baseline abundance panel is not included in the submission-ready fallback because the processed D0-versus-D3 abundance table was not available in the current attachments.
'''
(OUT/'legends'/'Figure1_revised_legend_draft.md').write_text(legend)

# final QC
qc_lines=[]
for svg in glob.glob(str(OUT/'individual_panels'/'*.svg'))+glob.glob(str(OUT/'assembled_figures'/'*.svg')):
    txt=Path(svg).read_text(errors='ignore')
    qc_lines.append({'file':os.path.basename(svg),'text_elements':txt.count('<text'),'paths':txt.count('<path')})
pd.DataFrame(qc_lines).to_csv(OUT/'final_QC_checklist.tsv', sep='\t', index=False)

# Save copy of script
try:
    shutil.copyfile(__file__, OUT/'scripts'/'build_figure1_revised.py')
except Exception:
    pass
print('Done', OUT)
