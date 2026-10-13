import os
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt

INFILE = '/mnt/data/figure_build/cleaned_data/gw_screen_hits_244_gene_level.tsv'
OUTDIR = '/mnt/data/figure1_annotation_literature_draft'
os.makedirs(OUTDIR, exist_ok=True)

df = pd.read_csv(INFILE, sep='\t')
genes = set(df['display_gene'])

annotations = {}

def assign(primary, secondary, names, confidence='high'):
    for g in names:
        if g not in genes:
            raise ValueError(f'Gene {g} is not in the 244-gene set')
        if g in annotations:
            raise ValueError(f'Duplicate primary assignment for {g}: {annotations[g]} vs {(primary, secondary)}')
        annotations[g] = (primary, secondary, confidence)

# 1. Developmental/lineage signaling
assign('Developmental signaling & lineage regulation', 'TGF-beta/BMP signaling', {
    'ACVR1B','ACVR2A','ACVR2B','BMPR1A','LEMD3','SMAD1','SMAD3','SMAD4','SMAD5'
})
assign('Developmental signaling & lineage regulation', 'Hippo/YAP signaling', {
    'AMOTL2','LATS2','NF2','SAV1','TAOK1','WWC3'
})
assign('Developmental signaling & lineage regulation', 'VEGF/endothelial signaling', {
    'CDH5','ETV2','KDR','PLCG1'
})
assign('Developmental signaling & lineage regulation', 'Mesoderm/cardiovascular transcription factors', {
    'EOMES','HAND1','MIXL1','SOX4','TBX3','TBXT'
})
assign('Developmental signaling & lineage regulation', 'Other developmental or signal transduction', {
    'ARTN','CSNK2B','KCTD10','MEMO1','PITPNC1','PLXNA2','PPP2R1A','PRKCD','PRKD1','PTEN','PYGO2','RFX3','STK40','VHL'
}, confidence='medium')

# 2. Transcription/chromatin
assign('Transcription & chromatin regulation', 'Chromatin modification/remodeling', {
    'ASH2L','BRD3','DPY30','EHMT2','HDAC2','KDM1A','KDM3B','MBD1','PRMT5','RBBP5','RESF1','RNF20','RNF40','SETD1B','SETDB1','SGF29','SMARCA4','SMARCA5','SMARCB1','SMARCC1','SMARCD3','SMARCE1','SPIN1','SSRP1','TADA2B','TAF6L','TASOR','TRIM24','TRRAP','WAC','WDR77','WIZ'
})
assign('Transcription & chromatin regulation', 'Mediator/general transcription machinery', {
    'CCNC','CCNT2','CTR9','ELOA','GTF2A1','GTF2A2','GTF2H4','GTF3C2','MED1','MED12','MED15','MED19','POLR2M','POLR3H','POLR3K','SPTY2D1','SUPT20H','TAF10','TAF2','TAF3','TAF4','TBP','UXT'
})
assign('Transcription & chromatin regulation', 'Transcriptional cofactors/repressors', {
    'PAXBP1','PPHLN1','RLF','ZBTB22','ZBTB48'
}, confidence='medium')

# 3. RNA processing
assign('RNA processing & surveillance', 'Splicing and RNA-binding', {
    'ATXN2L','BCAS2','CD2BP2','DDX42','DNAJC17','EWSR1','HNRNPM','HNRNPU','LARP4','LENG8','LSM12','PPIL4','PTBP1','RBM14-RBM4','RBM33','RBM4','SCAF4','SCAF8','SF1','SNRNP25','SNRNP27','SNUPN','SRSF6','SRRM2','SSB'
})
assign('RNA processing & surveillance', 'miRNA/RNA decay and surveillance', {
    'CNOT2','DBR1','DGCR8','DROSHA','EXOSC8','G3BP1','INTS10','INTS14','NCBP1','RBM7','SKIC8','SMG5','XPO5','XRN1','ZC3H8'
})

# 4. Translation/ribosome
assign('Translation, ribosome & protein targeting', 'Translation initiation/elongation', {
    'EEF2','EIF3F','EIF4E','EIF4G1','EIF4G2','EIF5A','RACK1'
})
assign('Translation, ribosome & protein targeting', 'Ribosome biogenesis/ribosomal proteins', {
    'CCT6B','LAGE3','MRPL32','NUFIP2','RPS15A','WDR12'
})
assign('Translation, ribosome & protein targeting', 'Signal recognition particle/protein targeting', {
    'BAG6','SRP19','SRP54','SRP68','SRPRA','SRPRB'
})

# 5. Proteostasis/PTM
assign('Proteostasis & post-translational modification', 'Ubiquitin/NEDD8/proteasome', {
    'ARIH1','CBLL1','COPS4','COPS6','CUL3','FBXW11','HUWE1','NEDD8','PSMB3','PSMD4','PSMG4','PTAR1','RBX1','RNF115','UFM1'
})
assign('Proteostasis & post-translational modification', 'Protein modification/chaperone/proteolysis', {
    'CAPN5','CHORDC1','NAA20','NAA25','OGT','RABGGTB','SPPL2B'
})

# 6. Cell cycle/DNA
assign('Cell cycle, chromosome & DNA repair', 'Chromosome segregation/cohesin/condensin', {
    'CENPC','MAU2','RAD21','SKA3','SMC2','STAG2','WAPL'
})
assign('Cell cycle, chromosome & DNA repair', 'DNA replication/repair and genome stability', {
    'ASCC3','ATMIN','ERCC2','ITPRID2','POLE4','TOP1','TLE6'
})

# 7. Trafficking/cytoskeleton
assign('Membrane trafficking, cytoskeleton & adhesion', 'Vesicle trafficking and membrane organization', {
    'ATP6AP1','CACFD1','CHMP3','CNEP1R1','EMC2','GPR107','LHFPL2','NUP93','RAB18','RAB6A','RGP1','RIC1','TBC1D9','TMTC3','VPS54'
})
assign('Membrane trafficking, cytoskeleton & adhesion', 'Cytoskeleton, cilia and cell adhesion', {
    'ACTL7B','CACNG1','CD99L2','CFAP119','IPO8','RIBC1'
}, confidence='medium')

# 8. Metabolism
assign('Metabolism, mitochondria & glycosylation', 'Mitochondrial/peroxisomal metabolism', {
    'ACAA1','ATP5PB','HINT2','PAM16'
})
assign('Metabolism, mitochondria & glycosylation', 'Lipid/carbohydrate metabolism', {
    'GCK','GPI','HSD17B12','NCEH1','OSBP','PBLD','SELENOI','TSTD2'
})
assign('Metabolism, mitochondria & glycosylation', 'Glycosylation and sugar transport', {
    'DPM1','DPM2','DPM3','ERO1A','FUT7','MANBAL','SLC29A1','SLC35C1'
})

# 9. Apoptosis/stress
assign('Apoptosis & stress response', 'Intrinsic apoptosis', {
    'APAF1','BAD','BAX','CASP3','CYCS'
})

# 10. Other / low-confidence
assign('Other / insufficiently characterized', 'Poorly characterized or multifunctional', {
    'C5orf22','KIAA0232','PEG10','PRR30'
}, confidence='low')

if set(annotations) != genes:
    missing = sorted(genes - set(annotations))
    extra = sorted(set(annotations) - genes)
    raise RuntimeError(f'Annotation coverage mismatch. Missing={missing}; extra={extra}')

# Literature-supported regulator evidence table.
# Tiers are deliberately conservative and distinguish direct endothelial evidence from broader developmental/pathway evidence.
known_rows = [
    # Tier 1: direct endothelial specification/function
    ('ETV2','Tier 1 - direct endothelial specification','Master/pioneer factor for endothelial lineage specification and reprogramming','20235229;32832668','10.1038/s41556-022-00901-3'),
    ('KDR','Tier 1 - direct endothelial specification','VEGFR2 signaling specifies endothelial fate from vascular progenitors','18391074;19706681',''),
    ('PLCG1','Tier 1 - direct endothelial specification','VEGFR2-PLCgamma1 relay is required for endothelial specification and vasculogenesis','19706681;11744703',''),
    ('CDH5','Tier 1 - endothelial identity/function','VE-cadherin is an endothelial-selective junctional regulator and manuscript-defined positive control','8770858',''),
    ('SMARCA4','Tier 1 - direct endothelial co-regulator','BRG1/SMARCA4 is recruited by ETV2 to open endothelial chromatin','', '10.1038/s41556-022-00901-3'),
    ('AMOTL2','Tier 1 - endothelial morphogenesis/angiogenesis','Required for endothelial migration, proliferation, polarity and aortic lumen expansion','21937427;24806444',''),
    ('LATS2','Tier 1 - endothelial fate maintenance','Endothelial LATS2 restrains endothelial-to-mesenchymal transition and fibrosis','39155965','10.1038/s44161-024-00508-x'),
    ('KCTD10','Tier 1 - vascular development','Required for embryonic angiogenesis and blood-vessel development; regulates Notch','25401743;24705121',''),
    # Tier 2: pathway/mesoderm/vascular development
    ('BMPR1A','Tier 2 - endothelial differentiation pathway','BMP receptor pathway supports mesoderm and endothelial differentiation from pluripotent stem cells','16846035','10.1038/ncb3205'),
    ('SMAD1','Tier 2 - endothelial differentiation pathway','BMP-SMAD1/5 signaling is activated during vascular/endothelial differentiation','16846035','10.1038/nbt.1605'),
    ('SMAD4','Tier 2 - endothelial differentiation pathway','Common SMAD required for TGF-beta/BMP developmental signaling relevant to vascular differentiation','16846035','10.1038/nbt.1605'),
    ('SMAD5','Tier 2 - endothelial differentiation pathway','BMP-SMAD1/5 signaling is activated during vascular/endothelial differentiation','16846035','10.1038/nbt.1605'),
    ('ACVR1B','Tier 2 - mesoderm/vascular pathway','Activin/TGF-beta receptor signaling controls early lineage patterning and endothelial differentiation context','16846035','10.1038/nbt.1605'),
    ('ACVR2A','Tier 2 - mesoderm/vascular pathway','Activin/TGF-beta receptor signaling controls early lineage patterning and endothelial differentiation context','16846035','10.1038/nbt.1605'),
    ('ACVR2B','Tier 2 - mesoderm/vascular pathway','Activin/TGF-beta receptor signaling controls early lineage patterning and endothelial differentiation context','16846035','10.1038/nbt.1605'),
    ('SMAD3','Tier 2 - mesoderm/vascular pathway','TGF-beta/Activin signaling regulates pluripotency, mesoderm patterning and endothelial expansion','16846035','10.1038/nbt.1605'),
    ('EOMES','Tier 2 - early mesoderm regulator','Primitive-streak/mesendoderm transcription factor upstream of cardiovascular/endothelial competence','27934856;27906170',''),
    ('MIXL1','Tier 2 - early mesoderm regulator','Mesendoderm/primitive-streak regulator in pluripotent-stem-cell differentiation','25004077;27906170',''),
    ('TBXT','Tier 2 - early mesoderm regulator','BRACHYURY is required for BMP-driven mesoderm induction preceding endothelial differentiation','21816365',''),
    ('HAND1','Tier 2 - cardiovascular mesoderm regulator','Regulates cardiovascular/extraembryonic mesoderm progenitor specification','40164946',''),
    ('VHL','Tier 2 - vascular oxygen-sensing regulator','VHL-HIF-VEGF axis is required for normal vascular development and patterning','22627278;42003386',''),
    ('PTEN','Tier 2 - endothelial angiogenesis regulator','Controls endothelial angiogenesis, permeability and Notch-dependent stalk-cell arrest','20705603','10.1038/ncomms8935'),
    # Tier 3: contextual endothelial/angiogenesis evidence or pathway-core genes
    ('NF2','Tier 3 - contextual Hippo/endothelial','Core Hippo/YAP regulator; endothelial relevance inferred from Hippo/YAP angiogenesis and mechanotransduction literature','',''),
    ('SAV1','Tier 3 - contextual Hippo/endothelial','Core Hippo scaffold; pathway-level endothelial relevance rather than direct lineage evidence','',''),
    ('TAOK1','Tier 3 - contextual Hippo/endothelial','Hippo pathway kinase; direct endothelial evidence remains limited','',''),
    ('WWC3','Tier 3 - contextual Hippo/endothelial','WWC family Hippo regulator; endothelial evidence is pathway-contextual','',''),
    ('SOX4','Tier 3 - contextual angiogenesis','SOX4 can promote angiogenic transcriptional programs in tumor-associated endothelium','30507376',''),
    ('LEMD3','Tier 3 - contextual BMP/TGF regulation','Nuclear-envelope antagonist of BMP/TGF-beta signaling; pathway-level relevance','',''),
]
known = pd.DataFrame(known_rows, columns=['gene','known_EC_evidence_tier','literature_basis','PMID','DOI'])
known['known_EC_regulator'] = known['known_EC_evidence_tier'].str.startswith(('Tier 1','Tier 2'))
known['literature_supported_vascular_or_lineage_regulator'] = True

# Label recommendations: strict manuscript labels are priority 1; direct literature and central narrative are priority 2.
priority1 = {'CDH5','KDR','ETV2'}
priority2 = {'PLCG1','SMARCA4','AMOTL2','LATS2','KCTD10','BMPR1A','SMAD1','SMAD4','SMAD5','ACVR2A','NF2','SAV1','TAOK1'}
known['must_label'] = known['gene'].isin(priority1)
known['label_priority'] = np.select([known['gene'].isin(priority1), known['gene'].isin(priority2)], [1,2], default=3)
known['recommended_use'] = np.select(
    [known['label_priority']==1, known['label_priority']==2],
    ['Main Figure 1b/1c label','Optional label if space or manuscript narrative requires'],
    default='Review/secondary annotation'
)

# Merge functional annotations into screen table.
ann = df.copy()
ann[['primary_category','secondary_category','annotation_confidence']] = ann['display_gene'].map(annotations).apply(pd.Series)
ann = ann.merge(known[['gene','known_EC_regulator','literature_supported_vascular_or_lineage_regulator','known_EC_evidence_tier','literature_basis','PMID','DOI','must_label','label_priority','recommended_use']],
                left_on='display_gene', right_on='gene', how='left', suffixes=('','_known'))
ann.drop(columns=['gene_known'], inplace=True)
ann['known_EC_regulator'] = ann['known_EC_regulator'].fillna(False)
ann['literature_supported_vascular_or_lineage_regulator'] = ann['literature_supported_vascular_or_lineage_regulator'].fillna(False)
ann['known_EC_evidence_tier'] = ann['known_EC_evidence_tier'].fillna('No specific endothelial/vascular evidence curated')
ann['literature_basis'] = ann['literature_basis'].fillna('Broad molecular-function assignment; no gene-specific endothelial literature curated in this draft')
ann['PMID'] = ann['PMID'].fillna('')
ann['DOI'] = ann['DOI'].fillna('')
ann['must_label'] = ann['must_label'].fillna(False)
ann['label_priority'] = ann['label_priority'].fillna(4).astype(int)
ann['recommended_use'] = ann['recommended_use'].fillna('No main-figure label recommended')

# Broad essential/viability flag is intentionally cautious.
likely_essential = {
    'RNA processing & surveillance','Translation, ribosome & protein targeting',
    'Proteostasis & post-translational modification','Cell cycle, chromosome & DNA repair',
    'Apoptosis & stress response'
}
possible_essential = {'Transcription & chromatin regulation','Metabolism, mitochondria & glycosylation','Membrane trafficking, cytoskeleton & adhesion'}
ann['broad_essential_or_viability'] = np.select(
    [ann['primary_category'].isin(likely_essential), ann['primary_category'].isin(possible_essential)],
    ['likely - core cellular machinery','possible - context dependent'],
    default='not primary interpretation'
)
ann['predefined_positive_control'] = ann['display_gene'].isin(priority1)
ann['annotation_status'] = np.where(ann['annotation_confidence']=='low','draft - low confidence / author review','draft - literature-informed / author review')
ann['notes'] = np.where(
    ann['annotation_confidence']=='low',
    'Primary category is provisional because the gene is poorly characterized or multifunctional.',
    'Primary category is a broad literature-informed molecular-function assignment; secondary category should be reviewed by the authors.'
)

# Select and order columns for deliverable.
out_cols = [
    'display_gene','original_gene','original_target_ids','screen_direction',
    'primary_category','secondary_category','annotation_confidence',
    'known_EC_regulator','literature_supported_vascular_or_lineage_regulator','known_EC_evidence_tier',
    'predefined_positive_control','broad_essential_or_viability','must_label','label_priority','recommended_use',
    'literature_basis','PMID','DOI','annotation_status','notes',
    'u_test_effect_log2_CDH5pos_vs_CDH5neg','u_test_fdr','rank_by_u_test_effect_ascending'
]
out = ann[out_cols].rename(columns={'display_gene':'gene'})
out_path = os.path.join(OUTDIR,'screen_regulators_244_literature_annotation_draft.tsv')
out.to_csv(out_path, sep='\t', index=False)

known_path = os.path.join(OUTDIR,'known_endothelial_regulators_literature_review.tsv')
known.sort_values(['label_priority','known_EC_evidence_tier','gene']).to_csv(known_path, sep='\t', index=False)

# Category summary for Fig 1d.
summary = (out.groupby(['primary_category','screen_direction'])
           .size().rename('n_unique_genes').reset_index())
totals = summary.groupby('primary_category')['n_unique_genes'].sum().sort_values(ascending=False)
summary['category_total'] = summary['primary_category'].map(totals)
summary = summary.sort_values(['category_total','primary_category','screen_direction'], ascending=[False,True,True])
summary_path = os.path.join(OUTDIR,'figure1d_category_summary.tsv')
summary.to_csv(summary_path, sep='\t', index=False)

# Wide summary.
wide = summary.pivot_table(index='primary_category',columns='screen_direction',values='n_unique_genes',fill_value=0).astype(int)
wide['total'] = wide.sum(axis=1)
wide = wide.sort_values('total',ascending=True)
wide.to_csv(os.path.join(OUTDIR,'figure1d_category_summary_wide.tsv'), sep='\t')

# Draft diverging bar panel.
required_label = 'Required for EC differentiation'
restrict_label = 'Restricts EC differentiation'
req = wide.get(required_label, pd.Series(0,index=wide.index))
res = wide.get(restrict_label, pd.Series(0,index=wide.index))

fig, ax = plt.subplots(figsize=(8.6,5.8), dpi=220)
y = np.arange(len(wide))
ax.barh(y, -req.values, color='#2F76B7', edgecolor='white', linewidth=0.5, label='Required for EC differentiation')
ax.barh(y, res.values, color='#D95F4A', edgecolor='white', linewidth=0.5, label='Restricts EC differentiation')
ax.axvline(0,color='#555555',linewidth=0.8)
ax.set_yticks(y)
ax.set_yticklabels(wide.index, fontsize=10)
ax.set_xlabel('Number of unique screen-derived regulators', fontsize=11)
ax.set_title('Functional organization of the 244 screen-derived regulators', fontsize=14, fontweight='bold', pad=12)
# symmetric-ish labels; left values displayed positive
max_left = int(req.max())
max_right = max(int(res.max()), 2)
ax.set_xlim(-(max_left+7), max_right+5)
xticks = ax.get_xticks()
ax.set_xticklabels([str(abs(int(x))) if float(x).is_integer() else f'{abs(x):.1f}' for x in xticks])
for i,(a,b) in enumerate(zip(req.values,res.values)):
    if a>0:
        ax.text(-a-0.8,i,str(int(a)),va='center',ha='right',fontsize=9,color='#204F7A')
    if b>0:
        ax.text(b+0.35,i,str(int(b)),va='center',ha='left',fontsize=9,color='#9E3427')
ax.legend(frameon=False,loc='lower left',bbox_to_anchor=(0,-0.18),ncol=2,fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.grid(axis='x',color='#E8E8E8',linewidth=0.6,zorder=0)
ax.set_axisbelow(True)
fig.tight_layout()
for ext in ['png','pdf','svg']:
    fig.savefig(os.path.join(OUTDIR,f'Fig1d_functional_categories_literature_draft.{ext}'), bbox_inches='tight', dpi=350 if ext=='png' else None)
plt.close(fig)

# Reference table used for literature review.
refs = pd.DataFrame([
    ['ETV2 endothelial specification','20235229;32832668','10.1038/s41556-022-00901-3'],
    ['VEGFR2/KDR and PLCG1 endothelial specification','18391074;19706681;11744703',''],
    ['BMP/TGF-beta signaling in PSC-derived endothelial differentiation','16846035','10.1038/nbt.1605;10.1038/ncb3205'],
    ['Mesoderm/primitive-streak patterning before endothelial differentiation','21816365;27906170;27934856',''],
    ['AMOTL2 endothelial angiogenesis and lumen morphogenesis','21937427;24806444',''],
    ['Endothelial LATS2 and endothelial-to-mesenchymal transition','39155965','10.1038/s44161-024-00508-x'],
    ['KCTD10 embryonic angiogenesis and blood-vessel development','25401743;24705121',''],
    ['PTEN regulation of endothelial angiogenesis','20705603','10.1038/ncomms8935'],
    ['VHL/HIF regulation of vascular development','22627278;42003386',''],
    ['SOX4 contextual angiogenesis evidence','30507376',''],
], columns=['evidence_group','PMID','DOI'])
refs.to_csv(os.path.join(OUTDIR,'literature_reference_index.tsv'),sep='\t',index=False)

# README
readme = f'''# Figure 1d literature-informed annotation draft\n\n## Scope\n\nThis directory contains a draft functional annotation of all {len(out)} unique screen-derived regulators. The source set is the cleaned 244-gene table, with the duplicate NEDD8 target collapsed at the biological-gene level.\n\n## Annotation approach\n\n- Primary categories are broad molecular-function/process categories designed for a compact Figure 1d summary.\n- Categories are literature-informed and organized using the general logic of curated resources such as Gene Ontology and Reactome, but this is not a formal GO enrichment analysis.\n- Direct gene-specific literature review was performed for a conservative subset of established endothelial, vascular, mesodermal and pathway regulators.\n- Tier 1 indicates direct endothelial specification, identity or endothelial morphogenesis evidence.\n- Tier 2 indicates established vascular-development, mesoderm or endothelial-differentiation pathway evidence.\n- Tier 3 indicates contextual endothelial/angiogenesis or pathway-core evidence without strong direct endothelial-lineage evidence.\n- Genes without a curated endothelial literature entry are not assumed to lack endothelial relevance.\n\n## Author review required\n\nThe table is a draft. Secondary categories, broad-essential flags and Tier 3 assignments should be reviewed by the authors before the panel is finalized. The four low-confidence genes are C5orf22, KIAA0232, PEG10 and PRR30.\n\n## Main files\n\n- screen_regulators_244_literature_annotation_draft.tsv\n- known_endothelial_regulators_literature_review.tsv\n- figure1d_category_summary.tsv\n- Fig1d_functional_categories_literature_draft.svg/pdf/png\n- literature_reference_index.tsv\n'''
with open(os.path.join(OUTDIR,'README_annotation_method.md'),'w') as f:
    f.write(readme)

print('Wrote',OUTDIR)
print('Category totals:')
print(wide.sort_values('total',ascending=False).to_string())
print('\nKnown regulator tiers:')
print(known.groupby('known_EC_evidence_tier').size().to_string())
