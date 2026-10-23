# Optional single-cell signature distributions for selected perturbations
# Requires adata_int loaded with per-cell program usages and target carrier labels.
import scanpy as sc
SELECTED_TARGETS = ['NF2','LATS2','SAV1','TAOK1','AMOTL2','PRKCD','SMAD1','SMAD4','SMAD5','KDR','ETV2']
CONTROL_LABEL = 'NT_control'
carrier_col = 'target_carrier'  # replace with your actual per-cell target/guide-carrier column
EC_PROGRAMS = ['Program 2','Program 3','Program 8','Program 11']
ALT_PROGRAMS = ['Program 34','Program 27']
adata_int.obs['EC_endpoint_score'] = adata_int.obs[EC_PROGRAMS].mean(axis=1)
adata_int.obs['ALT_endpoint_score'] = adata_int.obs[ALT_PROGRAMS].mean(axis=1)
adata_int.obs['P34_score'] = adata_int.obs['Program 34']
adata_int.obs['P17_score'] = adata_int.obs['Program 17']
adata_int.obs['P15_score'] = adata_int.obs['Program 15']
for score in ['EC_endpoint_score','ALT_endpoint_score','P34_score','P17_score','P15_score']:
    sc.pl.violin(adata_int[adata_int.obs[carrier_col].isin(SELECTED_TARGETS+[CONTROL_LABEL])], keys=score, groupby=carrier_col, rotation=90, stripplot=False)
