#!/usr/bin/env python3
import argparse, os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HIPPO_TARGETS_DEFAULT = ['NF2','LATS2','AMOTL2','TAOK1','SAV1']
MAIN_PROGRAMS = [34,27,24,2]
SUPPLEMENT_PROGRAMS = [34,27,24,35,2,3,8,11]
SECONDARY_ALT_PROGRAMS = [27,24]
ALT_CONTEXT_PROGRAMS = [34,27,24]
HIPPO_COLOR = '#D55E00'; POINT_GREY='#BDBDBD'; BLUE='#2166AC'; RED='#B2182B'
PROGRAM_LABELS_MAIN = {34:'P34\nYAP/TAZ\nmesenchyme',27:'P27\nalt. context',24:'P24\nmesenchymal\ncontext',2:'P2\narterial EC'}
PROGRAM_LABELS_SUPP = {34:'P34',27:'P27',24:'P24',35:'P35\nexcluded',2:'P2',3:'P3',8:'P8',11:'P11'}

def savefig(fig,outdir,stem,dpi=350):
    outdir=Path(outdir); fig.savefig(outdir/f'{stem}.png',dpi=dpi,bbox_inches='tight'); fig.savefig(outdir/f'{stem}.pdf',bbox_inches='tight'); plt.close(fig)

def clean_axes(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.tick_params(axis='both', labelsize=10); return ax

def get_nt_baseline(ec):
    control=ec[ec['row_type'].eq('control_group')].copy(); nt=ec[ec['row_type'].eq('individual_NT_guide_excl_Hippo')].copy()
    if len(control):
        row=control.iloc[0]; return float(row['pct_EC_D3']), int(row['n_D3_cells']), int(row['n_D3_EC_cells']), nt, control
    n=int(nt['n_D3_cells'].sum()); ec_n=int(nt['n_D3_EC_cells'].sum()); pct=100*ec_n/n if n else np.nan; return pct,n,ec_n,nt,control

def summarize_ec(ec, targets):
    guide=ec[ec['row_type'].eq('targeting_guide') & ec['target'].isin(targets)].copy()
    rows=[]
    for t in targets:
        sub=guide[guide['target'].eq(t)].copy()
        if len(sub)==0: continue
        n=int(sub['n_D3_cells'].sum()); e=int(sub['n_D3_EC_cells'].sum()); ne=int(sub['n_D3_nonEC_cells'].sum())
        rows.append({'target':t,'n_guides':int(sub['guide_name'].nunique()),'n_D3_cells':n,'n_D3_EC_cells':e,'n_D3_nonEC_cells':ne,'pct_EC_D3_pooled':100*e/n if n else np.nan,'mean_guide_pct_EC_D3':sub['pct_EC_D3'].mean(),'median_guide_pct_EC_D3':sub['pct_EC_D3'].median()})
    return pd.DataFrame(rows), guide

def pivot_effects(pert, targets, programs):
    d3=pert[pert['day'].str.upper().eq('D3') & pert['target'].isin(targets) & pert['program'].isin(programs)].copy()
    effect=d3.pivot_table(index='target',columns='program',values='approx_log2FC',aggfunc='mean').reindex(index=targets,columns=programs).dropna(how='all')
    fdr=d3.pivot_table(index='target',columns='program',values='adj_pval',aggfunc='min').reindex(index=effect.index,columns=programs)
    return effect,fdr,d3

def compute_comp(effect):
    comp=pd.DataFrame(index=effect.index)
    comp['P34_YAP_TAZ_mesenchyme']=effect[34] if 34 in effect.columns else np.nan
    secondary=[p for p in SECONDARY_ALT_PROGRAMS if p in effect.columns]
    comp['P27_P24_secondary_alt_mean']=effect[secondary].mean(axis=1) if secondary else np.nan
    comp['P2_arterial_EC']=effect[2] if 2 in effect.columns else np.nan
    comp['P34_minus_P2']=comp['P34_YAP_TAZ_mesenchyme']-comp['P2_arterial_EC']
    alt=[p for p in ALT_CONTEXT_PROGRAMS if p in effect.columns]
    comp['mean_P34_P27_P24_minus_P2']=effect[alt].mean(axis=1)-comp['P2_arterial_EC'] if alt else np.nan
    return comp

def plot_ec_guides(guide,nt,target_summary,nt_pct,outdir,targets):
    cats=['NT guides']+[t for t in targets if t in set(guide['target'])]; xmap={c:i for i,c in enumerate(cats)}
    fig,ax=plt.subplots(figsize=(8.8,5.2)); rng=np.random.default_rng(1)
    if len(nt):
        ntp=nt.dropna(subset=['pct_EC_D3']); x=np.full(len(ntp),xmap['NT guides'],float)+rng.uniform(-.18,.18,len(ntp)); sizes=15+70*(ntp['n_D3_cells'].fillna(0)/max(ntp['n_D3_cells'].max(),1)); ax.scatter(x,ntp['pct_EC_D3'],s=sizes,color=POINT_GREY,alpha=.55,linewidths=0)
    for t in targets:
        if t not in xmap: continue
        sub=guide[guide['target'].eq(t)].dropna(subset=['pct_EC_D3'])
        if len(sub)==0: continue
        x=np.full(len(sub),xmap[t],float)+rng.uniform(-.16,.16,len(sub)); sizes=25+105*(sub['n_D3_cells'].fillna(0)/max(guide['n_D3_cells'].max(),1)); ax.scatter(x,sub['pct_EC_D3'],s=sizes,color=HIPPO_COLOR,alpha=.82,linewidths=0)
    ax.scatter([xmap['NT guides']],[nt_pct],marker='D',s=80,color='black',zorder=5)
    for _,row in target_summary.iterrows():
        if row['target'] in xmap: ax.scatter([xmap[row['target']]],[row['pct_EC_D3_pooled']],marker='D',s=90,color='black',zorder=5)
    ax.axhline(nt_pct,ls='--',color='black',lw=1,alpha=.65); ax.set_xticks(range(len(cats))); ax.set_xticklabels(cats,rotation=35,ha='right'); ax.set_ylabel('D3 EC cells (%)'); ax.set_title('Hippo perturbations reduce D3 endothelial output',fontweight='bold'); ax.set_ylim(0,max(100,np.nanmax(guide['pct_EC_D3'])+5)); clean_axes(ax); fig.tight_layout(); savefig(fig,outdir,'01_D3_EC_percent_by_individual_guide_Hippo_vs_NT')

def plot_ec_delta(summary,nt_pct,outdir,targets):
    df=summary.copy(); df['delta_EC_pct_vs_NT']=df['pct_EC_D3_pooled']-nt_pct; df['EC_loss_pp_vs_NT']=nt_pct-df['pct_EC_D3_pooled']; df['target']=pd.Categorical(df['target'],categories=[t for t in targets if t in set(df['target'])],ordered=True); df=df.sort_values('target')
    fig,ax=plt.subplots(figsize=(6.3,4.6)); x=np.arange(len(df)); ax.bar(x,df['delta_EC_pct_vs_NT'],color=HIPPO_COLOR,alpha=.85); ax.axhline(0,color='black',lw=1)
    for i,v in enumerate(df['delta_EC_pct_vs_NT']): ax.text(i,v-1.5 if v<0 else v+1.5,f'{v:.1f}',ha='center',va='top' if v<0 else 'bottom',fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels(df['target'],rotation=35,ha='right'); ax.set_ylabel('D3 EC percentage-point change vs NT'); ax.set_title('D3 endothelial output change by Hippo target',fontweight='bold'); clean_axes(ax); fig.tight_layout(); savefig(fig,outdir,'02_D3_EC_percent_delta_vs_NT_by_Hippo_target'); return df

def heatmap_matrix(mat,title,outdir,stem,fdr=None,col_labels=None,cbar_label='D3 approx_log2FC'):
    mat=mat.copy(); vals=mat.values.astype(float); row_labels=mat.index.astype(str).tolist(); col_labels=col_labels or [str(c) for c in mat.columns]
    maxabs=np.nanmax(np.abs(vals)); maxabs=1 if not np.isfinite(maxabs) or maxabs==0 else maxabs; norm=TwoSlopeNorm(vmin=-maxabs,vcenter=0,vmax=maxabs)
    fig,ax=plt.subplots(figsize=(max(5,0.95*len(col_labels)+2.4),max(2.8,0.45*len(row_labels)+1.4))); im=ax.imshow(vals,aspect='auto',cmap='RdBu_r',norm=norm)
    ax.set_xticks(np.arange(len(col_labels))); ax.set_xticklabels(col_labels,rotation=35,ha='right'); ax.set_yticks(np.arange(len(row_labels))); ax.set_yticklabels(row_labels); ax.set_title(title,fontweight='bold')
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            v=vals[i,j]
            if np.isfinite(v):
                lab=f'{v:.2f}'
                if fdr is not None:
                    try:
                        q=fdr.loc[mat.index[i],mat.columns[j]]
                        if np.isfinite(q) and q<0.05: lab+='*'
                    except Exception: pass
                ax.text(j,i,lab,ha='center',va='center',fontsize=8,color='black')
    cb=fig.colorbar(im,ax=ax,fraction=.046,pad=.02); cb.set_label(cbar_label)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.tick_params(length=0); fig.tight_layout(); savefig(fig,outdir,stem)

def plot_scatter(ec_cnmf,x_col,x_label,stem,title,outdir):
    df=ec_cnmf.dropna(subset=[x_col,'EC_loss_pp_vs_NT']).copy(); fig,ax=plt.subplots(figsize=(5.5,4.8)); ax.scatter(df[x_col],df['EC_loss_pp_vs_NT'],s=90,color=HIPPO_COLOR,alpha=.88,linewidths=0)
    for _,row in df.iterrows(): ax.text(row[x_col]+0.04,row['EC_loss_pp_vs_NT']+0.3,row['target'],fontsize=10,fontweight='bold')
    if len(df)>=3:
        r,p=stats.pearsonr(df[x_col],df['EC_loss_pp_vs_NT']); xs=np.linspace(df[x_col].min(),df[x_col].max(),100); slope,intercept,*_=stats.linregress(df[x_col],df['EC_loss_pp_vs_NT']); ax.plot(xs,intercept+slope*xs,color='black',lw=1); ax.text(.04,.96,f'r = {r:.2f}\np = {p:.2g}',transform=ax.transAxes,va='top',fontsize=10)
    ax.set_xlabel(x_label); ax.set_ylabel('D3 EC loss vs NT\npercentage points'); ax.set_title(title,fontweight='bold'); clean_axes(ax); fig.tight_layout(); savefig(fig,outdir,stem)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--ec-guide-csv',default='/mnt/data/d3_EC_percent_by_guide_with_individual_NT_excl_Hippo.csv')
    parser.add_argument('--cnmf-xlsx',default='/mnt/data/cNMF_50_0_2_032526_6co.xlsx')
    parser.add_argument('--outdir',default='/mnt/data/d3_EC_percent_cnmf_panel_updated_programs')
    parser.add_argument('--targets',default=','.join(HIPPO_TARGETS_DEFAULT))
    args=parser.parse_args(); outdir=Path(args.outdir); outdir.mkdir(parents=True,exist_ok=True); targets=[x.strip() for x in args.targets.split(',') if x.strip()]
    ec=pd.read_csv(args.ec_guide_csv); nt_pct,nt_n,nt_ec,nt,control=get_nt_baseline(ec); summary,guide=summarize_ec(ec,targets); summary['nt_pct_EC_D3']=nt_pct; summary['delta_EC_pct_vs_NT']=summary['pct_EC_D3_pooled']-nt_pct; summary['EC_loss_pp_vs_NT']=nt_pct-summary['pct_EC_D3_pooled']; summary['target']=pd.Categorical(summary['target'],categories=targets,ordered=True); summary=summary.sort_values('target'); guide['target']=pd.Categorical(guide['target'],categories=targets,ordered=True); guide=guide.sort_values(['target','guide_name'])
    summary.to_csv(outdir/'source_Hippo_D3_EC_percent_target_summary.csv',index=False); guide.to_csv(outdir/'source_D3_EC_percent_by_guide_for_plot.csv',index=False); nt.to_csv(outdir/'source_individual_NT_guides_excl_Hippo.csv',index=False); control.to_csv(outdir/'source_control_group_NT_excl_Hippo.csv',index=False)
    print('reading cnmf...')
    pert=pd.read_excel(args.cnmf_xlsx,sheet_name='Perturbation Association 1')
    pert=pert.rename(columns={'target_name':'target','program_name':'program','Sample':'day'}); pert['target']=pert['target'].astype(str); pert['day']=pert['day'].astype(str); pert['program']=pd.to_numeric(pert['program'],errors='coerce').astype('Int64'); pert['approx_log2FC']=pd.to_numeric(pert['approx_log2FC'],errors='coerce'); pert['adj_pval']=pd.to_numeric(pert['adj_pval'],errors='coerce')
    effect_main,fdr_main,long_main=pivot_effects(pert,targets,MAIN_PROGRAMS); effect_main.to_csv(outdir/'source_D3_cNMF_Hippo_main_program_approx_log2FC_P34_P27_P24_P2.csv'); fdr_main.to_csv(outdir/'source_D3_cNMF_Hippo_main_program_adj_pval_P34_P27_P24_P2.csv'); long_main.to_csv(outdir/'source_D3_cNMF_Hippo_main_program_long_P34_P27_P24_P2.csv',index=False)
    effect_supp,fdr_supp,long_supp=pivot_effects(pert,targets,SUPPLEMENT_PROGRAMS); effect_supp.to_csv(outdir/'source_D3_cNMF_Hippo_supplement_program_approx_log2FC.csv'); fdr_supp.to_csv(outdir/'source_D3_cNMF_Hippo_supplement_program_adj_pval.csv'); long_supp.to_csv(outdir/'source_D3_cNMF_Hippo_supplement_program_long.csv',index=False)
    comp=compute_comp(effect_main); comp.to_csv(outdir/'source_D3_cNMF_Hippo_composite_program_scores_updated.csv')
    ec_cnmf=summary.merge(comp.reset_index().rename(columns={'index':'target'}),on='target',how='left'); ec_cnmf.to_csv(outdir/'source_D3_EC_loss_vs_updated_cNMF_effects.csv',index=False)
    print('plotting...')
    plot_ec_guides(guide,nt,summary,nt_pct,outdir,targets); plot_ec_delta(summary,nt_pct,outdir,targets)
    heatmap_matrix(effect_main,'D3 program effects of Hippo perturbations',outdir,'03_D3_cNMF_Hippo_program_effect_heatmap_main_P34_P27_P24_P2',fdr=fdr_main,col_labels=[PROGRAM_LABELS_MAIN.get(c,f'P{c}') for c in effect_main.columns])
    labels={'P34_YAP_TAZ_mesenchyme':'P34\nYAP/TAZ\nmesenchyme','P27_P24_secondary_alt_mean':'P27/P24\nsecondary\nalt mean','P2_arterial_EC':'P2\narterial EC','P34_minus_P2':'P34 − P2','mean_P34_P27_P24_minus_P2':'mean(P34,P27,P24)\n− P2'}
    comp_disp=comp[[c for c in labels if c in comp.columns]].copy(); comp_disp.columns=[labels[c] for c in comp_disp.columns]; heatmap_matrix(comp_disp,'Composite D3 program effects',outdir,'04_D3_cNMF_Hippo_composite_effect_heatmap_updated',cbar_label='effect')
    plot_scatter(ec_cnmf,'P34_YAP_TAZ_mesenchyme','D3 P34 YAP/TAZ-mesenchyme effect\napprox_log2FC','05_D3_EC_loss_vs_P34_YAP_TAZ_effect_scatter','P34 induction tracks D3 EC loss',outdir)
    plot_scatter(ec_cnmf,'P34_minus_P2','D3 P34 − P2 program effect','06_D3_EC_loss_vs_P34_minus_P2_scatter','P34-over-P2 effect tracks D3 EC loss',outdir)
    heatmap_matrix(effect_supp,'Supplemental D3 program effects including excluded/EC programs',outdir,'S1_D3_cNMF_Hippo_program_effect_heatmap_supplement_all_context',fdr=fdr_supp,col_labels=[PROGRAM_LABELS_SUPP.get(c,f'P{c}') for c in effect_supp.columns])
    with open(outdir/'README.txt','w') as f: f.write(f"Updated D3 EC% + cNMF Hippo panels based on program annotation notes.\n\nMain program set: P34, P27, P24, P2.\nP34 = primary YAP/TAZ-mesenchyme / Hippo program.\nP27/P24 = secondary alternative-context programs.\nP2 = arterial endothelial reference program.\nP35 excluded from main panel; included only in supplement context heatmap.\n\nNT pooled D3 EC%: {nt_pct:.3f}\nNT pooled n_D3_cells: {nt_n}\nNT pooled n_D3_EC_cells: {nt_ec}\n")
    print('Done',outdir); print(summary[['target','n_D3_cells','pct_EC_D3_pooled','delta_EC_pct_vs_NT','EC_loss_pp_vs_NT']].to_string(index=False)); print(effect_main.round(3).to_string()); print(comp.round(3).to_string())

if __name__=='__main__': main()
