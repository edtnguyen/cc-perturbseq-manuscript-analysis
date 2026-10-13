#!/usr/bin/env python3
"""Control/reference-defined coarse developmental-state atlas for CC-Perturb-seq.

No target-to-program effects or perturbation-driven clustering are used. Feature selection,
latent fitting, UMAP fitting, clustering and classifier training are restricted to either an
independent unperturbed reference or pure safe/non-targeting controls.
"""
from __future__ import annotations

import json, math, os, platform, shutil, sys, textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import joblib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from matplotlib.colors import ListedColormap
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler

try:
    import anndata as ad
    import scanpy as sc
except ImportError as e:  # pragma: no cover
    raise SystemExit("Install scanpy/anndata from requirements.txt") from e
try:
    import umap
except ImportError as e:  # pragma: no cover
    raise SystemExit("Install umap-learn from requirements.txt") from e
try:
    from adjustText import adjust_text
except Exception:  # pragma: no cover
    adjust_text = None

DAY_ORDER=["D0","D1","D2","D3"]
DEFAULT_STATE_ORDER=[
    "hiPSC","Primitive streak","Posterior / alternative primitive streak",
    "Lateral mesoderm","Vascular / pre-EC progenitor","Endothelial cell",
    "Contractile / mesenchymal alternative","Other / unresolved",
]
DEFAULT_STATE_COLORS={
    "hiPSC":"#6E6E6E","Primitive streak":"#4C9A5B",
    "Posterior / alternative primitive streak":"#7EAA65",
    "Lateral mesoderm":"#8A63B8","Vascular / pre-EC progenitor":"#D78B35",
    "Endothelial cell":"#C94949","Contractile / mesenchymal alternative":"#3F78B5",
    "Other / unresolved":"#BDBDBD",
}
DEFAULT_DAY_COLORS={"D0":"#6E6E6E","D1":"#4C9A5B","D2":"#D78B35","D3":"#C94949"}
EXPECTED_DAY_BY_STATE={
    "hiPSC":{"D0"},"Primitive streak":{"D1"},
    "Posterior / alternative primitive streak":{"D1"},
    "Lateral mesoderm":{"D2"},"Vascular / pre-EC progenitor":{"D2"},
    "Endothelial cell":{"D3"},"Contractile / mesenchymal alternative":{"D3"},
    "Other / unresolved":set(DAY_ORDER),
}

@dataclass
class AtlasModel:
    source: str
    genes: List[str]
    scaler: StandardScaler
    svd: TruncatedSVD
    umap_model: Any
    classifier: KNeighborsClassifier
    state_classes: List[str]
    training_latent: np.ndarray
    training_embedding: np.ndarray
    training_labels: np.ndarray
    training_clusters: np.ndarray
    training_cell_ids: np.ndarray
    training_days: np.ndarray
    training_replicates: np.ndarray


def normalize_day(x:Any)->str:
    s=str(x).strip().upper().replace("DAY", "D").replace(" ","")
    return {"0":"D0","0.0":"D0","D0":"D0","1":"D1","1.0":"D1","D1":"D1",
            "2":"D2","2.0":"D2","D2":"D2","3":"D3","3.0":"D3","D3":"D3"}.get(s,s)

def resolve_path(base:Path, value:Optional[str])->Optional[Path]:
    if value in (None,"","null"): return None
    p=Path(value).expanduser()
    return p if p.is_absolute() else (base/p).resolve()

def load_config(path:str|Path)->Tuple[Dict[str,Any],Path]:
    p=Path(path).expanduser().resolve()
    cfg=yaml.safe_load(p.read_text())
    base=p.parent.parent if p.parent.name=="config" else p.parent
    for key in ["query_h5ad","reference_h5ad","guide_metadata_tsv","approved_cluster_state_map_tsv",
                "approved_state_color_map_tsv","approved_day_color_map_tsv","marker_sets_tsv"]:
        if key in cfg: cfg[key]=str(resolve_path(base,cfg.get(key))) if resolve_path(base,cfg.get(key)) else None
    out=resolve_path(base,cfg.get("output_dir","control_defined_state_atlas_results"))
    cfg["output_dir"]=str(out)
    return cfg,base

def ensure_csr(x:Any)->sp.csr_matrix:
    if sp.issparse(x): return x.tocsr().astype(np.float32)
    return sp.csr_matrix(np.asarray(x,dtype=np.float32))

def get_matrix(adata:ad.AnnData, layer:Optional[str])->Any:
    if layer in (None,"","X"): return adata.X
    if layer=="raw":
        if adata.raw is None: raise ValueError("Requested raw but adata.raw is absent")
        return adata.raw.X
    if layer not in adata.layers: raise KeyError(f"Layer {layer!r} not found; available={list(adata.layers.keys())}")
    return adata.layers[layer]

def gene_symbols(adata:ad.AnnData, col:Optional[str])->pd.Index:
    vals=adata.var[col].astype(str) if col and col in adata.var.columns else adata.var_names.astype(str)
    return pd.Index([str(x).strip() for x in vals],dtype="object")

def first_symbol_index(symbols:pd.Index)->Dict[str,int]:
    out={}
    for i,g in enumerate(symbols):
        if g and g.lower()!="nan" and g not in out: out[g]=i
    return out

def normalize_log1p(x:Any, raw_counts:bool, target_sum:float=1e4)->sp.csr_matrix:
    x=ensure_csr(x).copy()
    if not raw_counts: return x
    totals=np.asarray(x.sum(axis=1)).ravel()
    totals=np.where(totals>0,totals,1.0)
    x=sp.diags(target_sum/totals).dot(x).tocsr()
    x.data=np.log1p(x.data)
    return x

def detection_fraction(x:sp.csr_matrix)->np.ndarray:
    return np.asarray((x>0).mean(axis=0)).ravel()

def sparse_variance(x:sp.csr_matrix)->np.ndarray:
    m=np.asarray(x.mean(axis=0)).ravel(); m2=np.asarray(x.multiply(x).mean(axis=0)).ravel()
    return np.maximum(m2-m*m,0)

def audit_adata(obj:ad.AnnData,name:str,path:Path)->Dict[str,Any]:
    return {"dataset":name,"path":str(path),"n_obs":int(obj.n_obs),"n_vars":int(obj.n_vars),
            "obs_columns":json.dumps(list(obj.obs.columns)),"var_columns":json.dumps(list(obj.var.columns)),
            "obsm_keys":json.dumps(list(obj.obsm.keys())),"uns_keys":json.dumps(list(obj.uns.keys())),
            "layers":json.dumps(list(obj.layers.keys())),"obsp_keys":json.dumps(list(obj.obsp.keys())),
            "raw_present":bool(obj.raw is not None),"file_size_bytes":path.stat().st_size}

def read_color_map(path:Optional[str], key_col:str)->Dict[str,str]:
    if not path: return {}
    p=Path(path); sep="\t" if p.suffix.lower() in {".tsv",".txt"} else ","
    d=pd.read_csv(p,sep=sep)
    c="color" if "color" in d.columns else "hex"
    return dict(zip(d[key_col].astype(str),d[c].astype(str)))

def adata_palette(adata:ad.AnnData,col:str)->Dict[str,str]:
    key=f"{col}_colors"
    if col not in adata.obs or key not in adata.uns: return {}
    s=adata.obs[col]
    if not isinstance(s.dtype,pd.CategoricalDtype): return {}
    cats=list(map(str,s.cat.categories)); colors=list(map(str,adata.uns[key]))
    return dict(zip(cats,colors)) if len(cats)==len(colors) else {}

def select_pure_controls(query:ad.AnnData,cfg:Mapping[str,Any],outdir:Path)->Tuple[np.ndarray,str,pd.DataFrame]:
    key=cfg.get("query_guide_assignment_obsm","guide_assignment")
    names_key=cfg.get("query_guide_names_uns","guide_names")
    if key not in query.obsm: raise KeyError(f"{key} absent from query.obsm; available={list(query.obsm.keys())}")
    if names_key not in query.uns: raise KeyError(f"{names_key} absent from query.uns; available={list(query.uns.keys())}")
    assign=ensure_csr(query.obsm[key]); guides=list(map(str,query.uns[names_key]))
    if assign.shape[1]!=len(guides): raise ValueError("guide_assignment columns do not match guide_names length")
    meta=pd.read_csv(cfg["guide_metadata_tsv"],sep="\t")
    if "guide_id" not in meta.columns: raise ValueError("Guide metadata requires guide_id")
    meta=meta.drop_duplicates("guide_id").set_index(meta["guide_id"].astype(str),drop=False)
    missing=[g for g in guides if g not in meta.index]
    if missing:
        pd.Series(missing,name="guide_id").to_csv(outdir/"guides_missing_from_metadata.tsv",sep="\t",index=False)
        raise ValueError(f"{len(missing)} query guides missing from metadata; see guides_missing_from_metadata.tsv")
    m=meta.loc[guides]
    if "type" in m.columns:
        types=m["type"].astype(str).str.lower().to_numpy()
        control=np.isin(types,["safe-targeting","non-targeting","ntc"])
        real=~control
    else:
        real=m["targeting"].fillna(False).astype(bool).to_numpy(); control=~real
    c_count=np.asarray(assign[:,control].sum(axis=1)).ravel()
    r_count=np.asarray(assign[:,real].sum(axis=1)).ravel()
    mask=(c_count>0)&(r_count==0)
    qc=pd.DataFrame({"cell_id":query.obs_names.astype(str),"n_control_guides":c_count,
                     "n_real_targeting_guides":r_count,"is_pure_control":mask,
                     "mixed_control_and_real":(c_count>0)&(r_count>0)})
    definition=("At least one safe-targeting or non-targeting guide and zero real-targeting guides. "
                "Cells carrying control and real-targeting guides together are excluded from controls.")
    return mask,definition,qc

def load_markers(path:str)->pd.DataFrame:
    p=Path(path); sep="\t" if p.suffix.lower() in {".tsv",".txt"} else ","
    d=pd.read_csv(p,sep=sep)
    need={"state","gene"}
    if not need.issubset(d.columns): raise ValueError(f"Marker file needs {need}")
    return d

def select_features(train_x:sp.csr_matrix,train_symbols:pd.Index,marker_genes:Sequence[str],cfg:Mapping[str,Any])->List[str]:
    idx=first_symbol_index(train_symbols)
    valid=np.array([i for g,i in idx.items()],dtype=int)
    unique_genes=np.array(list(idx.keys()),dtype=object)
    x=train_x[:,valid]
    det=detection_fraction(x); var=sparse_variance(x)
    keep=det>=float(cfg.get("minimum_gene_detection_fraction",0.01))
    order=np.argsort(var[keep])[::-1]
    candidates=unique_genes[keep][order][:int(cfg.get("max_genes",3500))].tolist()
    for g in marker_genes:
        if g in idx and g not in candidates: candidates.append(g)
    return candidates

def aligned_matrix(adata:ad.AnnData, layer:Optional[str], raw_counts:bool, symbols_col:Optional[str], genes:Sequence[str])->sp.csr_matrix:
    syms=gene_symbols(adata,symbols_col); idx=first_symbol_index(syms)
    missing=[g for g in genes if g not in idx]
    if missing: raise ValueError(f"{len(missing)} selected genes missing from dataset; examples={missing[:10]}")
    cols=[idx[g] for g in genes]
    return normalize_log1p(get_matrix(adata,layer),raw_counts)[:,cols]

def fit_latent(train_x:sp.csr_matrix, all_x:sp.csr_matrix, cfg:Mapping[str,Any]):
    scaler=StandardScaler(with_mean=False,copy=True)
    train_scaled=scaler.fit_transform(train_x); all_scaled=scaler.transform(all_x)
    ncomp=min(int(cfg.get("n_components",50)),train_scaled.shape[1]-1,train_scaled.shape[0]-1)
    svd=TruncatedSVD(n_components=max(2,ncomp),random_state=int(cfg.get("random_seed",20260727)))
    train_latent=svd.fit_transform(train_scaled); all_latent=svd.transform(all_scaled)
    reducer=umap.UMAP(n_neighbors=int(cfg.get("n_neighbors",30)),min_dist=float(cfg.get("umap_min_dist",0.3)),
                      metric="cosine",random_state=int(cfg.get("random_seed",20260727)),transform_seed=int(cfg.get("random_seed",20260727)))
    train_emb=reducer.fit_transform(train_latent); all_emb=reducer.transform(all_latent)
    return scaler,svd,reducer,train_latent,all_latent,train_emb,all_emb

def cluster_training(latent:np.ndarray,cfg:Mapping[str,Any])->np.ndarray:
    a=ad.AnnData(X=np.asarray(latent,dtype=np.float32))
    sc.pp.neighbors(a,n_neighbors=int(cfg.get("n_neighbors",30)),use_rep="X",random_state=int(cfg.get("random_seed",20260727)))
    sc.tl.leiden(a,resolution=float(cfg.get("leiden_resolution",0.65)),key_added="source_cluster",
                 random_state=int(cfg.get("random_seed",20260727)))
    return a.obs["source_cluster"].astype(str).to_numpy()

def marker_expression(adata:ad.AnnData,rows:np.ndarray,cfg:Mapping[str,Any],markers:pd.DataFrame,is_reference:bool=False):
    prefix="reference" if is_reference else "query"
    layer=cfg.get(f"{prefix}_layer")
    raw=bool(cfg.get(f"{prefix}_matrix_is_raw_counts",True))
    symcol=cfg.get(f"{prefix}_gene_symbol_col")
    all_genes=markers["gene"].drop_duplicates().astype(str).tolist()
    syms=gene_symbols(adata,symcol); idx=first_symbol_index(syms)
    avail=[g for g in all_genes if g in idx]
    cols=[idx[g] for g in avail]
    x=normalize_log1p(get_matrix(adata,layer),raw)[rows][:,cols]
    return x,avail

def auto_annotate_clusters(adata:ad.AnnData,rows:np.ndarray,clusters:np.ndarray,days:np.ndarray,reps:np.ndarray,
                           markers:pd.DataFrame,cfg:Mapping[str,Any],is_reference:bool=False):
    x,genes=marker_expression(adata,rows,cfg,markers,is_reference=is_reference)
    mean=np.asarray(x.mean(axis=0)).ravel(); mean2=np.asarray(x.multiply(x).mean(axis=0)).ravel()
    sd=np.sqrt(np.maximum(mean2-mean*mean,1e-6))
    gidx={g:i for i,g in enumerate(genes)}
    rows_out=[]; score_rows=[]
    for cl in sorted(pd.unique(clusters),key=lambda s:(len(str(s)),str(s))):
        mask=clusters==cl; n=int(mask.sum())
        dayvc=pd.Series(days[mask]).value_counts(normalize=True)
        repvc=pd.Series(reps[mask]).value_counts(normalize=True)
        dom_day=str(dayvc.index[0]) if len(dayvc) else "NA"
        cluster_mean=np.asarray(x[mask].mean(axis=0)).ravel()
        cluster_z=(cluster_mean-mean)/sd
        state_scores={}; state_counts={}; state_present={}
        for state,mdf in markers.groupby("state"):
            avail=[g for g in mdf["gene"].astype(str) if g in gidx]
            state_present[state]=avail
            state_counts[state]=len(avail)
            if avail:
                vals=np.asarray([cluster_z[gidx[g]] for g in avail],dtype=float)
                state_scores[state]=float(np.nanmean(vals))
            else: state_scores[state]=np.nan
        valid=[(s,v) for s,v in state_scores.items() if np.isfinite(v)]
        valid=sorted(valid,key=lambda kv:kv[1],reverse=True)
        top_state,top_score=valid[0] if valid else ("Other / unresolved",np.nan)
        second_score=valid[1][1] if len(valid)>1 else -np.inf
        margin=top_score-second_score if np.isfinite(top_score) else np.nan
        enough=state_counts.get(top_state,0)>=int(cfg.get("minimum_markers_present_for_named_state",2))
        score_ok=np.isfinite(top_score) and top_score>=float(cfg.get("minimum_marker_score",0.15))
        margin_ok=np.isfinite(margin) and margin>=float(cfg.get("minimum_score_margin",0.10))
        day_ok=dom_day in EXPECTED_DAY_BY_STATE.get(top_state,set(DAY_ORDER))
        provisional=top_state if enough and score_ok and margin_ok and day_ok else "Other / unresolved"
        if provisional=="Other / unresolved": confidence="low"
        elif margin>=0.35 and top_score>=0.50: confidence="high"
        else: confidence="medium"
        rows_out.append({"source_cluster":cl,"n_cells":n,"dominant_day":dom_day,
                         "dominant_day_fraction":float(dayvc.iloc[0]) if len(dayvc) else np.nan,
                         "dominant_replicate":str(repvc.index[0]) if len(repvc) else "NA",
                         "day_composition_json":json.dumps({str(k):float(v) for k,v in dayvc.items()}),
                         "replicate_composition_json":json.dumps({str(k):float(v) for k,v in repvc.items()}),
                         "provisional_state":provisional,"annotation_confidence":confidence,
                         "top_marker_state":top_state,"top_marker_score":top_score,
                         "second_marker_score":second_score,"score_margin":margin,
                         "markers_present_for_top_state":state_counts.get(top_state,0),
                         "top_state_markers_present":";".join(state_present.get(top_state,[])),
                         "annotation_source":"marker scores plus observed day; source cells only"})
        for state,v in state_scores.items():
            score_rows.append({"source_cluster":cl,"state":state,"marker_score":v,
                               "n_markers_available":state_counts[state],"markers_available":";".join(state_present[state])})
    return pd.DataFrame(rows_out),pd.DataFrame(score_rows),x,genes

def load_approved_map(path:Optional[str])->Optional[pd.DataFrame]:
    if not path: return None
    p=Path(path); sep="\t" if p.suffix.lower() in {".tsv",".txt"} else ","
    d=pd.read_csv(p,sep=sep,comment="#")
    if not {"source_cluster","approved_state"}.issubset(d.columns): raise ValueError("Approved map requires source_cluster and approved_state")
    return d

def resolve_cluster_labels(review:pd.DataFrame,cfg:Mapping[str,Any])->Tuple[Dict[str,str],str]:
    approved=load_approved_map(cfg.get("approved_cluster_state_map_tsv"))
    if approved is not None and len(approved):
        mapping=dict(zip(approved["source_cluster"].astype(str),approved["approved_state"].astype(str)))
        missing=set(review["source_cluster"].astype(str))-set(mapping)
        if missing: raise ValueError(f"Approved cluster map lacks clusters: {sorted(missing)}")
        return mapping,"author-approved cluster-to-state map"
    if "annotation_confidence" in review.columns and len(review) and review["annotation_confidence"].astype(str).eq("approved").all():
        return dict(zip(review["source_cluster"].astype(str),review["source_cluster"].astype(str))),"author-approved independent reference states"
    if bool(cfg.get("allow_provisional_state_labels_in_final_panel",False)):
        return dict(zip(review["source_cluster"].astype(str),review["provisional_state"].astype(str))),"provisional automatic annotation"
    raise RuntimeError("Named final state labels require an approved_cluster_state_map_tsv. Review control_state_annotation_review.tsv first.")

def top_markers_from_selected(x:sp.csr_matrix,genes:Sequence[str],clusters:np.ndarray,labels:Mapping[str,str],top_n:int=25)->pd.DataFrame:
    rows=[]; x=ensure_csr(x)
    for cl in sorted(pd.unique(clusters),key=str):
        m=clusters==cl; n=int(m.sum()); rest=~m
        cm=np.asarray(x[m].mean(axis=0)).ravel(); rm=np.asarray(x[rest].mean(axis=0)).ravel()
        score=cm-rm; order=np.argsort(score)[::-1][:top_n]
        frac=np.asarray((x[m]>0).mean(axis=0)).ravel()
        for rank,j in enumerate(order,1):
            rows.append({"source_cluster":cl,"state":labels.get(str(cl),"Other / unresolved"),"rank":rank,
                         "gene":genes[j],"mean_log_expression":cm[j],"mean_difference_vs_other_clusters":score[j],
                         "fraction_expressing":frac[j],"n_source_cells":n})
    return pd.DataFrame(rows)

def state_counts(assign:pd.DataFrame)->pd.DataFrame:
    group_cols=["population","observed_day","biological_replicate","assigned_state"]
    d=assign.groupby(group_cols,dropna=False).size().reset_index(name="n_cells")
    denom=d.groupby(["population","observed_day","biological_replicate"],dropna=False)["n_cells"].transform("sum")
    d["fraction_within_day_replicate"]=d["n_cells"]/denom
    return d

def sample_idx(n:int,max_n:int,seed:int)->np.ndarray:
    if n<=max_n: return np.arange(n)
    rng=np.random.default_rng(seed); return np.sort(rng.choice(n,size=max_n,replace=False))

def plot_atlas(embedding:np.ndarray,states:np.ndarray,days:np.ndarray,state_colors:Mapping[str,str],day_colors:Mapping[str,str],
               outbase:Path,title:str,note:str,composition:Optional[pd.DataFrame],cfg:Mapping[str,Any]):
    maxp=int(cfg.get("max_points_display",120000)); idx=sample_idx(len(states),maxp,int(cfg.get("random_seed",20260727)))
    emb=embedding[idx]; st=states[idx]
    if composition is None:
        fig,ax=plt.subplots(figsize=(float(cfg.get("figure_width_inches",8.2)),float(cfg.get("figure_height_inches",5.8))))
        axes=[ax]
    else:
        fig=plt.figure(figsize=(float(cfg.get("figure_width_inches",8.2)),float(cfg.get("figure_height_inches",5.8))))
        gs=fig.add_gridspec(1,2,width_ratios=[4.3,1.5],wspace=0.18)
        ax=fig.add_subplot(gs[0,0]); axc=fig.add_subplot(gs[0,1]); axes=[ax,axc]
    order=[s for s in DEFAULT_STATE_ORDER if s in set(states)] + [s for s in pd.unique(states) if s not in DEFAULT_STATE_ORDER]
    for s in order:
        m=st==s
        if not np.any(m): continue
        ax.scatter(emb[m,0],emb[m,1],s=float(cfg.get("point_size_all_cells",2.0)),c=state_colors.get(s,"#BDBDBD"),
                   alpha=float(cfg.get("point_alpha",0.7)),linewidths=0,rasterized=bool(cfg.get("rasterize_points",True)),label=s)
    texts=[]
    for s in order:
        m=states==s
        if not np.any(m): continue
        x,y=np.median(embedding[m,0]),np.median(embedding[m,1])
        texts.append(ax.text(x,y,s,fontsize=8.7,fontweight="bold",ha="center",va="center",
                             bbox=dict(boxstyle="round,pad=0.22",fc="white",ec=state_colors.get(s,"#999999"),lw=0.7,alpha=0.90)))
    if adjust_text is not None and texts:
        adjust_text(texts,ax=ax,arrowprops=dict(arrowstyle="-",color="#777777",lw=0.5),time_lim=1.5)
    ax.set_title(title,loc="left",fontsize=13,fontweight="bold",pad=10)
    ax.text(0.0,1.01,note,transform=ax.transAxes,ha="left",va="bottom",fontsize=8.5,color="#555555")
    ax.set_xlabel("atlas embedding 1"); ax.set_ylabel("atlas embedding 2")
    ax.set_xticks([]); ax.set_yticks([])
    for spn in ax.spines.values(): spn.set_visible(False)
    if composition is not None:
        comp=composition.copy(); comp["observed_day"]=pd.Categorical(comp["observed_day"],DAY_ORDER,ordered=True)
        comp=comp.groupby(["observed_day","assigned_state"],observed=True)["n_cells"].sum().reset_index()
        comp["fraction"]=comp.groupby("observed_day",observed=True)["n_cells"].transform(lambda x:x/x.sum())
        left=np.zeros(len(DAY_ORDER))
        y=np.arange(len(DAY_ORDER))
        for s in order:
            vals=[]
            for d in DAY_ORDER:
                q=comp[(comp["observed_day"]==d)&(comp["assigned_state"]==s)]
                vals.append(float(q["fraction"].iloc[0]) if len(q) else 0.0)
            axc.barh(y,vals,left=left,color=state_colors.get(s,"#BDBDBD"),height=0.68,label=s)
            left+=np.asarray(vals)
        axc.set_yticks(y,DAY_ORDER); axc.invert_yaxis(); axc.set_xlim(0,1)
        axc.set_xlabel("fraction of cells"); axc.set_title("State composition\nby observed day",fontsize=9.5,fontweight="bold")
        axc.spines[["top","right","left"]].set_visible(False); axc.tick_params(axis="y",length=0)
    handles=[mpl.lines.Line2D([0],[0],marker="o",color="none",markerfacecolor=state_colors.get(s,"#BDBDBD"),markersize=6,label=s) for s in order]
    fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(0.5,-0.01),ncol=min(4,len(handles)),frameon=False,fontsize=7.3)
    fig.tight_layout(rect=[0,0.07,1,1])
    for ext in ["png","pdf","svg"]:
        fig.savefig(outbase.with_suffix("."+ext),dpi=350 if ext=="png" else None,bbox_inches="tight")
    plt.close(fig)

def build_source_training(query:ad.AnnData,reference:Optional[ad.AnnData],control_mask:np.ndarray,cfg:Mapping[str,Any],markers:pd.DataFrame,outdir:Path):
    source_pref=str(cfg.get("state_definition_source","auto"))
    use_ref=reference is not None and source_pref in {"auto","reference_approved","reference_transfer"}
    if use_ref:
        source=reference; source_name="independent_reference"
        day_col=cfg["reference_day_col"]; rep_col=cfg["reference_replicate_col"]
        rows=np.arange(reference.n_obs); is_ref=True
        syms=gene_symbols(reference,cfg.get("reference_gene_symbol_col")); x=normalize_log1p(get_matrix(reference,cfg.get("reference_layer")),bool(cfg.get("reference_matrix_is_raw_counts",True)))
    else:
        source=query; source_name="pure_controls"
        day_col=cfg["query_day_col"]; rep_col=cfg["query_replicate_col"]
        rows=np.where(control_mask)[0]; is_ref=False
        syms=gene_symbols(query,cfg.get("query_gene_symbol_col")); x=normalize_log1p(get_matrix(query,cfg.get("query_layer")),bool(cfg.get("query_matrix_is_raw_counts",True)))[rows]
    days=source.obs.iloc[rows][day_col].map(normalize_day).to_numpy(); reps=source.obs.iloc[rows][rep_col].astype(str).to_numpy()
    marker_genes=markers["gene"].astype(str).drop_duplicates().tolist()
    genes=select_features(x,syms,marker_genes,cfg)
    query_gene_set=set(gene_symbols(query,cfg.get("query_gene_symbol_col")))
    genes=[g for g in genes if g in query_gene_set]
    if len(genes)<100:
        raise RuntimeError(f"Only {len(genes)} shared source/query genes remain after filtering")
    tx=aligned_matrix(source,cfg.get("reference_layer" if is_ref else "query_layer"),bool(cfg.get("reference_matrix_is_raw_counts" if is_ref else "query_matrix_is_raw_counts",True)),cfg.get("reference_gene_symbol_col" if is_ref else "query_gene_symbol_col"),genes)[rows]
    qx=aligned_matrix(query,cfg.get("query_layer"),bool(cfg.get("query_matrix_is_raw_counts",True)),cfg.get("query_gene_symbol_col"),genes)
    scaler,svd,reducer,train_latent,all_latent,train_emb,all_emb=fit_latent(tx,qx,cfg)
    if use_ref and bool(cfg.get("reference_states_approved",False)) and cfg.get("reference_state_col") in reference.obs.columns:
        labels=reference.obs.iloc[rows][cfg["reference_state_col"]].astype(str).to_numpy(); clusters=labels.copy()
        review=pd.DataFrame({"source_cluster":pd.unique(clusters)})
        review["provisional_state"]=review["source_cluster"]; review["annotation_confidence"]="approved"
        review["annotation_source"]="author-approved independent reference states"
        score_table=pd.DataFrame()
    else:
        clusters=cluster_training(train_latent,cfg)
        review,score_table,_,_=auto_annotate_clusters(source,rows,clusters,days,reps,markers,cfg,is_reference=is_ref)
        review.to_csv(outdir/"control_state_annotation_review.tsv",sep="\t",index=False)
        score_table.to_csv(outdir/"control_cluster_marker_scores.tsv",sep="\t",index=False)
        labels=np.array([dict(zip(review["source_cluster"].astype(str),review["provisional_state"]))[str(c)] for c in clusters])
    return source_name,source,rows,days,reps,genes,tx,qx,scaler,svd,reducer,train_latent,all_latent,train_emb,all_emb,clusters,labels,review,score_table,is_ref

def run_review(config_path:str|Path)->Path:
    cfg,_=load_config(config_path); outdir=Path(cfg["output_dir"]); outdir.mkdir(parents=True,exist_ok=True)
    qpath=Path(cfg["query_h5ad"])
    if not qpath.exists(): raise FileNotFoundError(f"Query H5AD not found: {qpath}")
    query=sc.read_h5ad(qpath)
    rpath=Path(cfg["reference_h5ad"]) if cfg.get("reference_h5ad") else None
    reference=sc.read_h5ad(rpath) if rpath and rpath.exists() else None
    audits=[audit_adata(query,"query",qpath)]
    if reference is not None: audits.append(audit_adata(reference,"reference",rpath))
    pd.DataFrame(audits).to_csv(outdir/"data_audit.tsv",sep="\t",index=False)
    for c in [cfg["query_day_col"],cfg["query_replicate_col"]]:
        if c not in query.obs: raise KeyError(f"Query obs lacks {c}")
    control_mask,control_def,control_qc=select_pure_controls(query,cfg,outdir)
    control_qc.to_csv(outdir/"control_cell_definition_qc.tsv",sep="\t",index=False)
    control_counts=(query.obs.assign(is_pure_control=control_mask,observed_day=query.obs[cfg["query_day_col"]].map(normalize_day))
                    .groupby(["observed_day",cfg["query_replicate_col"],"is_pure_control"],dropna=False).size().reset_index(name="n_cells"))
    control_counts.to_csv(outdir/"control_counts_by_day_and_replicate.tsv",sep="\t",index=False)
    if int(control_mask.sum())<int(cfg.get("minimum_controls_total",500)) and reference is None:
        raise RuntimeError(f"Only {control_mask.sum()} pure controls and no reference. Control definition was not broadened.")
    byday=pd.Series(query.obs.loc[control_mask,cfg["query_day_col"]].map(normalize_day)).value_counts()
    low={d:int(byday.get(d,0)) for d in DAY_ORDER if int(byday.get(d,0))<int(cfg.get("minimum_controls_per_day",75))}
    if low and reference is None: raise RuntimeError(f"Pure controls below minimum by day: {low}. Control definition was not broadened.")
    markers=load_markers(cfg["marker_sets_tsv"])
    bundle=build_source_training(query,reference,control_mask,cfg,markers,outdir)
    (source_name,source,rows,days,reps,genes,tx,qx,scaler,svd,reducer,train_latent,all_latent,train_emb,all_emb,clusters,labels,review,score_table,is_ref)=bundle
    pd.Series(genes,name="gene").to_csv(outdir/"atlas_feature_genes.tsv",sep="\t",index=False)
    if "control_state_annotation_review.tsv" not in [p.name for p in outdir.iterdir()]:
        review.to_csv(outdir/"control_state_annotation_review.tsv",sep="\t",index=False)
    # provisional candidate A for review only
    state_colors=dict(DEFAULT_STATE_COLORS); state_colors.update(read_color_map(cfg.get("approved_state_color_map_tsv"),"state"))
    plot_atlas(train_emb,labels,days,state_colors,DEFAULT_DAY_COLORS,outdir/"Fig2_control_defined_state_atlas_candidateA_review",
               "Control/reference-defined developmental states",
               f"State boundaries learned from {source_name}; n={len(labels):,} source cells. Labels are provisional until reviewed.",None,cfg)
    joblib.dump({"source_name":source_name,"genes":genes,"scaler":scaler,"svd":svd,"umap_model":reducer,
                 "train_latent":train_latent,"all_latent":all_latent,"train_embedding":train_emb,"all_embedding":all_emb,
                 "clusters":clusters,"provisional_labels":labels,"training_rows":rows,"training_days":days,"training_replicates":reps,
                 "control_mask":control_mask,"control_definition":control_def,"is_reference":is_ref},outdir/"review_model_bundle.joblib",compress=3)
    (outdir/"REVIEW_NEXT_STEPS.md").write_text(
        "# Review required\n\nInspect `control_state_annotation_review.tsv`, `control_cluster_marker_scores.tsv`, "
        "and `Fig2_control_defined_state_atlas_candidateA_review.pdf`. Copy the template "
        "`approved_cluster_state_map.tsv` into the output directory, fill one row per source cluster, "
        "set its path in the YAML config, then run finalize mode.\n")
    return outdir

def run_finalize(config_path:str|Path)->Path:
    cfg,_=load_config(config_path); outdir=Path(cfg["output_dir"]); outdir.mkdir(parents=True,exist_ok=True)
    qpath=Path(cfg["query_h5ad"]); query=sc.read_h5ad(qpath)
    rpath=Path(cfg["reference_h5ad"]) if cfg.get("reference_h5ad") else None
    reference=sc.read_h5ad(rpath) if rpath and rpath.exists() else None
    control_mask,control_def,control_qc=select_pure_controls(query,cfg,outdir)
    markers=load_markers(cfg["marker_sets_tsv"])
    bundle=build_source_training(query,reference,control_mask,cfg,markers,outdir)
    (source_name,source,rows,days,reps,genes,tx,qx,scaler,svd,reducer,train_latent,all_latent,train_emb,all_emb,clusters,provisional,review,score_table,is_ref)=bundle
    mapping,mapping_source=resolve_cluster_labels(review,cfg)
    train_states=np.array([mapping[str(c)] for c in clusters])
    clf=KNeighborsClassifier(n_neighbors=min(int(cfg.get("knn_classifier_neighbors",35)),len(train_states)-1),weights="distance",metric="euclidean")
    clf.fit(train_latent,train_states); prob=clf.predict_proba(all_latent); assigned=clf.classes_[np.argmax(prob,axis=1)]
    conf=prob.max(axis=1); ent=-np.sum(np.clip(prob,1e-12,1)*np.log(np.clip(prob,1e-12,1)),axis=1)/np.log(prob.shape[1]) if prob.shape[1]>1 else np.zeros(len(prob))
    obs_day=query.obs[cfg["query_day_col"]].map(normalize_day).to_numpy(); obs_rep=query.obs[cfg["query_replicate_col"]].astype(str).to_numpy()
    assignments=pd.DataFrame({"cell_id":query.obs_names.astype(str),"observed_day":obs_day,"biological_replicate":obs_rep,
                              "is_pure_control":control_mask,"assigned_state":assigned,"assignment_confidence":conf,
                              "assignment_entropy":ent,"atlas_embedding_1":all_emb[:,0],"atlas_embedding_2":all_emb[:,1],
                              "state_definition_source":source_name,"state_mapping_source":mapping_source})
    for j,s in enumerate(clf.classes_): assignments[f"probability__{s}"]=prob[:,j]
    assignments.to_csv(outdir/"control_state_assignments.tsv",sep="\t",index=False)
    # top markers from selected feature matrix, source cells only
    marker_summary=top_markers_from_selected(tx,genes,clusters,mapping,top_n=25)
    marker_summary.to_csv(outdir/"state_marker_summary.tsv",sep="\t",index=False)
    # annotation review final, including marker and developmental correspondence summaries
    top_marker_strings=(marker_summary.sort_values(["source_cluster","rank"])
                        .groupby("source_cluster")["gene"].apply(lambda x:";".join(x.head(10))).to_dict())
    review=review.copy()
    review["approved_state"]=review["source_cluster"].astype(str).map(mapping)
    review["state_mapping_source"]=mapping_source
    review["top_marker_genes"]=review["source_cluster"].astype(str).map(top_marker_strings)
    review["reference_stage_correspondence"]=review.get("dominant_day",pd.Series(["NA"]*len(review))).astype(str)
    review.to_csv(outdir/"control_state_annotation_review.tsv",sep="\t",index=False)
    count_rows=[]
    for population,mask in [("all_CC-Perturb-seq_cells",np.ones(query.n_obs,dtype=bool)),("pure_controls",control_mask)]:
        d=assignments.loc[mask,["observed_day","biological_replicate","assigned_state"]].copy(); d["population"]=population; count_rows.append(d)
    counts=state_counts(pd.concat(count_rows,ignore_index=True)); counts.to_csv(outdir/"state_counts_by_day_and_replicate.tsv",sep="\t",index=False)
    state_colors=dict(DEFAULT_STATE_COLORS); state_colors.update(read_color_map(cfg.get("approved_state_color_map_tsv"),"state"))
    used_states=[s for s in DEFAULT_STATE_ORDER if s in set(assigned)] + [s for s in pd.unique(assigned) if s not in DEFAULT_STATE_ORDER]
    pd.DataFrame({"state":used_states,"color":[state_colors.get(s,"#BDBDBD") for s in used_states],
                  "source":["approved/fallback config" for _ in used_states]}).to_csv(outdir/"state_color_map.tsv",sep="\t",index=False)
    day_colors=dict(DEFAULT_DAY_COLORS); day_colors.update(read_color_map(cfg.get("approved_day_color_map_tsv"),"day"))
    # candidate A source/reference only
    plot_atlas(train_emb,train_states,days,state_colors,day_colors,outdir/"Fig2_control_defined_state_atlas_candidateA",
               "Control/reference-defined developmental states",
               f"State boundaries learned from {source_name}; n={len(train_states):,} source cells.",None,cfg)
    # candidate B all cells with day composition strip
    plot_atlas(all_emb,assigned,obs_day,state_colors,day_colors,outdir/"Fig2_control_defined_state_atlas_candidateB",
               "Control-defined states resolve the endothelial differentiation landscape",
               f"All CC-Perturb-seq cells projected into a fixed {source_name}-defined atlas; n={query.n_obs:,} cells.",counts[counts["population"]=="all_CC-Perturb-seq_cells"],cfg)
    selected=str(cfg.get("selected_main_candidate","all_cells_projected"))
    srcbase=outdir/("Fig2_control_defined_state_atlas_candidateB" if selected=="all_cells_projected" else "Fig2_control_defined_state_atlas_candidateA")
    for ext in ["png","pdf","svg"]: shutil.copyfile(srcbase.with_suffix("."+ext),outdir/("Fig2_control_defined_state_atlas."+ext))
    # QC
    maxsil=min(20000,len(train_states)); rng=np.random.default_rng(int(cfg.get("random_seed",20260727))); si=rng.choice(len(train_states),maxsil,replace=False) if len(train_states)>maxsil else np.arange(len(train_states))
    sil=float(silhouette_score(train_latent[si],train_states[si])) if len(pd.unique(train_states))>1 else np.nan
    qc=pd.DataFrame([{"state_definition_source":source_name,"state_mapping_source":mapping_source,"n_source_cells":len(train_states),
                      "n_query_cells":query.n_obs,"n_pure_controls":int(control_mask.sum()),"n_states":len(pd.unique(assigned)),
                      "source_state_silhouette_qc":sil,"median_assignment_confidence":float(np.median(conf)),
                      "median_assignment_entropy":float(np.median(ent))}])
    qc.to_csv(outdir/"state_atlas_qc.tsv",sep="\t",index=False)
    # notes and legend
    legend=("**Figure 2. Control-defined states resolve the endothelial differentiation landscape.** "
            f"A fixed coarse-state atlas was defined using {source_name}. CC-Perturb-seq cells were projected into this atlas without allowing perturbed cells to influence feature selection, latent-space fitting, UMAP geometry, clustering or state-classifier training. "
            "Colors denote conservative developmental or alternative states supported by multiple markers and observed timing. The adjacent composition strip shows the fraction of projected cells assigned to each state within observed D0-D3 collection days. The display is cross-sectional and does not represent lineage tracing.")
    (outdir/"draft_panel_legend.md").write_text(legend+"\n")
    methods=("# Methods note\n\n"
             f"State-definition source: **{source_name}**. Mapping source: **{mapping_source}**.\n\n"
             f"Control definition: {control_def}\n\n"
             "Only the independent reference or pure controls contributed to gene selection, gene scaling, latent-space fitting, UMAP fitting, source clustering and k-nearest-neighbor classifier training. All other CC-Perturb-seq cells were transformed into the fixed model. Cells carrying both control and real-targeting guides were excluded from the pure-control training set but retained in the all-cell projection. Marker evidence was evaluated using multiple genes per state; a named state was not assigned from one marker alone. Silhouette was used only as a separation QC metric, not as biological validation.\n")
    (outdir/"methods_note.md").write_text(methods)
    (outdir/"software_versions.json").write_text(json.dumps({"python":sys.version,"platform":platform.platform(),"scanpy":sc.__version__,"anndata":ad.__version__,"numpy":np.__version__,"pandas":pd.__version__},indent=2))
    joblib.dump({"source_name":source_name,"genes":genes,"scaler":scaler,"svd":svd,"umap_model":reducer,"classifier":clf,"state_classes":list(clf.classes_),"mapping":mapping},outdir/"control_defined_state_atlas_model.joblib",compress=3)
    return outdir
