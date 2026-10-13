#!/usr/bin/env python3
"""Read-only discovery of candidate single-cell objects and state-atlas assets on Sherlock/OAK."""
from __future__ import annotations
import argparse, json, os
from pathlib import Path
from typing import Any, Dict, List
import pandas as pd

H5_EXT={'.h5ad','.h5mu','.loom','.h5','.hdf5'}
SCRIPT_EXT={'.py','.ipynb','.r','.rmd','.sh'}
TABLE_EXT={'.csv','.tsv','.txt','.parquet','.feather','.xlsx','.xls','.yaml','.yml'}
KEYWORDS=[
 'perturb','cc-perturb','cc_perturb','300gene','multiome','reference','d0','d1','d2','d3',
 'state','cell_type','annotation','marker','umap','latent','projection','label_transfer',
 'guide_assignment','safe-target','non-target','palette','color','h5ad','adata'
]

def args():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--roots',nargs='+',required=True)
 p.add_argument('--outdir',required=True)
 p.add_argument('--max-depth',type=int,default=8)
 p.add_argument('--max-files',type=int,default=250000)
 p.add_argument('--inspect-h5ad',action='store_true')
 return p.parse_args()

def score(path:Path)->int:
 s=str(path).lower(); return sum(int(k in s) for k in KEYWORDS)

def role(path:Path)->str:
 s=str(path).lower(); ext=path.suffix.lower()
 if ext in H5_EXT:
  if 'multiome' in s or 'reference' in s: return 'candidate_reference_object'
  if 'perturb' in s or '300gene' in s: return 'candidate_query_object'
  return 'candidate_single_cell_object'
 if ext in SCRIPT_EXT:
  if any(k in s for k in ['project','transfer','label','state','umap']): return 'candidate_projection_or_state_script'
  return 'candidate_analysis_script'
 if ext in TABLE_EXT:
  if any(k in s for k in ['palette','color']): return 'candidate_palette'
  if any(k in s for k in ['state','annotation','marker']): return 'candidate_state_annotation'
  if 'guide' in s: return 'candidate_guide_metadata'
  return 'candidate_table'
 return 'other'

def inspect_h5ad(path:Path)->Dict[str,Any]:
 out={}
 try:
  import anndata as ad
  obj=ad.read_h5ad(path,backed='r')
  out.update({
   'n_obs':int(obj.n_obs),'n_vars':int(obj.n_vars),
   'obs_columns':json.dumps(list(obj.obs.columns)),
   'var_columns':json.dumps(list(obj.var.columns)),
   'obsm_keys':json.dumps(list(obj.obsm.keys())),
   'uns_keys':json.dumps(list(obj.uns.keys())),
   'layers':json.dumps(list(obj.layers.keys())),
   'obsp_keys':json.dumps(list(obj.obsp.keys())),
   'raw_present':bool(obj.raw is not None),
  })
  summaries={}
  for c in ['day','stage','timepoint','biological_sample','replicate','sample_id','cell_type']:
   if c in obj.obs.columns:
    summaries[c]=obj.obs[c].astype(str).value_counts(dropna=False).head(40).to_dict()
  out['obs_value_summaries']=json.dumps(summaries)
  try: obj.file.close()
  except Exception: pass
 except Exception as e: out['inspection_error']=repr(e)
 return out

def main():
 a=args(); outdir=Path(a.outdir); outdir.mkdir(parents=True,exist_ok=True)
 rows:List[Dict[str,Any]]=[]; seen=0
 for root_s in a.roots:
  root=Path(root_s).expanduser().resolve()
  if not root.exists():
   rows.append({'root':str(root),'path':str(root),'exists':False,'role_guess':'missing_root'}); continue
  for dirpath,dirnames,filenames in os.walk(root):
   d=Path(dirpath)
   try: depth=len(d.relative_to(root).parts)
   except Exception: depth=999
   if depth>a.max_depth:
    dirnames[:]=[]; continue
   dirnames[:]=[x for x in dirnames if x not in {'.git','__pycache__','.snakemake','node_modules'}]
   for name in filenames:
    seen+=1
    if seen>a.max_files: break
    p=d/name; ext=p.suffix.lower(); sc=score(p)
    if ext not in H5_EXT|SCRIPT_EXT|TABLE_EXT and sc==0: continue
    try:
     st=p.stat(); row={'root':str(root),'path':str(p),'name':p.name,'suffix':ext,
      'size_bytes':st.st_size,'mtime':pd.to_datetime(st.st_mtime,unit='s',utc=True).isoformat(),
      'keyword_score':sc,'role_guess':role(p),'exists':True}
     if a.inspect_h5ad and ext=='.h5ad': row.update(inspect_h5ad(p))
     rows.append(row)
    except OSError as e: rows.append({'root':str(root),'path':str(p),'exists':True,'inspection_error':repr(e)})
   if seen>a.max_files: break
 inv=pd.DataFrame(rows)
 if not inv.empty: inv=inv.sort_values(['keyword_score','size_bytes'],ascending=[False,False],na_position='last')
 inv.to_csv(outdir/'single_cell_input_inventory.tsv',sep='\t',index=False)
 if not inv.empty:
  likely=inv[(inv['keyword_score'].fillna(0)>0)|inv['suffix'].isin(H5_EXT)].copy()
 else: likely=inv
 likely.to_csv(outdir/'single_cell_likely_inputs.tsv',sep='\t',index=False)
 print(f'Scanned {seen:,} files')
 print(outdir/'single_cell_input_inventory.tsv')
 print(outdir/'single_cell_likely_inputs.tsv')
if __name__=='__main__': main()
