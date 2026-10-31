#!/usr/bin/env python3
from pathlib import Path
import argparse, re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Rectangle

DEFAULT_BASE = Path('/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50')
PROGRAMS = [34, 27, 24, 2]
PROGRAM_LABELS = {34:'P34\nYAP/TAZ mesenchyme',27:'P27\nalternative context',24:'P24\nalternative context',2:'P2\narterial EC'}
MARKERS = {34:['ANKRD1','KRT8','KRT18'],27:[],24:['GATA4','RSPO2'],2:['MECOM','NRP1','EFNB2']}

mpl.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.transparent':False,'figure.facecolor':'white','axes.facecolor':'white','font.family':'DejaVu Sans'})

def read_matrix(path):
    for kwargs in [dict(sep='\t'), dict(sep=r'\s+', engine='python'), dict(sep=',')]:
        try:
            df = pd.read_csv(path, index_col=0, **kwargs)
            df = df.loc[:, ~df.columns.astype(str).str.startswith('Unnamed')]
            if df.shape[0] > 1 and df.shape[1] > 1:
                break
        except Exception:
            df = None
    if df is None:
        raise ValueError(f'Could not read {path}')
    # orient as programs x genes
    if df.shape[0] > 200 and df.shape[1] <= 200:
        df = df.T
    elif not (df.shape[0] <= 200 and df.shape[1] > 200):
        idx_num = np.mean([bool(re.fullmatch(r'\s*(?:P|Program)?\s*\d+\s*', str(x), re.I)) for x in df.index])
        col_num = np.mean([bool(re.fullmatch(r'\s*(?:P|Program)?\s*\d+\s*', str(x), re.I)) for x in df.columns])
        if col_num > idx_num:
            df = df.T
    df.columns = df.columns.astype(str).str.strip()
    return df

def parse_program(x):
    m = re.search(r'(\d+)', str(x))
    return int(m.group(1)) if m else None

def raw_lookup(df):
    out = {}
    for label in df.index:
        p = parse_program(label)
        if p is not None and p not in out:
            out[p] = label
    return out

def choose_offset(df):
    lookup = raw_lookup(df)
    upper = {g.upper():g for g in df.columns}
    scores = {}
    for off in [0,1,-1]:
        vals=[]
        for p, genes in MARKERS.items():
            if p+off not in lookup: continue
            row=lookup[p+off]
            for g in genes:
                gg=upper.get(g.upper())
                if gg is not None:
                    v=pd.to_numeric(pd.Series([df.loc[row,gg]]),errors='coerce').iloc[0]
                    if np.isfinite(v): vals.append(float(v))
        scores[off]=np.mean(vals) if vals else -np.inf
    print('Offset scores:', scores)
    return max(scores,key=scores.get)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--base-dir', type=Path, default=DEFAULT_BASE)
    ap.add_argument('--matrix-type', choices=['consensus','gene_score'], default='consensus')
    ap.add_argument('--outdir', type=Path, default=None)
    ap.add_argument('--p27-fallback-top-n', type=int, default=0)
    args=ap.parse_args()

    pattern='*.spectra.k_50.dt_0_2.consensus.txt' if args.matrix_type=='consensus' else '*.gene_spectra_score.k_50.dt_0_2.txt'
    matches=sorted(args.base_dir.glob(pattern))
    if not matches: raise FileNotFoundError(f'No match for {pattern} in {args.base_dir}')
    path=matches[0]
    print('Reading:', path)
    df=read_matrix(path)
    offset=choose_offset(df)
    lookup=raw_lookup(df)
    rows=[]
    for p in PROGRAMS:
        if p+offset not in lookup: raise KeyError(f'Missing P{p}; expected raw row {p+offset}')
        s=df.loc[lookup[p+offset]].copy(); s.name=p; rows.append(s)
    pm=pd.DataFrame(rows).apply(pd.to_numeric, errors='coerce')

    markers={p:list(g) for p,g in MARKERS.items()}
    if args.p27_fallback_top_n>0:
        markers[27]=pm.loc[27].dropna().sort_values(ascending=False).head(args.p27_fallback_top_n).index.tolist()
        print('P27 fallback genes:', markers[27])

    genes=[]
    for p in PROGRAMS:
        for g in markers[p]:
            if g not in genes: genes.append(g)
    upper={g.upper():g for g in pm.columns}
    mat=pd.DataFrame(index=PROGRAMS, columns=genes, dtype=float)
    missing=[]
    for p in PROGRAMS:
        for g in genes:
            actual=upper.get(g.upper())
            if actual is None: missing.append(g)
            else: mat.loc[p,g]=pm.loc[p,actual]

    out=args.outdir or (args.base_dir/'fig4_program_marker_loading_heatmap')
    out.mkdir(parents=True, exist_ok=True)
    mat.to_csv(out/'source_marker_loading_matrix.csv')
    pd.DataFrame([{'program':p,'gene':g,'from_notes':not(p==27 and args.p27_fallback_top_n>0)} for p in PROGRAMS for g in markers[p]]).to_csv(out/'source_marker_metadata.csv',index=False)
    pd.DataFrame({'missing_gene':sorted(set(missing))}).to_csv(out/'source_missing_genes.csv',index=False)

    vals=mat.to_numpy(float); finite=vals[np.isfinite(vals)]
    if finite.size==0: raise ValueError('No selected genes found in matrix')
    if np.nanmin(finite)<0:
        vmax=max(abs(np.nanmin(finite)),abs(np.nanmax(finite))); norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax); cmap=mpl.colormaps['RdBu_r']
    else:
        norm=Normalize(vmin=0,vmax=max(np.nanmax(finite),1e-12)); cmap=mpl.colormaps['Reds']

    nrow,ncol=vals.shape
    fig=plt.figure(figsize=(max(8,0.62*ncol+3.2),max(5.2,0.62*nrow+2.6)),facecolor='white')
    ax=fig.add_axes([0.20,0.24,min(0.66,0.62*ncol/fig.get_figwidth()),min(0.58,0.62*nrow/fig.get_figheight())])
    for i in range(nrow):
        for j in range(ncol):
            v=vals[i,j]; color=(0.94,0.94,0.94,1) if not np.isfinite(v) else (*cmap(norm(v))[:3],1)
            ax.add_patch(Rectangle((j,i),1,1,facecolor=color,edgecolor='white',linewidth=1.4))
    ax.set_xlim(0,ncol); ax.set_ylim(nrow,0); ax.set_aspect('equal')
    ax.set_xticks(np.arange(ncol)+0.5); ax.set_xticklabels(genes,rotation=58,ha='right',fontsize=10.5,fontstyle='italic')
    ax.set_yticks(np.arange(nrow)+0.5); ax.set_yticklabels([PROGRAM_LABELS[p] for p in PROGRAMS],fontsize=12)
    ax.tick_params(length=0,pad=7)
    for s in ax.spines.values(): s.set_visible(False)

    gene_source={g:p for p in PROGRAMS for g in markers[p]}
    start=0
    while start<ncol:
        source=gene_source[genes[start]]; end=start+1
        while end<ncol and gene_source[genes[end]]==source: end+=1
        ax.text((start+end)/2,-0.34,f'P{source} markers',ha='center',va='bottom',fontsize=9.5,fontweight='bold',clip_on=False)
        if end<ncol: ax.plot([end,end],[0,nrow],color='white',linewidth=3.0)
        start=end

    title='consensus loading' if args.matrix_type=='consensus' else 'gene spectra score'
    fig.suptitle(f'Program marker genes\n{title}',x=0.48,y=0.95,fontsize=17)
    cax=fig.add_axes([ax.get_position().x1+0.025,ax.get_position().y0,0.02,ax.get_position().height])
    sm=mpl.cm.ScalarMappable(norm=norm,cmap=cmap); sm.set_array([])
    cb=fig.colorbar(sm,cax=cax); cb.set_label('cNMF gene loading' if args.matrix_type=='consensus' else 'cNMF gene score',fontsize=11,labelpad=10); cb.ax.tick_params(labelsize=9)
    fig.text(0.20,0.065,f'Marker genes from notes table; program offset {offset:+d}. P27 has no notes markers unless fallback enabled.',fontsize=8.5,color='#444444')

    stem=out/f'Program_marker_gene_{args.matrix_type}_dt_0_2'
    fig.savefig(stem.with_suffix('.png'),dpi=600,bbox_inches='tight',facecolor='white',transparent=False)
    fig.savefig(stem.with_suffix('.pdf'),bbox_inches='tight',facecolor='white',transparent=False)
    fig.savefig(stem.with_suffix('.svg'),bbox_inches='tight',facecolor='white',transparent=False)
    plt.close(fig)
    print('Done:', out)

if __name__=='__main__':
    main()
