import os, zipfile, shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'DejaVu Sans'

outdir = Path('/mnt/data/Figure2_new_panel_mock_layouts_v2')
(outdir / 'combined_mock_layouts').mkdir(parents=True, exist_ok=True)
(outdir / 'individual_mock_panels').mkdir(parents=True, exist_ok=True)
(outdir / 'scripts').mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(44)

D1 = '#8A5A44'; D2 = '#2F6B8F'; D3 = '#6A4C93'; grey = '#9B9B9B'; light = '#E7E7E7'; dark='#202020'

def mock_usage_pair(n=3600, rho=-0.45, center=(0.48,0.48), sd=0.19):
    cov = np.array([[sd**2, rho*sd**2],[rho*sd**2, sd**2]])
    xy = rng.multivariate_normal(center, cov, n)
    x = np.clip(xy[:,0], 0, 1)
    y = np.clip(xy[:,1], 0, 1)
    # add mild tails to make hexes informative
    tails = rng.choice([0,1,2], n, p=[0.82,0.09,0.09])
    x[tails==1] = np.clip(rng.normal(0.20, 0.07, (tails==1).sum()),0,1)
    y[tails==1] = np.clip(rng.normal(0.73, 0.09, (tails==1).sum()),0,1)
    x[tails==2] = np.clip(rng.normal(0.73, 0.09, (tails==2).sum()),0,1)
    y[tails==2] = np.clip(rng.normal(0.20, 0.07, (tails==2).sum()),0,1)
    return x,y

pair_specs = [
    ('D1', 'P30 usage', 'P13 usage', -0.36, D1),
    ('D2', 'P15 usage', 'P17 usage', -0.41, D2),
    ('D3', 'P3 usage', 'P34 usage', -0.48, D3),
]
programs = [
    ('P30','D1',[0.05,0.74,0.20,0.05]),
    ('P13','D1',[0.04,0.68,0.16,0.05]),
    ('P15','D2',[0.03,0.11,0.75,0.15]),
    ('P17','D2',[0.03,0.09,0.70,0.20]),
    ('P3','D3',[0.02,0.04,0.13,0.70]),
    ('P34','D3',[0.02,0.03,0.12,0.58]),
]
days=['D0','D1','D2','D3']

def add_header(fig, title, subtitle=None):
    fig.text(0.055, 0.965, title, fontsize=15, fontweight='bold', ha='left', va='top')
    if subtitle:
        fig.text(0.055, 0.935, subtitle, fontsize=8.5, ha='left', va='top', color='#555555')

def draw_panelF(fig, bbox, label='F'):
    x0,y0,w,h = bbox
    fig.text(x0, y0+h+0.045, label, fontsize=15, fontweight='bold', va='bottom')
    fig.text(x0+0.035, y0+h+0.047, 'Program-pair usage structure', fontsize=11.5, fontweight='bold', va='bottom')
    fig.text(x0+0.035, y0+h+0.025, 'Raw cell-level usage for two individual cNMF programs', fontsize=7.8, color='#555555', va='bottom')
    gap = 0.025
    panel_w = (w - 2*gap)/3
    for i,(day,xlab,ylab,rho,col) in enumerate(pair_specs):
        ax = fig.add_axes([x0+i*(panel_w+gap), y0, panel_w, h])
        x,y = mock_usage_pair(rho=rho)
        ax.hexbin(x,y, gridsize=33, mincnt=1, cmap='Greys', linewidths=0, alpha=0.96, extent=(0,1,0,1))
        ax.plot([0,1],[1,0], color='#CFCFCF', lw=0.8, ls='--')
        ax.text(0.04,0.94, f'{day}\nr={rho:.2f}', transform=ax.transAxes, ha='left', va='top', fontsize=8.4, fontweight='bold', color=col)
        ax.set_xlim(0,1); ax.set_ylim(0,1)
        ax.set_xticks([0,0.5,1]); ax.set_yticks([0,0.5,1])
        ax.tick_params(labelsize=6.5, length=2)
        ax.set_xlabel(xlab, fontsize=8)
        ax.set_ylabel(ylab, fontsize=8)
        for sp in ax.spines.values(): sp.set_color('#555555'); sp.set_linewidth(0.7)

def draw_panelG_distribution(fig, bbox, label='G'):
    x0,y0,w,h = bbox
    fig.text(x0, y0+h+0.045, label, fontsize=15, fontweight='bold', va='bottom')
    fig.text(x0+0.035, y0+h+0.047, 'Temporal prevalence of selected programs', fontsize=11.5, fontweight='bold', va='bottom')
    fig.text(x0+0.035, y0+h+0.025, 'Cell-level usage distributions with replicate medians', fontsize=7.8, color='#555555', va='bottom')
    ax = fig.add_axes([x0,y0,w,h])
    ax.set_xlim(-0.55,3.55); ax.set_ylim(-0.55,5.55)
    for j in range(4):
        ax.add_patch(Rectangle((j-0.45,-0.5),0.9,6,facecolor='#FAFAFA' if j%2==0 else 'white', edgecolor='none', zorder=0))
    for row,(prog,day,vals) in enumerate(programs):
        yy = 5-row
        col = D1 if day=='D1' else D2 if day=='D2' else D3
        ax.hlines(yy, -0.4, 3.4, color=light, lw=0.7, zorder=0)
        ax.text(-0.48, yy, prog, ha='right', va='center', fontsize=9, fontweight='bold', color=col)
        for j,mu in enumerate(vals):
            u = np.clip(rng.normal(mu, 0.07, 180),0,1)
            xs = rng.normal(j,0.035,180)
            ys = yy + (u-0.5)*0.56
            active = mu > 0.45
            ax.scatter(xs,ys,s=4,color=col if active else grey,alpha=0.18 if active else 0.08, edgecolor='none')
            rep = np.clip(rng.normal(mu,0.04,3),0,1)
            ax.scatter([j-0.10,j,j+0.10], yy+(rep-0.5)*0.56, s=22, color=col if active else grey, edgecolor='white', lw=0.4, zorder=4, alpha=0.95)
    ax.set_xticks(range(4)); ax.set_xticklabels(days, fontsize=8, fontweight='bold')
    ax.set_yticks([]); ax.set_xlabel('Differentiation day', fontsize=8)
    ax.text(3.58,5.1,'cell usage\nreplicate median', ha='left', va='center', fontsize=7, color='#555555')
    for sp in ax.spines.values(): sp.set_visible(False)

def draw_panelG_tile(fig, bbox, label='G'):
    x0,y0,w,h = bbox
    fig.text(x0, y0+h+0.045, label, fontsize=15, fontweight='bold', va='bottom')
    fig.text(x0+0.035, y0+h+0.047, 'Stage enrichment', fontsize=11.5, fontweight='bold', va='bottom')
    fig.text(x0+0.035, y0+h+0.025, 'Mean usage by day', fontsize=7.8, color='#555555', va='bottom')
    ax = fig.add_axes([x0,y0,w,h])
    vals = np.array([p[2] for p in programs])
    ax.imshow(vals, cmap='Greys', vmin=0, vmax=0.80, aspect='auto')
    ax.set_xticks(range(4)); ax.set_xticklabels(days, fontsize=8, fontweight='bold')
    ax.set_yticks(range(6)); ax.set_yticklabels([p[0] for p in programs], fontsize=8, fontweight='bold')
    for i,(p,day,_) in enumerate(programs):
        ax.get_yticklabels()[i].set_color(D1 if day=='D1' else D2 if day=='D2' else D3)
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            if vals[i,j] > 0.45:
                ax.text(j,i,'●', ha='center',va='center',color='white',fontsize=8)
    for sp in ax.spines.values(): sp.set_color('#555555'); sp.set_linewidth(0.7)

def option1():
    fig = plt.figure(figsize=(11,8.4), facecolor='white')
    add_header(fig, 'New Figure 2 panels: descriptive program atlas', 'Option 1: density triptych over replicate-aware temporal prevalence')
    draw_panelF(fig,(0.075,0.54,0.85,0.285),'F')
    draw_panelG_distribution(fig,(0.14,0.13,0.72,0.27),'G')
    fig.text(0.94,0.025,'MOCK LAYOUT - simulated values', ha='right',fontsize=7.5,color='#777777')
    return fig

def option2():
    fig = plt.figure(figsize=(12.0,5.7), facecolor='white')
    add_header(fig, 'New Figure 2 panels: compact side-by-side version', 'Option 2: more main-figure friendly when vertical space is limited')
    draw_panelF(fig,(0.055,0.22,0.60,0.46),'F')
    draw_panelG_tile(fig,(0.72,0.28,0.22,0.33),'G')
    fig.text(0.94,0.04,'MOCK LAYOUT - simulated values', ha='right',fontsize=7.5,color='#777777')
    return fig

def option3():
    fig = plt.figure(figsize=(11.5,6.5), facecolor='white')
    add_header(fig, 'New Figure 2 panels: integrated atlas-strip version', 'Option 3: triptych plus a compact temporal strip')
    draw_panelF(fig,(0.08,0.38,0.84,0.34),'F')
    draw_panelG_tile(fig,(0.24,0.10,0.52,0.16),'G')
    fig.text(0.94,0.04,'MOCK LAYOUT - simulated values', ha='right',fontsize=7.5,color='#777777')
    return fig

outputs=[]
for name, func in [('option1_stacked',option1),('option2_side_by_side',option2),('option3_integrated_strip',option3)]:
    fig=func()
    base=outdir/'combined_mock_layouts'/f'Figure2_new_panels_mock_{name}'
    for ext in ['png','pdf','svg']:
        p=base.with_suffix('.'+ext)
        fig.savefig(p, dpi=300 if ext=='png' else None, bbox_inches='tight')
        outputs.append(p)
    plt.close(fig)
# individual F/G
fig=plt.figure(figsize=(9.8,3.7), facecolor='white'); draw_panelF(fig,(0.08,0.18,0.84,0.58),'F'); fig.text(0.94,0.03,'MOCK LAYOUT - simulated values',ha='right',fontsize=7.5,color='#777777')
base=outdir/'individual_mock_panels'/'Fig2F_program_pair_density_triptych_mock'
for ext in ['png','pdf','svg']:
    p=base.with_suffix('.'+ext); fig.savefig(p,dpi=300 if ext=='png' else None,bbox_inches='tight'); outputs.append(p)
plt.close(fig)
fig=plt.figure(figsize=(8.2,3.7), facecolor='white'); draw_panelG_distribution(fig,(0.16,0.18,0.70,0.52),'G'); fig.text(0.94,0.03,'MOCK LAYOUT - simulated values',ha='right',fontsize=7.5,color='#777777')
base=outdir/'individual_mock_panels'/'Fig2G_temporal_prevalence_mock'
for ext in ['png','pdf','svg']:
    p=base.with_suffix('.'+ext); fig.savefig(p,dpi=300 if ext=='png' else None,bbox_inches='tight'); outputs.append(p)
plt.close(fig)
# readme/script/package
(outdir/'README.md').write_text('''Figure 2 mock layout previews

These are mock layout previews with simulated values, not data-driven plots.

Panel F: raw program-pair usage density triptych (P30/P13, P15/P17, P3/P34).
Panel G: temporal prevalence of selected programs across D0-D3.

To make data-driven versions, upload the AnnData-derived export or re-upload the AnnData object.
''')
shutil.copyfile(__file__, outdir/'scripts'/'create_figure2_mock_layouts_v2.py')
zip_path=Path('/mnt/data/Figure2_new_panel_mock_layouts_v2_package.zip')
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for p in outdir.rglob('*'):
        z.write(p, arcname=str(p.relative_to(outdir.parent)))
print('Created v2 mock files')
for p in outputs:
    print(p)
print('ZIP:',zip_path)
