import os, zipfile
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patches
plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42; plt.rcParams['font.family']='DejaVu Sans'
OUT='/mnt/data/program_first_hippo_discovery_figures'
BASE1='/mnt/data/outputs_program_level_temporal_map_v4_unzipped'
BASE2='/mnt/data/outputs_program_first_advanced_unzipped'
POS='#E68613'; NEG='#2F76B7'
ann=pd.read_csv(os.path.join(BASE1,'curated_program_annotations_used.csv'))
corr=pd.read_csv(os.path.join(BASE1,'curated_program_program_pearson_correlation.csv')).set_index('program_id')
d2d3=pd.read_csv(os.path.join(BASE2,'D2_to_D3_raw_program_effect_correlations.csv'))
piv=d2d3.pivot(index='D2_program',columns='D3_program',values='pearson_r')

def get_corr(a,b):
    try: return float(corr.loc[a,b])
    except Exception: return float(corr.loc[b,a])

def draw_network(stem='lagged_program_coupling_network_minimal_clean', include_p27=False):
    pos={'iPSC':(0,0),'P30':(1.35,0.8),'P13':(1.35,-0.8),'P15':(3.0,0.8),'P17':(3.0,-0.8),'P3':(4.7,0.8),'P34':(4.7,-0.8)}
    if include_p27: pos['P27']=(4.7,-1.55)
    labels={'iPSC':'iPSC','P30':'P30\nanterior PS','P13':'P13\nposterior PS','P15':'P15\nvascular\nprogenitor','P17':'P17\nYAP/TAZ\nLM','P3':'P3\nEC endpoint','P34':'P34\nYAP/TAZ\nmesenchyme','P27':'P27\nECM\nmesenchyme'}
    fill={'iPSC':'#F4F4F4','P30':'#FFF1E5','P13':'#FFF1E5','P15':'#EAF3FB','P17':'#EAF3FB','P3':'#EAF3FB','P34':'#FCE8E4','P27':'#FCE8E4'}
    edges=[('P30','P15',get_corr('P30','P15'),0.08),('P30','P17',get_corr('P30','P17'),-0.10),('P13','P17',get_corr('P13','P17'),-0.08),('P13','P15',get_corr('P13','P15'),0.12),('P15','P17',get_corr('P15','P17'),0.0),('P15','P3',float(piv.loc['P15','P3']),0.06),('P15','P34',float(piv.loc['P15','P34']),-0.08),('P17','P34',float(piv.loc['P17','P34']),0.04)]
    if include_p27:
        edges += [('P15','P27',float(piv.loc['P15','P27']),-0.18),('P17','P27',float(piv.loc['P17','P27']),0.05)]
    fig, ax=plt.subplots(figsize=(8.4,4.8 if include_p27 else 4.0))
    ax.set_xlim(-0.45,5.25); ax.set_ylim(-2.05 if include_p27 else -1.45,1.35); ax.axis('off')
    # day labels
    for day,x in [('D0',0),('D1',1.35),('D2',3.0),('D3',4.7)]:
        ax.text(x,1.18,day,ha='center',va='center',fontsize=13,weight='bold')
        ax.plot([x-0.35,x+0.35],[1.06,1.06],color='#222222',lw=0.8)
    # start arrows
    for node in ['P30','P13']:
        ax.annotate('',xy=(pos[node][0]-0.48,pos[node][1]),xytext=(pos['iPSC'][0]+0.30,0),arrowprops=dict(arrowstyle='->',lw=1.0,color='#606060'))
    # nodes first (behind edges? no, edges need behind nodes)
    # edges
    for a,b,r,rad in edges:
        (x1,y1),(x2,y2)=pos[a],pos[b]
        color=POS if r>0 else NEG
        lw=0.7+4.5*abs(r)
        ls='--' if abs(r)<0.30 else '-'
        con=patches.FancyArrowPatch((x1+0.38,y1),(x2-0.42,y2),arrowstyle='-|>',mutation_scale=10,
                                     linewidth=lw,linestyle=ls,color=color,alpha=0.90,connectionstyle=f'arc3,rad={rad}',zorder=1)
        ax.add_patch(con)
        lx,ly=(x1+x2)/2,(y1+y2)/2
        # nudges
        if (a,b)==('P15','P17'): lx+=0.38; ly+=0.00
        elif y1==y2: ly+=0.15
        else: ly+=0.13 if r>0 else -0.05
        if b=='P27': ly-=0.04
        ax.text(lx,ly,f'{r:+.2f}',color=color,fontsize=8.5,weight='bold',ha='center',va='center',bbox=dict(fc='white',ec='none',alpha=0.7,pad=0.6),zorder=3)
    # nodes
    for node,(x,y) in pos.items():
        w=0.74 if node=='iPSC' else 0.84
        h=0.38 if node=='iPSC' else 0.55
        rect=patches.FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.035,rounding_size=0.06',linewidth=1.0,edgecolor='#333333',facecolor=fill[node],zorder=4)
        ax.add_patch(rect)
        ax.text(x,y,labels[node],ha='center',va='center',fontsize=8.1 if node!='iPSC' else 10,weight='bold',zorder=5)
    # legend bottom left
    yleg=-1.28 if not include_p27 else -1.92
    ax.plot([0.0,0.35],[yleg,yleg],color=POS,lw=3)
    ax.text(0.42,yleg,'positive correlation',va='center',fontsize=8)
    ax.plot([1.55,1.90],[yleg,yleg],color=NEG,lw=3)
    ax.text(1.97,yleg,'negative correlation',va='center',fontsize=8)
    ax.text(3.55,yleg,'line width ~ |r|; dashed = weaker link',va='center',fontsize=8,color='#444444')
    ax.set_title('Program-effect coupling across D1–D3',fontsize=14,weight='bold',pad=12)
    for ext in ['png','pdf','svg']:
        fig.savefig(os.path.join(OUT,f'{stem}.{ext}'),dpi=350 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)

draw_network(stem='lagged_program_coupling_network_minimal_clean', include_p27=False)
draw_network(stem='lagged_program_coupling_network_with_P27_clean', include_p27=True)

# update zip
zip_path='/mnt/data/program_first_hippo_discovery_figures_package.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for root,_,files in os.walk(OUT):
        for fn in files:
            fp=os.path.join(root,fn)
            z.write(fp,arcname=os.path.relpath(fp,OUT))
print('updated',zip_path)
