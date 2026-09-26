"""Shared publication geometry, colours and vector export."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

ROOT=Path(__file__).resolve().parents[1]
SD=ROOT/'source_data'
TEAL='#167C80'; OCHRE='#BE7434'; INK='#263442'; GRAY='#80909A'; LIGHT='#E5EAEB'; PURPLE='#79618C'
CMAP=LinearSegmentedColormap.from_list('zeb_balance',[TEAL,'#F5F4EF',OCHRE])
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':.55,'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2.5,'ytick.major.size':2.5,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white','figure.facecolor':'white','lines.linewidth':1.1,'legend.frameon':False,'legend.fontsize':7})
def read(name): return pd.read_csv(SD/(name+'.tsv'),sep='\t')
def original(rel): return pd.read_csv(SD/'inputs'/rel,sep='\t')
def canvas(number,title,height=175):
    f=plt.figure(figsize=(180/25.4,height/25.4)); f._mm_height=height
    f.text(9/180,1-7/height,f'{number:02d}',fontsize=15,fontweight='bold',color=TEAL,va='top')
    f.text(23/180,1-8/height,title,fontsize=11,fontweight='bold',va='top')
    return f
def ax(f,x,y,w,h,label=None,title=None):
    a=f.add_axes([x/180,1-(y+h)/f._mm_height,w/180,h/f._mm_height])
    if label: f.text((x-7)/180,1-(y-3)/f._mm_height,label,fontweight='bold',fontsize=10,va='top')
    if title: a.set_title(title,loc='left',pad=7,fontweight='bold')
    return a
def zero(a,vertical=False):
    (a.axvline if vertical else a.axhline)(0,color='#ADB6BC',lw=.65,ls=(0,(3,3)),zorder=0)
def tidy(a,grid='y'):
    if grid:a.grid(axis=grid,color='#E9EDEE',lw=.45,zorder=0)
    a.set_axisbelow(True)
def footer(f,text):
    # Keep detailed qualifications in the separate scientific legend, not in
    # tiny, crowded text below the final-size axes.
    f._legend_note=text
def save(f,name,section='main'):
    d=ROOT/'figures'/section;d.mkdir(parents=True,exist_ok=True)
    for ext in ['pdf','svg','png','tiff']:
        kw={'dpi':600 if ext=='tiff' else 300}
        if ext=='tiff':kw['pil_kwargs']={'compression':'tiff_lzw'}
        f.savefig(d/(name+'.'+ext),**kw)
    # Explicit print-scale review raster: 120 dpi at the saved physical size.
    f.savefig(ROOT/'review'/(name+'_120dpi.png'),dpi=120)
    if hasattr(f,'_legend_note'):
        (ROOT/'logs'/(name+'_legend_note.txt')).write_text(f._legend_note,encoding='utf-8')
    plt.close(f)
    print('SAVED',name,flush=True)
def short(s):
    return str(s).replace('TAL1 αβ-like','TAL1 αβ').replace('TAL1 DP-like','TAL1 DP').replace('LMO2 γδ-like','LMO2 γδ').replace('TME-enriched','TME enriched')
def colour(v,limit=2): return CMAP(Normalize(-limit,limit,clip=True)(v))
def strip(a,groups,values,labels,colors,vertical=True,seed=11,size=12):
    rng=np.random.default_rng(seed)
    for i,(key,label,c) in enumerate(zip(groups,labels,colors)):
        v=np.asarray(values[key],float); jitter=rng.uniform(-.16,.16,len(v))
        if vertical:
            a.scatter(i+jitter,v,s=size,c=c,edgecolors='white',linewidths=.3,zorder=3)
            a.plot([i-.22,i+.22],[np.median(v)]*2,c=INK,lw=1.3,zorder=4)
        else:
            a.scatter(v,i+jitter,s=size,c=c,edgecolors='white',linewidths=.3,zorder=3)
            a.plot([np.median(v)]*2,[i-.22,i+.22],c=INK,lw=1.3,zorder=4)
    if vertical:a.set_xticks(range(len(groups)),labels)
    else:a.set_yticks(range(len(groups)),labels)
