"""Publication display only. Read frozen plot tables; never fit/reselect/test."""
from pathlib import Path
import os
import json,datetime
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch,Rectangle
from PIL import ImageCms,Image
ICC=ImageCms.ImageCmsProfile(ImageCms.createProfile('sRGB')).tobytes()
ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parents[1]
DATA=Path(os.environ.get('ZEB1_DATA_ROOT',str(PROJECT/'data/processed')))/'total/rebuild_2026/data'; D=DATA/'deepen_2026'; SD=DATA/'source_data_rebuilt'; V=DATA/'editorial_visual_sources'
OUT=PROJECT/'results/figures'; OUT.mkdir(parents=True,exist_ok=True)
Z1='#E64B35'; Z2='#3C5488'; ETP='#00A087'; GREY='#A9A9A9'; LIGHT='#BDBDBD'; INK='#252525'; NA='#D5D5D5'
CMAP=LinearSegmentedColormap.from_list('locked',[Z2,'#FFFFFF',Z1]);CMAP.set_bad(NA)
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':8.5,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'axes.linewidth':.55,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white'})
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def color(s):
    s=str(s).lower();return Z2 if 'bcl11b' in s or 'zeb2' in s else ETP if 'etp' in s else Z1 if 'tlx3' in s or 'zeb1' in s else GREY
def short(s):
    if str(s)=='ETP-like':return 'ETP-like'
    for a,b in [('<U+03B3>','γ'),('<U+03B4>','δ'),('<U+03B1>','α'),('<U+03B2>','β')]:s=str(s).replace(a,b)
    return s.replace('-like','').replace('immature','imm.')
def fig(height,letter_only=True):
    f=plt.figure(figsize=(180/25.4,height/25.4));f._letter_only=letter_only;return f
def ax(f,rect,letter,title):
    a=f.add_axes(rect);head(a,letter,title);return a
def box(f,x,top,width,height):
    """Axes geometry in mm, measured from the top of the composite."""
    total=f.get_figheight()*25.4
    return [x/180,(total-top-height)/total,width/180,height/total]
def anchor(f,x,top):
    return (x/180,1-top/(f.get_figheight()*25.4))
def head(a,l,t):
    p=a.get_position()
    if l:a.figure.text(p.x0-.04,p.y1+.005,l,fontsize=11.5,weight='bold',gid='panel_letter')
    a._source_panel_description=t;a.tick_params(length=2,pad=2);return a
def grid(a,axis='y'):a.grid(axis=axis,color='#ECECEC',lw=.4);a.set_axisbelow(True)
def zeros(a):a.axhline(0,color=LIGHT,lw=.5);a.axvline(0,color=LIGHT,lw=.5)
def sidekey(f,y=.015):
    f.legend(handles=[Line2D([],[],marker='o',ls='',color=Z2,label='ZEB2 side'),Line2D([],[],marker='o',ls='',color=Z1,label='ZEB1 side')],loc='lower center',bbox_to_anchor=(.5,y),ncol=2,frameon=False)
def swatchkey(f,items,xy,ncol):
    key=f.legend(handles=[Patch(facecolor=co,label=label) for co,label in items],loc='center',bbox_to_anchor=xy,ncol=ncol,frameon=False,handlelength=.9,handletextpad=.35,columnspacing=.85,borderpad=0)
    f.add_artist(key)
def scale(f,im,rect,label,ticks):
    c=f.add_axes(rect);b=f.colorbar(im,cax=c,orientation='horizontal',ticks=ticks);b.ax.tick_params(length=2,pad=1);b.set_label(label,labelpad=1);return b
def checkpoint(stage,detail):
    return None

def save(f,n):
    from render_figure import capture
    capture(f,str(n))
