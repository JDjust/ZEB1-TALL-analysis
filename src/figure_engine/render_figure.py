"""Render locked source plots and export editable panel PDFs; no data download."""
from pathlib import Path
import argparse, hashlib, importlib, io, json, os, re, sys, tempfile, time
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.fonttype']='none'
import matplotlib.pyplot as plt
import pandas as pd
import fitz

PROJECT=Path(__file__).resolve().parents[2]
TARGET=None
SELECTED=None
JOURNAL='canonical'
DEST=None
FILES=[]
CAPTURED=False
REGISTRY=[]

def isolated_pdf(f,axes=(),texts=(),legends=()):
    keep=set(axes)|set(texts)|set(legends)
    artists=set(f.axes)|set(f.legends)|set(f.texts)|set(f.patches)|set(f.artists)
    hidden=[a for a in artists if a not in keep and a.get_visible()]
    for artist in hidden:artist.set_visible(False)
    buf=io.BytesIO()
    try:f.savefig(buf,format='pdf',dpi=600)
    finally:
        for artist in hidden:artist.set_visible(True)
    return fitz.open(stream=buf.getvalue(),filetype='pdf')

def register(a,letter,title):
    a._panel_letter=letter
    a._panel_description=title
    if letter:REGISTRY.append((letter,a))

def _render(f,number):
    global CAPTURED
    if str(number)!=TARGET:return
    CAPTURED=True
    is_s=TARGET.startswith('S');folder='SFigure' if is_s else 'Figure'
    name=('SFigure'+TARGET[1:]) if is_s else ('Figure'+TARGET)
    temp=None
    if SELECTED and not DEST and JOURNAL=='canonical':
        temp=tempfile.TemporaryDirectory(prefix='zeb1_panel_')
        output=Path(temp.name).resolve()
        if not output.is_relative_to(Path(tempfile.gettempdir()).resolve()):raise RuntimeError('Unexpected temporary directory')
    else:
        output=DEST or (PROJECT/'results/figures'/name/JOURNAL if JOURNAL!='canonical' else PROJECT/'results/figures'/name)
    output.mkdir(parents=True,exist_ok=True)
    # Main scripts are locked at 180 mm; journal exports preserve point size.
    if JOURNAL!='canonical':
        f.set_size_inches((177 if JOURNAL=='Blood Advances' else 180)/25.4,f.get_figheight())
        for t in f.findobj(matplotlib.text.Text):
            if t.get_gid()=='panel_letter' and not is_s:t.set_text(t.get_text().upper() if JOURNAL=='Blood Advances' else t.get_text().lower())
    for t in f.findobj(matplotlib.text.Text):
        t.set_text(t.get_text().translate(str.maketrans('₀₁₂₃₄₅₆₇₈₉','0123456789')))
    f.canvas.draw()
    texts=[t for t in f.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    canvas=f.bbox
    outside=[]
    for t in texts:
        bb=t.get_window_extent(f.canvas.get_renderer())
        if not t.get_clip_on() and (bb.x0 < -1 or bb.y0 < -1 or bb.x1 > canvas.x1+1 or bb.y1 > canvas.y1+1):
            outside.append(t.get_text())
    base=f'Figure{TARGET}'
    pdf=output/(base+'.pdf');png=output/(base+'.png')
    f.savefig(pdf,dpi=600);f.savefig(png,dpi=300)
    f.savefig(output/(base+'.svg'))
    # Save panel groups by source letter. Blank continuation axes join the last
    # labelled group; explicit legacy subgrids are handled below.
    groups={};last=None
    for a in f.axes:
        label=getattr(a,'_panel_letter',None)
        if label:last=label
        if label is not None and last:groups.setdefault(last,[]).append(a)
    if TARGET=='S10':
        groups={'A':f.axes[:5],'C':[f.axes[5]],'B':f.axes[6:8],'D':[f.axes[8]]}
    if TARGET=='S8':
        # Source sc8(): a, b1, b2, c, d, e, footer colourbar.
        groups={'A':[f.axes[0]],'B':f.axes[1:3],'C':[f.axes[3]],'D':[f.axes[4]],'E':[f.axes[5]]}
    if TARGET=='S12':
        # Source function creates an unlettered companion to its E panel.
        pass
    if not groups:raise RuntimeError(f'No panel labels in Figure{TARGET}')
    renderer=f.canvas.get_renderer();height=f.get_figheight()*72;width=f.get_figwidth()*72
    panels=[]
    for letter,axes in groups.items():
        if SELECTED and SELECTED!=letter:continue
        # Tight extents retain axis ticks, legends attached to axes and titles.
        boxes=[a.get_tightbbox(renderer).transformed(f.dpi_scale_trans.inverted()) for a in axes]
        xmin=min(b.x0 for b in boxes)*72-8;xmax=max(b.x1 for b in boxes)*72+8
        ymin=min(b.y0 for b in boxes)*72-12;ymax=max(b.y1 for b in boxes)*72+14
        own_letters=[t for t in f.texts if t.get_gid()=='panel_letter' and t.get_text().upper()==letter]
        for t in own_letters:
            box=t.get_window_extent(renderer).transformed(f.dpi_scale_trans.inverted())
            xmin=min(xmin,box.x0*72-3);xmax=max(xmax,box.x1*72+3);ymax=max(ymax,box.y1*72+3)
        rect=fitz.Rect(max(0,xmin),max(0,height-ymax),min(width,xmax),min(height,height-ymin))
        # Main figure panel letters and global keys are figure-level artists;
        # keep them in the plate and document shared keys for standalone views.
        paneldir=(output/'panels'/letter) if DEST or JOURNAL!='canonical' else PROJECT/'scripts/06_figures'/folder/name/letter
        paneldir.mkdir(parents=True,exist_ok=True)
        specific_keys=(TARGET=='7' and letter=='F')
        legends=[] if specific_keys else f.legends
        doc=isolated_pdf(f,legends=legends);sub=fitz.open()
        keys=[]
        for legend in legends:
            box=legend.get_window_extent(renderer).transformed(f.dpi_scale_trans.inverted())
            key=fitz.Rect(max(0,box.x0*72-3),max(0,height-box.y1*72-3),min(width,box.x1*72+3),min(height,height-box.y0*72+3))
            if not key.is_empty:keys.append(key)
        pagewidth=max([rect.width]+[k.width for k in keys])
        keyheight=(18+sum(k.height+4 for k in keys)) if keys else 0
        if specific_keys:
            pagewidth=max(pagewidth,220)
            keyheight=76
        page=sub.new_page(width=pagewidth,height=rect.height+keyheight)
        # Export the actual plot artists in isolation. Hide neighbouring axes and
        # other panel letters before clipping, so nearby panels cannot leak in.
        owntexts=[]
        for t in f.texts:
            gid=t.get_gid() or ''
            if gid=='panel_letter':
                if t.get_text().upper()==letter:owntexts.append(t)
            elif str(gid).startswith('panel_title_'):
                if gid=='panel_title_'+letter:owntexts.append(t)
            else:
                x,y=t.get_position()
                def distance(a):
                    q=a.get_position()
                    return max(q.x0-x,0,x-q.x1)**2+max(q.y0-y,0,y-q.y1)**2
                labelled=[a for g in groups.values() for a in g]
                if labelled and min(labelled,key=distance) in axes:owntexts.append(t)
        scene=isolated_pdf(f,axes=axes,texts=owntexts)
        page.show_pdf_page(fitz.Rect(0,0,rect.width,rect.height),scene,0,clip=rect)
        scene.close()
        if keys:
            page.insert_text((4,rect.height+11),'Shared figure keys',fontsize=7)
            y=rect.height+18
            for key in keys:
                page.show_pdf_page(fitz.Rect(0,y,key.width,y+key.height),doc,0,clip=key)
                y+=key.height+4
        if specific_keys:
            # F encodes cohort AND adult subtype; the original shared footer
            # explains only cohort. Add the actual subtype colours here without
            # changing the complete plate or its underlying observations.
            common=importlib.import_module('plot_common')
            from matplotlib.colors import to_rgb
            page.insert_text((4,rect.height+11),'Panel F keys',fontsize=8)
            entries=[('Discovery',common.GREY,False),('Adult: BCL11B',common.color('BCL11B'),False),
                     ('Adult: ETP-like',common.color('ETP-like'),False),('Adult: TAL1 DP',common.color('TAL1 DP'),False),
                     ('Adult: TLX3 (n=2)',common.color('TLX3'),True)]
            for i,(label,colour,open_marker) in enumerate(entries):
                x=8+(i%2)*110;y=rect.height+26+(i//2)*15;rgb=to_rgb(colour)
                page.draw_circle((x,y-2.5),2.5,color=rgb,fill=None if open_marker else rgb,width=.8)
                page.insert_text((x+7,y),label,fontsize=7)
        sub.save(paneldir/f'panel{letter}.pdf',garbage=4,deflate=True)
        page.get_pixmap(dpi=300,alpha=False).save(paneldir/f'panel{letter}.png')
        (paneldir/f'panel{letter}.svg').write_text(page.get_svg_image(text_as_path=False),encoding='utf8')
        sub.close();doc.close()
        desc='; '.join(dict.fromkeys(getattr(a,'_panel_description','') for a in axes if getattr(a,'_panel_description','')))
        panels.append({'panel':letter,'description':desc,'crop_pt':list(rect),'source_axes':len(axes)})
    if SELECTED and not panels:raise RuntimeError(f'Figure{TARGET} has no panel {SELECTED}')
    inventory={'figure':TARGET,'engine':engine_name(TARGET),'inputs':sorted(set(FILES)),'input_trace_scope':'read_csv paths; see data manifest for other input formats','panels':panels,'width_mm':f.get_figwidth()*25.4,'height_mm':f.get_figheight()*25.4,'journal':JOURNAL,'minimum_font_pt':min(t.get_fontsize() for t in texts),'text_outside_canvas':outside,'palette_name':'zeb_semantic_locked','palette_rationale':'Preserves the manuscript-defined ZEB1/ZEB2/ETP roles; custom semantic mapping documented rather than attributed to Nature.','rendered_at':time.strftime('%Y-%m-%dT%H:%M:%S')}
    inventory['input_sha256']={p:hashlib.sha256((PROJECT/p).read_bytes()).hexdigest() for p in inventory['inputs'] if (PROJECT/p).is_file()}
    receipt_dir=paneldir if SELECTED and not DEST else output
    (receipt_dir/'render_receipt.json').write_text(json.dumps(inventory,indent=2,ensure_ascii=False),encoding='utf8')
    print(json.dumps({'figure':TARGET,'panels':[p['panel'] for p in panels],'inputs':len(set(FILES)),'output':str(receipt_dir)},ensure_ascii=False),flush=True)
    plt.close(f)
    if temp:temp.cleanup()

def capture(f,n):_render(f,n)

def engine_name(target):
    if not target.startswith('S'):return f'figure{target}.py'
    n=int(target[1:])
    if n>=11:return f'supplement{n}.py'
    if n==6:return 'locked_supplement6.py'
    if n in [3,9]:return 'editorial_upgrade_supplement.py'
    if n==10:return 'editorial_upgrade_10.py'
    return 'editorial_upgrade_supplement_rest.py'

def run(target,panel=None,journal='canonical',out=None):
    global TARGET,SELECTED,JOURNAL,DEST,CAPTURED,FILES,REGISTRY
    TARGET=target.removeprefix('Figure').removeprefix('SFigure');SELECTED=panel;JOURNAL=journal;DEST=Path(out) if out else None
    CAPTURED=False;FILES=[];REGISTRY=[]
    original_read=pd.read_csv
    def read(path,*args,**kwargs):
        if isinstance(path,(str,Path)):
            p=Path(path)
            if p.exists():
                resolved=p.resolve()
                FILES.append(str(resolved.relative_to(PROJECT)) if resolved.is_relative_to(PROJECT) else str(resolved))
        return original_read(path,*args,**kwargs)
    pd.read_csv=read
    try:
        if not TARGET.startswith('S') or int(TARGET[1:])>=11:
            common=importlib.import_module('plot_common');head=common.head
            def tagged(a,l,t):
                register(a,l,t)
                result=head(a,l,t)
                if TARGET.startswith('S'):
                    pos=a.get_position()
                    if l:
                        label=a.figure.texts[-1]
                        label.set_position((pos.x0-.025,pos.y1+.013))
                    a.figure.text(pos.x0,pos.y1+10/(a.figure.get_figheight()*72),t,fontsize=8.5,weight='bold',gid='panel_title_'+l)
                return result
            common.head=tagged
            original_scale=common.scale
            def scale(*args,**kwargs):
                colorbar=original_scale(*args,**kwargs)
                register(colorbar.ax,'','Color scale')
                return colorbar
            common.scale=scale
            importlib.import_module(engine_name(TARGET)[:-3]).main()
            common.head=head
            common.scale=original_scale
        else:
            b=importlib.import_module('editorial_rebuild_figures')
            style=importlib.import_module('final_visual_style_lock')
            original_panel=style.panel
            def tagged(a,l,t,subtitle=None):register(a,l,t);return original_panel(a,l,t,subtitle)
            style.panel=tagged;b.panel=tagged
            module=importlib.import_module(engine_name(TARGET)[:-3]);module.panel=tagged
            # Existing derived detection summary avoids reopening raw h5ad.
            if TARGET=='S1':
                def detection():
                    q=read(b.DATA/'editorial_visual_sources/S1D_Visium_library_detection.tsv',sep='\t')
                    assert len(q)==16 and (q.ZEB2>=.1).sum()==2
                    return q,['ZEB1','ZEB2','CD1A','CD34','LYL1']
                module.detection_matrix=detection
            if TARGET in ['S6','S10']:
                if TARGET=='S6':
                    original_heading=module.heading
                    def heading(a,l,t):register(a,l,t);return original_heading(a,l,t)
                    module.heading=heading
                original_save=matplotlib.figure.Figure.savefig
                # Capture needs real savefig, so restore just for our export.
                def intercept(f,path,*args,**kwargs):
                    if str(path).endswith('Figure'+TARGET+'.pdf'):
                        matplotlib.figure.Figure.savefig=original_save
                        _render(f,TARGET)
                        matplotlib.figure.Figure.savefig=intercept
                matplotlib.figure.Figure.savefig=intercept
                try:
                    if TARGET=='S6':module.figS6()
                    else:module.sc10()
                finally:matplotlib.figure.Figure.savefig=original_save
            else:
                def save(f,n):
                    style.polish_text(f);_render(f,'S'+str(n))
                module.save=save;module._save=save
                getattr(module,'sc'+TARGET[1:])()
            style.panel=original_panel
        if not CAPTURED:raise RuntimeError('No figure captured')
    finally:pd.read_csv=original_read

if __name__=='__main__':
    # Make imports see this running script as the capture helper.
    sys.modules['render_figure']=sys.modules[__name__]
    p=argparse.ArgumentParser();p.add_argument('figure');p.add_argument('--panel');p.add_argument('--journal',default='canonical',choices=['canonical','Blood Advances','Scientific Reports']);p.add_argument('--out')
    a=p.parse_args();run(a.figure,a.panel,a.journal,a.out)
