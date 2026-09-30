"""Render frozen Figure 3 data on the approved 180 x 200 mm six-panel layout."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

import final_visual_style_lock as style
import editorial_rebuild_figures as base
style.apply_style()
style.lock_legacy_module(base)

D, DATA, SD, OUT = base.D, base.DATA, base.SD, base.OUT
Z1, Z2, ETP = style.ZEB1, style.ZEB2, style.ETP
GREY, LIGHT, DARK = style.GREY, style.REF, style.INK
CMAP = style.SIGNED_CMAP
AXES = {
    'A': [.08, .72, .35, .20],
    'B': [.50, .72, .43, .20],
    'C': [.08, .37, .85, .27],
    'D': [.08, .065, .20, .19],
    'E': [.34, .065, .31, .19],
    'F': [.70, .065, .23, .19],
}


def read(path):
    return pd.read_csv(path, sep='\t', low_memory=False)


def panel(ax, letter, title):
    x0, y0, width, height = ax.get_position().bounds
    title_y = y0 + height + (.004 if letter in 'DEF' else .010)
    ax.figure.text(x0 - .004, title_y, letter, fontsize=10,
                   weight='bold', color=DARK, va='bottom')
    ax.figure.text(x0 + .027, title_y, title, fontsize=7.8,
                   weight='semibold', color=DARK, va='bottom')
    ax.tick_params(labelsize=6.7, length=2, width=.45, pad=2)
    ax.xaxis.label.set_size(7.2)
    ax.yaxis.label.set_size(7.2)
    ax.spines['left'].set_color(style.REF)
    ax.spines['bottom'].set_color(style.REF)


def panel_a(ax):
    summary = read(D / 'gate2_program/freeze_summary.tsv').iloc[0]
    assert (summary.n_eligible, summary.n_pass, summary.n_pass_ZEB1_side,
            summary.n_pass_ZEB2_side, summary.n_frozen_ZEB1,
            summary.n_frozen_ZEB2) == (24356, 3560, 499, 3061, 50, 50)
    passed = read(D / 'gate2_program/effect_threshold_pass.tsv')
    assert len(passed) == 3560
    for sign, color, count in [(1, Z1, 499), (-1, Z2, 3061)]:
        q = passed[passed.beta_R * sign > 0].sort_values('partial_R2', ascending=False)
        assert len(q) == count
        ranks = np.arange(1, count + 1)
        ax.plot(ranks, q.partial_R2, color=color, alpha=.30, lw=.85)
        ax.plot(ranks[:50], q.partial_R2.iloc[:50], color=color, lw=1.7)
        ax.scatter([50], [q.partial_R2.iloc[49]], color=color, s=11, zorder=3)
    ax.axvline(50, color=GREY, ls=':', lw=.65)
    ax.set(xscale='log', xlim=(.85, 3400), xlabel='Partial R² rank', ylabel='Partial R²')
    ax.set_xticks([1, 10, 50, 500, 3000], ['1', '10', '50', '500', '3k'])
    ax.text(.98, .96, '24,356 eligible → 3,560 pass\n499 / 3,061 by side',
            transform=ax.transAxes, ha='right', va='top', fontsize=6.5, color=GREY)
    base.clean_grid(ax, 'y')
    panel(ax, 'A', 'Ranked residual effects')


def panel_b(ax):
    eligible = read(D / 'gate2_program/all_genes_residual_effects.tsv')
    frozen = read(D / 'gate2_program/frozen_ZEB_side50.tsv')
    assert len(eligible) == 24361 and len(frozen) == 100
    excluded = {'ZEB1', 'ZEB2', 'CD34', 'LYL1', 'CD1A'}
    eligible = eligible[~eligible.symbol.str.upper().isin(excluded)]
    assert len(eligible) == 24356
    ax.scatter(eligible.beta_R, eligible.partial_R2, s=1.6,
               color=style.BG_POINT, alpha=.5, rasterized=True)
    for side, color in [('ZEB1-side50', Z1), ('ZEB2-side50', Z2)]:
        q = frozen[frozen.side == side]
        assert len(q) == 50
        ax.scatter(q.beta_R, q.partial_R2, s=9, color=color,
                   alpha=.85, label=side, rasterized=True)
    ax.axvline(0, color=LIGHT, lw=.65)
    ax.set(xlabel='Leukemia residual coefficient βR', ylabel='Partial R²')
    ax.legend(frameon=False, fontsize=6.5, loc='upper right',
              handletextpad=.35, labelspacing=.2)
    base.clean_grid(ax)
    panel(ax, 'B', 'Genome-wide residual effects')


def panel_c(ax, fig):
    x = read(DATA / 'editorial_visual_sources/F3C_frozen100_residual_deciles.tsv')
    ann = read(DATA / 'editorial_visual_sources/F3C_residual_decile_annotations.tsv')
    assert len(x) == 100 and len(ann) == 10 and ann.n_patients.sum() == 1309
    x['order_side'] = np.where(x.side.eq('ZEB2-side50'), 0, 1)
    x = x.sort_values(['order_side', 'rank']).reset_index(drop=True)
    assert x.side.iloc[:50].eq('ZEB2-side50').all()
    assert x.side.iloc[50:].eq('ZEB1-side50').all()
    z = x[[f'decile_{i}' for i in range(1, 11)]].to_numpy(float)
    im = ax.imshow(z, aspect='auto', interpolation='nearest', origin='upper',
                   extent=(0, 10, 100, 0), cmap=CMAP, vmin=-1.25,
                   vmax=1.25, rasterized=True)
    ax.set_xlim(-3.15, 10.15)
    ax.set_ylim(101, -18.2)
    ax.axhline(50, color='white', lw=1.3)
    ax.set_xticks(np.arange(.5, 10.5), [f'D{i}' for i in range(1, 11)])
    ax.set_yticks([25, 75], ['ZEB2 · 50', 'ZEB1 · 50'])
    ax.tick_params(axis='y', length=0)
    for i, row in ann.sort_values('decile').reset_index(drop=True).iterrows():
        left, total = float(i), 0.
        for pct, color in [(row.pct_BCL11B, Z2),
                           (row.pct_ETP_like, ETP), (row.pct_TLX3, Z1)]:
            frac = min(max(float(pct) / 100, 0), 1 - total)
            if frac > 0:
                ax.add_patch(Rectangle((left, -16), frac, 3,
                                       facecolor=color, edgecolor='none'))
                left += frac
                total += frac
        if total < 1:
            ax.add_patch(Rectangle((left, -16), 1 - total, 3,
                                   facecolor=style.NA, edgecolor='none'))
        for y, value in [(-11, row.median_R), (-6, row.median_D)]:
            color = CMAP((np.clip(value, -2, 2) + 2) / 4)
            ax.add_patch(Rectangle((i, y), 1, 3, facecolor=color,
                                   edgecolor='white', lw=.25))
    for y, label in [(-14.5, 'Subtype'), (-9.5, 'R'), (-4.5, 'D')]:
        ax.text(-1.66, y, label, va='center', fontsize=6.5, color=DARK)
    for y, label, color in [(25, 'ZEB2', Z2), (75, 'ZEB1', Z1)]:
        ax.text(10.05, y, label, rotation=90, ha='center', va='center',
                fontsize=6.5, color=color)
    a6 = read(D / 'a6_yayon/a6_2_gene_effects.tsv').set_index('symbol')['class']
    a7 = read(D / 'a7c_decompose/a7c_zeb2_split.tsv').set_index('symbol')['subset']
    labels = set(x[x.side.eq('ZEB1-side50')].nsmallest(5,'rank').symbol)
    labels |= set(x[x.side.eq('ZEB2-side50')].nsmallest(5,'rank').symbol)
    labels |= {'AOAH','MAP3K5','ADRB2'} & set(x.symbol)
    requested = []
    for i, row in x.iterrows():
        for xpos, color in [(-.20, Z2 if i < 50 else Z1),
                            (-.43, {'concordant':Z1,'discordant':Z2,
                                    'neutral':GREY}.get(a6.get(row.symbol), style.NA)),
                            (-.66, {'myeloid_concordant':ETP,
                                    'non_myeloid_concordant':DARK}.get(a7.get(row.symbol),style.NA))]:
            ax.add_patch(Rectangle((xpos,i),.16,1,facecolor=color,edgecolor='none'))
        if row.symbol in labels: requested.append((i,row.symbol))
    # Deterministic label spacing changes text positions only; gene rows stay fixed.
    placed=[]
    for i, symbol in requested:
        y=max(float(i)+.5, placed[-1][0]+4.6) if placed else float(i)+.5
        placed.append([y,i,symbol])
    if placed and placed[-1][0]>98:
        shift=placed[-1][0]-98
        for item in placed: item[0]-=shift
    for y,i,symbol in placed:
        ax.plot([-.76,-.92,-1.05],[i+.5,i+.5,y],color=GREY,lw=.35)
        ax.text(-1.12,y,symbol,ha='right',va='center',fontsize=6.5,color=DARK)
    ax.text(-.18,-.9,'Side  A6  A7',ha='center',fontsize=6.5,color=DARK)
    ax.set_xlabel('Patient residual decile (low → high)')
    cb = fig.add_axes([.77, .336, .14, .008])
    fig.colorbar(im, cax=cb, orientation='horizontal', ticks=[-1, 0, 1])
    cb.set_xlabel('Median expression z', fontsize=6.5, labelpad=1)
    cb.tick_params(labelsize=6.5, length=1, pad=1)
    panel(ax, 'C', 'Frozen 100-gene expression across residual deciles')


def panel_d(ax):
    normal = read(D / 'a6_yayon/a6_1_donor_program.tsv')
    assert len(normal) == 5
    for side, color in [('ZEB1', Z1), ('ZEB2', Z2)]:
        keys = [f'{stage}_{side}_side50' for stage in ('early', 'DP', 'SP')]
        for _, row in normal.iterrows():
            ax.plot(range(3), [row[k] for k in keys], color=color,
                    alpha=.28, lw=.65)
        ax.plot(range(3), [normal[k].median() for k in keys],
                color=color, lw=1.8, marker='o', ms=2.5)
    ax.set_xticks(range(3), ['Early', 'DP', 'SP'])
    ax.set_xlim(-.1, 2.1)
    ax.set_ylabel('Program score')
    base.clean_grid(ax)
    panel(ax, 'D', 'Normal donors (n=5)')


def panel_e(ax):
    q = read(D / 'a6_yayon/a6_2_gene_effects.tsv')
    assert len(q) == 90
    assert q['class'].value_counts().to_dict() == {
        'concordant': 13, 'discordant': 32, 'neutral': 45}
    for cls, marker, alpha in [('concordant', 'o', .9),
                                ('discordant', 'x', .9), ('neutral', '.', .55)]:
        t = q[q['class'] == cls]
        ax.scatter(t.delta_normal, t.beta_R,
                   s=20 if marker != '.' else 10,
                   color=[Z1 if s.startswith('ZEB1') else Z2 for s in t.side],
                   marker=marker, alpha=alpha,
                   linewidths=.7 if marker == 'x' else .2, rasterized=True)
    ax.axhline(0, color=LIGHT, lw=.65)
    ax.axvline(0, color=LIGHT, lw=.65)
    ax.set(xlabel='Normal thymus early → DP effect', ylabel='Leukemia βR')
    ax.text(.98, .96, '13 concordant\n32 discordant\n45 weak\nrho = -0.50',
            transform=ax.transAxes, ha='right', va='top', fontsize=6.5,
            color=DARK, linespacing=1.2)
    base.clean_grid(ax)
    panel(ax, 'E', 'Normal versus leukemia association')
    return q


def panel_f(ax, q):
    z = q[q.side.eq('ZEB2-side50')].copy()
    matched = read(SD / 'revision4_targeted/matched_BCL11B_ETP_gene_results.tsv')
    z = z.merge(matched[['symbol', 'logFC']], on='symbol', validate='one_to_one')
    assert len(z) == 48
    assert (z['class'].eq('concordant').sum(),
            z['class'].isin(['neutral', 'discordant']).sum()) == (10, 38)
    z = z.sort_values('logFC').reset_index(drop=True)
    y = np.arange(len(z))
    for cls, color in [('concordant', GREY), ('neutral', Z2), ('discordant', Z2)]:
        mask = z['class'].eq(cls)
        ax.hlines(y[mask], 0, z.loc[mask, 'logFC'], color=color,
                  lw=.55, alpha=.7)
        ax.scatter(z.loc[mask, 'logFC'], y[mask], s=5.5,
                   color=color, zorder=2, rasterized=True)
    strip_colors = {'concordant': GREY, 'discordant': Z2, 'neutral': style.BG_POINT}
    for i, cls in enumerate(z['class']):
        ax.add_patch(Rectangle((-.86, i - .4), .13, .8,
                               facecolor=strip_colors[cls], edgecolor='none'))
    ax.axvline(0, color=LIGHT, lw=.65)
    ax.set_xlim(-1.05, 4.2)
    ax.set_ylim(-1, len(z))
    ax.set_yticks([])
    ax.set_xlabel('BCL11B − matched ETP\nlog2 effect')
    base.clean_grid(ax, 'x')
    panel(ax, 'F', 'Matched gene effects')


def main():
    fig = base.mmfig(180, 200)
    axs = {letter: fig.add_axes(bounds) for letter, bounds in AXES.items()}
    panel_a(axs['A'])
    panel_b(axs['B'])
    panel_c(axs['C'], fig)
    # Three fixed annotation strips at C's left margin: a compact key in the
    # existing C-to-D/E/F whitespace, outside all heatmap cells and gene names.
    def swatch(x, y, color):
        fig.patches.append(Rectangle((x,y-.004),.010,.008,
                                     transform=fig.transFigure,
                                     facecolor=color,edgecolor='none'))
    fig.text(.08,.322,'Side',fontsize=6.5,color=DARK,va='center')
    for x,label,color in [(.125,'ZEB2',Z2),(.205,'ZEB1',Z1)]:
        swatch(x,.322,color)
        fig.text(x+.013,.322,label,fontsize=6.5,color=DARK,va='center')
    fig.text(.295,.322,'A6',fontsize=6.5,color=DARK,va='center')
    for x,label,color in [(.337,'concordant',Z1),(.483,'discordant',Z2),
                          (.634,'weak',GREY)]:
        swatch(x,.322,color)
        fig.text(x+.013,.322,label,fontsize=6.5,color=DARK,va='center')
    fig.text(.08,.304,'A7',fontsize=6.5,color=DARK,va='center')
    for x,label,color in [(.125,'myeloid-concordant',ETP),
                          (.365,'non-myeloid-concordant',DARK),
                          (.672,'unavailable',style.NA)]:
        swatch(x,.304,color)
        fig.text(x+.013,.304,label,fontsize=6.5,color=DARK,va='center')
    panel_d(axs['D'])
    q = panel_e(axs['E'])
    panel_f(axs['F'], q)
    assert all(np.allclose(axs[k].get_position().bounds, bounds)
               for k, bounds in AXES.items())
    assert axs['D'].get_position().x1 < axs['E'].get_position().x0
    assert axs['D'].get_position().y0 == axs['E'].get_position().y0 == axs['F'].get_position().y0
    assert axs['D'].get_position().y1 == axs['E'].get_position().y1 == axs['F'].get_position().y1
    style.polish_text(fig)
    fig.savefig(OUT / 'Figure3.pdf', bbox_inches=None)
    fig.savefig(OUT / 'Figure3.png', dpi=300, bbox_inches=None)
    plt.close(fig)
    print(OUT / 'Figure3.pdf')
    print(OUT / 'Figure3.png')


if __name__ == '__main__':
    main()
