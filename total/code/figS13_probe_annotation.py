"""
Figure S13 | Probe-annotation validation (imported from the companion GPL570
oligonucleotide audit, D:/ZEB1).

Why it matters for this paper: a widely reused Affymetrix probeset, 208078_s_at,
is deposited in GEO with a ZEB1 representative sequence (NM_030751) but its 11
physical 25-mers match SIK1 (chr21), not ZEB1 (chr10). In the inspected GSE13159
cohort its signal is lower in T-ALL than in several comparison classes. This does
not establish which probeset any earlier mechanistic study used. Our study uses
the four ZEB1-specific probesets
(210875_s_at, 212758_s_at, 212764_at, 239952_at), which map only to ZEB1.

Panels
  A  Four-layer identity of the candidate probesets (annotation trap)
  B  Individual 208078_s_at oligonucleotide hits
  C  Sample distributions of the non-ZEB1 probe signal
  D  Complete correlation matrix of the five candidate probes
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import _style as S

S.set_style()
C = S.C

AUDIT_DIR = os.path.join(S.DATA_DIR, "probe_audit")
psum = pd.read_csv(os.path.join(AUDIT_DIR, "gpl570_zeb1_probe_summary.tsv"), sep="\t")
byclass = pd.read_csv(os.path.join(AUDIT_DIR, "gse13159_208078_by_geo_class.tsv"), sep="\t")
pw = pd.read_csv(os.path.join(AUDIT_DIR, "gse26713_probe_pairwise_spearman.tsv"), sep="\t")
oligo = pd.read_csv(os.path.join(AUDIT_DIR, "probe_25mer_transcriptome_hits.tsv"), sep="\t")
platform_audit = pd.read_csv(os.path.join(AUDIT_DIR, "gpl570_probe_audit.tsv"), sep="\t")
trap = oligo[oligo["probe"] == "208078_s_at"]
assert len(trap) == 11 and trap["hits_SIK1"].sum() == 11 and trap["hits_ZEB1"].sum() == 0

# ---- write reproducible source-data copies --------------------------------
psum.to_csv(os.path.join(S.DATA_DIR, "SourceData_FigS13_probe_identity.tsv"), sep="\t", index=False)
cls = (byclass.groupby("leukemia_class")["probe_208078_s_at"]
       .agg(["mean", "median", "count"]).reset_index()
       .sort_values("median"))
cls.to_csv(os.path.join(S.DATA_DIR, "SourceData_FigS13_208078_by_class.tsv"), sep="\t", index=False)

audit = (platform_audit["audit_class"].value_counts()
         .rename_axis("class").reset_index(name="n"))
audit.to_csv(os.path.join(S.DATA_DIR, "SourceData_FigS13_gpl570_audit.tsv"), sep="\t", index=False)

fig = S.new_fig(183, 150)
axd = fig.subplot_mosaic([["A", "A", "B"], ["C", "C", "D"]],
                         gridspec_kw=dict(width_ratios=[1, 1, 1.15]))

# ---- A: identity table -----------------------------------------------------
axA = axd["A"]
axA.axis("off")
axA.set_title("Four-layer identity of candidate ZEB1 probesets (GPL570)", loc="left", fontsize=7.5)
probes = ["208078_s_at", "210875_s_at", "212758_s_at", "212764_at", "239952_at", "232470_at"]
info = psum.set_index("ID")
rows = []
for p in probes:
    r = info.loc[p]
    rep = "ZEB1" if "ZEB1" in str(r["Nucleotide Title"]) else ("—" if pd.isna(r["Nucleotide Title"]) else "other")
    consensus = platform_audit.loc[platform_audit["probeset"] == p, "consensus_gene"]
    rows.append([p, r["Gene symbol"], str(r["Chromosome location"]), rep,
                 str(consensus.iloc[0]) if len(consensus) else "not audited"])
cols_h = ["Probeset", "Current\nsymbol", "Chr", "GEO rep.\nseq.", "Oligo\nconsensus"]
ncol = len(cols_h); nrow = len(rows)
x0, y0, dw, dh = 0.0, 0.86, 1.0 / ncol, 0.135
for j, h in enumerate(cols_h):
    axA.text(x0 + (j + 0.5) * dw, y0 + 0.09, h, ha="center", va="center",
             fontsize=5.6, fontweight="bold", transform=axA.transAxes)
for i, row in enumerate(rows):
    yy = y0 - i * dh
    is_zeb1_probe = row[0] in {"210875_s_at", "212758_s_at", "212764_at", "239952_at"}
    is_trap = row[0] == "208078_s_at"
    band = "#FDECEC" if is_trap else ("#EAF3EC" if is_zeb1_probe else "white")
    axA.add_patch(Rectangle((0, yy - dh / 2), 1, dh, transform=axA.transAxes,
                            fc=band, ec="none", zorder=0))
    for j, val in enumerate(row):
        col = "#222"
        if j == 4:
            col = C["tall"] if val == "SIK1" else (C["teal"] if val == "ZEB1" else "#222")
        if j == 3 and val == "ZEB1" and is_trap:
            col = C["tall"]
        axA.text(x0 + (j + 0.5) * dw, yy, str(val), ha="center", va="center",
                 fontsize=5.5, color=col, transform=axA.transAxes,
                 fontweight="bold" if j in (1, 4) else "normal")
axA.text(0.0, y0 - nrow * dh + 0.02,
         "208078_s_at: 11/11 oligos hit SIK1 + an overlapping model;\n"
         "0/11 hit ZEB1. SIK1 gene-level specificity is unresolved.\n"
         "Green: the four ZEB1 probesets used in this study.",
         fontsize=5.3, va="top", transform=axA.transAxes)
S.plabel(axA, "A")

# ---- B: inspect every physical oligonucleotide -----------------------------
axB = axd["B"]
targets = ["ZEB1", "SIK1", "ENSG00000275993"]
trap = trap.sort_values("oligo")
hits = np.array([[int(g in str(v).split(",")) for g in targets] for v in trap.genes])
from matplotlib.colors import ListedColormap
axB.imshow(hits, cmap=ListedColormap(["#F1F3F5", C["purple"]]), vmin=0, vmax=1, aspect="auto")
axB.set_yticks(range(len(trap)), [f"25-mer {v:02d}" for v in trap.oligo], fontsize=5.8)
axB.set_xticks(range(3), ["ZEB1", "SIK1", "Overlapping\nENSG model"], fontsize=6)
for i in range(len(trap)):
    for j in range(3):
        axB.text(j, i, "hit" if hits[i, j] else "—", ha="center", va="center",
                 fontsize=5.4, color="white" if hits[i, j] else "#777")
axB.set_title("208078_s_at: exact 25-mer matches", loc="left", fontsize=7.3)
axB.set_xlabel("GENCODE v47; both strands", fontsize=6)
S.plabel(axB, "B")
pd.DataFrame(hits, index=trap.oligo, columns=targets).to_csv(
    os.path.join(S.DATA_DIR, "SourceData_FigS13_oligo_hits.tsv"), sep="\t")

# ---- C: 208078 by leukemia class ------------------------------------------
axC = axd["C"]
keep = ["T-ALL", "B-ALL", "AML", "CML", "CLL", "MDS", "healthy BM", "Healthy BM", "Normal BM"]
cc = cls[cls["leukemia_class"].isin(keep) | cls["leukemia_class"].str.contains("BM|healthy|Normal", case=False, na=False)]
cc = cls.copy()
cc = cc[cc["count"] >= 30].sort_values("median")
colsC = [C["tall"] if "T-ALL" in str(s) else C["grey"] for s in cc["leukemia_class"]]
S.add_grid(axC, "x")
yy = np.arange(len(cc))
distributions = [byclass.loc[byclass.leukemia_class.eq(g), "probe_208078_s_at"].dropna().to_numpy()
                 for g in cc.leukemia_class]
bp = axC.boxplot(distributions, positions=yy, vert=False, widths=.58,
                 patch_artist=True, showfliers=False,
                 medianprops=dict(color="#222", linewidth=.8),
                 whiskerprops=dict(color="#777", linewidth=.5),
                 capprops=dict(color="#777", linewidth=.5))
for patch, col in zip(bp["boxes"], colsC):
    patch.set(facecolor=col, edgecolor=col, alpha=.7)
axC.set_yticks(yy); axC.set_yticklabels(
    [f"{g} (n={n})" for g, n in zip(cc.leukemia_class, cc["count"])], fontsize=5.6)
axC.set_xlabel("208078_s_at signal, GSE13159")
axC.set_title("208078_s_at signal by leukemia class", loc="left", fontsize=7.3)
S.style_ax(axC); S.plabel(axC, "C")

# ---- D: probe concordance --------------------------------------------------
axD = axd["D"]
order = ["208078_s_at", "210875_s_at", "212758_s_at", "212764_at", "239952_at"]
matrix = pw.pivot(index="probe_a", columns="probe_b", values="rho").reindex(index=order, columns=order)
matrix = matrix.combine_first(matrix.T)
assert matrix.notna().all().all()
im = axD.imshow(matrix, cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
labels = [v.split("_")[0] for v in order]
axD.set_xticks(range(5), labels, rotation=45, ha="right", fontsize=5.6)
axD.set_yticks(range(5), labels, fontsize=5.6)
for i in range(5):
    for j in range(5):
        value = matrix.iloc[i, j]
        axD.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=5.4,
                 color="white" if abs(value) > .65 else "#222")
fig.colorbar(im, ax=axD, fraction=.045, pad=.03, label="Spearman ρ")
axD.set_title("Probe concordance: GSE26713 (n=117)", loc="left", fontsize=7.3)
S.plabel(axD, "D")

supp = os.path.join(S.FIG_DIR, "Supplementary")
os.makedirs(os.path.join(supp, "png"), exist_ok=True)
fig.savefig(os.path.join(supp, "FigureS13_probe_annotation.pdf"))
fig.savefig(os.path.join(supp, "png", "FigureS13_probe_annotation.png"), dpi=400)
plt.close(fig)
print("saved Supplementary/FigureS13_probe_annotation.pdf + png")
