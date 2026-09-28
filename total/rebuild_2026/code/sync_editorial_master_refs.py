"""Correct panel references after the frozen visual redesign."""
from pathlib import Path

p=Path(__file__).resolve().parents[3]/'submission/_build/editorial_rebuild/Editorial_Master.md'
s=p.read_text(encoding='utf-8')
repls={
 'eight prespecified early, DP and SP poles are shown in developmental order in Figure 1C, with the complete state audit in Supplementary Figure S1':
 'all 33 eligible author states are displayed in Figure 1C, with the five-donor comparison as an inset and an aligned state audit in Supplementary Figure S1',
 'both frozen program scores were tracked from early through DP to SP states (Fig. 3D)':
 'both frozen program scores were tracked from early through DP to SP states (Fig. 3D, inset within 3E)',
 '(Fig. 5A–D; Supplementary Fig. S6)':
 '(Fig. 5A,B; Supplementary Fig. S6)',
 'effect correlations 0.60, 0.64, 0.77, 0.62 and 0.86 (Fig. 5E)':
 'effect correlations 0.60, 0.64, 0.77, 0.62 and 0.86 (Fig. 5C)',
 '(ρ = −0.21, interval including zero; Fig. 5F; Supplementary Fig. S8)':
 '(ρ = −0.21, interval including zero; Fig. 5D,E; Supplementary Fig. S8)',
 'The direction remained under both adult age restrictions (Fig. 7E).':
 'The score directions remained under both adult age restrictions (Fig. 7C inset).',
 'within the adult cohort (Fig. 7F; Supplementary Fig. S9)':
 'within the adult cohort (Fig. 7E; Supplementary Fig. S9)',
 'partial portability (Fig. 7G)':
 'partial portability (Fig. 7F)',
 'and pooled loops and seven scATAC peak sets did not reproduce the polarity (Supplementary Fig. S10B–D)':
 'and pooled loops and seven scATAC peak sets did not reproduce the polarity (Supplementary Fig. S10B,C; Supplementary Table S10)',
}
for old,new in repls.items():
    if old not in s:
        if new not in s: raise AssertionError(old)
        continue
    s=s.replace(old,new)
s=s.replace('However, the preservation of directional structure in blast-enriched T-ALL and in author-defined malignant-cell pseudobulks argues against passive myeloid admixture as its sole explanation (Fig. 6E). A 45-gene normal-myeloid-concordant slice and three non-myeloid-concordant genes (AOAH, MAP3K5 and ADRB2) are an audit decomposition, not a new three-gene program (Supplementary Fig. S7).',
'''A measured 48-gene by four-cohort heatmap shows that the 45 normal-myeloid-concordant genes and three non-myeloid-concordant genes (AOAH, MAP3K5 and ADRB2) both retain a leukemia-side directional component (Fig. 6E). Together with blast-enriched and author-labelled malignant-cell observations, this argues against passive myeloid admixture as the sole explanation without proving a fully cell-intrinsic mechanism. The three genes are an audit decomposition, not a new three-gene program (Supplementary Fig. S7).''')
p.write_text(s,encoding='utf-8')
