"""Typeset the Revision 2 within-cohort residual table."""
from pathlib import Path
import csv
from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape, A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph

base = Path(__file__).resolve().parents[1]
src = base / 'data/source_data_rebuilt/supplementary_tables/Table_S4_full_subtype_residuals.tsv'
dst = base / 'manuscript/Supplementary_Table_S4.pdf'
with src.open(encoding='utf-8-sig', newline='') as fh:
    rows = list(csv.DictReader(fh, delimiter='\t'))

def f(x): return f'{float(x):.2f}'
def p(x):
    x=float(x)
    return f'{x:.1e}' if x < .001 else f'{x:.3f}'

data = [["Subtype", "n", "Cohort residual", "Bootstrap 95% CI",
         "Adjusted effect", "Bootstrap 95% CI", "BH FDR", "Normal sensitivity",
         "Outside range"]]
for r in rows:
    label=(r['subtype'].replace('γδ','gamma-delta').replace('αβ','alpha-beta')
           .replace('<U+03B3><U+03B4>','gamma-delta')
           .replace('<U+03B1><U+03B2>','alpha-beta'))
    data.append([label, r['n'], f(r['residual_within_median']),
      f(r['residual_within_ci_low'])+' to '+f(r['residual_within_ci_high']),
      f(r['adjusted_effect']),
      f(r['adjusted_effect_ci_low'])+' to '+f(r['adjusted_effect_ci_high']),
      p(r['adjusted_effect_bh_fdr']),
      f(r['normal_reference_median_sensitivity']),
      r['n_outside_normal_dev_range']])

doc=SimpleDocTemplate(str(dst),pagesize=landscape(A4),leftMargin=28,
                      rightMargin=28,topMargin=26,bottomMargin=25)
title=ParagraphStyle('title',fontName='Helvetica-Bold',fontSize=13,
                     leading=16,spaceAfter=5)
note=ParagraphStyle('note',fontName='Helvetica',fontSize=7.3,
                    leading=10,spaceBefore=8)
table=Table(data,colWidths=[115,24,66,93,68,93,57,82,64],
            rowHeights=[25]+[22]*len(rows),repeatRows=1)
table.setStyle(TableStyle([
  ('BACKGROUND',(0,0),(-1,0),colors.HexColor('#21384b')),
  ('TEXTCOLOR',(0,0),(-1,0),colors.white),
  ('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
  ('FONTNAME',(0,1),(-1,-1),'Helvetica'),
  ('FONTSIZE',(0,0),(-1,-1),8),
  ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
  ('ALIGN',(1,0),(-1,-1),'CENTER'),
  ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f1f6f8')]),
  ('LINEBELOW',(0,0),(-1,0),.7,colors.HexColor('#21384b')),
  ('BOTTOMPADDING',(0,0),(-1,0),7),
  ('TOPPADDING',(0,0),(-1,0),7)]))
text=('Cohort residuals are medians of observed minus within-cohort spline-expected '
      'balance. Adjusted effects come from the joint developmental-spline plus '
      'subtype model, using sum-to-zero contrasts relative to the equal-weighted '
      'average of subtype intercepts. Both intervals are percentile 95% intervals '
      'from 3,000 within-subtype patient bootstrap replicates with model refitting. '
      'BH FDR corrects HC3 robust coefficient P values across 17 subtypes. '
      'The normal-stage projection is a separately standardized cross-cohort '
      'sensitivity; its magnitude is not an absolute biological offset. Outside '
      'range counts patients beyond the observed normal-stage coordinate domain. '
      'MLLT10 has a negative residual median but no independent adjusted subtype '
      'effect after BH correction. Full-precision values are in the companion TSV.')
doc.build([Paragraph('Supplementary Table S4. Within-cohort developmental '
                     'adjustment across 17 T-ALL subtypes',title),table,
           Paragraph(text,note)])
print(dst)
