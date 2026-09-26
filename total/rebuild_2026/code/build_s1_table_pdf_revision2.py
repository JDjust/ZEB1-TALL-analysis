"""Typeset clinical attenuation with endpoint multiplicity shown explicitly."""
from pathlib import Path
import csv
from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape, A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph

base=Path(__file__).resolve().parents[1]
src=base/'data/source_data_rebuilt/supplementary_tables/Table_S1_clinical_models.tsv'
dst=base/'manuscript/Supplementary_Table_S1.pdf'
with src.open(encoding='utf-8-sig',newline='') as fh:
    rows=list(csv.DictReader(fh,delimiter='\t'))
def pv(x):
    x=float(x)
    return f'{x:.1e}' if x<.001 else f'{x:.3f}'
data=[['Endpoint','Measure','Adjustment','n','Events','Estimate','95% CI','P','BH FDR']]
for r in rows:
    data.append([r['endpoint'],r['measure'],r['adjustment'],r['n'],r['events'],
       f"{float(r['estimate']):.2f}",
       f"{float(r['ci_low']):.2f} to {float(r['ci_high']):.2f}",
       pv(r['p']),pv(r['endpoint_bh_fdr'])])
doc=SimpleDocTemplate(str(dst),pagesize=landscape(A4),leftMargin=28,
                      rightMargin=28,topMargin=30,bottomMargin=25)
title=ParagraphStyle('title',fontName='Helvetica-Bold',fontSize=13,leading=16,
                     spaceAfter=8)
note=ParagraphStyle('note',fontName='Helvetica',fontSize=7.6,leading=10,
                    spaceBefore=9)
table=Table(data,colWidths=[136,42,106,38,47,60,100,55,58],
            rowHeights=[28]+[28]*len(rows),repeatRows=1)
table.setStyle(TableStyle([
 ('BACKGROUND',(0,0),(-1,0),colors.HexColor('#21384b')),
 ('TEXTCOLOR',(0,0),(-1,0),colors.white),
 ('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
 ('FONTNAME',(0,1),(-1,-1),'Helvetica'),
 ('FONTSIZE',(0,0),(-1,-1),8.7),
 ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
 ('ALIGN',(1,0),(-1,-1),'CENTER'),
 ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f1f6f8')]),
 ('LINEBELOW',(0,0),(-1,0),.7,colors.HexColor('#21384b'))]))
note_text=('Frozen odds ratios (OR) or hazard ratios (HR) per unit higher balance. '
  'BH FDR is calculated across six endpoints separately within each adjustment '
  'level. Induction failure and MRD >=0.1% attenuate after molecular-subtype '
  'adjustment. M2/M3 morphology retains an adjusted association after endpoint '
  'correction; MRD >=0.01% does not. These exploratory models do not establish '
  'independent clinical utility. Full model identifiers and exact estimates '
  'are in the companion TSV.')
doc.build([Paragraph('Supplementary Table S1. Clinical associations and subtype '
                     'attenuation',title),table,Paragraph(note_text,note)])
print(dst)
