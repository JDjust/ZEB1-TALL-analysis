"""Build four journal-specific submission packages from the frozen master."""
from pathlib import Path
import re, shutil, subprocess, tempfile
import fitz, openpyxl
from openpyxl.styles import Font, PatternFill, Alignment

ROOT=Path(__file__).resolve().parents[3]
SUB=ROOT/'submission'; OUT=SUB/'article'; OUT.mkdir(exist_ok=True)
SRC=ROOT/'total/rebuild_2026/manuscript'
DATA=ROOT/'total/rebuild_2026/data/source_data_rebuilt/revision4_targeted'
MASTER=(SRC/'ZEB1_ZEB2_TALL_manuscript_polished.md').read_text(encoding='utf-8')
SI=(SRC/'Supplementary_Information.md').read_text(encoding='utf-8')
PANDOC=shutil.which('pandoc'); CHROME=Path(r'C:/Program Files/Google/Chrome/Application/chrome.exe')
REF=SUB/'_build/reference.docx'
assert PANDOC and CHROME.exists() and REF.exists()

TITLES={
 'npj Systems Biology and Applications':'Developmental context distinguishes molecular subtype ZEB expression states in T cell acute lymphoblastic leukemia',
 'JCMM':'BCL11B-associated ZEB expression and transcriptomic programs in developmentally matched T-cell acute lymphoblastic leukemia',
 'British Journal of Haematology':'Subtype-associated ZEB1–ZEB2 expression distinguishes BCL11B from developmentally matched ETP-like T-ALL',
 'Scientific Reports':'Subtype-associated ZEB1–ZEB2 expression patterns in T-cell acute lymphoblastic leukemia after developmental-expression adjustment'}

ABSTRACTS={
 'npj Systems Biology and Applications':'''Normal thymopoiesis patterns ZEB1 and ZEB2 expression, complicating interpretation of their balance in T-cell acute lymphoblastic leukemia (T-ALL). We integrated four normal-thymus resources with diagnostic RNA from 1,309 pediatric and young-adult patients across 17 molecular subtypes. A ZEB-independent three-gene developmental-expression proxy explained 25.6% of balance variance; subtype added 13.9 percentage points. Cohort-internal residuals placed BCL11B at the ZEB2-skewed pole and TLX3 at the opposite pole. The polarity persisted with subtype-held-out fitting, alternative expression scores and spline degrees of freedom 2–4. Eighteen BCL11B cases differed from developmental-score-matched ETP-like cases in residual balance (mean paired difference −1.92; bootstrap 95% confidence interval −2.98 to −0.67). All 12 lesion-defined BCL11B-rearranged cases were ZEB2-dominant, including seven without a ZEB2 rearrangement partner. Independent bulk and patient-paired single-cell data supported narrower contrasts. These findings distinguish subtype-associated expression from a developmental proxy without establishing developmental identity, a causal mechanism or clinical utility.''',
 'JCMM':'''The relative expression of ZEB1 and ZEB2 varies during human thymopoiesis, but its relationship to molecular T-cell acute lymphoblastic leukemia (T-ALL) subtypes remains unclear. We integrated four normal-thymus resources with 1,309 pediatric and young-adult diagnostic cases across 17 subtypes. A ZEB-independent three-gene developmental-expression proxy explained 25.6% of balance variance, and subtype added 13.9 percentage points. BCL11B showed the lowest cohort-internal residual (median −1.76; 95% confidence interval −2.47 to −0.98), while TLX3 showed the highest (+0.64). Eighteen BCL11B cases differed from developmental-score-matched ETP-like cases (mean paired residual difference −1.92; bootstrap interval −2.98 to −0.67). In an exploratory pair-adjusted transcriptomic comparison, interferon-γ response and IL6–JAK–STAT3 Hallmark programs were higher in BCL11B cases (competitive-test FDR 0.0025 and 0.0028). All 12 external lesion-defined BCL11B-rearranged cases were ZEB2-dominant, including seven without a ZEB2 partner. Subtype-held-out and alternative-score analyses supported the major polarity; independent bulk and patient-paired pseudobulks supported narrower expression contrasts. These observations describe a developmentally contextualized, subtype-associated expression phenotype, without establishing direct regulation or clinical utility.''',
 'British Journal of Haematology':'''T-cell acute lymphoblastic leukaemia (T-ALL) subtypes may resemble different stages of thymopoiesis, complicating interpretation of ZEB1 and ZEB2 expression. We compared four normal-thymus resources with diagnostic RNA from 1,309 paediatric and young-adult patients across 17 molecular subtypes. A three-gene developmental-expression proxy explained 25.6% of the variance in standardized ZEB1-minus-ZEB2 balance; subtype added 13.9 percentage points. BCL11B had the lowest cohort-internal residual (median −1.76; 95% confidence interval −2.47 to −0.98), whereas TLX3 had the highest (+0.64). BCL11B cases retained lower balance than 18 developmental-score-matched ETP-like cases (mean paired residual difference −1.92; bootstrap interval −2.98 to −0.67). All 12 lesion-defined BCL11B-rearranged cases were ZEB2-dominant, including seven without a ZEB2 rearrangement partner. Subtype-held-out fitting, alternative proxy scores and spline sensitivity preserved the main polarity; independent bulk and patient-paired single-cell data supported narrower contrasts. The findings distinguish BCL11B-associated ZEB2 skewing from developmental immaturity alone, while leaving mechanism, complete external subtype replication and clinical utility unresolved.''',
 'Scientific Reports':'''ZEB1 and ZEB2 expression changes during thymopoiesis, yet molecular T-cell acute lymphoblastic leukemia (T-ALL) subtypes may have additional expression differences. We integrated four normal-thymus resources with diagnostic RNA from 1,309 pediatric and young-adult patients assigned to 17 subtypes. A three-gene developmental-expression proxy excluding ZEB genes explained 25.6% of standardized ZEB1-minus-ZEB2 balance variance; subtype added 13.9 percentage points. BCL11B had the lowest within-cohort median residual (−1.76; 95% confidence interval −2.47 to −0.98), and TLX3 the highest (+0.64). Eighteen BCL11B patients differed from developmental-score-matched ETP-like patients (mean paired residual difference −1.92; bootstrap interval −2.98 to −0.67). All 12 external lesion-defined BCL11B-rearranged cases were ZEB2-dominant. Subtype-held-out fitting, alternative expression scores, spline degrees of freedom 2–4 and cross-validated model comparison supported the main subtype polarity. Independent bulk and patient-paired single-cell analyses tested narrower contrasts. No external cohort reproduced the complete 17-subtype residual ranking. These findings describe subtype-associated expression conditional on a developmental proxy; causal mechanism and clinical utility remain unresolved.'''}

def wc(s):return len(re.findall(r'\b[\w]+(?:[–\'-][\w]+)*\b',s,re.U))

ADDED_METHODS='''
### Prespecified spline-degree and normal-stage sensitivity

The same 1,309 diagnostic patients and frozen three-gene score were used to refit the cohort-internal balance model with natural spline degrees of freedom 2, 3 and 4. For each fit, we recomputed BCL11B, ETP-like and TLX3 median residuals, the fixed 18-pair BCL11B-minus-ETP-like mean residual difference, and the incremental raw R² after adding molecular subtype. No degree of freedom was selected by significance. In normal thymus, the existing 20 measured reference units were summarized by annotated stage. Eight GSE195812 pooled stage libraries were treated as stage summaries. The 12 GSE206710 sample-stage aggregates were first averaged within each of three donors at DN-early, DP-immature and DP-mature stages. Donor-level values are descriptive; there are three donors, not 12 independent people. The proxy is not expected to be monotonic throughout DP development.

### Ridge Cox sensitivity for unstable survival fits

For event-free and overall survival among complete cases, Cox ridge regression penalized one-hot subtype coefficients only; balance, age, sex and log10 white-cell count were unpenalized. Five-fold cross-validation selected the minimum-deviance penalty separately per endpoint with seed 20260926; Breslow ties were specified. The fixed selected penalty was refitted in 500 patient-level bootstrap samples to obtain percentile intervals conditional on penalty selection. One empty-model bootstrap fit per endpoint was excluded (499/500 valid for each endpoint). The one-standard-error penalty and penalty path are retained as source data because the one-standard-error rule can almost erase subtype effects. These intervals exclude tuning uncertainty and are not evidence of calibrated prediction or independent clinical utility. The original unpenalized Efron-ties Cox output and its convergence warnings remain in Table S4.

### Exploratory matched transcriptomic programs

The 18 previously fixed one-to-one BCL11B–ETP-like developmental-score matches were used without rematching or changing the outcome. Source patient counts were filtered by expression, TMM normalized and analysed with limma–voom in a pair-blocked linear model. The coefficient is BCL11B minus ETP-like within matched pairs. Hallmark gene sets were obtained with msigdbr [51] and tested by limma camera [52], accounting for intergene correlation; 50 eligible sets received Benjamini–Hochberg correction. Gene-level significance in Figure S11A is used for display only (BH FDR <0.05 and absolute log₂ fold change >0.5). This one-cohort program analysis is exploratory and does not establish direct ZEB targets or regulation.
'''
ADDED_LEGENDS='''
## Supplementary Figure S10. Normal-proxy and spline-degree sensitivity

(**A**) Median cohort-internal residuals for BCL11B (18 patients), ETP-like (235) and TLX3 (212) after fitting natural splines with two, three or four degrees of freedom to the same 1,309 diagnostic patients. Each line joins one subtype under the three fixed specifications; zero is the cohort-fitted expectation. No confidence intervals are attached to these descriptive estimates. (**B**) Raw incremental R², in percentage points, when 17 molecular-subtype terms are added after each developmental spline. This is a within-specification comparison, not a causal variance partition; adjusted and held-out comparisons are in main Figure 7C. (**C**) Three-gene proxy values for DN-early, DP-immature and DP-mature donor-stage means in GSE206710. Each line is one donor; D1 averages repeated sample-stage aggregates. All three donors have higher proxy values at both DP stages than at DN-early, but DP-immature can exceed DP-mature, so the score is not a monotonic clock. GSE195812 contributes eight pooled stage libraries and is summarized in Table S9B. Exact values are in Table S9A–B.

## Supplementary Figure S11. Exploratory programs in matched BCL11B versus ETP-like patients

(**A**) Pair-adjusted limma–voom BCL11B-minus-ETP-like contrast from the fixed 18 matched pairs (36 distinct patients). Each point is one of 22,359 expression-filtered genes. Orange and blue denote genes with BH FDR <0.05 and absolute log₂ fold change >0.5 in positive and negative directions; grey genes do not meet both display thresholds. Vertical dashed lines mark ±0.5 log₂ fold change. The y axis uses nominal gene-level P values for visual resolution; color uses corrected values. (**B**) Ten Hallmark gene sets with the lowest BH FDR among 50 tested using pair-blocked limma camera. Position is −log10 BH FDR, dot size is tested gene count, and the dashed line is FDR 0.05. Up denotes higher aggregate expression in BCL11B. Correlated pathways can reflect broader subtype composition; no ZEB-driven mechanism or independent replication is inferred. Full results are in Table S9D.

## Supplementary Table S9. Targeted sensitivity results

S9A contains spline-degree results; S9B normal-stage proxy summaries; S9C ridge Cox estimates and penalty selection; S9D all 50 matched Hallmark tests. Full-precision source tables and executable scripts are in the versioned code archive.
'''

def main_variant(journal):
    md=re.sub(r'^# .*$', '# '+TITLES[journal], MASTER,count=1,flags=re.M)
    md=re.sub(r'(?s)(## Abstract\n\n).*?(\n\n---\n\n## Introduction)',
              lambda m:m.group(1)+ABSTRACTS[journal]+m.group(2),md,count=1)
    common=('The fixed spline-degree sensitivity gave BCL11B median residuals of −1.78, −1.76 and −1.79 for df = 2, 3 and 4; TLX3 remained positive (+0.63 to +0.64). The 18-pair BCL11B–ETP-like difference remained −1.92 under each specification (Supplementary Fig. S10A,B; Table S9A). ')
    proxy=('Within three GSE206710 donors, the three-gene proxy was higher at both DP stages than at DN-early, although DP-immature exceeded DP-mature. This supports early-to-cortical discrimination, not a monotonic clock (Supplementary Fig. S10C; Table S9B). ')
    paragraph=(proxy+common if journal.startswith('npj') else common if journal=='British Journal of Haematology' else common+proxy)
    assert 'HiChIP is retained as exploratory chromatin context' in md
    md=md.replace('HiChIP is retained as exploratory chromatin context',paragraph+'\n\nHiChIP is retained as exploratory chromatin context',1)
    if journal=='JCMM':
        program=('In an exploratory comparison of the fixed 18 developmental-score-matched BCL11B–ETP-like pairs, pair-blocked limma–voom and competitive Hallmark tests found higher interferon-γ response (FDR = 0.0025) and IL6–JAK–STAT3 signaling (FDR = 0.0028) in BCL11B cases (Supplementary Fig. S11; Table S9D). These correlated programs provide within-cohort molecular-phenotype context; they do not identify direct ZEB targets or a causal mechanism.\n\n')
        md=md.replace('### Independent bulk expression reproduces',program+'### Independent bulk expression reproduces',1)
    clinical=('Ridge Cox sensitivity penalizing subtype coefficients while leaving balance and clinical covariates unpenalized gave directionally similar balance estimates. Conditional bootstrap intervals used 499 valid fits of 500 attempts per endpoint (Table S9C); penalty choice changes subtype shrinkage, so this does not establish independent prognostic utility. ')
    md=md.replace('These endpoint-specific results do not establish calibrated or transportable clinical utility.',clinical+'These endpoint-specific results do not establish calibrated or transportable clinical utility.',1)
    if journal=='British Journal of Haematology':
        md=re.sub(r'(?s)(### Clinical associations are largely explained by molecular subtype composition\n\n).*?(\n\n---\n\n## Discussion)',
          r'\1Exploratory clinical associations attenuated after subtype adjustment. Sparse-subtype Cox fits were unstable; ridge sensitivity regularized subtype coefficients but did not establish a transportable prognostic marker (Supplementary Fig. S8; Tables S1, S4 and S9C).\2',md,count=1)
    methods=('Spline degrees of freedom 2–4 were compared in the same patients. Normal-stage proxy summaries respected donor nesting. Survival sensitivity used ridge Cox with subtype indicators penalized and balance, age, sex and log10 white-cell count unpenalized; minimum-deviance five-fold CV selected the penalty, Breslow ties were specified, and 500 patient bootstraps provided conditional intervals. The matched program analysis used pair-blocked limma–voom and limma camera Hallmark tests [51,52]. Details are in the Supplementary Methods.\n\n')
    md=md.replace('### Reproducibility\n\n','### Targeted sensitivity analyses\n\n'+methods+'### Reproducibility\n\n',1)
    md=md.replace('S1–S8','S1–S9').replace('version tag: submission-locked-v2-2026-09-26','version tag: v1.0-submission')
    if journal in ('JCMM','British Journal of Haematology'):
        md=re.sub(r'(?s)\n---\n\n\*\*Table 1\..*?\n---\n\n(?=## Figure legends)','\n---\n\n',md,count=1)
        md=md.replace('Table 1 and Supplementary Table S2','Supplementary Table S2').replace('Table 1','Supplementary Table S2')
    if journal.startswith('npj'):
        methods_from_si=SI.split('## Supplementary methods and interpretation notes',1)[1].split('## Supplementary Figure S1.',1)[0]+ADDED_METHODS
        md=md.replace('Post-review cross-fitting, alternative scores, repeated patient-level cross-validation and focused BCL11B contrasts are detailed in the Supplementary Methods.','Cross-fitting, alternative scores, repeated patient-level cross-validation and focused BCL11B contrasts are detailed below.')
        md=md.replace('Details are in the Supplementary Methods.','Additional procedure details follow below.')
        md=md.replace('Detailed definitions and full precision are in Supplementary Methods and Table S8','Detailed definitions and full precision are in Methods and Table S8')
        md=md.replace('## Data and code availability','### Additional computational methods\n'+methods_from_si+'\n## Data and code availability',1)
        funding=re.search(r'(?s)## Funding\n\n(.*?)\n\n---',md).group(1)
        md=re.sub(r'(?s)\n## Funding\n\n.*?\n\n---\n','\n',md,count=1)
        md=md.replace('## Acknowledgements\n\n','## Acknowledgements\n\n'+funding+'\n\n',1)
    if journal in ('Scientific Reports','npj Systems Biology and Applications'):
        ai=re.search(r'\*\*Artificial Intelligence:\*\* .*?\n',md).group(0)
        md=md.replace(ai,'',1)
        md=md.replace('### Reproducibility\n\n','### Use of artificial intelligence tools\n\n'+ai.replace('**Artificial Intelligence:** ','')+'\n### Reproducibility\n\n',1)
    md=md.rstrip()+'''\n\n51. Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database Hallmark Gene Set Collection. *Cell Syst*. 2015;1(6):417–425. doi:10.1016/j.cels.2015.12.004\n\n52. Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res*. 2012;40(17):e133. doi:10.1093/nar/gks461\n'''
    assert wc(ABSTRACTS[journal]) <= (150 if journal.startswith('npj') else 200),(journal,wc(ABSTRACTS[journal]))
    assert not journal.startswith('npj') or wc(TITLES[journal])<=15
    return md

def si_variant(journal):
    si=SI.replace('## Supplementary Tables S1–S8','## Supplementary Tables S1–S9')
    si=si.replace('## Supplementary Figure S1.',ADDED_METHODS+'\n## Supplementary Figure S1.',1)
    si+='\n'+ADDED_LEGENDS
    if journal.startswith('npj'):
        a=si.index('## Supplementary methods and interpretation notes');b=si.index('## Supplementary Figure S1.')
        si=si[:a]+'## Supplementary figures and tables\n\n'+si[b:]
    return si

def append_figures(md,kind):
    if kind=='main':
        for i in range(1,8):
            next_marker=f'**Figure {i+1}.' if i<7 else '## References'
            image=f'![](figures/main/Figure{i}.png){{width=6.5in}}\n\n'
            assert next_marker in md
            md=md.replace(next_marker,image+next_marker,1)
    else:
        for i in range(1,12):
            next_marker=f'## Supplementary Figure S{i+1}.' if i<11 else '## Supplementary Table S9.'
            image='\n\n'.join(f'![](figures/supplementary/FigureS{i}_page{p}.png){{width=6.5in}}' for p in (1,2)) if i<=3 else f'![](figures/supplementary/FigureS{i}.png){{width=6.5in}}'
            assert next_marker in md
            md=md.replace(next_marker,image+'\n\n'+next_marker,1)
    return md

def docx(md,target,scratch):
    md=re.sub(r'<sup>(.*?)</sup>',lambda m:'^'+m.group(1).replace(' ','\\ ')+'^',md)
    src=scratch/(target.stem+'.md');src.write_text(md,encoding='utf-8')
    subprocess.run([PANDOC,str(src),'-f','markdown+pipe_tables+superscript-raw_html',
                    '-t','docx','--reference-doc',str(REF),'--resource-path',str(SUB),
                    '-o',str(target)],check=True)

CSS="""@page{size:A4;margin:19mm 20mm 18mm 20mm}body{font-family:'Times New Roman';font-size:11pt;line-height:1.45;color:#17252d}h1,h2,h3{font-family:Arial;color:#203744;break-after:avoid}h1{font-size:17pt}h2{font-size:13pt;margin-top:1.2em}h3{font-size:11pt}p{margin:.22em 0 .55em;orphans:3;widows:3}table{width:100%;border-collapse:collapse;font-family:Arial;font-size:7.7pt}thead{display:table-header-group;background:#244455;color:white}tr{break-inside:avoid}th,td{padding:4px 5px;border-bottom:.4px solid #cbd5d8}math[display='block']{display:block;text-align:center;margin:.7em 0;break-inside:avoid}"""

def pdf(md,target,scratch):
    src=scratch/(target.stem+'.md');html=scratch/(target.stem+'.html')
    src.write_text(md,encoding='utf-8')
    subprocess.run([PANDOC,str(src),'-f','markdown+pipe_tables+superscript-raw_html',
                    '-t','html5','-s','--mathml','-o',str(html)],check=True)
    html.write_text(html.read_text(encoding='utf-8').replace('</head>','<style>'+CSS+'</style></head>',1),encoding='utf-8')
    with tempfile.TemporaryDirectory(prefix='zeb_journal_',dir=scratch) as profile:
        subprocess.run([str(CHROME),'--headless=new','--disable-gpu','--no-first-run',
                        '--no-pdf-header-footer',f'--user-data-dir={profile}',
                        f'--print-to-pdf={target}',html.as_uri()],check=True,timeout=180,
                        stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    assert target.exists() and target.stat().st_size>10000

def combine(parts,target):
    dest=fitz.open()
    for p in parts:
        with fitz.open(p) as src:dest.insert_pdf(src)
    dest.save(target,garbage=4,deflate=True);dest.close()

def workbook(target):
    wb=openpyxl.load_workbook(SUB/'artical/Supplementary_Tables_S1-S8.xlsx')
    sheets={'S9A_Spline_df':'spline_df_2_3_4.tsv',
            'S9B1_Normal_sources':'normal_proxy_stage_summary.tsv',
            'S9B2_Donor_summary':'normal_proxy_donor_summary.tsv',
            'S9B3_Donor_stage':'normal_proxy_donor_stage.tsv',
            'S9C_Ridge_Cox':'ridge_cox_sensitivity.tsv',
            'S9D_Matched_Hallmark':'matched_BCL11B_ETP_hallmark_camera.tsv'}
    labels={
      'S9A_Spline_df':['Spline df','Subtype','Patients','Median residual','Matched BCL11B–ETP-like mean difference','Development-only R²','Added subtype ΔR²'],
      'S9B1_Normal_sources':['Resource','Measured units','Donors','Stage-order Spearman ρ','Early score median','DP score median'],
      'S9B2_Donor_summary':['Donor','Stages','Stage-order Spearman ρ','DN→immature DP change','DN→mature DP change'],
      'S9B3_Donor_stage':['Donor','Annotated stage','Stage order','Mean proxy score','Aggregates'],
      'S9C_Ridge_Cox':['Endpoint','Patients','Events','Model','CV λ minimum','CV λ one-SE','Balance HR','Bootstrap 95% CI','Valid / attempted bootstraps','Bootstrap fraction HR<1'],
      'S9D_Matched_Hallmark':['Hallmark program','Tested genes','Direction','Competitive P','BH FDR']}
    for name,file in sheets.items():
        ws=wb.create_sheet(name)
        lines=(DATA/file).read_text(encoding='utf-8').splitlines()
        ws.append(labels[name])
        for line in lines[1:]:
            vals=[]
            for v in line.split('\t'):
                try:vals.append(float(v) if '.' in v or 'e' in v.lower() else int(v))
                except ValueError:vals.append(v)
            if name=='S9C_Ridge_Cox':
                vals=[vals[0],vals[1],vals[2],vals[3],vals[4],vals[5],vals[6],
                      f'{vals[7]:.3f}–{vals[8]:.3f}',f'{vals[9]} / {vals[10]}',vals[11]]
            ws.append(vals)
        ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions;ws.sheet_view.showGridLines=False
        ws.row_dimensions[1].height=32
        for c in ws[1]:
            c.font=Font(name='Arial',size=9,bold=True,color='FFFFFF')
            c.fill=PatternFill('solid',fgColor='28495A')
            c.alignment=Alignment(wrap_text=True,vertical='center')
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width=min(35,max(12,max(len(str(c.value or '')) for c in col[:30])+2))
        for row in ws.iter_rows(min_row=2):
            for c in row:
                c.font=Font(name='Arial',size=9)
                c.alignment=Alignment(vertical='center')
                if isinstance(c.value,float):
                    c.number_format='0.00E+00' if (name=='S9D_Matched_Hallmark' and c.column>=4 and c.value<.001) else '0.000'
    wb.save(target)

def render_s9_table_pdf(book,target,scratch):
    wb=openpyxl.load_workbook(book,read_only=True,data_only=True)
    names=[s for s in wb.sheetnames if s.startswith('S9')]
    import html as html_lib
    parts=['<html><head><meta charset="utf-8"><style>@page{size:A4 landscape;margin:14mm}body{font-family:Arial,sans-serif;color:#213744;font-size:9pt}h2{break-before:page;font-size:14pt;margin:0 0 8mm}h2:first-of-type{break-before:auto}table{border-collapse:collapse;width:100%;table-layout:auto}th{background:#28495a;color:white}th,td{border-bottom:.4px solid #ccd6db;padding:4px 5px;text-align:left;vertical-align:top;overflow-wrap:anywhere}tr{break-inside:avoid}thead{display:table-header-group}</style></head><body>']
    for name in names:
        ws=wb[name]
        parts.append('<h2>Supplementary Table '+html_lib.escape(name.replace('_',' — '))+'</h2><table>')
        for j,row in enumerate(ws.values):
            tag='th' if j==0 else 'td'
            parts.append('<tr>'+''.join(f'<{tag}>'+html_lib.escape(f'{v:.3g}' if isinstance(v,float) else str(v if v is not None else ''))+f'</{tag}>' for v in row)+'</tr>')
        parts.append('</table>')
    parts.append('</body></html>')
    page=scratch/'s9_tables.html';page.write_text(''.join(parts),encoding='utf-8')
    with tempfile.TemporaryDirectory(prefix='zeb_s9_',dir=scratch) as profile:
        subprocess.run([str(CHROME),'--headless=new','--disable-gpu','--no-first-run','--no-pdf-header-footer',f'--user-data-dir={profile}',f'--print-to-pdf={target}',page.as_uri()],check=True,timeout=180,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    assert target.exists() and target.stat().st_size>10000

def cover(journal):
    angle={
      'npj Systems Biology and Applications':'a developmental-context model separating a ZEB-independent expression proxy from subtype-associated ZEB states, tested with held-out subtype fitting, alternative scores and matched contrasts',
      'JCMM':'the molecular phenotype of BCL11B-associated leukemia relative to developmentally matched ETP-like disease, with an exploratory pair-adjusted transcriptomic comparison',
      'British Journal of Haematology':'a hematologically focused distinction between ETP-like immaturity and additional ZEB2-skewed expression in BCL11B-associated disease',
      'Scientific Reports':'a transparent public-data analysis with sensitivity tests, explicit external-replication boundaries and versioned source code'}[journal]
    extra=''
    if journal=='JCMM':
        extra=('\n**Recent relevant publications by corresponding author Chao Wu (past 3–5 years):**\n\n'
          '1. Wu C, Yu X, Li X, et al. Aberrant METTL14 gene expression contributes to malignant transformation of benzene-exposed myeloid cells. *Ecotoxicol Environ Saf*. 2024;276:116302. doi:10.1016/j.ecoenv.2024.116302.\n'
          '2. Wu C, Liu W, Hu X, et al. Targeting TNK2/ACK1 reverses the immunosuppressive tumor microenvironment and synergizes with immunochemotherapy in pancreatic cancer. *Nat Commun*. 2025;17:512. doi:10.1038/s41467-025-67197-3.\n'
          '3. Yang L, Tang C, Cheng X, et al. Syrosingopine enhances immune checkpoint blockade efficacy by inhibition of MCTs and metabolic reprogramming in triple-negative breast cancer. *Biochim Biophys Acta Mol Basis Dis*. 2026;1872(7):168308. doi:10.1016/j.bbadis.2026.168308.\n'
          '4. Li B, Li X, Ma M, et al. Analysis of long non-coding RNAs associated with disulfidptosis for prognostic signature and immunotherapy response in uterine corpus endometrial carcinoma. *Sci Rep*. 2023;13:22220. doi:10.1038/s41598-023-49750-6.\n')
    return f'''# Cover letter to {journal}

Dear Editor,

Please consider our Original Article, “{TITLES[journal]}”, for publication in *{journal}*.

The study asks whether relative ZEB1–ZEB2 expression across modern T-cell acute lymphoblastic leukemia subtypes is explained by a developmental-expression proxy. Its contribution for this journal is {angle}. The primary analysis includes 1,309 pediatric and young-adult patients across 17 molecular subtypes. BCL11B and TLX3 occupy opposite residual poles, and matched BCL11B and ETP-like patients retain a substantial contrast. Independent datasets are interpreted according to their actual scope; exploratory clinical and chromatin analyses are not presented as causal or predictive validation.

We have limited the claims to observational subtype-associated expression. The underlying datasets are public and de-identified. The analysis code and aggregate outputs are versioned at the repository identified in the manuscript. The submitting authors will complete the journal's originality, author-approval and disclosure declarations in the portal.

{extra}
Sincerely,\nChao Wu and Limei Li\nCorresponding authors\nchaowutjmuch@163.com; lilimei116@126.com
'''

GUIDES={
 'npj Systems Biology and Applications':'Official: https://www.nature.com/npjsba/content-types and https://www.nature.com/npjsba/for-authors-and-referees/submission-guidelines. Title ≤15 words; abstract ≤150; all Methods in main; funding in Acknowledgements; one SI PDF; figure legends ≤350 words.',
 'JCMM':'Author-provided limits: Original Article ≤4,500 words, abstract ≤200, references ≤60, ≤7 combined figures/tables. The Wiley journal page was inaccessible during build; recheck the submission portal. Cover letter requires verified relevant senior-author publications.',
 'British Journal of Haematology':'Author-provided limits: Original Paper ≤3,000 main words, abstract ≤200, references ≤60, ≤7 combined figures/tables. The Wiley journal page was inaccessible during build; recheck the submission portal.',
 'Scientific Reports':'Official: https://www.nature.com/srep/author-instructions/submission-guidelines. Recommended main text ≤4,500, title ≤20, abstract ≤200, legends ≤350; consult its submission checklist for mandatory limits.'}

def readme(journal,md):
    body=md.split('## Introduction',1)[1].split('## Data and code availability',1)[0]
    display='7 figures + Table 1' if journal not in ('JCMM','British Journal of Haematology') else '7 figures; dataset table in S2'
    return f'''# {journal} submission package

Title: {TITLES[journal]}
Abstract: {wc(ABSTRACTS[journal])} words; Introduction–Methods–Results–Discussion: approximately {wc(body)} words by build token count.
Main display: {display}; 52 references; 11 supplementary figures; workbook Tables S1–S9.

## Ready-to-read files

- `Manuscript.docx`: editable article and full legends.
- `Manuscript_with_figures.docx`: reading copy with all seven main figures embedded.
- `Main_Article_with_Figures.pdf`: reading PDF with native vector figures.
- `Supplementary_Information.docx` and `.pdf`: complete SI; npj Methods are moved into main.
- `Supplementary_Tables_S1-S9.xlsx`: reader workbook; S9 contains targeted analyses.
- `Figures/`: seven separate final-size PDF and TIFF plates.
- `Cover_Letter.docx`: journal-specific cover-letter draft.
- `Author_Information.xlsx`: ordered author list and affiliation mapping; unknown emails/ORCIDs are blank for author completion.

## Requirements used

{GUIDES[journal]}

## Before submitting

Confirm author approval of this journal variant, corresponding-author ORCID, funding/conflicts, and any preprint overlap. JCMM citations were verified by DOI against Crossref; confirm final author selections in the portal. Check final portal PDF conversion for formula and image pagination. Resolve the public `v1.0-submission` code tag; no Zenodo DOI is claimed until an archive is deposited.

Ridge Cox bootstrap intervals condition on penalty selection, and there is no external clinical prediction validation. Matched programs are exploratory and do not identify direct ZEB targets.
'''

def author_sheet(target):
    authors=[
      ('Qian Zhou','1','co-first',''),('Qingyun Ni','2','co-first',''),
      ('Fuhua Liang','3','data collection; analysis',''),
      ('Haonan Wang','3','',''),('Xu Liang','3','',''),('Yanmin Li','4','',''),
      ('Huixia Xu','1','',''),('Chao Wu','3','corresponding','chaowutjmuch@163.com'),
      ('Limei Li','5; 1','corresponding','lilimei116@126.com')]
    wb=openpyxl.Workbook();ws=wb.active;ws.title='Authors'
    ws.append(['Order','Full name','Affiliation number(s)','Role note','Email','ORCID'])
    for i,(name,aff,role,email) in enumerate(authors,1):ws.append([i,name,aff,role,email,''])
    ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions;ws.sheet_view.showGridLines=False
    widths=[9,24,24,32,34,26]
    for i,w in enumerate(widths,1):ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width=w
    for c in ws[1]:
        c.font=Font(name='Arial',size=10,bold=True,color='FFFFFF')
        c.fill=PatternFill('solid',fgColor='28495A')
    for row in ws.iter_rows(min_row=2):
        for c in row:c.font=Font(name='Arial',size=10)
    aff=wb.create_sheet('Affiliations')
    for line in MASTER.splitlines()[4:10]:
        if line.startswith('<sup>') and 'Department' in line:aff.append([re.sub('<[^>]+>','',line)])
    aff.column_dimensions['A'].width=140
    wb.save(target)

def build(journal):
    folder=OUT/journal;folder.mkdir(parents=True,exist_ok=True)
    scratch=folder/'_working';scratch.mkdir(exist_ok=True)
    md=main_variant(journal);si=si_variant(journal)
    (folder/'Manuscript.md').write_text(md,encoding='utf-8')
    (folder/'Supplementary_Information.md').write_text(si,encoding='utf-8')
    docx(md,folder/'Manuscript.docx',scratch)
    docx(append_figures(md,'main'),folder/'Manuscript_with_figures.docx',scratch)
    docx(append_figures(si,'si'),folder/'Supplementary_Information.docx',scratch)
    pdf(md,scratch/'main_text.pdf',scratch)
    combine([scratch/'main_text.pdf']+[SUB/f'figures/main/Figure{i}.pdf' for i in range(1,8)],folder/'Main_Article_with_Figures.pdf')
    pdf(si,scratch/'si_text.pdf',scratch)
    workbook(folder/'Supplementary_Tables_S1-S9.xlsx')
    render_s9_table_pdf(folder/'Supplementary_Tables_S1-S9.xlsx',scratch/'s9_tables.pdf',scratch)
    combine([scratch/'si_text.pdf']+[SUB/f'figures/supplementary/FigureS{i}.pdf' for i in range(1,12)]+[SUB/'first_submission/Supplementary_Tables_S1-S8_readable.pdf',scratch/'s9_tables.pdf'],folder/'Supplementary_Information.pdf')
    figfolder=folder/'Figures';figfolder.mkdir(exist_ok=True)
    for i in range(1,8):
        for ext in ('pdf','tiff'):shutil.copy2(SUB/f'figures/main/Figure{i}.{ext}',figfolder/f'Figure{i}.{ext}')
    letter=cover(journal);(folder/'Cover_Letter.md').write_text(letter,encoding='utf-8')
    docx(letter,folder/'Cover_Letter.docx',scratch)
    author_sheet(folder/'Author_Information.xlsx')
    (folder/'README.md').write_text(readme(journal,md),encoding='utf-8')
    if journal=='British Journal of Haematology':
        enquiry=f'''# Presubmission enquiry — British Journal of Haematology

Dear Editors,

Would an Original Paper reporting a public-data developmental-context analysis of ZEB1–ZEB2 expression across 1,309 pediatric and young-adult T-ALL patients be of interest? Its central observation is that BCL11B-associated ZEB2 skewing persists relative to developmental-score-matched ETP-like disease. We supply the abstract below.

## Abstract

{ABSTRACTS[journal]}

Sincerely,\nChao Wu and Limei Li
'''
        docx(enquiry,folder/'Presubmission_Enquiry.docx',scratch)
    shutil.rmtree(scratch)
    print(journal,'abstract',wc(ABSTRACTS[journal]),'body',wc(md.split('## Introduction',1)[1].split('## Data and code availability',1)[0]))

if __name__=='__main__':
    for journal in TITLES:build(journal)
