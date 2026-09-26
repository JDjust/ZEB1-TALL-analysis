"""Build self-contained main/supplement Markdown and DOCX from reviewed text."""
from pathlib import Path
import csv, re, shutil, subprocess

base=Path(__file__).resolve().parents[1]
man=base/'manuscript'
main=man/'ZEB1_ZEB2_TALL_manuscript.md'
supp=man/'Supplementary_Information.md'
refs=(man/'References.md').read_text(encoding='utf-8')
main_text=main.read_text(encoding='utf-8')
legend_file=man/'Figure_Legends.md'
if not legend_file.exists():
    legend_file.write_text('# Figure legends\n\n'+ '\n\n'.join(
        (man/f'Figure_{i}_legend.md').read_text(encoding='utf-8').strip()
        for i in range(1,8))+'\n', encoding='utf-8')
if '## Figure legends' not in main_text:
    heading='## References\n\n'
    if heading not in main_text:
        raise SystemExit('Cannot locate references insertion point for figure legends')
    legend_text=legend_file.read_text(encoding='utf-8').replace('# Figure legends','## Figure legends',1)
    main_text=main_text.replace(heading,legend_text+'\n'+heading,1)
    main.write_text(main_text,encoding='utf-8')
if 'The numbered, DOI-identity-checked bibliography' in main_text:
    main_text=main_text.replace(
        '## References\n\nThe numbered, DOI-identity-checked bibliography is in the accompanying `References.md` and is appended to the submission DOCX/PDF.',
        '## References\n\n'+refs.replace('# References\n\n','',1).replace(
            'Bibliographic DOI/title identity was checked against Crossref on 2026-09-25. Citation-to-claim suitability remains a manuscript-level judgment.\n\n','',1).rstrip())
    main.write_text(main_text+'\n',encoding='utf-8')
if main_text.count('## References')!=1:
    raise SystemExit('Main manuscript references heading missing or duplicated')

supp_text=supp.read_text(encoding='utf-8')
if '## Supplementary Figure S1.' not in supp_text:
    supp_text+='\n'+(man/'Supplementary_Figure_Legends_consolidated.md').read_text(encoding='utf-8')
    table_dir=base/'data/source_data_rebuilt/supplementary_tables'
    for number,title,filename in [
        ('S1','Frozen clinical models and subtype attenuation','Table_S1_clinical_models.tsv'),
        ('S2','Datasets, biological units and evidentiary roles','Table_S2_dataset_units.tsv'),
        ('S3','Deprioritized exploratory lines','Table_S3_deprioritized_lines.tsv'),
        ('S4','Full 17-subtype developmental residual summary','Table_S4_full_subtype_residuals.tsv')]:
        with (table_dir/filename).open(encoding='utf-8-sig',newline='') as fh:
            rows=list(csv.reader(fh,delimiter='\t'))
        supp_text+=f'\n## Supplementary Table {number}. {title}\n\n'
        supp_text+='| '+' | '.join(rows[0])+' |\n'
        supp_text+='| '+' | '.join(['---']*len(rows[0]))+' |\n'
        for row in rows[1:]:
            vals=[]
            for val in row:
                val=val.replace('|','\\|').replace('\n',' ')
                if re.fullmatch(r'-?\d+\.\d+(?:[eE][+-]?\d+)?',val):
                    try: val=f'{float(val):.4g}'
                    except ValueError: pass
                vals.append(val)
            supp_text+='| '+' | '.join(vals)+' |\n'
        supp_text+='\n'
    supp.write_text(supp_text,encoding='utf-8')

pandoc=shutil.which('pandoc')
if not pandoc:
    raise SystemExit('Pandoc unavailable; self-contained Markdown was built')
for source,target in [(main,man/'ZEB1_ZEB2_TALL_manuscript.docx'),
                      (supp,man/'ZEB1_ZEB2_TALL_Supplementary.docx')]:
    subprocess.run([pandoc,str(source),'-f','gfm','-t','docx','-o',str(target),
                    '--standalone'],check=True)
    print(target.name,target.stat().st_size)
print('Main words',len(main.read_text(encoding='utf-8').split()),
      'Supplement words',len(supp.read_text(encoding='utf-8').split()))
