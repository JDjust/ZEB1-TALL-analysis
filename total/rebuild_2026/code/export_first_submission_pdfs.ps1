$ErrorActionPreference = 'Stop'
$root = 'D:\_bioinformation\ZEB1'
$first = Join-Path $root 'submission\first_submission'
python (Join-Path $root 'total\rebuild_2026\code\render_supplement_tables_pdf.py')
if ($LASTEXITCODE -ne 0) { throw 'Supplementary table PDF rendering failed' }
python (Join-Path $root 'total\rebuild_2026\code\build_vector_preview.py')
if ($LASTEXITCODE -ne 0) { throw 'Vector manuscript PDF rendering failed' }
Get-ChildItem -LiteralPath $first -Filter '*.pdf' |
    Where-Object { $_.Name -in @('Main_Article_with_Figures.pdf','Supplementary_Information_complete.pdf') } |
    Select-Object Name, Length
