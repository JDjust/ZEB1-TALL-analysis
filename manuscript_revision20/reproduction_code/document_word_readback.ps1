$ErrorActionPreference='Stop'
$base=Split-Path -Parent $PSScriptRoot
$app=New-Object -ComObject Word.Application
$app.Visible=$false
$app.DisplayAlerts=0
$reports=@()
try {
 foreach($name in @('MANUSCRIPT','SUPPLEMENTARY_INFORMATION')) {
  $path=Join-Path $base "manuscript/$name.docx"
  $before=(Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash
  Write-Output "$name WORD_OPENING"
  $doc=$app.Documents.Open($path,$false,$true,$false)
  try {
   $equations=$doc.OMaths.Count
   $paragraphs=$doc.Paragraphs.Count
   $reports+=@{document=$name;office_version=$app.Version;opened=$true;equation_objects=$equations;paragraphs=$paragraphs;pdf_exported=$false;pdf_engine='Microsoft Edge documented fallback';original_sha256=$before;original_hash_preserved=($before -eq (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash)}
   Write-Output "$name WORD_READBACK $paragraphs paragraphs $equations equations"
  } finally {$doc.Close(0);[void][Runtime.InteropServices.Marshal]::ReleaseComObject($doc)}
 }
 $reports | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $base 'QA/DOCUMENT_OFFICE_QA.json') -Encoding utf8
} finally {$app.Quit();[void][Runtime.InteropServices.Marshal]::ReleaseComObject($app)}
