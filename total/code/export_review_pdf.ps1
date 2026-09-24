param([string]$DocumentName = 'ZEB1_TALL_manuscript', [switch]$LayoutOnly)
$ErrorActionPreference = 'Stop'
$root = Split-Path $PSScriptRoot -Parent
$outDir = Join-Path $root 'manuscript/review_pdf'
New-Item -ItemType Directory -Path $outDir -Force | Out-Null
$log = Join-Path $outDir "$DocumentName.export.log"
function Stage([string]$message) {
    "$(Get-Date -Format o) $message" | Add-Content -LiteralPath $log -Encoding utf8
}
Stage 'Creating dedicated Word automation instance'
$wordBefore = @(Get-Process WINWORD -ErrorAction SilentlyContinue | ForEach-Object { $_.Id })
$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
$wordNew = @(Get-CimInstance Win32_Process -Filter "Name='WINWORD.EXE'" | Where-Object { $_.ProcessId -notin $wordBefore -and $_.CommandLine -match '/Automation' })
if ($wordNew.Count -eq 1) {
    $wordNew[0].ProcessId | Set-Content -LiteralPath (Join-Path $outDir "$DocumentName.wordpid")
}
try {
    foreach ($name in @($DocumentName)) {
        $source = Join-Path $outDir "$name.proof.docx"
        Copy-Item -LiteralPath (Join-Path $root "manuscript/$name.docx") -Destination $source -Force
        $target = Join-Path $outDir "$name.pdf"
        $document = $null
        try {
            Stage 'Opening read-only proof copy'
            $document = $word.Documents.Open($source, $false, $true, $false)
            if ($LayoutOnly) {
                Stage 'Opened; inspecting Word pagination without PDF export'
                $pages = $document.ComputeStatistics(2)
                $records = @()
                foreach ($paragraph in $document.Paragraphs) {
                    $range = $paragraph.Range
                    $text = $range.Text.Trim()
                    $imageCount = $range.InlineShapes.Count
                    if ($imageCount -gt 0 -or $text.Length -gt 0) {
                        $first = $range.Duplicate
                        $first.Collapse(1)
                        $last = $range.Duplicate
                        if ($last.End -gt $last.Start) { $last.End = $last.End - 1 }
                        $last.Collapse(0)
                        $records += [pscustomobject]@{text=$text;images=$imageCount;page_start=$first.Information(3);page_end=$last.Information(3)}
                    }
                }
                [pscustomobject]@{pages=$pages;records=$records} | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $outDir "$name.word_layout.json") -Encoding utf8
                Stage "Layout inspection complete: $pages pages"
                continue
            }
            Stage 'Opened; exporting PDF'
            $document.ExportAsFixedFormat($target, 17)
            Stage 'Exported; counting pages'
            $pages = $document.ComputeStatistics(2)
            Stage "Complete: $pages pages"
            Write-Output "$name : $pages pages; $target"
        } finally {
            if ($null -ne $document) { $document.Close(0) }
        }
    }
} finally {
    $word.Quit()
    [void][System.Runtime.InteropServices.Marshal]::ReleaseComObject($word)
}
