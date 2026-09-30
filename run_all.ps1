param([string[]]$Figures=@())
$ErrorActionPreference='Stop'
$r=if ($env:ZEB1_RSCRIPT) { $env:ZEB1_RSCRIPT } else { 'D:/R/R-4.6.0/bin/Rscript.exe' }
& $r --vanilla (Join-Path $PSScriptRoot 'scripts/06_figures/run_all.R') @Figures
if ($LASTEXITCODE -ne 0) { throw "Rendering failed: $LASTEXITCODE" }
