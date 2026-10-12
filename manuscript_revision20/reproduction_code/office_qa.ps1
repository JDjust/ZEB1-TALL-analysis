param([string[]]$Names=@('Fig1','Fig2','Fig3','Fig4','Fig5','Fig6','SFig1','SFig2','SFig3','SFig4','SFig5','SFig6'))
$ErrorActionPreference='Stop'
$base=Split-Path -Parent $PSScriptRoot
$qa=Join-Path $base 'QA'
$app=New-Object -ComObject PowerPoint.Application
$reports=@()
try {
 foreach($name in $Names) {
  $path=Join-Path $base "figures/$name/$name.pptx"
  $before=(Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash
  $pres=$app.Presentations.Open($path,-1,0,0)
  try {
   $height=[int](2400*$pres.PageSetup.SlideHeight/$pres.PageSetup.SlideWidth)
   $pres.Slides.Item(1).Export((Join-Path $qa ($name+'_office_original.png')),'PNG',2400,$height)
  } finally {$pres.Close();[void][Runtime.InteropServices.Marshal]::ReleaseComObject($pres)}
  $probe=Join-Path $qa ($name+'_edit_probe.pptx')
  $pres=$app.Presentations.Open($probe,0,0,0)
  $changed=$false
  $shapeEdited=$false
  $vectorId=0
  $vectorBefore=0.0
  $vectorAfter=0.0
  try {
   foreach($shape in $pres.Slides.Item(1).Shapes) {
    if($shape.HasTextFrame -eq -1 -and $shape.TextFrame.TextRange.Text -eq 'EDIT') {$shape.TextFrame.TextRange.Text='EDIT-OFFICE';$changed=$true;break}
   }
   if(-not $changed){throw "Office editable text not found: $name"}
   foreach($shape in $pres.Slides.Item(1).Shapes) {
    if($shape.Name -eq 'Editable plot path' -and $shape.Width -gt 2 -and $shape.Height -gt 2 -and $shape.Width -lt 50 -and $shape.Height -lt 50) {
     $vectorId=$shape.Id
     $vectorBefore=[double]$shape.Left
     $shape.Left=$shape.Left+8
     $vectorAfter=[double]$shape.Left
     $shapeEdited=$true
     break
    }
   }
   if(-not $shapeEdited){throw "Office editable native drawing not found: $name"}
   $saved=Join-Path $qa ($name+'_office_saved.pptx')
   $pres.SaveCopyAs($saved,24)
  } finally {$pres.Close();[void][Runtime.InteropServices.Marshal]::ReleaseComObject($pres)}
  $pres=$app.Presentations.Open($saved,-1,0,0)
  try {
   $found=$false
   $shapePersisted=$false
   foreach($shape in $pres.Slides.Item(1).Shapes){
    if($shape.HasTextFrame -eq -1 -and $shape.TextFrame.TextRange.Text -eq 'EDIT-OFFICE'){$found=$true}
    if($shape.Id -eq $vectorId -and [math]::Abs($shape.Left-$vectorAfter) -lt .02){$shapePersisted=$true}
   }
   if(-not $found){throw "Saved Office edit did not persist: $name"}
   if(-not $shapePersisted){throw "Saved Office native drawing edit did not persist: $name"}
   $pres.Slides.Item(1).Export((Join-Path $qa ($name+'_office_edit.png')),'PNG',2400,$height)
  } finally {$pres.Close();[void][Runtime.InteropServices.Marshal]::ReleaseComObject($pres)}
  $after=(Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash
  $reports+=@{figure=$name;office_open_export=$true;native_office_edit_save_reopen=$true;native_office_shape_edit_save_reopen=$shapePersisted;shape_edit=@{id=$vectorId;before_left_pt=$vectorBefore;after_left_pt=$vectorAfter;delta_pt=8};original_hash_preserved=($before -eq $after);original_sha256=$after;office_version=$app.Version}
  Write-Output "$name OFFICE_CHECKED"
 }
 $target=Join-Path $qa 'OFFICE_QA.json'
 if(Test-Path -LiteralPath $target){
  $active=@('Fig1','Fig2','Fig3','Fig4','Fig5','Fig6','SFig1','SFig2','SFig3','SFig4','SFig5','SFig6')
  foreach($record in (Get-Content -LiteralPath $target -Raw | ConvertFrom-Json)){
   if($record.figure -notin $Names -and $record.figure -in $active){$reports+=$record}
  }
 }
 $reports | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $target -Encoding utf8
} finally {[void][Runtime.InteropServices.Marshal]::ReleaseComObject($app)}
