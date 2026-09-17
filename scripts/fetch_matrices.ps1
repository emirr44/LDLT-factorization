# fetch_matrices.ps1 - download test matrices from the SuiteSparse Matrix
# Collection (Harwell-Boeing group) into data\suitesparse\<name>.mtx.
#
#   .\scripts\fetch_matrices.ps1                 default bcsstk set
#   .\scripts\fetch_matrices.ps1 bcsstk01 nos4   specific matrices (HB group)
#   .\scripts\fetch_matrices.ps1 GHS_indef/laser other groups as group/name
#   .\scripts\fetch_matrices.ps1 -Indefinite     default symmetric indefinite set
#
# Files are ~0.1-2 MB each; the directory is git-ignored.  Requires tar.exe
# (present on Windows 10 1803+).

param([string[]]$Names, [switch]$Indefinite)

$ErrorActionPreference = "Stop"
if ($Indefinite) {
    $Names = @("GHS_indef/laser","GHS_indef/qpband","GHS_indef/aug3d","GHS_indef/tuma2",
               "GHS_indef/sit100","GHS_indef/stokes64")
} elseif (-not $Names) {
    $Names = @("bcsstk01","bcsstk02","bcsstk03","bcsstk04","bcsstk05","bcsstk06",
               "bcsstk08","bcsstk09","bcsstk10","bcsstk11","bcsstk13","bcsstk14")
}

$root = Split-Path $PSScriptRoot -Parent
$dest = Join-Path $root "data\suitesparse"
$tmp  = Join-Path $dest "_tmp"
New-Item -ItemType Directory -Force $dest | Out-Null
New-Item -ItemType Directory -Force $tmp  | Out-Null

foreach ($spec in $Names) {
    $group = "HB"; $name = $spec
    if ($spec -match "^([^/]+)/(.+)$") { $group = $matches[1]; $name = $matches[2] }
    $out = Join-Path $dest "$name.mtx"
    if (Test-Path $out) { Write-Host "have    $name"; continue }
    $url = "https://suitesparse-collection-website.herokuapp.com/MM/$group/$name.tar.gz"
    $tgz = Join-Path $tmp "$name.tar.gz"
    Write-Host "fetch   $name"
    Invoke-WebRequest -Uri $url -OutFile $tgz -UseBasicParsing
    & tar -xzf $tgz -C $tmp
    if ($LASTEXITCODE -ne 0) { throw "tar failed for $name" }
    Move-Item (Join-Path $tmp "$name\$name.mtx") $out -Force
}
Remove-Item -Recurse -Force $tmp
Write-Host "done -> $dest"
