
param(
    [switch]$Debug,
    [switch]$Clean
)

$ErrorActionPreference = "Stop"
$root = $PSScriptRoot
$out  = Join-Path $root "build"
$gen  = Join-Path $out "generated"

if ($Clean -and (Test-Path $out)) { Remove-Item -Recurse -Force $out }
New-Item -ItemType Directory -Force $gen | Out-Null

$rootFwd = $root -replace '\\', '/'
$config  = "#pragma once`n#define LDLT_PROJECT_ROOT `"$rootFwd`"`n"
[System.IO.File]::WriteAllText((Join-Path $gen "ldlt_config.hpp"), $config,
                               (New-Object System.Text.UTF8Encoding $false))

function Invoke-Compiler([string]$exe, [string[]]$arguments) {
    $saved = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    & $exe @arguments 2>&1 | ForEach-Object { "$_" }
    $code = $LASTEXITCODE
    $ErrorActionPreference = $saved
    if ($code -ne 0) { throw "$exe failed with exit code $code" }
}

$sqliteObj = Join-Path $out "sqlite3.o"
if (-not (Test-Path $sqliteObj)) {
    Write-Host "compiling sqlite (one time, about half a minute)..."
    Invoke-Compiler "gcc" @("-O2", "-DSQLITE_THREADSAFE=0", "-DSQLITE_OMIT_LOAD_EXTENSION",
                            "-c", (Join-Path $root "third_party\sqlite\sqlite3.c"), "-o", $sqliteObj)
}

$flags = @("-std=c++17", "-Wall", "-Wextra", "-Wpedantic", "-Wshadow", "-Wconversion",
           "-I$root\src\core", "-I$root\third_party\sqlite", "-I$gen")
if ($Debug) { $flags += @("-g", "-O0", "-D_GLIBCXX_ASSERTIONS") }
else        { $flags += @("-O2", "-DNDEBUG") }

$lib = Get-ChildItem (Join-Path $root "src\core\*.cpp") | ForEach-Object { $_.FullName }
$cli = Get-ChildItem (Join-Path $root "src\cli\*.cpp")  | ForEach-Object { $_.FullName }

Write-Host "compiling ldlt..."
Invoke-Compiler "g++" ($flags + $lib + $cli + @($sqliteObj, "-o", (Join-Path $out "ldlt.exe")))

Write-Host "build ok -> $out\ldlt.exe"
Write-Host "run it with:  .\build\ldlt.exe"
