# build/build.ps1 -- AOT-compile Kora into a native shared library via juliac.
#
# Usage (from the Kora.jl root):
#   .\build\build.ps1 [-OutputDir <dir>]
#
# Overridable env vars:
#   KORA_LIB_DIR   output directory (default: julia_lib/ inside Kora.jl root)

param(
    [string]$OutputDir = ""
)

$ScriptDir   = Split-Path -Parent $MyInvocation.MyCommand.Path
$ProjectRoot = Split-Path -Parent $ScriptDir
$EntryFile   = Join-Path $ScriptDir "bridge_aot.jl"

if ($OutputDir -eq "") {
    $OutputDir = if ($env:KORA_LIB_DIR) { $env:KORA_LIB_DIR } `
                 else { Join-Path $ProjectRoot "julia_lib" }
}

New-Item -ItemType Directory -Force $OutputDir | Out-Null

$LibStem = Join-Path $OutputDir "kora_bridge"

juliac --project=$ProjectRoot --output-lib $LibStem --trim=safe --compile-ccallable --experimental $EntryFile
