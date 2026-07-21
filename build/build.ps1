# build/build.ps1 -- Compile Kora.jl into a deployable library via juliac.
#
# Usage (from the Kora.jl root):
#   .\build\build.ps1 [-Mode <mode>] [-OutputDir <dir>]
#
# Modes:
#   native    (default) juliac --trim=safe --bundle -> trimmed native shared lib + bundled
#                       runtime; no Julia on target, same-architecture machines only
#   bundled             juliac --bundle      -> shared lib + libjulia, stdlibs, and
#                       artifacts bundled automatically; no Julia on target required
#   sysimage            juliac --output-sysimage -> sysimage loaded by Julia at startup,
#                       requires Julia on target, portable to any machine with Julia
#
# Overridable env vars:
#   KORA_LIB_DIR   output directory (default: build/dist/<mode> inside Kora.jl root)

param(
    [ValidateSet("native", "bundled", "sysimage")]
    [string]$Mode      = "native",
    [string]$OutputDir = ""
)

$ScriptDir   = Split-Path -Parent $MyInvocation.MyCommand.Path
$ProjectRoot = Split-Path -Parent $ScriptDir
$EntryFile   = Join-Path $ScriptDir "bridge_aot.jl"

if ($OutputDir -eq "") {
    $OutputDir = if ($env:KORA_LIB_DIR) { $env:KORA_LIB_DIR } `
                 else { Join-Path $ProjectRoot "build/dist/$Mode" }
}

# Clean any stale output from a previous run -- juliac's exact output
# layout (e.g. whether kora_bridge.dll lands flat or under bin\) has
# changed across versions, and leftover artifacts from an older layout
# or an older source revision can silently shadow the fresh build
# (see: kf_new_dhw_trajectory missing from a stale bin\kora_bridge.dll
# while the flat copy at $OutputDir root was current).
if (Test-Path $OutputDir) { Remove-Item -Recurse -Force $OutputDir }
New-Item -ItemType Directory -Force $OutputDir | Out-Null

$LibStem = Join-Path $OutputDir "kora_bridge"

# Multi-target CPU dispatch — makes the sysimage (and other build outputs)
# usable across different x86_64 microarchitectures without recompilation.
$env:JULIA_CPU_TARGET = "generic;x86_64,sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)"

switch ($Mode) {
    "native" {
        Write-Host "Mode: native (trimmed, no Julia runtime required on target with reduced compatibility)"
        $BuildLog = Join-Path $OutputDir "build.log"
        Write-Host "Build log: $BuildLog"
        Measure-Command {
            juliac --verbose --project=$ProjectRoot --output-lib $LibStem `
                --bundle $OutputDir --trim=safe --compile-ccallable --experimental $EntryFile `
                *>&1 | Tee-Object -FilePath $BuildLog
        }
    }

    "bundled" {
        Write-Host "Mode: bundled (libjulia, stdlibs, and artifacts bundled via --bundle)"
        $BuildLog = Join-Path $OutputDir "build.log"
        Write-Host "Build log: $BuildLog"
        Measure-Command {
            juliac --verbose --project=$ProjectRoot --output-lib $LibStem `
                --bundle $OutputDir --compile-ccallable --experimental $EntryFile `
                *>&1 | Tee-Object -FilePath $BuildLog
        }
    }

    "sysimage" {
        Write-Host "Mode: sysimage (requires Julia on target, portable across Julia-supported platforms)"
        Measure-Command {
            juliac --project=$ProjectRoot --output-sysimage $LibStem `
                --compile-ccallable --experimental $EntryFile
        }
    }
}

# Remove import library (.dll.a) — not needed for distribution
Get-ChildItem (Join-Path $OutputDir "*.dll.a") -ErrorAction SilentlyContinue | ForEach-Object {
    Remove-Item $_.FullName -Force
    Write-Host "Removed $($_.Name)"
}
