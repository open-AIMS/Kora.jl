#!/usr/bin/env bash
# build/build.sh -- Compile Kora.jl into a deployable library via juliac.
#
# Usage (from the Kora.jl root):
#   ./build/build.sh [--mode <mode>] [--output-dir <dir>]
#
# Modes:
#   native    (default) juliac --trim=safe --bundle -> trimmed native shared lib + bundled
#                       runtime; no Julia on target, same-architecture machines only
#   bundled             juliac --bundle      -> shared lib + libjulia, stdlibs, and
#                       artifacts bundled automatically; no Julia on target required
#   sysimage            juliac --output-sysimage -> sysimage loaded by Julia at startup,
#                       requires Julia on target, portable to any machine with Julia
#   worker              juliac --output-exe --bundle -> standalone kora-worker executable
#                       for the web server backend; reads WireSimParams from stdin,
#                       writes WireEnsembleResult to stdout, loops until stdin closes.
#                       Build and test on Linux (WSL or CI) — not Windows native.
#   server              juliac --output-exe --bundle, untrimmed -> standalone kora-server
#                       executable (Oxygen.jl HTTP app). Compiled against
#                       build/server/Project.toml (Kora path dep + real Oxygen/HTTP
#                       deps) rather than the top-level Project.toml, so Oxygen stays
#                       out of every other consumer's dependency tree. No --trim:
#                       Oxygen's macro-based routing needs the dynamic compiler kept
#                       around (same rationale as ReefGuideWorker.jl's build).
#
# Overridable env vars:
#   KORA_LIB_DIR   output directory (default: build/dist/<mode> inside Kora.jl root)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
ENTRY_FILE="$SCRIPT_DIR/bridge_aot.jl"
WORKER_ENTRY_FILE="$SCRIPT_DIR/worker_main.jl"
SERVER_ENTRY_FILE="$SCRIPT_DIR/kora_server_main.jl"
SERVER_PROJECT_DIR="$SCRIPT_DIR/server"

MODE="native"
OUTPUT_DIR="${KORA_LIB_DIR:-}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)       MODE="$2";       shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

OUTPUT_DIR="${OUTPUT_DIR:-$PROJECT_ROOT/build/dist/$MODE}"
mkdir -p "$OUTPUT_DIR"

# Multi-target CPU dispatch — makes the sysimage (and other build outputs)
# usable across different x86_64 microarchitectures without recompilation.
export JULIA_CPU_TARGET="generic;x86_64,sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)"

# JuliaC only adds -lm on i686; on x86_64, floorf (used by Kora's numeric
# code) becomes a libcall to libm and the link fails regardless of --trim
# mode. Point JULIA_CC at the returned wrapper script to inject -lm.
_make_lm_wrapper() {
    local wrapper
    wrapper="$(mktemp /tmp/gcc-wrapper-XXXXXX.sh)"
    printf '#!/bin/sh\nexec gcc "$@" -lm\n' > "$wrapper"
    chmod +x "$wrapper"
    echo "$wrapper"
}

case "$MODE" in
    native)
        echo "Mode: native (trimmed, no Julia runtime required on target with reduced compatibility)"
        BUILD_LOG="$OUTPUT_DIR/build.log"
        echo "Build log: $BUILD_LOG"
        time juliac --verbose --project="$PROJECT_ROOT" --output-lib "$OUTPUT_DIR/kora_bridge" \
            --bundle "$OUTPUT_DIR" --trim=safe --compile-ccallable --experimental "$ENTRY_FILE" \
            2>&1 | tee "$BUILD_LOG"
        ;;

    bundled)
        echo "Mode: bundled (libjulia, stdlibs, and artifacts bundled via --bundle)"
        BUILD_LOG="$OUTPUT_DIR/build.log"
        echo "Build log: $BUILD_LOG"
        time juliac --verbose --project="$PROJECT_ROOT" --output-lib "$OUTPUT_DIR/kora_bridge" \
            --bundle "$OUTPUT_DIR" --compile-ccallable --experimental "$ENTRY_FILE" \
            2>&1 | tee "$BUILD_LOG"
        ;;

    sysimage)
        echo "Mode: sysimage (requires Julia on target, portable across Julia-supported platforms)"
        juliac --project="$PROJECT_ROOT" --output-sysimage "$OUTPUT_DIR/kora_bridge" \
            --compile-ccallable --experimental "$ENTRY_FILE"
        ;;

    worker)
        echo "Mode: worker (standalone exe, no Julia runtime required on target)"
        # PackageCompiler's artifact bundler errors if a destination artifact dir
        # from a previous build already exists (it doesn't pass force=true) —
        # start from a clean output dir each time.
        rm -rf "$OUTPUT_DIR"
        mkdir -p "$OUTPUT_DIR"
        BUILD_LOG="$OUTPUT_DIR/build.log"
        echo "Build log: $BUILD_LOG"
        _GCC_WRAPPER="$(_make_lm_wrapper)"
        JULIA_CC="$_GCC_WRAPPER" \
        time juliac --verbose --project="$PROJECT_ROOT" --output-exe kora-worker \
            --bundle "$OUTPUT_DIR" --trim=safe --experimental "$WORKER_ENTRY_FILE" \
            2>&1 | tee "$BUILD_LOG"
        rm -f "$_GCC_WRAPPER"
        ;;

    server)
        echo "Mode: server (standalone exe, no Julia runtime required on target, untrimmed)"
        rm -rf "$OUTPUT_DIR"
        mkdir -p "$OUTPUT_DIR"
        BUILD_LOG="$OUTPUT_DIR/build.log"
        echo "Build log: $BUILD_LOG"
        julia --project="$SERVER_PROJECT_DIR" -e 'using Pkg; Pkg.instantiate()'
        _GCC_WRAPPER="$(_make_lm_wrapper)"
        JULIA_CC="$_GCC_WRAPPER" \
        time juliac --verbose --project="$SERVER_PROJECT_DIR" --output-exe kora-server \
            --bundle "$OUTPUT_DIR" --experimental "$SERVER_ENTRY_FILE" \
            2>&1 | tee "$BUILD_LOG"
        rm -f "$_GCC_WRAPPER"
        ;;

    *)
        echo "ERROR: Unknown mode '$MODE'. Valid modes: native, bundled, sysimage, worker, server" >&2
        exit 1
        ;;
esac

# Remove import library (.dll.a) — not needed for distribution
for f in "$OUTPUT_DIR"/*.dll.a; do
    [[ -e "$f" ]] || continue
    rm -f "$f"
    echo "Removed $(basename "$f")"
done
