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
#
# Overridable env vars:
#   KORA_LIB_DIR   output directory (default: build/dist/<mode> inside Kora.jl root)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
ENTRY_FILE="$SCRIPT_DIR/bridge_aot.jl"

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

    *)
        echo "ERROR: Unknown mode '$MODE'. Valid modes: native, bundled, sysimage" >&2
        exit 1
        ;;
esac

# Remove import library (.dll.a) — not needed for distribution
for f in "$OUTPUT_DIR"/*.dll.a; do
    [[ -e "$f" ]] || continue
    rm -f "$f"
    echo "Removed $(basename "$f")"
done
