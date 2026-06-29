#!/usr/bin/env bash
# build/build.sh -- AOT-compile Kora into a native shared library via juliac.
#
# Usage (from the Kora.jl root):
#   ./build/build.sh [--output-dir <dir>]
#
# Overridable env vars:
#   KORA_LIB_DIR   output directory (default: julia_lib/ inside Kora.jl root)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
ENTRY_FILE="$SCRIPT_DIR/bridge_aot.jl"

OUTPUT_DIR="${KORA_LIB_DIR:-}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

OUTPUT_DIR="${OUTPUT_DIR:-$PROJECT_ROOT/build/dist}"

mkdir -p "$OUTPUT_DIR"

juliac --project="$PROJECT_ROOT" --output-lib "$OUTPUT_DIR/kora_bridge" \
    --trim=safe --compile-ccallable --experimental "$ENTRY_FILE"
