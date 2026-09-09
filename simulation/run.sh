#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
BUILD_DIR="$SCRIPT_DIR/build"

# 只需要在这里选择配置文件。
CONFIG_FILE="$SCRIPT_DIR/config/apv25_1.jsonc"

cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR"
cmake --build "$BUILD_DIR" -j

cd "$BUILD_DIR"
./microTPCsimulation "$CONFIG_FILE"
./microTPCsample "$CONFIG_FILE"

# JSON 允许整行注释，因此这里直接取 outputDir 的字符串值。
OUTPUT_DIR=$(sed -nE \
    's/^[[:space:]]*"outputDir"[[:space:]]*:[[:space:]]*"([^"]+)".*/\1/p' \
    "$CONFIG_FILE" | head -n 1)

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: outputDir is missing in $CONFIG_FILE" >&2
    exit 1
fi

TIER_DIR="$OUTPUT_DIR/tiers"
MEAN_DIR="$OUTPUT_DIR/mean"
mkdir -p "$MEAN_DIR"
tier=0
while [[ -f "$TIER_DIR/tier${tier}.root" ]]; do
    root -l -b -q \
        "$SCRIPT_DIR/analysis/mean.C(${tier},\"${TIER_DIR}\",\"${MEAN_DIR}\")"
    tier=$((tier + 1))
done

if [[ "$tier" -eq 0 ]]; then
    echo "Error: no tier files were generated in $TIER_DIR" >&2
    exit 1
fi

root -l -b -q \
    -e ".L $SCRIPT_DIR/analysis/matix.C" \
    -e "matrix(\"$MEAN_DIR\", \"$OUTPUT_DIR/matrix.root\");"

echo "Done: generated $tier tiers in $OUTPUT_DIR"
