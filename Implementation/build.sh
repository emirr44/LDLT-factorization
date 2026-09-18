#!/usr/bin/env bash
# build.sh - build the ldlt command line tool (macOS, Linux)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
OUT="$ROOT/build"
GEN="$OUT/generated"

if [ "${1:-}" = "--clean" ]; then
    rm -rf "$OUT"
fi
mkdir -p "$GEN"

CC="${CC:-cc}"
CXX="${CXX:-c++}"

printf '#pragma once\n#define LDLT_PROJECT_ROOT "%s"\n' "$ROOT" > "$GEN/ldlt_config.hpp"

if [ ! -f "$OUT/sqlite3.o" ]; then
    echo "compiling sqlite (one time, about half a minute)..."
    "$CC" -O2 -DSQLITE_THREADSAFE=0 -DSQLITE_OMIT_LOAD_EXTENSION \
        -c "$ROOT/third_party/sqlite/sqlite3.c" -o "$OUT/sqlite3.o"
fi

echo "compiling ldlt..."
"$CXX" -std=c++17 -O2 -DNDEBUG -Wall -Wextra \
    -I"$ROOT/src/core" -I"$ROOT/third_party/sqlite" -I"$GEN" \
    "$ROOT"/src/core/*.cpp "$ROOT"/src/cli/*.cpp "$OUT/sqlite3.o" \
    -o "$OUT/ldlt"

echo "build ok -> $OUT/ldlt"
echo "run it with:  ./build/ldlt"
