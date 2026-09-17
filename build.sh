#!/usr/bin/env bash
# build.sh - build the ldlt command line tool (macOS, Linux).
#
#   ./build.sh            release build   ->  build/ldlt
#   ./build.sh --clean    rebuild everything, including SQLite
#
# Needs a C and a C++17 compiler. On macOS that is Xcode's command line tools:
#   xcode-select --install
# Windows: use .\build.ps1 instead.

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

# Record where the project lives, so ldlt finds data/ from any directory.
printf '#pragma once\n#define LDLT_PROJECT_ROOT "%s"\n' "$ROOT" > "$GEN/ldlt_config.hpp"

# SQLite is plain C and large: compile it once and reuse the object file.
if [ ! -f "$OUT/sqlite3.o" ]; then
    echo "compiling sqlite (one time, about half a minute)..."
    "$CC" -O2 -DSQLITE_THREADSAFE=0 -DSQLITE_OMIT_LOAD_EXTENSION \
        -c "$ROOT/third_party/sqlite/sqlite3.c" -o "$OUT/sqlite3.o"
fi

echo "compiling ldlt..."
"$CXX" -std=c++17 -O2 -DNDEBUG -Wall -Wextra \
    -I"$ROOT/include" -I"$ROOT/third_party/sqlite" -I"$GEN" \
    "$ROOT"/src/*.cpp "$ROOT"/apps/cli/*.cpp "$OUT/sqlite3.o" \
    -o "$OUT/ldlt"

echo "build ok -> $OUT/ldlt"
echo "run it with:  ./build/ldlt"
