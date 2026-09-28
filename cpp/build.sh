#!/usr/bin/env bash
# Build the C++ kernel (qef_kernel.cpp) as the shared library that
# src/qef/kernel.py loads through ctypes:
#   cpp/build/libqef_kernel.dylib on macOS, cpp/build/libqef_kernel.so elsewhere.
#
# The compiler is $CXX if set, otherwise c++ (on macOS /usr/bin/c++, Apple
# clang). -ffp-contract=off keeps every product and sum separately rounded,
# as in the Python reference; -ffast-math must never be added. The library is
# written to a temporary name and renamed, so a concurrent build never leaves
# a partial file behind.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CXX="${CXX:-c++}"
case "$(uname -s)" in
    Darwin) EXT=dylib ;;
    *) EXT=so ;;
esac
OUT="$HERE/build/libqef_kernel.$EXT"
TMP="$HERE/build/.libqef_kernel.$$.$EXT"

mkdir -p "$HERE/build"
FLAGS=(-O2 -std=c++17 -fPIC -shared -ffp-contract=off -Wall -Wextra)
echo "compiler: $("$CXX" --version 2>/dev/null | head -n 1)"
echo "command: $CXX ${FLAGS[*]} -o cpp/build/libqef_kernel.$EXT cpp/qef_kernel.cpp"
"$CXX" "${FLAGS[@]}" -o "$TMP" "$HERE/qef_kernel.cpp"
mv -f "$TMP" "$OUT"
echo "built: $OUT"
