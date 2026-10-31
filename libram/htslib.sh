#!/bin/sh
# Builds htslib and samtools 1.22.1 that read RAM files, with libram and
# htslib-ram.patch, so that samtools view, flagstat and the rest take a .ram
# file, and a region of one, as they take a BAM.
# Usage: htslib.sh <empty directory>; set CPPFLAGS if xxhash.h is elsewhere.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$1"
cd "$1"
git clone -q --depth 1 --branch 1.22.1 --recurse-submodules --shallow-submodules https://github.com/samtools/htslib.git
git clone -q --depth 1 --branch 1.22.1 https://github.com/samtools/samtools.git
cd htslib
git apply "$here/htslib-ram.patch"
mkdir ram
cp "$here/ram.c" "$here/ram.h" "$here/rntuple.c" "$here/rntuple.h" ram/
make -j4 lib-static htslib_static.mk CPPFLAGS="${CPPFLAGS:-}"
cd ../samtools
autoheader
autoconf -Wno-syntax
./configure --with-htslib=../htslib --without-curses LIBS=-lzstd
make -j4 samtools
echo "built $(pwd)/samtools"
