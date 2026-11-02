#!/usr/bin/env bash
# Build a relocatable Linux BioSpring distribution (MDDriver + FreeSASA + OpenMP) for the
# architecture of the machine it runs on (x86_64 or aarch64), bundle every non-glibc shared
# library, verify it, and pack it into dist/.
#
# Needs (Debian/Ubuntu): build-essential cmake git autoconf automake libtool pkg-config
#                        libnetcdf-dev libnetcdf-c++4-dev patchelf file
# Build on the oldest distribution you want to support: the glibc of the build machine is the
# minimum glibc of the result.
#
# Environment (all optional):
#   MDDRIVER_REF   MDDriver commit to build   (default: pinned below)
#   FREESASA_TAG   FreeSASA tag to build      (default: 2.1.2)
#   WORK           scratch directory          (default: ./build-release)
#   SKIP_TESTS=1   skip ctest
set -euo pipefail

MDDRIVER_REPO="https://github.com/LBT-CNRS/MDDriver.git"
MDDRIVER_REF="${MDDRIVER_REF:-67e1f56696c67b1bd2f87e9f789dd930f5b5ed16}"
FREESASA_REPO="https://github.com/mittinatten/freesasa.git"
FREESASA_TAG="${FREESASA_TAG:-2.1.2}"

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
WORK="${WORK:-$PWD/build-release}"
DEPS="$WORK/deps"
ARCH="$(uname -m)"
VERSION="$(git -C "$SRC" describe --tags --always --dirty)"
NAME="biospring-${VERSION}-linux-${ARCH}"
STAGE="$WORK/$NAME"
LIBDIR="$STAGE/lib/bundled"
JOBS="$(nproc)"

rm -rf "$WORK"
mkdir -p "$DEPS" "$SRC/dist"

git clone --depth 1 --branch "$FREESASA_TAG" --recurse-submodules --shallow-submodules \
    "$FREESASA_REPO" "$WORK/freesasa"
(
    cd "$WORK/freesasa"
    autoreconf -i
    ./configure --prefix="$DEPS" --disable-json --disable-xml
    # FreeSASA hard-codes -lc++ for its command-line tool, which BioSpring does not use.
    make -j"$JOBS" freesasa_LDADD=libfreesasa.a
    make install freesasa_LDADD=libfreesasa.a
)

git init -q "$WORK/mddriver"
git -C "$WORK/mddriver" fetch -q --depth 1 "$MDDRIVER_REPO" "$MDDRIVER_REF"
git -C "$WORK/mddriver" checkout -q FETCH_HEAD
cmake -S "$WORK/mddriver" -B "$WORK/mddriver-build" \
    -DCMAKE_INSTALL_PREFIX="$DEPS" -DCMAKE_BUILD_TYPE=Release
cmake --build "$WORK/mddriver-build" --parallel "$JOBS"
cmake --install "$WORK/mddriver-build"

cmake -S "$SRC" -B "$WORK/build" \
    -DCMAKE_INSTALL_PREFIX="$STAGE" -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON -DMDDRIVER_SUPPORT=ON -DFREESASA_SUPPORT=ON -DOPENMP_SUPPORT=ON \
    -DMDDriver_DIR="$DEPS/share/cmake" -DFREESASA_PREFIX="$DEPS"
cmake --build "$WORK/build" --parallel "$JOBS"
if [[ -z "${SKIP_TESTS:-}" ]]; then
    ctest --test-dir "$WORK/build" --output-on-failure
fi
cmake --install "$WORK/build"

# Executables only: the *-config.sh helpers are scripts.
exes=()
for f in "$STAGE"/bin/*; do
    if file "$f" | grep -q "ELF.*executable"; then
        exes+=("$f")
    fi
done

# Libraries that must come from the host: the glibc family and the dynamic loader.
is_host_lib() {
    case "$1" in
        linux-vdso.so.*|ld-linux*|libc.so.*|libm.so.*|libdl.so.*|libpthread.so.*|librt.so.*|libutil.so.*|libresolv.so.*) return 0 ;;
    esac
    return 1
}

mkdir -p "$LIBDIR"
pending=("${exes[@]}")
: > "$WORK/bundled.txt"
while [[ ${#pending[@]} -gt 0 ]]; do
    current="${pending[0]}"
    pending=("${pending[@]:1}")
    while read -r name arrow path _; do
        [[ "$arrow" == "=>" && "$path" == /* ]] || continue
        is_host_lib "$name" && continue
        grep -qxF "$name" "$WORK/bundled.txt" && continue
        echo "$name" >> "$WORK/bundled.txt"
        cp -L "$path" "$LIBDIR/$name"
        chmod u+w "$LIBDIR/$name"
        pending+=("$LIBDIR/$name")
    done < <(LD_LIBRARY_PATH="$DEPS/lib" ldd "$current")
done

# shellcheck disable=SC2016 # $ORIGIN is for the dynamic loader, not for the shell
for exe in "${exes[@]}"; do
    patchelf --set-rpath '$ORIGIN/../lib/bundled' "$exe"
done
# shellcheck disable=SC2016
for lib in "$LIBDIR"/*; do
    patchelf --set-rpath '$ORIGIN' "$lib"
done

# Verification: every dependency must resolve inside the bundle or be a host (glibc) library.
fail=0
for f in "${exes[@]}" "$LIBDIR"/*; do
    while read -r name arrow path _; do
        [[ "$arrow" == "=>" ]] || continue
        if [[ "$path" == "not" ]]; then
            echo "ERROR: $f: $name not found"
            fail=1
        elif [[ "$path" == /* && "$path" != "$STAGE"/* ]] && ! is_host_lib "$name"; then
            echo "ERROR: $f: $name resolves to $path, outside the bundle"
            fail=1
        fi
    done < <(env -i ldd "$f")
done
[[ "$fail" -eq 0 ]] || exit 1

help_text="$(env -i "$STAGE/bin/biospring" --help 2>&1 || true)"
for flag in --port --wait --sasa; do
    grep -q -- "$flag" <<<"$help_text" || { echo "ERROR: biospring --help lacks $flag"; exit 1; }
done

glibc_floor="$(objdump -T "${exes[@]}" "$LIBDIR"/* 2>/dev/null | grep -o 'GLIBC_[0-9.]*' | sort -uV | tail -n 1)"
echo "minimum glibc required: ${glibc_floor#GLIBC_}"

tar -C "$WORK" -czf "$SRC/dist/$NAME.tar.gz" "$NAME"
(cd "$SRC/dist" && sha256sum "$NAME.tar.gz" > "$NAME.tar.gz.sha256")
echo "built $SRC/dist/$NAME.tar.gz ($(cat "$SRC/dist/$NAME.tar.gz.sha256"))"
