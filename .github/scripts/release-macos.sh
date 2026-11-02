#!/usr/bin/env bash
# Build a relocatable macOS BioSpring distribution (MDDriver + FreeSASA + OpenMP) for the
# architecture of the machine it runs on (arm64 or x86_64), bundle every non-system dylib,
# sign it, verify it, and pack it into dist/.
#
# Needs Homebrew: netcdf netcdf-cxx libomp automake autoconf libtool dylibbundler.
#
# Environment (all optional):
#   MDDRIVER_REF         MDDriver commit to build            (default: pinned below)
#   FREESASA_TAG         FreeSASA tag to build               (default: 2.1.2)
#   WORK                 scratch directory                   (default: ./build-release)
#   MACOS_SIGN_IDENTITY  Developer ID Application identity   (default: ad-hoc signature)
#   SKIP_TESTS=1         skip ctest
set -euo pipefail

MDDRIVER_REPO="https://github.com/LBT-CNRS/MDDriver.git"
MDDRIVER_REF="${MDDRIVER_REF:-67e1f56696c67b1bd2f87e9f789dd930f5b5ed16}"
FREESASA_REPO="https://github.com/mittinatten/freesasa.git"
FREESASA_TAG="${FREESASA_TAG:-2.1.2}"

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
WORK="${WORK:-$PWD/build-release}"
DEPS="$WORK/deps"
BREW="$(brew --prefix)"
LIBOMP="$(brew --prefix libomp)"
ARCH="$(uname -m)"
VERSION="$(git -C "$SRC" describe --tags --always --dirty)"
NAME="biospring-${VERSION}-macos-${ARCH}"
STAGE="$WORK/$NAME"
JOBS="$(sysctl -n hw.ncpu)"

rm -rf "$WORK"
mkdir -p "$DEPS" "$SRC/dist"

git clone --depth 1 --branch "$FREESASA_TAG" --recurse-submodules --shallow-submodules \
    "$FREESASA_REPO" "$WORK/freesasa"
(
    cd "$WORK/freesasa"
    autoreconf -i
    ./configure --prefix="$DEPS" --disable-json --disable-xml
    make -j"$JOBS"
    make install
)

git init -q "$WORK/mddriver"
git -C "$WORK/mddriver" fetch -q --depth 1 "$MDDRIVER_REPO" "$MDDRIVER_REF"
git -C "$WORK/mddriver" checkout -q FETCH_HEAD
cmake -S "$WORK/mddriver" -B "$WORK/mddriver-build" \
    -DCMAKE_INSTALL_PREFIX="$DEPS" -DCMAKE_BUILD_TYPE=Release
cmake --build "$WORK/mddriver-build" --parallel "$JOBS"
cmake --install "$WORK/mddriver-build"

# Apple clang has no built-in OpenMP: point it at Homebrew's libomp.
omp_flags="-Xpreprocessor -fopenmp -I$LIBOMP/include"
# Machines that also have MacPorts or an old /opt/lib carry x86_64 libraries there (NetCDF, gtest)
# that the Find modules would otherwise pick up and fail to link on arm64.
cmake -S "$SRC" -B "$WORK/build" \
    -DCMAKE_IGNORE_PREFIX_PATH=/opt/local -DCMAKE_IGNORE_PATH="/opt/lib;/opt/include" \
    -DCMAKE_INSTALL_PREFIX="$STAGE" -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON -DMDDRIVER_SUPPORT=ON -DFREESASA_SUPPORT=ON -DOPENMP_SUPPORT=ON \
    -DMDDriver_DIR="$DEPS/share/cmake" -DFREESASA_PREFIX="$DEPS" \
    -DCMAKE_CXX_FLAGS="$omp_flags" -DCMAKE_C_FLAGS="$omp_flags" \
    -DOpenMP_CXX_FLAGS="$omp_flags" -DOpenMP_C_FLAGS="$omp_flags" \
    -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_C_LIB_NAMES=omp \
    -DOpenMP_omp_LIBRARY="$LIBOMP/lib/libomp.dylib" \
    -DNetCDF_C_INCLUDE_DIR="$BREW/include" -DNetCDF_C_LIBRARY="$BREW/lib/libnetcdf.dylib" \
    -DNetCDF_CXX_INCLUDE_DIR="$BREW/include" -DNetCDF_CXX_LIBRARY="$BREW/lib/libnetcdf-cxx4.dylib" \
    -DNetCDF_INCLUDE_DIR="$BREW/include" -DNetCDF_LIBRARY="$BREW/lib/libnetcdf.dylib" \
    -DNetCDFCXX_INCLUDE_DIR="$BREW/include" -DNetCDFCXX_LIBRARY="$BREW/lib/libnetcdf-cxx4.dylib"
cmake --build "$WORK/build" --parallel "$JOBS"
if [[ -z "${SKIP_TESTS:-}" ]]; then
    ctest --test-dir "$WORK/build" --output-on-failure
fi
cmake --install "$WORK/build"

# Executables only: the *-config.sh helpers are scripts.
exes=()
fix_args=()
for f in "$STAGE"/bin/*; do
    if file "$f" | grep -q "Mach-O.*executable"; then
        exes+=("$f")
        fix_args+=(-x "$f")
    fi
done
dylibbundler -b "${fix_args[@]}" -d "$STAGE/lib/bundled" -p @executable_path/../lib/bundled \
    -of -cd -ns -s "$DEPS/lib" -s "$LIBOMP/lib"

sign_identity="${MACOS_SIGN_IDENTITY:--}"
sign_flags=(--force --sign "$sign_identity")
if [[ "$sign_identity" != "-" ]]; then
    sign_flags+=(--options runtime --timestamp)
fi
for f in "$STAGE"/lib/bundled/*.dylib "${exes[@]}"; do
    codesign "${sign_flags[@]}" "$f"
done

# Verification: the bundle must not reference anything outside the system or itself.
fail=0
for f in "$STAGE"/lib/bundled/*.dylib "${exes[@]}"; do
    codesign --verify --strict "$f" || fail=1
    archs="$(lipo -archs "$f")"
    [[ "$archs" == "$ARCH" ]] || { echo "ERROR: $f is $archs, expected $ARCH"; fail=1; }
    if otool -L "$f" | tail -n +2 | awk '{print $1}' \
        | grep -Ev '^(/usr/lib/|/System/|@executable_path/|@loader_path/|@rpath/)'; then
        echo "ERROR: $f references libraries outside the bundle (listed above)"
        fail=1
    fi
done
[[ "$fail" -eq 0 ]] || exit 1

help_text="$(env -i "$STAGE/bin/biospring" --help 2>&1 || true)"
for flag in --port --wait --sasa; do
    grep -q -- "$flag" <<<"$help_text" || { echo "ERROR: biospring --help lacks $flag"; exit 1; }
done

mkdir -p "$SRC/dist"
tar -C "$WORK" -czf "$SRC/dist/$NAME.tar.gz" "$NAME"
(cd "$SRC/dist" && shasum -a 256 "$NAME.tar.gz" > "$NAME.tar.gz.sha256")
echo "built $SRC/dist/$NAME.tar.gz ($(cat "$SRC/dist/$NAME.tar.gz.sha256"))"
