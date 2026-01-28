#!/bin/bash
set -e

# Determine if sudo is needed
if [ "$(id -u)" -eq 0 ]; then
    SUDO=""
else
    SUDO="sudo"
fi

QUANTLIB_VERSION="1.41"
QUANTLIB_FILENAME="QuantLib-${QUANTLIB_VERSION}.tar.gz"
QUANTLIB_URL="https://github.com/lballabio/QuantLib/releases/download/v${QUANTLIB_VERSION}/${QUANTLIB_FILENAME}"

CACHE_DIR="/host/quantlib-cache"

# Check Cache
if [ -f "${CACHE_DIR}/lib/libQuantLib.a" ]; then
    echo "Restoring QuantLib from cache..."
    $SUDO cp -r ${CACHE_DIR}/* /usr/local/
    echo "QuantLib restored from cache."
    exit 0
fi

echo "QuantLib cache not found."
echo "Building QuantLib ${QUANTLIB_VERSION} from source..."

# Download
wget -q -O "${QUANTLIB_FILENAME}" "${QUANTLIB_URL}"
tar xzf "${QUANTLIB_FILENAME}"

# Build
cd "QuantLib-${QUANTLIB_VERSION}"
./configure --disable-shared --prefix=/usr/local --enable-std-classes --with-pic CXXFLAGS="-O3 -fPIC"
make -j$(nproc)
$SUDO make install

# Save to Cache
echo "Saving QuantLib to cache..."
$SUDO mkdir -p ${CACHE_DIR}/include ${CACHE_DIR}/lib
$SUDO cp -r /usr/local/include/ql ${CACHE_DIR}/include/
$SUDO cp -r /usr/local/lib/libQuantLib* ${CACHE_DIR}/lib/

# Cleanup
cd ..
rm -rf "QuantLib-${QUANTLIB_VERSION}" "${QUANTLIB_FILENAME}"
echo "QuantLib installation complete."
