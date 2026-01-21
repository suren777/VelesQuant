#!/bin/bash
set -e

# Define Boost version and paths
BOOST_VERSION="1.83.0"
BOOST_VERSION_UNDERSCORE="1_83_0"
BOOST_FILENAME="boost_${BOOST_VERSION_UNDERSCORE}.tar.gz"
BOOST_URL="https://boostorg.jfrog.io/artifactory/main/release/${BOOST_VERSION}/source/${BOOST_FILENAME}"
BOOST_URL_BACKUP="https://archives.boost.io/release/${BOOST_VERSION}/source/${BOOST_FILENAME}"

CACHE_DIR="/host/boost-cache"

# Check Cache
if [ -f "${CACHE_DIR}/include/boost/version.hpp" ]; then
    echo "Restoring Boost from cache..."
    cp -r ${CACHE_DIR}/* /usr/local/
    exit 0
fi

echo "Boost cache not found or incomplete."
echo "Building Boost ${BOOST_VERSION} from source..."

# robust download with retries
MAX_ATTEMPTS=3
attempts=0
success=0

while [ $attempts -lt $MAX_ATTEMPTS ]; do
    echo "Downloading Boost (Attempt $((attempts + 1))/${MAX_ATTEMPTS})..."
    
    # Try primary URL
    wget -q -O "${BOOST_FILENAME}" "${BOOST_URL}"

    # Check if primary is valid
    if ! gzip -t "${BOOST_FILENAME}" 2>/dev/null; then
        echo "Primary URL download invalid or failed. Trying backup..."
        rm -f "${BOOST_FILENAME}"
        wget -q -O "${BOOST_FILENAME}" "${BOOST_URL_BACKUP}"
    fi
    
    # Verify the download is a valid gzip file
    if gzip -t "${BOOST_FILENAME}" 2>/dev/null; then
        echo "Download verified successfully."
        success=1
        break
    else
        echo "Download invalid (gzip check failed)."
        rm -f "${BOOST_FILENAME}"
    fi
    
    attempts=$((attempts + 1))
    echo "Retrying in 5 seconds..."
    sleep 5
done

if [ $success -ne 1 ]; then
    echo "Failed to download and verify ${BOOST_FILENAME} after ${MAX_ATTEMPTS} attempts."
    echo "Dumping file content (head) for debugging:"
    head -n 20 "${BOOST_FILENAME}" || true
    exit 1
fi

# Extract
echo "Extracting Boost..."
tar xzf "${BOOST_FILENAME}" || { echo "Failed to extract Boost archive"; exit 1; }

# Install
cd "boost_${BOOST_VERSION_UNDERSCORE}"
echo "Bootstrapping..."
./bootstrap.sh
echo "Installing..."
./b2 install --prefix=/usr/local

# Save to Cache
echo "Saving Boost to cache..."
mkdir -p ${CACHE_DIR}/include ${CACHE_DIR}/lib
cp -r /usr/local/include/boost ${CACHE_DIR}/include/
cp -r /usr/local/lib/libboost* ${CACHE_DIR}/lib/

# Cleanup
cd ..
rm -rf "boost_${BOOST_VERSION_UNDERSCORE}" "${BOOST_FILENAME}"
echo "Boost installation complete."
