#!/bin/bash
# scripts/configure_dev.sh
# Configures CMake for local development ensuring the active Poetry Python environment is used.

# 1. Get the Python executable from Poetry
POETRY_PYTHON=$(poetry run which python)
echo "Using Python: $POETRY_PYTHON"

# 2. Get the pybind11 CMake directory
PYBIND11_CMAKE=$(poetry run python -c 'import pybind11; print(pybind11.get_cmake_dir())')
echo "Using pybind11: $PYBIND11_CMAKE"

# 3. Configure CMake
# We explicitly set Python3_EXECUTABLE to ensure FindPython3 picks the right one.
# We set pybind11_DIR to ensure it is found even if not in standard paths.
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Debug \
    -DPython3_EXECUTABLE="$POETRY_PYTHON" \
    -Dpybind11_DIR="$PYBIND11_CMAKE"

echo "Configuration complete. Run 'cmake --build build' to build."
