#!/bin/bash
# scripts/clean_all.sh

# Remove CMake in-source build artifacts
rm -rf CMakeFiles
rm -f CMakeCache.txt
rm -f Makefile
rm -f cmake_install.cmake
rm -f CTestTestfile.cmake
rm -rf Testing
rm -rf _deps

# Remove build directories
rm -rf build
rm -rf build_core
rm -rf build_coverage
rm -rf build_debug
rm -rf build_tests
rm -rf bin
rm -rf lib

# Remove Python artifacts
rm -rf .mypy_cache
rm -rf .pytest_cache
rm -rf .ruff_cache
rm -rf dist
rm -rf *.egg-info
find . -name "*.so" -delete
find . -name "__pycache__" -type d -exec rm -rf {} +

echo "Clean complete."
