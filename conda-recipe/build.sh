#!/bin/bash
set -euo pipefail

# ECCsplorer_slim build script
# -----------------------------
# Copies standalone Python and R scripts from the bin/ directory
# into the conda environment's bin/ directory so they are on PATH.

mkdir -p ${PREFIX}/bin

# Copy all Python scripts
cp -r ${SRC_DIR}/*.py ${PREFIX}/bin/ 2>/dev/null || true

# Copy all R scripts
cp -r ${SRC_DIR}/*.R ${PREFIX}/bin/ 2>/dev/null || true

# Make them executable
chmod +x ${PREFIX}/bin/*.py ${PREFIX}/bin/*.R 2>/dev/null || true

echo "ECCsplorer_slim scripts installed to ${PREFIX}/bin"
