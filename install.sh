#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status.
set -e

export MOLPERT="$(pwd)/Molpert"
export MOLECULE_AUTO_CORRECT="$(pwd)"

# Clone the Molpert repository dependency
# git clone https://github.com/AlanKerstjens/Molpert.git

# Install molpert
mkdir ${MOLPERT}/build && cd ${MOLPERT}/build
cmake ..
make install

# Create the build directory and navigate into it
mkdir -p "${MOLECULE_AUTO_CORRECT}/build"
cd "${MOLECULE_AUTO_CORRECT}/build"

# Run CMake with the given include directory (Ensure $MOLPERT is set before running this script)
cmake -DMolpert_INCLUDE_DIRS="${MOLPERT}/source" ..

# Build and install MoleculeAutoCorrect
make install

echo "Adding repos to python path"
export PYTHONPATH="${PYTHONPATH}:${MOLECULE_AUTO_CORRECT}/lib"
export PYTHONPATH="${PYTHONPATH}:${MOLPERT}/lib"

echo "Installation complete"
