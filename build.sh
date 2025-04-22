#!/bin/bash

# Simple build script for AIHMM using CMake

# Exit on error
set -e

# Create and navigate to build directory
echo "Creating build directory..."
mkdir -p build
cd build

# Configure with CMake
echo "Configuring with CMake..."
cmake ..

# Build
echo "Building AIHMM..."
make -j$(nproc)

echo "Build completed successfully!"
echo "The AIHMM executable can be found at: $(pwd)/bin/AIHMM"
echo ""
echo "Run with: ./bin/AIHMM --numChromosomes 22 --numIndividuals 55 --startFile ceu8StateHMMParams.txt --fileNamePrefix ../data/ceuRevisedRatiosChrom --outputFile ../output/outputFile.txt" 