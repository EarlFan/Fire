#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

# Enter the directory of the script (i.e. the case directory)
cd "$(dirname "$0")"

# Parse parameters
BUILD_TYPE="Release"
JOBS=8
RUN_AFTER_BUILD=true

while [[ $# -gt 0 ]]; do
    case $1 in
        --debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        -j|--jobs)
            JOBS="$2"
            shift; shift
            ;;
        --no-run)
            RUN_AFTER_BUILD=false
            shift
            ;;
        *)
            echo "Usage: $0 [--debug] [-j JOBS] [--no-run]"
            exit 1
            ;;
    esac
done

# Clean old compilation files and executables
echo "[1/4] Cleaning..."
rm -rf build Fire_*Species

# Compile the project
echo "[2/4] Building ($BUILD_TYPE)..."
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE ..
make -j$JOBS
cd ..

source config.sh

# Locate the generated executable (according to CMake rules)
EXECUTABLE="Fire_${BUILD_TYPE}_${NS}Species"
SOURCE_PATH="build/${EXECUTABLE}"  # Correct path to bin directory

if [[ -f "$SOURCE_PATH" ]]; then
    echo "[3/4] Copying executable..."
    cp "$SOURCE_PATH" .  # Correct copy command format
else
    echo "Error: Executable not found at $SOURCE_PATH!"
    echo "Available files in build/bin:"
    ls -l build/bin || true
    exit 1
fi

# Run the case
if $RUN_AFTER_BUILD; then
    echo "[4/4] Running $EXECUTABLE..."
    ./$EXECUTABLE

    echo "Plotting 1D detonation results in the case output dir ..."
    cp plot_1D_deton.py results/Detonation_1D/
    cd results/Detonation_1D/
    python plot_1D_deton.py
    cd ../../
else
    echo "[4/4] Skip running."
    echo "The case can be run with \"./$EXECUTABLE\""
fi

echo "✅ All done."
