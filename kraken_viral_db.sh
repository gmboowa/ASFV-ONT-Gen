#!/bin/bash

# Activate your conda environment first (optional safety check)
# conda activate AFSV_ont

echo "🔍 Checking if snpEff is already installed..."
if command -v snpEff &> /dev/null; then
    echo "✅ snpEff is already installed and available in the environment."
else
    echo "⏳ Installing snpEff via conda..."
    conda install -y -c bioconda -c conda-forge snpeff

    if command -v snpEff &> /dev/null; then
        echo "✅ snpEff installation completed and is now available."
    else
        echo "❌ snpEff installation failed. Please check conda environment or install manually."
        exit 1
    fi
fi

# Confirm where snpEff.jar is located
echo "🔍 Searching for snpEff.jar..."
JAR_PATH=$(find "$CONDA_PREFIX" -type f -name "snpEff.jar" | head -n 1)

if [ -n "$JAR_PATH" ]; then
    echo "✅ snpEff.jar found at: $JAR_PATH"
else
    echo "❌ snpEff.jar not found. The pipeline may fail to locate it. Check your installation."
    exit 1
fi