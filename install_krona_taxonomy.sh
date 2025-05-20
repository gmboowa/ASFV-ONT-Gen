#!/bin/bash

# Detect the active Conda environment prefix
if [ -z "$CONDA_PREFIX" ]; then
  echo "❌ No active Conda environment detected. Please activate a Conda environment with Krona installed."
  exit 1
fi

# Locate Krona installation within the conda environment
KRONA_DIR="$CONDA_PREFIX/opt/krona"
UPDATE_SCRIPT="$KRONA_DIR/updateTaxonomy.sh"

# Ensure the Krona directory and update script exist
if [ ! -d "$KRONA_DIR" ]; then
  echo "❌ Krona directory not found at $KRONA_DIR"
  exit 1
fi

if [ ! -x "$UPDATE_SCRIPT" ]; then
  echo "❌ updateTaxonomy.sh not found or not executable at $UPDATE_SCRIPT"
  exit 1
fi

# Create the taxonomy directory if it doesn't exist
mkdir -p "$KRONA_DIR/taxonomy"

# Move into the Krona directory
cd "$KRONA_DIR" || {
  echo "❌ Failed to change directory to $KRONA_DIR"
  exit 1
}

# Run the update script
echo "🔄 Running updateTaxonomy.sh to download taxonomy data..."
./updateTaxonomy.sh

# Check if the update succeeded
if [ $? -eq 0 ]; then
  echo "✅ Krona taxonomy database installed successfully in: $KRONA_DIR/taxonomy"
else
  echo "❌ Failed to update Krona taxonomy database."
  exit 1
fi