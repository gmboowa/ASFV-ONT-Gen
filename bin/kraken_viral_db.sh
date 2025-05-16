#!/bin/bash
# setup_kraken_viral_db.sh - Install Kraken2 viral database to user's home

# Check if kraken2 is installed
if ! command -v kraken2 &> /dev/null; then
    echo "Error: kraken2 not found in PATH. Please install Kraken2 first."
    exit 1
fi

# Set default database location
DEFAULT_DB_DIR="$HOME/kraken2_viral_db"

# Allow custom install location through first argument
DB_DIR="${1:-$DEFAULT_DB_DIR}"

# Create database directory
echo "Installing Kraken2 viral database to: $DB_DIR"
mkdir -p "$DB_DIR"

# Check if database directory exists and is empty
if [ "$(ls -A "$DB_DIR")" ]; then
    echo "Warning: Database directory is not empty!"
    read -p "Delete existing contents and reinstall? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "${DB_DIR:?}/"*
    else
        echo "Aborting installation."
        exit 1
    fi
fi

# Get available threads (use 50% of total)
THREADS=$(($(nproc)/2))
THREADS=${THREADS:-4} # Fallback to 4 threads if nproc fails

echo "Starting database installation with $THREADS threads..."

# Build database
set -euo pipefail  # Enable strict error checking

kraken2-build --download-taxonomy --db "$DB_DIR"
kraken2-build --download-library viral --db "$DB_DIR"
kraken2-build --build --threads "$THREADS" --db "$DB_DIR"
kraken2-build --clean --db "$DB_DIR"

# Verify installation
if [ -f "$DB_DIR/hash.k2d" ] && [ -f "$DB_DIR/taxo.k2d" ]; then
    echo -e "\nSuccessfully installed Kraken2 viral database!"
    echo "Database location: $DB_DIR"
    echo "To use: kraken2 --db $DB_DIR [...]"
else
    echo -e "\nError: Database files not found! Installation may have failed."
    exit 1
fi