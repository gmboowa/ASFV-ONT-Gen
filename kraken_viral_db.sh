#!/bin/bash

# Kraken2 Viral Database Installer for Home Directory
# Downloads and installs the pre-built viral database to $HOME

set -e  # Exit on error

# Configuration
DB_NAME="kraken2_viral_db"
DB_DIR="$HOME/$DB_NAME"
DB_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20230605.tar.gz"
TAR_FILE="k2_viral_20230605.tar.gz"
EXPECTED_MD5="a9c8c3c5e5a6f5b5a9c8c3c5e5a6f5b5"  # Example MD5 - replace with actual

# Check if Kraken2 is installed
if ! command -v kraken2 &> /dev/null; then
    echo "Error: Kraken2 not found in PATH. Please install Kraken2 first."
    exit 1
fi

# Check if database already exists
if [ -d "$DB_DIR" ]; then
    echo "Database directory already exists at: $DB_DIR"
    echo "Would you like to overwrite it? [y/N]"
    read -r response
    if [[ ! "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
        echo "Installation cancelled."
        exit 0
    fi
    echo "Removing existing database..."
    rm -rf "$DB_DIR"
fi

# Create database directory in home
echo "Creating database directory: $DB_DIR"
mkdir -p "$DB_DIR"
cd "$DB_DIR"

# Download database
echo "Downloading Kraken2 viral database (this may take a while)..."
wget "$DB_URL"

# Verify download integrity
echo "Verifying download integrity..."
if ! md5sum -c <<<"$EXPECTED_MD5  $TAR_FILE"; then
    echo "Error: MD5 checksum verification failed!"
    exit 1
fi

# Extract database
echo "Extracting database..."
tar -xvzf "$TAR_FILE" --strip-components=1

# Clean up
echo "Cleaning up..."
rm "$TAR_FILE"

# Verify essential files
ESSENTIAL_FILES=("hash.k2d" "opts.k2d" "taxo.k2d")
for file in "${ESSENTIAL_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Error: Essential file $file missing after extraction!"
        exit 1
    fi
done

echo "✅ Kraken2 viral database successfully installed in: $DB_DIR"
echo ""
echo "To use this database, you can:"
echo "1. Add to your environment:"
echo "   export KRAKEN2_DB_PATH=\"$DB_DIR\""
echo "2. Then use with:"
echo "   kraken2 --db \$KRAKEN2_DB_PATH [YOUR_OPTIONS]"
echo ""
echo "Or reference directly:"
echo "   kraken2 --db $DB_DIR [YOUR_OPTIONS]"
