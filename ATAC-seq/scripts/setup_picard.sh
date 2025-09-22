#!/bin/bash
# setup_picard.sh - Download Picard tools

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

TOOLS_DIR="tools"
PICARD_VERSION="3.1.1"
PICARD_URL="https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar"

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

echo "Setting up Picard tools..."

# Create tools directory
mkdir -p "$TOOLS_DIR"

# Download Picard if not already present
if [[ ! -f "$TOOLS_DIR/picard.jar" ]]; then
    echo "Downloading Picard version $PICARD_VERSION..."
    wget "$PICARD_URL" -O "$TOOLS_DIR/picard.jar"
else
    echo "Picard already exists"
fi

echo "Picard setup complete!"