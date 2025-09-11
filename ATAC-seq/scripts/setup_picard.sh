#!/bin/bash
# Download Picard tools

# Configuration
TOOLS_DIR="tools"
PICARD_VERSION="3.1.1"  # Latest stable version
PICARD_URL="https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar"

echo "Setting up Picard tools..."

# Create tools directory
mkdir -p "$TOOLS_DIR"

# Download Picard if not already present
if [[ ! -f "$TOOLS_DIR/picard.jar" ]]; then
    echo "Downloading Picard version $PICARD_VERSION..."
    wget "$PICARD_URL" -O "$TOOLS_DIR/picard.jar"
    
    if [[ $? -eq 0 ]]; then
        echo "✓ Picard downloaded successfully"
    else
        echo "✗ Failed to download Picard"
        exit 1
    fi
else
    echo "Picard already exists at $TOOLS_DIR/picard.jar"
fi

echo "Picard setup complete!"