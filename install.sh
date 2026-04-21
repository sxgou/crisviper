#!/usr/bin/env bash
# Installation script for LineageTracer-Amplicon

set -e

echo "Installing LineageTracer-Amplicon..."

# Check Python version
python_version=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
echo "Python version: $python_version"

# Check if pip is available
if ! command -v pip3 &> /dev/null; then
    echo "Error: pip3 not found. Please install pip3."
    exit 1
fi

# Install in development mode
echo "Installing package in development mode..."
pip3 install -e .

echo ""
echo "Installation complete!"
echo ""
echo "To verify installation, run:"
echo "  lineage-tracer --help"
echo ""
echo "To generate a configuration file:"
echo "  lineage-tracer config-wizard --output my_config.json"
echo ""
echo "To analyze single-end data:"
echo "  lineage-tracer analyze --fastq_r1 sample_R1.fastq.gz --config my_config.json --output_dir results"
echo ""
echo "To analyze paired-end data (requires fastp):"
echo "  lineage-tracer analyze --fastq_r1 sample_R1.fastq.gz --fastq_r2 sample_R2.fastq.gz --config my_config.json --output_dir results"
echo ""
echo "For testing:"
echo "  python test/test_suite.py"
echo ""
echo "Note: For paired-end analysis, fastp must be installed separately:"
echo "  conda install -c bioconda fastp"
echo "  or visit: https://github.com/OpenGene/fastp"