"""
Simple test to verify imports work.
"""

import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

try:
    from lineage_tracer_amplicon.io.config_parser import load_config
    print("✓ Config parser import successful")
    
    from lineage_tracer_amplicon.core.primer_validation import validate_primers
    print("✓ Primer validation import successful")
    
    from lineage_tracer_amplicon.core.anchor_alignment import AnchorIndex
    print("✓ Anchor alignment import successful")
    
    from lineage_tracer_amplicon.core.mutation_classification import classify_deletion_with_insertion
    print("✓ Mutation classification import successful")
    
    from lineage_tracer_amplicon.analysis.engine import LineageTracerAnalyzer
    print("✓ Analysis engine import successful")
    
    from lineage_tracer_amplicon.cli import main as cli_main
    print("✓ CLI import successful")
    
    print("\n✓ All imports successful!")
    
except Exception as e:
    print(f"✗ Import error: {e}")
    import traceback
    traceback.print_exc()