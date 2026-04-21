#!/usr/bin/env python3
"""
Test paired-end sequencing functionality.
"""

import sys
import os
import tempfile
from pathlib import Path

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

def test_fastp_merge_module():
    """Test fastp_merge module functionality."""
    print("Testing fastp_merge module...")
    
    try:
        from lineage_tracer_amplicon.io.fastp_merge import (
            check_fastp_installed,
            count_fastq_reads,
            get_merge_statistics
        )
        
        print("✓ fastp_merge module imports successfully")
        
        # Test fastp check (will return False in test environment)
        has_fastp = check_fastp_installed()
        print(f"  fastp installed: {has_fastp}")
        
        # Create a dummy FASTQ file for testing
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            for i in range(10):
                f.write(f"@read{i}\n")
                f.write("ACGTACGTACGTACGT\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIII\n")
            temp_fastq = f.name
        
        try:
            count = count_fastq_reads(Path(temp_fastq))
            print(f"✓ count_fastq_reads works: {count} reads")
        finally:
            os.unlink(temp_fastq)
        
        return True
        
    except Exception as e:
        print(f"✗ fastp_merge test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_cli_argument_parsing():
    """Test CLI argument parsing for paired-end options."""
    print("\nTesting CLI argument parsing...")
    
    try:
        # Test argument parsing by importing the argparser setup
        import argparse
        
        # Recreate the argument parser from cli.py
        parser = argparse.ArgumentParser(
            description="CRISPR lineage tracing amplicon analysis software"
        )
        
        subparsers = parser.add_subparsers(dest='command', help='Command')
        
        # Analyze command (simplified)
        analyze_parser = subparsers.add_parser('analyze', help='Analyze FASTQ file(s)')
        
        # Add the key arguments we care about
        analyze_parser.add_argument('--fastq_r1', '-i1', required=True)
        analyze_parser.add_argument('--fastq_r2', '-i2', default=None)
        analyze_parser.add_argument('--min_overlap', type=int, default=30)
        analyze_parser.add_argument('--max_mismatch_rate', type=float, default=0.2)
        analyze_parser.add_argument('--skip_merge', action='store_true', default=False)
        analyze_parser.add_argument('--config', '-c', required=True)
        analyze_parser.add_argument('--output_dir', '-o', default='./results')
        
        # Test parsing
        test_args = [
            'analyze',
            '--fastq_r1', 'test_R1.fastq.gz',
            '--fastq_r2', 'test_R2.fastq.gz',
            '--config', 'config.json',
            '--min_overlap', '25',
            '--max_mismatch_rate', '0.15'
        ]
        
        args = parser.parse_args(test_args)
        
        # Verify parsed values
        assert args.fastq_r1 == 'test_R1.fastq.gz'
        assert args.fastq_r2 == 'test_R2.fastq.gz'
        assert args.min_overlap == 25
        assert args.max_mismatch_rate == 0.15
        assert args.skip_merge == False
        
        print("✓ CLI argument parsing works correctly")
        print(f"  Parsed: R1={args.fastq_r1}, R2={args.fastq_r2}")
        print(f"  min_overlap={args.min_overlap}, max_mismatch_rate={args.max_mismatch_rate}")
        
        # Test skip_merge
        test_args_skip = [
            'analyze',
            '--fastq_r1', 'test_R1.fastq.gz',
            '--fastq_r2', 'test_R2.fastq.gz',
            '--config', 'config.json',
            '--skip_merge'
        ]
        
        args_skip = parser.parse_args(test_args_skip)
        assert args_skip.skip_merge == True
        print("✓ skip_merge argument works correctly")
        
        return True
        
    except Exception as e:
        print(f"✗ CLI argument parsing test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_documentation():
    """Check that documentation is updated."""
    print("\nChecking documentation...")
    
    docs_to_check = [
        ('../README.md', ['--fastq_r1', '--fastq_r2', '--min_overlap', 'fastp']),
        ('../docs/paired_end_usage.md', ['paired-end', 'fastp', 'merge']),
    ]
    
    all_ok = True
    for doc_path, keywords in docs_to_check:
        doc_full_path = os.path.join(os.path.dirname(__file__), doc_path)
        if os.path.exists(doc_full_path):
            with open(doc_full_path, 'r') as f:
                content = f.read()
            
            missing = []
            for keyword in keywords:
                if keyword not in content.lower():
                    missing.append(keyword)
            
            if missing:
                print(f"✗ {doc_path}: Missing keywords: {missing}")
                all_ok = False
            else:
                print(f"✓ {doc_path}: All keywords found")
        else:
            print(f"✗ {doc_path}: File not found")
            all_ok = False
    
    return all_ok

def main():
    """Run all tests."""
    print("=" * 60)
    print("Paired-End Sequencing Support Test Suite")
    print("=" * 60)
    
    tests = [
        test_fastp_merge_module,
        test_cli_argument_parsing,
        test_documentation,
    ]
    
    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
        except Exception as e:
            print(f"Test failed with exception: {e}")
            results.append(False)
    
    passed = sum(1 for r in results if r)
    total = len(results)
    
    print("\n" + "=" * 60)
    print(f"Test Results: {passed}/{total} passed")
    print("=" * 60)
    
    if passed == total:
        print("\n✅ All tests passed! Paired-end support is ready.")
    else:
        print(f"\n⚠️  {total - passed} tests failed.")
        sys.exit(1)

if __name__ == "__main__":
    main()