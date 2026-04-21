"""
Command-line interface for lineage tracer amplicon analysis.
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

from .analysis.engine import LineageTracerAnalyzer


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="CRISPR lineage tracing amplicon analysis software",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  lineage-tracer analyze \\
    --fastq sample_R1.fastq.gz \\
    --config amplicon_structure.json \\
    --output_dir ./results \\
    --sample_name Embryo_Day5 \\
    --threads 8 \\
    --min_quality 20 \\
    --primer_max_mismatch_rate 0.15 \\
    --min_anchor_identity 0.8 \\
    --snv_window 3
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Command')
    
    # Analyze command
    analyze_parser = subparsers.add_parser('analyze', help='Analyze FASTQ file')
    analyze_parser.add_argument(
        '--fastq', '-i',
        type=str,
        required=True,
        help='Input FASTQ file (can be gzipped)'
    )
    analyze_parser.add_argument(
        '--config', '-c',
        type=str,
        required=True,
        help='Amplicon structure configuration file (JSON/YAML)'
    )
    analyze_parser.add_argument(
        '--output_dir', '-o',
        type=str,
        default='./results',
        help='Output directory (default: ./results)'
    )
    analyze_parser.add_argument(
        '--sample_name', '-s',
        type=str,
        default='sample',
        help='Sample name for output files (default: sample)'
    )
    analyze_parser.add_argument(
        '--threads', '-t',
        type=int,
        default=1,
        help='Number of threads (default: 1)'
    )
    analyze_parser.add_argument(
        '--min_quality', '-q',
        type=int,
        default=20,
        help='Minimum base quality threshold (default: 20)'
    )
    analyze_parser.add_argument(
        '--primer_max_mismatch_rate',
        type=float,
        default=0.15,
        help='Maximum mismatch rate for primer validation (default: 0.15)'
    )
    analyze_parser.add_argument(
        '--trim_primers',
        action='store_true',
        default=True,
        help='Trim primer sequences before analysis (default: True)'
    )
    analyze_parser.add_argument(
        '--no_trim_primers',
        dest='trim_primers',
        action='store_false',
        help='Do not trim primer sequences'
    )
    analyze_parser.add_argument(
        '--min_anchor_identity',
        type=float,
        default=0.8,
        help='Minimum identity for anchor matching (default: 0.8)'
    )
    analyze_parser.add_argument(
        '--kmer_size', '-k',
        type=int,
        default=9,
        help='k-mer size for anchor indexing (default: 9)'
    )
    analyze_parser.add_argument(
        '--min_deletion_size',
        type=int,
        default=20,
        help='Minimum size for reporting as large deletion (default: 20)'
    )
    analyze_parser.add_argument(
        '--snv_window',
        type=int,
        default=3,
        help='Window size around cutsite for valid SNVs (default: 3)'
    )
    analyze_parser.add_argument(
        '--max_reads',
        type=int,
        default=None,
        help='Maximum number of reads to process (for testing)'
    )
    analyze_parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output'
    )
    
    # Config wizard command
    config_parser = subparsers.add_parser('config-wizard', help='Interactive configuration wizard')
    config_parser.add_argument(
        '--output', '-o',
        type=str,
        default='amplicon_structure.json',
        help='Output configuration file (default: amplicon_structure.json)'
    )
    
    args = parser.parse_args()
    
    if args.command == 'analyze':
        analyze_command(args)
    elif args.command == 'config-wizard':
        config_wizard_command(args)
    else:
        parser.print_help()
        sys.exit(1)


def analyze_command(args) -> None:
    """Execute analyze command."""
    fastq_path = Path(args.fastq)
    config_path = Path(args.config)
    output_dir = Path(args.output_dir)
    
    if not fastq_path.exists():
        print(f"Error: FASTQ file not found: {fastq_path}")
        sys.exit(1)
    
    if not config_path.exists():
        print(f"Error: Configuration file not found: {config_path}")
        sys.exit(1)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Starting analysis...")
    print(f"  Input FASTQ: {fastq_path}")
    print(f"  Configuration: {config_path}")
    print(f"  Output directory: {output_dir}")
    print(f"  Sample name: {args.sample_name}")
    print(f"  Threads: {args.threads}")
    print(f"  Parameters:")
    print(f"    - Min quality: {args.min_quality}")
    print(f"    - Primer max mismatch rate: {args.primer_max_mismatch_rate}")
    print(f"    - Trim primers: {args.trim_primers}")
    print(f"    - Min anchor identity: {args.min_anchor_identity}")
    print(f"    - k-mer size: {args.kmer_size}")
    print(f"    - Min deletion size: {args.min_deletion_size}")
    print(f"    - SNV window: {args.snv_window}")
    
    try:
        # Initialize analyzer
        analyzer = LineageTracerAnalyzer(
            config_path=config_path,
            primer_max_mismatch_rate=args.primer_max_mismatch_rate,
            min_anchor_identity=args.min_anchor_identity,
            kmer_size=args.kmer_size,
            snv_window_size=args.snv_window,
            trim_primers=args.trim_primers,
            min_deletion_size=args.min_deletion_size
        )
        
        # Analyze FASTQ file
        stats = analyzer.analyze_fastq(
            fastq_path=fastq_path,
            output_dir=output_dir,
            max_reads=args.max_reads
        )
        
        # Print summary
        print("\nAnalysis complete!")
        print(f"  Total reads processed: {stats['total_reads']}")
        print(f"  Reads with both primers: {stats['both_primers_detected']}")
        print(f"  Valid reads for analysis: {stats['passed_reads']}")
        print(f"  Editing efficiency: {stats.get('editing_efficiency', 0):.3f}")
        
        # Output files generated
        print("\nOutput files generated:")
        output_files = [
            output_dir / "per_read_annotation.tsv",
            output_dir / "barcode_frequencies.tsv",
            output_dir / "statistics.json"
        ]
        
        for file_path in output_files:
            if file_path.exists():
                print(f"  - {file_path}")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def config_wizard_command(args) -> None:
    """Interactive configuration wizard."""
    output_path = Path(args.output)
    
    print("Amplicon Structure Configuration Wizard")
    print("=" * 50)
    print("\nThis wizard will guide you through creating an amplicon structure configuration.")
    print("Please have the following information ready:")
    print("1. Full reference sequence")
    print("2. Coordinates of primers, targets, conserved regions, cutsites, etc.")
    print("\nPress Ctrl+C to exit at any time.")
    
    try:
        # Get reference sequence
        reference = input("\nEnter the full reference sequence: ").strip().upper()
        
        if not reference:
            print("Error: Reference sequence cannot be empty")
            return
        
        # Validate sequence
        import re
        if not re.match(r'^[ACGTN]+$', reference):
            print("Warning: Sequence contains non-standard nucleotides")
        
        features = []
        
        # Add Primer5
        print("\n--- Primer5 ---")
        primer5_name = "Primer5"
        primer5_start = int(input("Enter Primer5 start position (0-based): "))
        primer5_end = int(input("Enter Primer5 end position (0-based exclusive): "))
        
        features.append({
            "name": primer5_name,
            "start": primer5_start,
            "end": primer5_end,
            "type": "primer"
        })
        
        # Add prefix
        print("\n--- Prefix ---")
        prefix_name = "prefix"
        prefix_start = int(input("Enter prefix start position: "))
        prefix_end = int(input("Enter prefix end position: "))
        
        features.append({
            "name": prefix_name,
            "start": prefix_start,
            "end": prefix_end,
            "type": "prefix"
        })
        
        # Add targets
        target_count = int(input("\nHow many targets? (e.g., 10): "))
        
        for i in range(1, target_count + 1):
            print(f"\n--- Target {i} ---")
            target_name = f"Target{i}"
            target_start = int(input(f"Enter Target{i} start position: "))
            target_end = int(input(f"Enter Target{i} end position: "))
            
            conserved_start = int(input(f"Enter conserved region start for Target{i}: "))
            conserved_end = int(input(f"Enter conserved region end for Target{i}: "))
            
            cutsite_start = int(input(f"Enter cutsite start for Target{i}: "))
            cutsite_end = int(input(f"Enter cutsite end for Target{i}: "))
            
            features.append({
                "name": target_name,
                "start": target_start,
                "end": target_end,
                "type": "target",
                "conserved": [conserved_start, conserved_end],
                "cutsite": [cutsite_start, cutsite_end]
            })
            
            # Add PAM_Linker for targets 1-9
            if i < target_count:
                print(f"\n--- PAM_Linker{i} ---")
                pam_linker_name = f"PAM_Linker{i}"
                pam_linker_start = int(input(f"Enter PAM_Linker{i} start position: "))
                pam_linker_end = int(input(f"Enter PAM_Linker{i} end position: "))
                
                pam_start = int(input(f"Enter PAM start position in PAM_Linker{i}: "))
                pam_end = int(input(f"Enter PAM end position in PAM_Linker{i}: "))
                
                features.append({
                    "name": pam_linker_name,
                    "start": pam_linker_start,
                    "end": pam_linker_end,
                    "type": "pam_linker",
                    "pam": [pam_start, pam_end],
                    "linker": [pam_end, pam_linker_end]
                })
        
        # Add postfix
        print("\n--- Postfix ---")
        postfix_name = "postfix"
        postfix_start = int(input("Enter postfix start position: "))
        postfix_end = int(input("Enter postfix end position: "))
        
        features.append({
            "name": postfix_name,
            "start": postfix_start,
            "end": postfix_end,
            "type": "postfix"
        })
        
        # Add Primer3
        print("\n--- Primer3 ---")
        primer3_name = "Primer3"
        primer3_start = int(input("Enter Primer3 start position: "))
        primer3_end = int(input("Enter Primer3 end position: "))
        
        features.append({
            "name": primer3_name,
            "start": primer3_start,
            "end": primer3_end,
            "type": "primer"
        })
        
        # Create configuration
        config = {
            "reference": reference,
            "features": features
        }
        
        # Write to file
        import json
        with open(output_path, 'w') as f:
            json.dump(config, f, indent=2)
        
        print(f"\nConfiguration saved to: {output_path}")
        print("\nYou can now use this configuration file with:")
        print(f"  lineage-tracer analyze --fastq your_data.fastq.gz --config {output_path} --output_dir ./results")
        
    except KeyboardInterrupt:
        print("\n\nWizard cancelled.")
    except Exception as e:
        print(f"\nError: {e}")


if __name__ == "__main__":
    main()