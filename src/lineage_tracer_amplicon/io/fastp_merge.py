"""
Fastp-based read merging for paired-end sequencing data.
"""

import subprocess
import tempfile
import os
from pathlib import Path
from typing import Tuple, Optional, List
import gzip
import shutil


def check_fastp_installed() -> bool:
    """
    Check if fastp is installed and available.
    
    Returns:
        True if fastp is available, False otherwise
    """
    try:
        result = subprocess.run(
            ["fastp", "--version"],
            capture_output=True,
            text=True,
            check=False
        )
        return result.returncode == 0
    except (subprocess.SubprocessError, FileNotFoundError):
        return False


def merge_paired_end_reads(
    fastq_r1: Path,
    fastq_r2: Path,
    output_dir: Path,
    sample_name: str = "sample",
    min_overlap: int = 10,
    max_mismatch_rate: float = 0.2,
    min_quality: int = 20,
    threads: int = 1,
    verbose: bool = False
) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    """
    Merge paired-end reads using fastp.
    
    Args:
        fastq_r1: Path to R1 FASTQ file
        fastq_r2: Path to R2 FASTQ file
        output_dir: Output directory
        sample_name: Sample name for output files
        min_overlap: Minimum overlap length required for merging (default: 10)
        max_mismatch_rate: Maximum mismatch rate allowed in overlap region (default: 0.2)
        min_quality: Minimum base quality (default: 20)
        threads: Number of threads for fastp (default: 1)
        verbose: Print verbose output (default: False)
        
    Returns:
        Tuple of (merged_fastq_path, fastp_json_report, fastp_html_report) or (None, None, None) on failure
    """
    # Check if fastp is installed
    if not check_fastp_installed():
        raise RuntimeError(
            "fastp is not installed or not in PATH. "
            "Please install fastp: conda install -c bioconda fastp"
        )
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Prepare output file paths
    merged_fastq = output_dir / f"{sample_name}_merged.fastq.gz"
    fastp_json = output_dir / f"{sample_name}_fastp.json"
    fastp_html = output_dir / f"{sample_name}_fastp.html"
    
    # Build fastp command
    cmd = [
        "fastp",
        "--in1", str(fastq_r1),
        "--in2", str(fastq_r2),
        "--stdout",  # Output merged reads to stdout
        "--merge",
        "--merged_out", str(merged_fastq),
        "--overlap_len_require", str(min_overlap),
        "--overlap_diff_limit", str(int(min_overlap * max_mismatch_rate)),
        "--correction",
        "--thread", str(threads),
        "--json", str(fastp_json),
        "--html", str(fastp_html),
        "--report_title", f"{sample_name} - fastp merge report"
    ]
    
    # Add quality filtering if specified
    if min_quality > 0:
        cmd.extend(["--qualified_quality_phred", str(min_quality)])
    
    # Add verbose flag if requested
    if verbose:
        cmd.append("--verbose")
    
    # Run fastp
    if verbose:
        print(f"Running fastp command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        if verbose:
            print(f"fastp stdout: {result.stdout[:500]}")
            if result.stderr:
                print(f"fastp stderr: {result.stderr[:500]}")
        
        # Check if output file was created
        if not merged_fastq.exists():
            raise RuntimeError(f"fastp did not create output file: {merged_fastq}")
        
        if verbose:
            print(f"Successfully merged reads: {merged_fastq}")
            print(f"fastp JSON report: {fastp_json}")
            print(f"fastp HTML report: {fastp_html}")
        
        return merged_fastq, fastp_json, fastp_html
        
    except subprocess.CalledProcessError as e:
        error_msg = f"fastp failed with exit code {e.returncode}: {e.stderr}"
        if verbose:
            print(f"fastp error: {error_msg}")
        raise RuntimeError(error_msg)
    except Exception as e:
        error_msg = f"Error running fastp: {str(e)}"
        if verbose:
            print(f"Error: {error_msg}")
        raise RuntimeError(error_msg)


def process_paired_end_input(
    fastq_r1: Path,
    fastq_r2: Optional[Path] = None,
    output_dir: Path = Path("./"),
    sample_name: str = "sample",
    min_overlap: int = 10,
    max_mismatch_rate: float = 0.2,
    min_quality: int = 20,
    threads: int = 1,
    skip_merge: bool = False,
    verbose: bool = False
) -> Path:
    """
    Process paired-end or single-end input, merging if needed.
    
    Args:
        fastq_r1: Path to R1 FASTQ file (required)
        fastq_r2: Path to R2 FASTQ file (optional, if None treat as single-end)
        output_dir: Output directory for intermediate files
        sample_name: Sample name
        min_overlap: Minimum overlap length for merging
        max_mismatch_rate: Maximum mismatch rate in overlap region
        min_quality: Minimum base quality
        threads: Number of threads
        skip_merge: Skip merging even for paired-end (use R1 only)
        verbose: Print verbose output
        
    Returns:
        Path to FASTQ file to analyze (merged or single-end)
    """
    # Create output directory for intermediate files
    intermediate_dir = output_dir / "intermediate"
    intermediate_dir.mkdir(parents=True, exist_ok=True)
    
    if fastq_r2 is None or skip_merge:
        # Single-end mode, just use R1
        if verbose:
            print(f"Single-end mode: using {fastq_r1}")
        return fastq_r1
    else:
        # Paired-end mode, merge with fastp
        if verbose:
            print(f"Paired-end mode: merging {fastq_r1} and {fastq_r2}")
            print(f"  min_overlap: {min_overlap}")
            print(f"  max_mismatch_rate: {max_mismatch_rate}")
            print(f"  min_quality: {min_quality}")
        
        merged_fastq, json_report, html_report = merge_paired_end_reads(
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            output_dir=intermediate_dir,
            sample_name=sample_name,
            min_overlap=min_overlap,
            max_mismatch_rate=max_mismatch_rate,
            min_quality=min_quality,
            threads=threads,
            verbose=verbose
        )
        
        if verbose:
            print(f"Merged reads saved to: {merged_fastq}")
        
        return merged_fastq


def count_fastq_reads(fastq_path: Path) -> int:
    """
    Count number of reads in a FASTQ file.
    
    Args:
        fastq_path: Path to FASTQ file (can be gzipped)
        
    Returns:
        Number of reads in the file
    """
    count = 0
    open_func = gzip.open if str(fastq_path).endswith('.gz') else open
    
    try:
        with open_func(fastq_path, 'rt') as f:
            for line in f:
                if line.startswith('@'):
                    count += 1
    except Exception as e:
        print(f"Warning: Could not count reads in {fastq_path}: {e}")
        return 0
    
    return count


def get_merge_statistics(json_report_path: Path) -> dict:
    """
    Parse fastp JSON report to get merge statistics.
    
    Args:
        json_report_path: Path to fastp JSON report
        
    Returns:
        Dictionary with merge statistics
    """
    import json
    
    try:
        with open(json_report_path, 'r') as f:
            report = json.load(f)
        
        # Extract relevant statistics
        summary = report.get('summary', {})
        merging = report.get('merging', {})
        
        stats = {
            'total_reads': summary.get('before_filtering', {}).get('total_reads', 0),
            'merged_reads': merging.get('merged_reads', 0),
            'merge_rate': merging.get('merge_rate', 0),
            'overlap_bases': merging.get('overlap_bases', {}),
            'adapter_trimming': report.get('adapter_cutting', {}),
            'quality_filtering': summary.get('after_filtering', {}),
        }
        
        return stats
        
    except Exception as e:
        print(f"Warning: Could not parse fastp report {json_report_path}: {e}")
        return {}