"""ltlib/metrics.py — diversity and heterogeneity metrics.

MATLAB equivalents from matlab_version/metrics/:
  - effective_alleles:     Shannon entropy → effective allele count
  - diversity_index:       effective_alleles normalized by total counts
  - alleles_per_cell:      average alleles per cell (bulk) or CB (SC)
  - singletons_per_cell:   singleton allele rate
  - bootstrap_*:           bootstrap CI for metrics
"""

import numpy as np
from typing import List, Optional


def effective_alleles(
    allele_freqs: List[int],
    is_template: Optional[List[bool]] = None,
) -> float:
    """Compute effective number of alleles (Shannon entropy-based).

    MATLAB equivalent: effective_alleles.m

    Excludes the template (wild-type) allele from the calculation.
    Returns 2^H where H is Shannon entropy of the allele frequency distribution,
    minus 1 if there are non-template alleles (MATLAB: 2^H - ~isempty(p)).

    Args:
        allele_freqs: List of allele frequencies (counts).
        is_template: Optional boolean mask, True for template allele.
                     If None, no template exclusion.

    Returns:
        Effective allele count (>= 0).
    """
    counts = np.array(allele_freqs, dtype=float)
    if is_template is not None:
        mask = np.array(is_template, dtype=bool)
        counts = counts[~mask]
    counts = counts[counts > 0]
    if len(counts) == 0:
        return 0.0
    p = counts / counts.sum()
    H = -np.sum(p * np.log2(p))
    return 2.0 ** H - 1.0 if len(p) > 0 else 0.0


def diversity_index(
    allele_freqs: List[int],
    is_template: Optional[List[bool]] = None,
    normalize_by_edited: bool = False,
) -> float:
    """Diversity index — effective alleles normalized by total calls.

    MATLAB equivalent: diversity_index.m

    Args:
        allele_freqs: List of allele frequencies.
        is_template: Boolean mask, True for template allele.
        normalize_by_edited: If True, normalize only by edited (non-template) calls.

    Returns:
        Diversity index value.
    """
    ec = effective_alleles(allele_freqs, is_template)
    total = sum(allele_freqs)
    if total == 0:
        return 0.0
    if normalize_by_edited and is_template is not None:
        edited_total = sum(f for f, t in zip(allele_freqs, is_template) if not t)
        return max(0.0, ec / edited_total) if edited_total > 0 else 0.0
    return ec / total


def alleles_per_cell(
    n_alleles: int,
    total_cells: int,
) -> float:
    """Average alleles per cell (or per UMI).

    MATLAB equivalent: alleles_per_cell.m

    Args:
        n_alleles: Number of distinct alleles.
        total_cells: Total number of cells (or UMIs).

    Returns:
        Average alleles per cell.
    """
    if total_cells == 0:
        return 0.0
    return n_alleles / total_cells


def singletons_per_cell(
    n_singletons: int,
    total_cells: int,
) -> float:
    """Singleton alleles per cell.

    MATLAB equivalent: singletons_per_cell.m

    Args:
        n_singletons: Number of singleton alleles (frequency == 1).
        total_cells: Total number of cells (or UMIs).

    Returns:
        Singleton rate per cell.
    """
    if total_cells == 0:
        return 0.0
    return n_singletons / total_cells


def carlin_potential(
    n_targets: int,
    n_modified_sites: int,
) -> float:
    """Mean CARLIN potential — average number of remaining editable targets.

    MATLAB equivalent: CARLIN_potential.m

    Args:
        n_targets: Total number of target sites in the CARLIN array.
        n_modified_sites: Number of sites that have been modified.

    Returns:
        Mean remaining editable targets.
    """
    return max(0.0, n_targets - n_modified_sites)
