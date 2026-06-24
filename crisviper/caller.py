"""Allele calling — coarse-grain and exact consensus.

Port of MATLAB's CallableCollection.call_alleles_coarse_grain and
call_alleles_exact.
"""

from typing import List, Tuple
from dataclasses import dataclass, field
from crisviper.models import AlignmentResult, MutationEvent, MutationType


# ═══════════════════════════════════════════════════════════════
# Data model
# ═══════════════════════════════════════════════════════════════

@dataclass
class CalledAllele:
    """A single called allele with its supporting reads."""
    sequence: str                     # consensus full-length sequence
    aligned_sequence: str             # consensus aligned sequence
    weight: int                       # total supporting read count
    n_sequences: int                  # number of unique sequences
    event_structure: str              # coarse-grain signature
    mutations: List[MutationEvent] = field(default_factory=list)


# ═══════════════════════════════════════════════════════════════
# Event structure signature
# ═══════════════════════════════════════════════════════════════

def _event_structure(mutations: List[MutationEvent]) -> str:
    """Generate a compact signature string from mutation events.

    Examples:
        []                     → "WT"
        [DEL(41,3)]            → "D3@41"
        [SUB(65)]              → "S@65"
        [DEL(41,3), SUB(68)]   → "D3@41+S@68"
    """
    if not mutations:
        return "WT"
    parts = []
    for m in sorted(mutations, key=lambda x: x.ref_pos):
        if m.type == MutationType.SUBSTITUTION:
            parts.append(f"S@{m.ref_pos}")
        elif m.type == MutationType.DELETION:
            parts.append(f"D{m.length}@{m.ref_pos}")
        elif m.type == MutationType.INSERTION:
            parts.append(f"I{m.length}@{m.ref_pos}")
        elif m.type == MutationType.INDEL:
            parts.append(f"C{m.length}@{m.ref_pos}")
    return "+".join(parts) if parts else "WT"


# ═══════════════════════════════════════════════════════════════
# Grouping helper
# ═══════════════════════════════════════════════════════════════

def _unique_by_freq(
    keys: List[str], weights: List[int]
) -> Tuple[List[str], List[int], List[int]]:
    """Like MATLAB unique_by_freq: unique keys sorted by descending frequency.

    Returns:
        (unique_keys, sorted_weights, reverse_index)
        where reverse_index[i] = index into unique_keys for original item i.
    """
    from collections import Counter
    weight_per_key = Counter()
    for k, w in zip(keys, weights):
        weight_per_key[k] += w
    sorted_keys = sorted(weight_per_key, key=weight_per_key.get, reverse=True)
    sorted_weights = [weight_per_key[k] for k in sorted_keys]
    key_to_idx = {k: i for i, k in enumerate(sorted_keys)}
    reverse_index = [key_to_idx[k] for k in keys]
    return sorted_keys, sorted_weights, reverse_index


# ═══════════════════════════════════════════════════════════════
# Majority-vote consensus per column
# ═══════════════════════════════════════════════════════════════

def _majority_consensus(sequences: List[str], weights: List[int]) -> str:
    """Build per-position majority-rule consensus.

    At each column, the base with the highest total weight wins.
    All sequences must have the same length.
    """
    if not sequences:
        return ""
    L = len(sequences[0])
    if not all(len(s) == L for s in sequences):
        raise ValueError(
            f"Majority consensus requires all sequences to have the same length, "
            f"got lengths {set(len(s) for s in sequences)}"
        )
    consensus = []
    for col in range(L):
        col_weights: dict = {}
        for seq, w in zip(sequences, weights):
            col_weights[seq[col]] = col_weights.get(seq[col], 0) + w
        consensus.append(max(col_weights, key=col_weights.get))
    return "".join(consensus)


# ═══════════════════════════════════════════════════════════════
# Shared allele calling helper
# ═══════════════════════════════════════════════════════════════

def _call_alleles_by_key(
    valid: List[AlignmentResult],
    keys: List[str],
    dominant_frac: float,
) -> List[CalledAllele]:
    """Group results by key, filter by dominant_frac, build consensus.

    Args:
        valid: Filtered list of successful alignment results.
        keys: Grouping key per result.
        dominant_frac: Minimum weight fraction to keep an allele.

    Returns:
        List of CalledAllele objects sorted by weight descending.
    """
    unique_keys, sorted_weights, reverse_index = _unique_by_freq(
        keys, [r.query.readCount for r in valid]
    )
    total_weight = sum(r.query.readCount for r in valid)
    if not unique_keys:
        return []

    alleles = []
    for i, (key, key_weight) in enumerate(zip(unique_keys, sorted_weights)):
        if key_weight / total_weight <= dominant_frac:
            break
        idx = [j for j, ri in enumerate(reverse_index) if ri == i]
        group = [valid[j] for j in idx]
        group_weights = [r.query.readCount for r in group]
        consensus_aq = _majority_consensus(
            [r.aligned_query for r in group], group_weights
        )
        alleles.append(CalledAllele(
            sequence=consensus_aq.replace("-", ""),
            aligned_sequence=consensus_aq,
            weight=key_weight,
            n_sequences=len(group),
            event_structure=key,
            mutations=group[0].mutations,
        ))
    return alleles


# ═══════════════════════════════════════════════════════════════
# Public callers
# ═══════════════════════════════════════════════════════════════

def call_alleles_coarse_grain(
    results: List[AlignmentResult],
    dominant_frac: float = 0.5,
) -> List[CalledAllele]:
    """Call alleles by coarse-grain event structure grouping.

    Groups alignment results by their mutation signature (event structure).
    If a group exceeds dominant_frac of total weight, it is kept as a
    called allele with majority-rule consensus sequence.

    Args:
        results: List of successful alignment results.
        dominant_frac: Minimum weight fraction to keep an allele (default 0.5).

    Returns:
        List of CalledAllele objects, sorted by weight descending.
    """
    valid = [r for r in results if r.success and r.stats]
    if not valid:
        return []
    keys = [_event_structure(r.mutations) for r in valid]
    return _call_alleles_by_key(valid, keys, dominant_frac)


def call_alleles_exact(
    results: List[AlignmentResult],
    dominant_frac: float = 0.5,
) -> List[CalledAllele]:
    """Call alleles by exact aligned-query sequence grouping.

    Groups alignment results by their exact aligned query sequence.
    If a group exceeds dominant_frac of total weight, it is kept.

    Args:
        results: List of successful alignment results.
        dominant_frac: Minimum weight fraction to keep an allele.

    Returns:
        List of CalledAllele objects, sorted by weight descending.
    """
    valid = [r for r in results if r.success and r.stats]
    if not valid:
        return []
    keys = [r.aligned_query for r in valid]
    if not all(keys):
        return []
    # Keys are aligned queries; event_structure derived from mutation events
    alleles = _call_alleles_by_key(valid, keys, dominant_frac)
    for a in alleles:
        a.event_structure = _event_structure(a.mutations)
    return alleles
