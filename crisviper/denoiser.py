"""UMI/CB denoising — directional adjacency top-down clustering.

Port of MATLAB's TaggedCollection.directional_adjacency_top_down_denoiser.
"""

import numpy as np
from typing import List, Optional


def directional_adjacency_top_down_denoiser(
    tags: List[str],
    weights: np.ndarray,
    exclude: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Cluster tags by Hamming-distance-1 adjacency, top-down by weight.

    For each tag processed in descending order of weight, unassigned tags
    with Hamming distance 1, the same length, weight within 2x, and not
    excluded are assigned to the same cluster.

    Args:
        tags: Tag sequences (UMI or barcode).
        weights: Weight (e.g. read count) per tag.
        exclude: Boolean mask of tags ineligible as children.

    Returns:
        tag_map[i] = index of the parent (representative) tag for tag i.
    """
    N = len(tags)
    assert len(weights) == N
    if exclude is None:
        exclude = np.zeros(N, dtype=bool)

    tag_map = np.full(N, -1, dtype=int)
    lengths = np.array([len(t) for t in tags])
    order = np.argsort(weights)[::-1]

    for i in order:
        if tag_map[i] == -1:
            tag_map[i] = i
        # unassigned, same-length, not-excluded, weight within 2x
        candidates = np.where(
            (weights[i] >= 2 * weights - 1)
            & (tag_map == -1)
            & (lengths == lengths[i])
            & ~exclude
        )[0]
        if len(candidates) == 0:
            continue
        tag_i = tags[i]
        for c in candidates:
            hd = sum(a != b for a, b in zip(tag_i, tags[c]))
            if hd == 1:
                tag_map[c] = tag_map[i]

    return tag_map
