"""Threshold computation for UMI/CB filtering.

Port of MATLAB's TaggedCollection.threshold_function.
"""

import numpy as np


def compute_threshold(
    freqs: np.ndarray,
    max_elem: int,
    n_reads: int,
    p: float = 0.01,
    read_floor: int = 2,
    read_override: float = np.nan,
    threshold_type: str = "UMI",
) -> dict:
    """Compute a statistical read-count threshold for UMI/cell filtering.

    Several heuristics are computed and the maximum is chosen (unless an
    explicit override is given).

    Args:
        freqs: Sorted read-count frequencies (descending).
        max_elem: Maximum expected molecules or cells.
        n_reads: Total read count.
        p: Estimated sequencing error rate.
        read_floor: Absolute minimum floor.
        read_override: Explicit override (NaN = use computed max).
        threshold_type: 'UMI' or 'CB'. Controls the field name in output.

    Returns:
        dict with keys: one_tenth_99_pctl, max_molecules|max_cells,
        equal_partition, err_floor, read_floor, chosen.
    """
    field = "max_molecules" if threshold_type == "UMI" else "max_cells"
    ind_99 = max(round(len(freqs) / 100), 1)

    thresholds = {
        "one_tenth_99_pctl": int(np.ceil(freqs[ind_99 - 1] / 10)),
        field: int(freqs[min(max_elem, len(freqs)) - 1] + 1),
        "equal_partition": int(np.ceil(n_reads / max_elem)),
        "err_floor": int(np.ceil(freqs[0] * p * (1 - p) ** 9)),
        "read_floor": read_floor,
    }

    if np.isnan(read_override):
        thresholds["chosen"] = max(thresholds.values())
    else:
        thresholds["override"] = int(read_override)
        thresholds["chosen"] = thresholds["override"]

    return thresholds
