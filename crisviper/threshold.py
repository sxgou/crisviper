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
    umi_length: int = 10,
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
        umi_length: UMI/barcode length for error floor calculation (default 10).

    Returns:
        dict with keys: one_tenth_99_pctl, max_molecules|max_cells,
        equal_partition, err_floor, read_floor, chosen.
    """
    field = "max_molecules" if threshold_type == "UMI" else "max_cells"
    ind_99 = max(round(len(freqs) / 100), 1)

    # Guard: empty frequency array → return minimal threshold
    if len(freqs) == 0:
        return {field: 0, "one_tenth_99_pctl": 0, "err_floor": 0,
                "equal_partition": int(np.ceil(n_reads / max_elem)) if max_elem else 0,
                "read_floor": read_floor, "chosen": read_floor}

    # Error floor: expected number of reads at which a UMI of length `umi_length`
    # has 0 errors given per-base error rate p:  freqs[0] * (1-p)^umi_length
    err_floor_val = int(np.ceil(freqs[0] * p * (1 - p) ** (umi_length - 1)))

    thresholds = {
        "one_tenth_99_pctl": int(np.ceil(freqs[ind_99 - 1] / 10)),
        field: int(freqs[min(max_elem, len(freqs)) - 1] + 1),
        "equal_partition": int(np.ceil(n_reads / max_elem)),
        "err_floor": err_floor_val,
        "read_floor": read_floor,
    }

    if np.isnan(read_override):
        thresholds["chosen"] = max(thresholds.values())
    else:
        thresholds["override"] = int(read_override)
        thresholds["chosen"] = thresholds["override"]

    return thresholds
