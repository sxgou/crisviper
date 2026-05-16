# Position-Dependent Gap Exit Bonus Profile

## Problem

Global alignment with gradient-based position-dependent gap penalties (lineage mode) can produce spurious isolated match segments at cutsite regions. These short match segments (1-4 bp) split what should be a single large deletion into multiple fragments, creating biologically incorrect allele calls.

In the CARLIN reference (332 bp, 10 targets with 27 bp period), 290 isolated match segments exist across 23,430 query sequences, occurring exclusively within 12 bp of a cutsite center.

## Root Cause

The root cause is a **frame shift artifact**: when a query with a large deletion (e.g., Target1→Target7) reaches the junction, the query sequence adjacent to the deletion boundary happens to contain, by chance or due to the repetitive amplicon structure, a short segment that perfectly matches a cutsite region near the deletion start. For example, `ref[71:85]` (Target2 cutsite) and `query[39:53]` (junction sequence) share 13/14 consecutive matches.

The existing DP mechanisms fail to prevent this:

| Feature | Problem |
|---------|---------|
| `isolated_base_penalty` | Uses sequence-level `run_len` — if the diagonal has consecutive matches (which is the case here, `run_len=13`), the penalty does not trigger |
| `short_match_discount` | Same `run_len` dependency — window=3 misses runs of 13 |
| `gap_exit_bonus` (scalar) | Applied uniformly; at cutsite (scale=1.0), -1.0 is too weak to overcome the frame shift benefit (~20 points) |

The match itself is not the error — the error is that *exiting a gap to take a match at a cutsite, then immediately re-entering a gap*, is too cheap relative to the benefit of frame-shifting the rest of the alignment.

## Solution: Position-Dependent gap_exit_bonus Profile

Replace the scalar `gap_exit_bonus` with a **position-dependent profile** (geb_profile) that applies inverted gradient scaling: the strongest gap exit penalty at cutsite centers (where matches are most suspicious), smoothly tapering to the weakest penalty at conservative inter-target regions.

### Formula

```
geb_profile[i] = base_geb × (max_scale + min_scale − effective_scale[i])
```

Where:
- `base_geb` = configurable base value (default -3.5)
- `max_scale`, `min_scale` = existing gradient bounds (default 6.0, 1.0)
- `effective_scale[i]` = combined gradient scale from all cutsites at position i

At key positions:

| Position | effective_scale | geb (base=-3.5) | Effect |
|----------|:---------------:|:----------------:|--------|
| Cutsite center | 1.0 | -3.5 × 6 = **-21.0** | Very hard to exit gap |
| Cutsite edge | 2.0 | -3.5 × 5 = **-17.5** | Strong |
| Mid-gradient (~3.5) | 3.5 | -3.5 × 3.5 = **-12.25** | Moderate |
| Conservative region | 6.0 | -3.5 × 1 = **-3.5** | Mild |

This ensures: at cutsites, exiting a gap costs ~21 extra points — enough to overcome any frame shift benefit (>20). At conservative regions, the exit penalty is mild (-3.5), preserving normal alignment behavior.

### Validation

Tested with the problematic query `seq2129` (Target1→Target7 deletion with isolated A at Target2 cutsite):

| Config | Score | Non-gap blocks | Non-gap block sizes | Isolated match? |
|--------|:-----:|:--------------:|:-------------------:|:---------------:|
| Current (scalar geb=-1.0) | 259.9 | 2 | [40, 125] | No |
| **With geb_profile (base=-3.5)** | **240.6** | **2** | **[40, 125]** | **No** |

Both produce a single continuous gap between positions 40 and 206 (167 gap chars), which is the correct Target1→Target7 deletion. The current scalar also prevents this case (because all features were enabled in the test), but the profile approach is more fundamental — it works at the DP transition level regardless of other feature settings.

The lower score with geb_profile reflects that the biologically correct path (continuous deletion) is taken instead of the higher-scoring but spurious frame-shifted path.

## Parameter Design

### New: `--geb-base` (CLI) / `geb_base` (PipelineConfig)

- Type: float
- Default: -3.5
- Meaning: controls the base strength of the position-dependent gap exit penalty
- The profile is auto-computed from this base + existing `--min-scale`/`--max-scale`
- Set to 0.0 to disable (equivalent to no gap exit penalty)

Users do **not** need to tune min/max separately — the existing gradient bounds already determine the profile shape. `--geb-base` is the single "how aggressive" knob.

### Removed / Superseded

- **`--gap-exit-bonus`**: Replaced by geb_profile. A uniform scalar cannot distinguish cutsite (suspicious) from conservative (normal) regions.

### Kept (unchanged)

- **`--isolated-base-penalty`**: Operates at the substitution score level (reduces match scores where `run_len == 1`). Still useful as a safety net for truly isolated matches, though geb_profile makes it less critical.
- **`--short-match-window` / `--short-match-discount`**: Orthogonal feature targeting short match runs.
- **`--dense-mismatch-penalty`**: Orthogonal, targets a different pathology.

## Implementation

### 1. `crisviper/lineage.py` — `build_gradient_profiles()`

Add geb_profile as a 4th return value:

```python
def build_gradient_profiles(
    ref_length: int,
    cutsites: List[CutsiteRegion],
    base_gap_open: float = -2.0,
    base_gap_extend: float = -0.1,
    mismatch_penalty: float = -3.0,
    min_scale: float = 1.0,
    max_scale: float = 6.0,
    cutsite_edge_scale: float = 2.0,
    gradient_radius: Optional[float] = None,
    base_gap_exit: float = -3.5,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
```

The profile only depends on the `effective_scale` array already computed:

```python
# geb_profile: inverted scaling — strongest at cutsite center
geb_profile = base_gap_exit * (max_scale + min_scale - effective_scale)
```

### 2. `crisviper/alignment.py` — `affine_gap_alignment_position_aware()`

Replace the scalar `gap_exit_bonus` parameter with `gap_exit_bonus_profile: Optional[np.ndarray] = None`.

In the DP recurrence (simplified, previously scalar):

```python
if gap_exit_profile_active:
    geb_i = gap_exit_bonus_profile[i-1]  # profile value at current row
    for j in range(1, n + 1):
        m_best = max(Mi_1[j-1],
                     Ixi_1[j-1] + geb_i,
                     Iyi_1[j-1] + geb_i)
        Mi[j] = s_row[j-1] + m_best
```

Same change in backtrace:

```python
if gap_exit_active:
    geb_i = gap_exit_bonus_profile[i-1]
    scores = [M[i-1, j-1],
              Ix[i-1, j-1] + geb_i,
              Iy[i-1, j-1] + geb_i]
```

The Ix/Iy recurrences are unchanged — only the M-state recurrence and M backtrace use the value.

### 3. `crisviper/lineage.py` — `lineage_tracer_align()`

Pass the fourth return value from `build_gradient_profiles()`:

```python
gap_open_profile, gap_extend_profile, mismatch_profile, geb_profile = \
    build_gradient_profiles(..., base_gap_exit=base_gap_exit)
```

### 4. `crisviper/models.py` — `PipelineConfig`

Add field:

```python
gap_exit_base: float = -3.5  # replaced gap_exit_bonus
```

### 5. CLI (`crisviper/cli.py`)

Replace `--gap-exit-bonus` with `--geb-base` (default -3.5).

### 6. Tests

- Add unit test: verify geb_profile values at known positions match expectation
- Add integration test: run seq2129 with geb_profile, verify no isolated match
- Update existing tests that reference `gap_exit_bonus`

## Scope

All 290 identified isolated segments are covered:
- 187 inside cutsite regions → directly at peak penalty (-21)
- 103 within 12 bp of cutsite center → well within gradient radius (~15 bp), penalty ≥ -12

No position outside the gradient zone has this artifact. No new false positives expected.

## Backward Compatibility

- `--gap-exit-bonus` was default 0.0 (disabled). `--geb-base` defaults to -3.5 (enabled). Users upgrading will see different alignment results specifically for sequences with large deletions near cutsites — this is the intended fix.
- `PipelineConfig.gap_exit_bonus` → `PipelineConfig.gap_exit_base`. Old field name kept as deprecated alias.
- `isolated_base_penalty`, `short_match_window`, `dense_mismatch_penalty` — all unchanged.
