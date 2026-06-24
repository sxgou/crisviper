# Changelog

## v1.2.0 (2026-06-24)

### Bug fixes
- **DP alignment**: Fixed Ix/Iy termination state — backtrace now starts from max(M, Ix, Iy) instead of M alone (#1)
- **Dense mismatch density**: Fixed diagonal end calculation — denom now correctly uses min(..., 2*half_w+1) (#2)
- **`_check_primer_quality`**: Fixed p3 match counting — range now correctly limited to primer region (N6)
- **Allele heatmap**: Fixed empty internal_key collision on short sequences (#3)
- **`_target_num`**: Fixed crash on target names without digits (#4)
- **`cutsites_from_list`**: Added explicit start/end key validation with clear ValueError (#7)
- **`_assemble_full_length`**: Fixed p3=0 edge case (#10)
- **BrokenProcessPool**: Fixed silent data loss — remaining chunks now reprocessed after pool crash (#14)
- **Empty QC string**: Fixed all() returning True on empty quality strings (#17)
- **`assert` → `ValueError`**: Replaced 3 `assert` statements with explicit `if` checks (caller.py, denoiser.py) to prevent suppression under `-O` (C1, D1)
- **`max_elem=0` guard**: Added ZeroDivisionError protection in `compute_threshold` (T1, T2)
- **`total_reads_success`**: Changed direct key access to `.get("readCount", 1)` in summary.py (#5)

### Improvements
- **`save_tsv`**: Dynamic fieldnames — now auto-derives from row keys, supports `original_read_names` (#6/N1)
- **CutsiteRegion**: Added integer type validation on start/end fields (N6)
- **GAGTCG threshold**: Relaxed from ≥3 to ≥1 motif for custom short amplicons (#24/N2)
- **`exclusive_end`**: Added INDEL handling alongside DELETION (M2)
- **`unmutated_reads`**: Added explicit override in HTML report (N3)
- **`run_len` condition**: Removed unnecessary gap_exit_bonus/dense_mismatch from needs_run_len (N1)
- **Profile broadcast protection**: Added length validation + np.pad for profile arrays (N3)
- **Read-to-allele mapping**: Added FASTQ direct input support
- **`merge_parent`**: Extracted shared argument parser, eliminated ~40 lines of duplicate CLI code (#35)
- **Code duplication**: Refactored `call_alleles_coarse_grain`/`call_alleles_exact` into shared `_call_alleles_by_key` (C2)

### Documentation
- ALGORITHM.md: Fixed --global flag reference, unified smoothstep description
- ALGORITHM.md: Unified cutsite JSON format example with CLI help
- USER_GUIDE.md: Updated chunk-size and allele-window-end default descriptions
- CLI help: Updated --chunk-size to "(PipelineConfig default: 500)"
- CLI help: Unified --allele-window-end to "inclusive, 0-indexed"

### Testing
- Added tests/test_summary.py: 22 new unit tests (allele_label, split_indel, mutation_type_label, mutation_overlaps_target)
- Added concurrency consistency test (parallel vs sequential)
- Enhanced assertions across 5 test files (precise assertion replacement)
- Renamed test_empty_edge → test_single_element_umi
- Total: 237 tests (+23)

### Lint fixes
- Removed 9 unused variables (F841): bpp_seen, diffs, ref_base, seg_len, w, alt_bg_rgb, avg_ins/del_len_reads, total_ins/del_reads
- Removed 17 unused imports (F401) across cli.py, config.py, io.py, pipeline.py, reporting.py
- Fixed f-string backslash compatibility for Python 3.8

## v1.1.0 (2026-05-15)

- Added lineage mode with gradient-based gap penalties
- Added HTML report generation with allele heatmaps
- Added FASTQ direct input and paired-end merging
- Added read-to-allele mapping
- Added summary tables (allele frequency, per-target editing, filter reasons, indel length distributions)
- Added Numba JIT acceleration for DP fill (28x speedup)
- Improved heatmap rendering performance (imshow, 24x speedup)
- i18n: Translated output to English
- Initial public release
