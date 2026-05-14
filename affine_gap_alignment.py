#!/usr/bin/env python3
"""
affine_gap_alignment — 兼容模块
===============================
此文件保留向后兼容。所有代码已迁移到 ltlib/ 包中。
"""
import sys
import warnings
from ltlib.logging_config import get_logger

log = get_logger(__name__)
warnings.warn(
    "affine_gap_alignment is deprecated; use ltlib instead",
    DeprecationWarning, stacklevel=2
)

from ltlib.alignment import (
    affine_gap_alignment, calculate_alignment_stats,
    count_gap_blocks, affine_gap_alignment_position_aware
)
from ltlib.lineage import (
    lineage_tracer_align, build_gap_penalty_profile,
    get_amplicon_structure,
)
from ltlib.corrections import (
    convert_dense_mismatch_to_indel, filter_point_mutations,
    correct_repetitive_misalignment, correct_target_misalignments,
    remove_isolated_matches
)
from ltlib.config import CutsiteRegion, AmpliconConfig


def print_alignment(aligned_ref, aligned_query, width=100):
    log.info("\n比对结果:")
    for i in range(0, len(aligned_ref), width):
        ref_slice = aligned_ref[i:i + width]
        query_slice = aligned_query[i:i + width]
        log.info("Ref  %4d: %s", i, ref_slice)
        log.info("Query%4d: %s", i, query_slice)
        match_line = " " * 9
        for r, q in zip(ref_slice, query_slice):
            if r == q and r != '-':
                match_line += '|'
            elif r != '-' and q != '-':
                match_line += '.'
            else:
                match_line += ' '
        log.info(match_line)
        log.info("")


def print_stats(stats):
    log.info("\n比对统计:")
    log.info("  比对得分: %.2f", stats.get('score', 0))
    log.info("  比对长度: %d", stats.get('alignment_length', 0))
    log.info("  匹配数: %d", stats.get('matches', 0))
    log.info("  错配数: %d", stats.get('mismatches', 0))
    log.info("  相似度: %.2f%%", stats.get('similarity', 0) * 100)
    log.info("  一致性: %.2f%%", stats.get('identity', 0) * 100)
    log.info("  Reference中gap数: %d", stats.get('gaps_in_ref', 0))
    log.info("  Query中gap数: %d", stats.get('gaps_in_query', 0))


def main():
    from ltlib.lineage import main as lineage_main
    lineage_main()


if __name__ == "__main__":
    main()
