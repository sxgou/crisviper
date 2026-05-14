#!/usr/bin/env python3
"""Quick test runner for core algorithms. Run with: python run_tests.py"""
import sys, os
sys.path.insert(0, '.')

passed = 0
failed = 0

def check(name, condition):
    global passed, failed
    if condition:
        passed += 1
        print(f"  ✓ {name}")
    else:
        failed += 1
        print(f"  ✗ {name} FAILED")

print("Testing affine_gap_alignment core...\n")

from ltlib import (
    count_gap_blocks, calculate_alignment_stats, get_amplicon_structure,
    CutsiteRegion, build_gap_penalty_profile, remove_isolated_matches,
    convert_dense_mismatch_to_indel, affine_gap_alignment,
    affine_gap_alignment_position_aware, lineage_tracer_align,
)
from tests.test_accuracy import run_accuracy_checks

# ── 1. count_gap_blocks ──
check("count_gap_blocks empty", count_gap_blocks("ACGT") == [])
check("count_gap_blocks single", count_gap_blocks("A--T") == [2])
check("count_gap_blocks double", count_gap_blocks("A--T--G") == [2, 2])
check("count_gap_blocks leading", count_gap_blocks("---A") == [3])
check("count_gap_blocks trailing", count_gap_blocks("A---") == [3])

# ── 2. calculate_alignment_stats ──
stats = calculate_alignment_stats("ACGT", "ACGT")
check("stats identical matches", stats["matches"] == 4)
check("stats identical mismatches", stats["mismatches"] == 0)
check("stats identical similarity", stats["similarity"] == 1.0)

stats = calculate_alignment_stats("ACGT", "A-GT")
check("stats with gap matches", stats["matches"] == 3)
check("stats with gap gaps_in_query", stats["gaps_in_query"] == 1)
check("stats with gap alignment_length", stats["alignment_length"] == 4)

stats = calculate_alignment_stats("ACGT-", "A-GTA")
check("stats complex matches", stats["matches"] == 3)
check("stats complex gaps_ref", stats["gaps_in_ref"] == 1)
check("stats complex gaps_qry", stats["gaps_in_query"] == 1)

# ── 3. Standard Gotoh ──
score, ar, aq, s = affine_gap_alignment("ACGT", "ACGT")
check("gotoh identical ar", ar == "ACGT")
check("gotoh identical aq", aq == "ACGT")
check("gotoh identical score", score == 8.0)

score, ar, aq, s = affine_gap_alignment("ACGT", "ACCT")
check("gotoh sub matches", s["matches"] == 3)
check("gotoh sub mismatches", s["mismatches"] == 1)

score, ar, aq, s = affine_gap_alignment("ACGT", "AGT")
check("gotoh del matches", s["matches"] == 3)
check("gotoh del gaps_qry", s["gaps_in_query"] >= 1)

score, ar, aq, s = affine_gap_alignment("AAAACGT", "CGT", semi_global=True)
check("gotoh semi_global", "CGT" in aq.replace("-", ""))

s_sg, _, _, _ = affine_gap_alignment("AAACCC", "CCC", semi_global=True)
s_g, _, _, _ = affine_gap_alignment("AAACCC", "CCC", semi_global=False)
check("gotoh sg vs global", s_sg >= s_g)

score, ar, aq, s = affine_gap_alignment("AAAAACGTGGGGG", "ACGT", fit_mode=True)
check("gotoh fit_mode", aq.replace("-", "") == "ACGT")

# ── 4. build_gap_penalty_profile ──
cs = [CutsiteRegion("T1", 5, 10)]
go, ge = build_gap_penalty_profile(20, cs,
    base_gap_open=-2.0, base_gap_extend=-0.1,
    cutsite_scale=1.0, flank_scale=2.0, far_scale=2.0, flank_width=3)
check("profile length go", len(go) == 20)
check("profile length ge", len(ge) == 20)
check("profile cutsite", all(abs(go[i] - (-2.0)) < 0.01 for i in range(5, 11)))
# far_scale=2.0 → far region = -2.0 * 2.0 = -4.0
check("profile far", abs(go[0] - (-4.0)) < 0.01 and abs(go[15] - (-4.0)) < 0.01)
check("profile flank", abs(go[2] - (-4.0)) < 0.01 and abs(go[13] - (-4.0)) < 0.01)

# ── 5. Position-aware alignment ──
go2, ge2 = build_gap_penalty_profile(8, [CutsiteRegion("T1", 2, 5)],
    base_gap_open=-2.0, base_gap_extend=-0.1)
score, ar, aq, s = affine_gap_alignment_position_aware("ACGTACGT", "ACGTACGT", go2, ge2)
check("pos-aware identical", ar.replace("-", "") == "ACGTACGT" and s["matches"] == 8)

# ── 6. remove_isolated_matches ──
ar2, aq2 = "ACGT", "A-G-"
_, result_qry, modified = remove_isolated_matches(ar2, aq2)
check("isolated matches modified", modified)
check("isolated matches result", result_qry == "A---")

# ── 7. convert_dense_mismatch_to_indel ──
cr, cq, mod = convert_dense_mismatch_to_indel("ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTACGT", threshold=0.34)
check("dense mismatch no change", not mod)

# Dense region should trigger: 8 mismatches in first 8 bases (100% density)
ar3 = "ACGTACGTACGT"
aq3 = "TTTTTTTTACGT"
cr, cq, mod = convert_dense_mismatch_to_indel(ar3, aq3, ar3, aq3, threshold=0.34)
check("dense mismatch triggers", mod and cr.count('-') > ar3.count('-'))

# ── 8. get_amplicon_structure ──
carlin_ref = (
    "TATGTGTGGGAGGGCTAAGAGGCCGCCGGACTGCACGACAGTCGACGATGGAGTCGACACGACTCGCGCATACGATGGAGTCGACTACAGTCGCTACGACGATGGAGTCGCGAGCGCTATGAGCGACTATGGAGTCGATACGATACGCGCACGCTATGGAGTCGAGAGCGCGCTCGTCAACGATGGAGTCGCGACTGTACGCACTCGCGATGGAGTCGATAGTATGCGTACACGCGATGGAGTCGACTGCACGACAGTCGACTATGGAGTCGATACGTAGCACGCACATGATGGGAGCTAGCTGTGCCTTCTAGTTGCCAGCCATCTGTTGT"
)
check("carlin ref length", len(carlin_ref) == 332)
cutsites = get_amplicon_structure(carlin_ref)
check("cutsite count", len(cutsites) == 10)
check("cutsite T1 start", cutsites[0].start == 41)
check("cutsite T1 end", cutsites[0].end == 47)
check("cutsite T10 start", cutsites[9].start == 284)
check("cutsite T10 end", cutsites[9].end == 290)

# ── 9. lineage_tracer_align smoke test ──
ref4 = "A" * 23 + "CGCCG" + "A" * 13 + "GAGTCGA" + "A" * 7 + "A" * 13 + "GAGTCGA" + "A" * 8 + "A" * 33
query4 = ref4
cs4 = [CutsiteRegion("Target1", 41, 47), CutsiteRegion("Target2", 68, 74)]
score, ar, aq, s = lineage_tracer_align(ref4, query4, cs4)
check("lineage identical matches", s["matches"] > 0)
check("lineage identical no mismatch", s["mismatches"] == 0)
check("lineage identical no gaps", s["gaps_in_query"] == 0 and s["gaps_in_ref"] == 0)

# ── 10. lineage_tracer with deletion ──
query5 = "A" * 23 + "CGCCG" + "A" * 13 + "GA" + "A" * 7 + "A" * 13 + "GAGTCGA" + "A" * 8 + "A" * 33
score, ar, aq, s = lineage_tracer_align(ref4, query5, cs4)
check("lineage deletion detected", s["gaps_in_query"] > 0)

# ── 11. 准确性验证测试 ──
run_accuracy_checks(check)

print(f"\n{'='*50}")
print(f"Results: {passed} passed, {failed} failed out of {passed+failed}")
if failed == 0:
    print("All tests passed ✓")
else:
    print(f"{failed} test(s) FAILED ✗")
    sys.exit(1)
