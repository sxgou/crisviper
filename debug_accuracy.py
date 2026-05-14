"""调试脚本：追踪 align_single 对已知突变的处理流程"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ltlib import (
    Pipeline, PipelineConfig, QueryRecord,
    align_single,
    get_amplicon_structure, CutsiteRegion,
    extract_mutations,
    MutationEvent, MutationType,
)

# CARLIN reference
CARLIN_REF = (
    "TATGTGTGGGAGGGCTAAGAGG"   # Primer5 (23bp)
    "CCGCC"                     # Prefix (5bp)
    "GACTGCACGACAGTCGA"         # Target1: 保守区13bp + cutsite 7bp
    "CGATGGAG"                  # Linker (7bp)
    "TCGACACGACTCGCGCA"         # Target2
    "TACGATGG"                  # Linker
    "AGTCGACTACAGTCGCTA"        # Target3
    "CGACGATG"                  # Linker
    "GAGTCGCGAGCGCTATG"         # Target4
    "AGCGACTA"                  # Linker
    "TGGAGTCGATACGATACG"        # Target5
    "CGCACGCT"                  # Linker
    "ATGGAGTCGAGAGCGCGC"        # Target6
    "TCGTCAAC"                  # Linker
    "GATGGAGTCGCGACTGTA"        # Target7
    "CGCACTCG"                  # Linker
    "CGATGGAGTCGATAGTAT"        # Target8
    "GCGTACAC"                  # Linker
    "GCGATGGAGTCGACTGCA"        # Target9
    "CGACAGTC"                  # Linker
    "GACTATGGAGTCGATACGTAGC"    # Target10
    "ACGCACATGATGGGAGCTAGCTGTGCCTTCTAGTTGCCAGCCATCTGTTGT"  # Postfix(8) + Primer3(33)
)

def _make_deletion(seq, pos, length):
    return seq[:pos] + seq[pos+length:]

def _make_insertion(seq, pos, insert):
    return seq[:pos] + insert + seq[pos:]

def _make_substitution(seq, pos, new_base):
    lst = list(seq)
    lst[pos] = new_base
    return ''.join(lst)

print("=" * 70)
print("CARLIN_REF length:", len(CARLIN_REF))
print("Region 38-52:", CARLIN_REF[38:52])
print("  index 42:", CARLIN_REF[42] if 42 < len(CARLIN_REF) else "N/A")
print("  index 44:", CARLIN_REF[44] if 44 < len(CARLIN_REF) else "N/A")
print()

# Test 1: 3bp deletion at position 42
del_query = _make_deletion(CARLIN_REF, 42, 3)
print("=== 3bp DEL at pos 42 ===")
print(f"  Original (40-50): {CARLIN_REF[40:50]}")
print(f"  Deleted  (40-50): {del_query[40:50]}")
print(f"  Length: orig={len(CARLIN_REF)}, del={len(del_query)}")

# Test 2: 5bp insertion at position 42
ins_query = _make_insertion(CARLIN_REF, 42, "ACGTA")
print("\n=== 5bp INS at pos 42 ===")
print(f"  Original (40-50): {CARLIN_REF[40:50]}")
print(f"  Inserted (40-55): {ins_query[40:55]}")
print(f"  Length: orig={len(CARLIN_REF)}, ins={len(ins_query)}")

# Test 3: substitution at position 44
sub_query = _make_substitution(CARLIN_REF, 44, "T")
print("\n=== SUB at pos 44 ===")
print(f"  Original (40-50): {CARLIN_REF[40:50]}")
print(f"  Subbed   (40-50): {sub_query[40:50]}")
print(f"  Length: orig={len(CARLIN_REF)}, sub={len(sub_query)}")

# Now run through the actual pipeline with detailed tracing
print("\n" + "=" * 70)
print("RUNNING PIPELINE - 3bp DELETION")
print("=" * 70)

config = PipelineConfig(
    lineage_mode=True,
    primer5_len=23,
    primer3_len=33,
    primer5_threshold=19,
    primer3_threshold=29,
    min_reads_snv=1,
    min_reads_indel=1,
)
query = QueryRecord(
    readName="del_3bp", cellBC="test", UMI="UMI",
    readCount=50, seq=del_query,
)
cutsites = get_amplicon_structure(CARLIN_REF)

# Manually trace through align_single (新方案流程)
q_seq = del_query
r_seq = CARLIN_REF
p5 = 23
p3 = 33

# Step 1: Full-length global alignment (332bp)
from ltlib.pipeline import _align_full_lineage
score, ar, aq, raw_stats = _align_full_lineage(r_seq, q_seq, cutsites, config)
print(f"\nFull-length alignment:")
print(f"  score={score}, len(ar)={len(ar)}, len(aq)={len(aq)}")
print(f"  raw_stats: matches={raw_stats['matches']}, mismatches={raw_stats['mismatches']}, "
      f"gaps_in_ref={raw_stats['gaps_in_ref']}, gaps_in_query={raw_stats['gaps_in_query']}")

# Step 2: Post-alignment primer quality check
from ltlib.pipeline import _check_primer_quality
p5_ok, p5_match, p3_ok, p3_match = _check_primer_quality(
    ar, aq, r_seq, p5, p3,
    config.primer5_threshold, config.primer3_threshold,
)
print(f"\nPrimer quality: p5={p5_match}/{p5} {'OK' if p5_ok else 'FAIL'}, "
      f"p3={p3_match}/{p3} {'OK' if p3_ok else 'FAIL'}")

# Step 3: Extract internal region
from ltlib.pipeline import _extract_internal_region
ar_int, aq_int, int_r = _extract_internal_region(ar, aq, r_seq, p5, p3)
print(f"\nInternal region:")
print(f"  int_r len={len(int_r)}: {int_r[:60]}...")
print(f"  ar_int len={len(ar_int)}: {ar_int[:60]}...")
print(f"  aq_int len={len(aq_int)}: {aq_int[:60]}...")
print(f"  int_r[19:30]: {int_r[19:30]}")
print(f"  aq_int[19:30]: {aq_int[19:30]}")

# Adjust cutsites to internal coordinates
from ltlib.pipeline import _adjust_cutsites_to_internal
int_r_len = len(r_seq) - p5 - p3
internal_cutsites = _adjust_cutsites_to_internal(cutsites, p5, int_r_len)
print(f"\nInternal cutsites:")
for cs in internal_cutsites:
    print(f"  {cs.name}: {cs.start}-{cs.end}")
print(f"Target1 in internal coords: ({internal_cutsites[0].start}, {internal_cutsites[0].end})")

# Show alignment around deletion area
if ar_int and aq_int:
    print(f"\n  Alignment around position 18-25 (Target1 cutsite):")
    print(f"  ar: {ar_int[13:30]}")
    print(f"  aq: {aq_int[13:30]}")

    # Try extract_mutations directly on internal alignment
    print("\n  extract_mutations directly on alignment:")
    muts = extract_mutations(ar_int, aq_int, cutsites=internal_cutsites, mutation_window=3)
    for m in muts:
        print(f"    {m.type}: ref_pos={m.ref_pos}, length={m.length}, in_window={m.in_cutsite_window}")

# Now do full align_single
print("\n" + "-" * 50)
print("Full align_single result:")
result = align_single(query, CARLIN_REF, config, cutsites)
print(f"  success: {result.success}")
if result.error:
    print(f"  error: {result.error}")
print(f"  score: {result.score}")
print(f"  mutations:")
for m in result.mutations:
    print(f"    {m.type}: ref_pos={m.ref_pos}, length={m.length}, in_window={m.in_cutsite_window}")

# Also try with lineage_mode=False to see if it's a lineage mode issue
print("\n" + "=" * 70)
print("WITH lineage_mode=False")
print("=" * 70)
config_std = PipelineConfig(
    lineage_mode=False,
    primer5_len=23,
    primer3_len=33,
    primer5_threshold=19,
    primer3_threshold=29,
    min_reads_snv=1,
    min_reads_indel=1,
)
result_std = align_single(query, CARLIN_REF, config_std, cutsites=None)
print(f"  success: {result_std.success}")
if result_std.error:
    print(f"  error: {result_std.error}")
print(f"  score: {result_std.score}")
print(f"  mode: {result_std.mode}")
print(f"  mutations:")
for m in result_std.mutations:
    print(f"    {m.type}: ref_pos={m.ref_pos}, length={m.length}, in_window={m.in_cutsite_window}")

# Show aligned sequences
if result_std.aligned_ref and result_std.aligned_query:
    print(f"\n  Aligned internal region:")
    # Find where the deletion might be
    for i, (ar, aq) in enumerate(zip(result_std.aligned_ref, result_std.aligned_query)):
        if aq == '-':
            print(f"    Gap at col {i}: ref={ar}")
            break
