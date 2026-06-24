"""Run full validation: align all examples data with geb_profile, then analyze."""
import sys, os, json, time

# Use available CPUs, cap at 12 for reasonable parallelism
CPU_COUNT = min(12, os.cpu_count() or 1)

# Must set fork BEFORE importing numpy or crisviper
import multiprocessing as mp
try:
    mp.set_start_method('fork')
except RuntimeError:
    pass

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from crisviper import (
    PipelineConfig, Pipeline, QueryRecord,
    read_reference_fasta, read_queries_tsv,
)

def find_isolated(aq, ar):
    if not aq or not ar:
        return []
    segs = []
    i = 0
    while i < len(aq):
        if aq[i] != '-':
            s = i
            while i < len(aq) and aq[i] != '-':
                i += 1
            l = i - s
            if l <= 4:
                left_gap = s > 0 and aq[s - 1] == '-'
                right_gap = i < len(aq) and aq[i] == '-'
                if left_gap and right_gap and any(ar[j] != '-' for j in range(s, i)):
                    segs.append((s, l))
        else:
            i += 1
    return segs


# ── Step 1: Align with geb_profile ──
ref_seq = read_reference_fasta("examples/reference.fa")
queries = read_queries_tsv("examples/output/all_queries.tsv")
records = [QueryRecord(**q) for q in queries]
print(f"Total: {len(records)} reads")

config = PipelineConfig(
    lineage_mode=True,
    gap_exit_strength=-3.5,
    threads=CPU_COUNT,
    min_reads_sub=0, min_reads_indel=0,
)

t0 = time.time()
pipeline = Pipeline(config, ref_seq)
result = pipeline.run(records)
t1 = time.time()

output = [r.to_dict() for r in result.results]
with open("examples/output/geb_alignment.json", "w") as f:
    json.dump(output, f, indent=2)

success = sum(1 for r in result.results if r.success)
failed = sum(1 for r in result.results if not r.success)
print(f"Done: {success} success, {failed} failed, {t1 - t0:.1f}s")

# ── Step 2: Analyze for isolated matches ──
cs_centers = [44, 71, 98, 125, 152, 179, 206, 233, 260, 287]
cs_names = ["Target1", "Target2", "Target3", "Target4", "Target5",
            "Target6", "Target7", "Target8", "Target9", "Target10"]

iso_count = 0
iso_seqs = []
total_iso_segments = 0
length_dist = {}
dist_to_cs = {}
by_cs = {}

for r in output:
    aq = r.get("aligned_query")
    ar = r.get("aligned_ref")
    segs = find_isolated(aq, ar)
    if segs:
        iso_count += 1
        total_iso_segments += len(segs)
        iso_seqs.append((r["readName"], segs))
        for pos, length in segs:
            length_dist[length] = length_dist.get(length, 0) + 1
            dists = [abs(pos - c) for c in cs_centers]
            min_d = min(dists)
            min_cs = cs_names[dists.index(min_d)]
            by_cs[min_cs] = by_cs.get(min_cs, 0) + 1
            if min_d == 0:
                dist_to_cs["0bp"] = dist_to_cs.get("0bp", 0) + 1
            elif min_d <= 3:
                dist_to_cs["1-3bp"] = dist_to_cs.get("1-3bp", 0) + 1
            elif min_d <= 6:
                dist_to_cs["4-6bp"] = dist_to_cs.get("4-6bp", 0) + 1
            elif min_d <= 12:
                dist_to_cs["7-12bp"] = dist_to_cs.get("7-12bp", 0) + 1
            else:
                dist_to_cs[">12bp"] = dist_to_cs.get(">12bp", 0) + 1

print(f"\n=== Analysis Results ===")
print(f"Sequences with isolated matches: {iso_count} / {success}")
print(f"Total isolated segments: {total_iso_segments}")
print(f"Length distribution: {length_dist}")
print(f"Distance to cutsite: {dist_to_cs}")
print(f"By cutsite: {by_cs}")

if iso_seqs:
    print(f"\n=== Top affected sequences ===")
    for name, segs in sorted(iso_seqs, key=lambda x: -len(x[1]))[:20]:
        seg_str = "; ".join([f"pos={p} len={l}" for p, l in segs])
        print(f"  {name}: {seg_str}")
else:
    print("\nNo isolated match segments found!")

# Save analysis
with open("examples/output/geb_analysis.json", "w") as f:
    json.dump({
        "total_sequences": success,
        "with_isolated_matches": iso_count,
        "total_isolated_segments": total_iso_segments,
        "length_distribution": length_dist,
        "distance_to_cutsite": dist_to_cs,
        "by_cutsite": by_cs,
    }, f, indent=2)

print(f"\nResults saved to examples/output/")
