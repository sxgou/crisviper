"""crisviper/summary.py — Summary statistics table generation module.

Produces 4 TSV summary tables saved to the output directory:
  - allele_frequency.tsv    Allele frequency table aggregated by mutation fingerprint
  - per_target_editing.tsv  Editing statistics per target site
  - filter_reason.tsv       Sequence discard reason breakdown
  - indel_length.tsv        Indel length distribution (grouped by actual length)
  - deletion_length.tsv     Deletion length distribution
  - insertion_length.tsv    Insertion length distribution
  - event_level_details.tsv Per-event statistics table
"""

import csv
import os
from typing import List, Dict, Optional
from collections import Counter

from crisviper.logging_config import get_logger

log = get_logger(__name__)


def save_summary_tables(results: List[Dict], output_dir: str,
                         ref_seq: str = "", cutsites: Optional[list] = None,
                         read_to_allele_path: Optional[str] = None,
                         target_region_left: int = 13,
                         target_region_right: int = 7) -> Dict:
    """Generate 5 summary TSV tables and return pre-computed summary statistics.

    Tables produced:
      1. allele_frequency.tsv  — Alleles sorted by frequency with read counts/percentages
      2. per_target_editing.tsv — Per-target mutation rates and type breakdowns
      3. filter_reason.tsv     — Sequence discard reason counts
      4a. deletion_length.tsv  — Deletion length distribution
      4b. insertion_length.tsv — Insertion length distribution
      5. event_level_details.tsv — Detailed per-event breakdown

    When original_read_names are present in results and read_to_allele_path
    is provided, also writes a read-to-allele mapping table.

    Args:
        results: List of alignment result dicts.
        output_dir: Output directory.
        ref_seq: Reference sequence (for cutsite detection).
        cutsites: CutsiteRegion list (per_target table is empty when absent).
        read_to_allele_path: Path for read-to-allele mapping TSV (optional).
        target_region_left: Bases left of cutsite start to extend target region
            (default 13 = CARLIN conserved region length).
        target_region_right: Bases right of cutsite end to extend target region
            (default 7 = CARLIN linker length).

    Returns:
        Dict with pre-computed statistics for HTML report.
    """
    os.makedirs(output_dir, exist_ok=True)
    successful = [r for r in results if "error" not in r]
    failed = [r for r in results if "error" in r]
    total_reads_success = sum(r["readCount"] for r in successful)

    # ── Helper: Parse INDEL events into separate del/ins descriptions ──
    def _split_indel(m: Dict) -> List[str]:
        """Parse an INDEL mutation into a list of independent del/ins descriptions."""
        ref_base = m.get("ref_base", "")
        query_base = m.get("query_base", "")
        ref_pos = m.get("ref_pos", 0)
        if not ref_base or not query_base:
            return [f"indel:{ref_pos}-{m.get('length', 1)}bp"]
        parts = []
        i = 0
        while i < len(ref_base):
            rb = ref_base[i]
            qb = query_base[i]
            if rb != '-' and qb == '-':
                del_start = ref_pos
                del_len = 0
                while i < len(ref_base) and ref_base[i] != '-' and query_base[i] == '-':
                    ref_pos += 1
                    del_len += 1
                    i += 1
                parts.append(f"del:{del_start}-{del_len}bp")
            elif rb == '-' and qb != '-':
                ins_start = ref_pos
                ins_start_idx = i
                while i < len(ref_base) and ref_base[i] == '-' and query_base[i] != '-':
                    i += 1
                ins_seq = query_base[ins_start_idx:i]
                parts.append(f"ins:{ins_start}+{ins_seq}")
            else:
                if rb != '-':
                    ref_pos += 1
                i += 1
        return parts if parts else [f"indel:{ref_pos}-{m.get('length', 1)}bp"]

    def _allele_label(mutations: List[Dict]) -> str:
        """Generate a compact allele label from a list of mutation dicts."""
        if not mutations:
            return "wt"
        parts = []
        for m in sorted(mutations, key=lambda x: (x.get("ref_pos", 0), x.get("type", ""), x.get("length", 0))):
            typ = m.get("type", "")
            pos = m.get("ref_pos", -1)
            if typ == "deletion":
                parts.append(f"del:{pos}-{m.get('length', 1)}bp")
            elif typ == "insertion":
                seq = m.get('query_base', '')
                if seq:
                    parts.append(f"ins:{pos}+{seq}")
                else:
                    parts.append(f"ins:{pos}+{m.get('length', 1)}bp")
            elif typ == "indel":
                parts.extend(_split_indel(m))
            elif typ == "substitution":
                parts.append(f"sub:{pos}{m.get('ref_base', '?')}>{m.get('query_base', '?')}")
            else:
                parts.append(f"{typ}:{pos}")
        return ";".join(parts) if parts else "wt"

    def _mutation_type_label(mutations: List[Dict]) -> str:
        """Classify a list of mutations into a type category label."""
        has_del = any(m.get("type") == "deletion" for m in mutations)
        has_ins = any(m.get("type") == "insertion" for m in mutations)
        has_indel = any(m.get("type") == "indel" for m in mutations)
        has_sub = any(m.get("type") == "substitution" for m in mutations)
        if not (has_del or has_ins or has_indel or has_sub):
            return "wt"
        parts = []
        if has_del: parts.append("del")
        if has_ins: parts.append("ins")
        if has_indel: parts.append("indel")
        if has_sub: parts.append("sub")
        return "+".join(parts)

    # ═══════════════════════════════════════════════════════════════
    # 1. allele_frequency.tsv
    # ═══════════════════════════════════════════════════════════════
    allele_groups: Dict[str, dict] = {}
    for r in successful:
        muts = r.get("mutations", [])
        label = _allele_label(muts)
        mtype = _mutation_type_label(muts)
        rc = r.get("readCount", 1)
        if label not in allele_groups:
            allele_groups[label] = {"allele": label, "mutation_type": mtype,
                                    "sequences": 0, "reads": 0}
        allele_groups[label]["sequences"] += 1
        allele_groups[label]["reads"] += rc

    sorted_alleles = sorted(allele_groups.values(), key=lambda x: -x["reads"])
    path_af = os.path.join(output_dir, "allele_frequency.tsv")
    with open(path_af, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Rank", "Allele", "Mutation_Type",
                                                "Sequences", "Reads", "Reads_Pct"],
                                delimiter='\t')
        writer.writeheader()
        for rank, a in enumerate(sorted_alleles, 1):
            pct = a["reads"] / total_reads_success * 100 if total_reads_success > 0 else 0
            writer.writerow({
                "Rank": rank, "Allele": a["allele"],
                "Mutation_Type": a["mutation_type"],
                "Sequences": a["sequences"], "Reads": a["reads"],
                "Reads_Pct": f"{pct:.2f}",
            })
    log.info("Allele frequency table saved to %s", path_af)

    # ═══════════════════════════════════════════════════════════════
    # 2. per_target_editing.tsv
    # ═══════════════════════════════════════════════════════════════
    def _mutation_overlaps_target(m: Dict, target_start: int, target_end: int) -> bool:
        """Check whether a mutation overlaps a target region on the full-length reference.

        Deletion/indel types use span overlap ([pos, pos+len-1] vs [ts, te]).
        Substitution/insertion types use point-in-region check.
        """
        mp = m.get("ref_pos", -1)
        if mp < 0:
            return False
        mt = m.get("type", "")
        if mt in ("deletion", "indel"):
            me = mp + m.get("length", 1) - 1
            return mp <= target_end and me >= target_start
        else:
            return target_start <= mp <= target_end

    def _is_intrasite(m: Dict, intervals: list) -> bool:
        """Determine if a del/indel mutation is entirely within a single 27bp target region.

        A mutation that spans two or more targets is classified as inter-site.
        Insertions and substitutions are always intra-site (single-base events).
        """
        mt = m.get("type", "")
        if mt in ("insertion", "substitution"):
            return True
        mp = m.get("ref_pos", -1)
        me = mp + m.get("length", 1) - 1
        for _, ts, te in intervals:
            if mp >= ts and me <= te:
                return True
        return False

    # Pre-build all 27bp target intervals for intra/inter classification
    all_target_intervals = []
    if cutsites:
        for cs in cutsites:
            if cs.name.startswith("Target"):
                all_target_intervals.append((cs.name, cs.start - target_region_left, cs.end + target_region_right))

    path_pt = os.path.join(output_dir, "per_target_editing.tsv")
    fieldnames_pt = ["Target", "Total", "Edited", "Rate_Pct",
                     "Del_intra", "Del_inter", "Ins",
                     "Indel_intra", "Indel_inter", "Sub",
                     "Avg_Mut_Length", "Wt"]
    per_target_data = {}
    with open(path_pt, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames_pt, delimiter='\t')
        writer.writeheader()
        if cutsites:
            for cs in cutsites:
                if not cs.name.startswith("Target"):
                    continue
                # Full target region: conserved + cutsite + linker
                target_start = cs.start - target_region_left
                target_end = cs.end + target_region_right
                covering = 0
                edited = 0
                del_intra_c = del_inter_c = 0
                ins_c = 0
                indel_intra_c = indel_inter_c = 0
                sub_c = 0
                total_len_weighted = 0.0
                total_len_weight = 0
                for r in successful:
                    ar = r.get("aligned_ref", "")
                    if not ar:
                        continue
                    ref_bases = sum(1 for c in ar if c != '-')
                    if ref_bases <= cs.start:
                        continue
                    rc = r.get("readCount", 1)
                    covering += rc
                    has_mut = False
                    for m in r.get("mutations", []):
                        if not _mutation_overlaps_target(m, target_start, target_end):
                            continue
                        has_mut = True
                        mt = m.get("type", "")
                        ml = m.get("length", 1)
                        intra = _is_intrasite(m, all_target_intervals)
                        if mt == "deletion":
                            if intra:
                                del_intra_c += rc
                            else:
                                del_inter_c += rc
                            total_len_weighted += ml * rc
                            total_len_weight += rc
                        elif mt == "insertion":
                            ins_c += rc
                            total_len_weighted += ml * rc
                            total_len_weight += rc
                        elif mt == "indel":
                            if intra:
                                indel_intra_c += rc
                            else:
                                indel_inter_c += rc
                            total_len_weighted += ml * rc
                            total_len_weight += rc
                        elif mt == "substitution":
                            sub_c += rc
                    if has_mut:
                        edited += rc
                rate = f"{edited/covering*100:.1f}" if covering else "NA"
                avg_len = f"{total_len_weighted/total_len_weight:.1f}" if total_len_weight > 0 else "NA"
                wt = covering - edited
                writer.writerow({
                    "Target": cs.name, "Total": covering,
                    "Edited": edited, "Rate_Pct": rate,
                    "Del_intra": del_intra_c, "Del_inter": del_inter_c,
                    "Ins": ins_c,
                    "Indel_intra": indel_intra_c, "Indel_inter": indel_inter_c,
                    "Sub": sub_c,
                    "Avg_Mut_Length": avg_len, "Wt": wt,
                })
                per_target_data[cs.name] = {"total": covering, "mutated": edited, "rate": rate}
        else:
            writer.writerow({"Target": "No cutsite data available",
                            "Total": "", "Edited": "", "Rate_Pct": "",
                            "Del_intra": "", "Del_inter": "",
                            "Ins": "", "Indel_intra": "", "Indel_inter": "",
                            "Sub": "", "Avg_Mut_Length": "", "Wt": ""})
    log.info("Per-target editing table saved to %s", path_pt)

    # ═══════════════════════════════════════════════════════════════
    # 3. filter_reason.tsv
    # ═══════════════════════════════════════════════════════════════
    reason_groups: Dict[str, dict] = {}
    for r in failed:
        err = r.get("error", "unknown")
        if err not in reason_groups:
            reason_groups[err] = {"reason": err, "sequences": 0, "reads": 0}
        reason_groups[err]["sequences"] += 1
        reason_groups[err]["reads"] += r.get("readCount", 1)

    sorted_reasons = sorted(reason_groups.values(), key=lambda x: -x["sequences"])
    path_fr = os.path.join(output_dir, "filter_reason.tsv")
    with open(path_fr, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Reason", "Sequences", "Reads"],
                                delimiter='\t')
        writer.writeheader()
        for item in sorted_reasons:
            writer.writerow({"Reason": item["reason"],
                             "Sequences": item["sequences"],
                             "Reads": item["reads"]})
        if not sorted_reasons:
            writer.writerow({"Reason": "none", "Sequences": 0, "Reads": 0})
    log.info("Filter reason table saved to %s", path_fr)

    # ═══════════════════════════════════════════════════════════════
    # 4a. deletion_length.tsv
    # ═══════════════════════════════════════════════════════════════
    del_len_counter: Dict[int, int] = Counter()
    del_len_reads: Dict[int, int] = Counter()
    ins_len_counter: Dict[int, int] = Counter()
    ins_len_reads: Dict[int, int] = Counter()

    for r in successful:
        rc = r.get("readCount", 1)
        for block in r.get("stats", {}).get("gap_blocks_query", []):
            del_len_counter[block] += 1
            del_len_reads[block] += rc
        for block in r.get("stats", {}).get("gap_blocks_ref", []):
            ins_len_counter[block] += 1
            ins_len_reads[block] += rc

    total_del_reads = sum(del_len_reads.values())
    total_ins_reads = sum(ins_len_reads.values())

    path_dl = os.path.join(output_dir, "deletion_length.tsv")
    with open(path_dl, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Length_bp", "Events", "Reads", "Reads_Pct"],
                                delimiter='\t')
        writer.writeheader()
        for length in sorted(del_len_counter):
            pct = del_len_reads[length] / total_del_reads * 100 if total_del_reads > 0 else 0
            writer.writerow({
                "Length_bp": length, "Events": del_len_counter[length],
                "Reads": del_len_reads[length], "Reads_Pct": f"{pct:.2f}",
            })
        if not del_len_counter:
            writer.writerow({"Length_bp": "", "Events": 0, "Reads": 0, "Reads_Pct": ""})
    log.info("Deletion length distribution table saved to %s", path_dl)

    # ═══════════════════════════════════════════════════════════════
    # 4b. insertion_length.tsv
    # ═══════════════════════════════════════════════════════════════
    path_il = os.path.join(output_dir, "insertion_length.tsv")
    with open(path_il, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Length_bp", "Events", "Reads", "Reads_Pct"],
                                delimiter='\t')
        writer.writeheader()
        for length in sorted(ins_len_counter):
            pct = ins_len_reads[length] / total_ins_reads * 100 if total_ins_reads > 0 else 0
            writer.writerow({
                "Length_bp": length, "Events": ins_len_counter[length],
                "Reads": ins_len_reads[length], "Reads_Pct": f"{pct:.2f}",
            })
        if not ins_len_counter:
            writer.writerow({"Length_bp": "", "Events": 0, "Reads": 0, "Reads_Pct": ""})
    log.info("Insertion length distribution table saved to %s", path_il)

    # ═══════════════════════════════════════════════════════════════
    # 5. event_level_details.tsv
    # ═══════════════════════════════════════════════════════════════
    _save_event_level_table(successful, output_dir, cutsites,
                            target_region_left=target_region_left,
                            target_region_right=target_region_right)

    # ═══════════════════════════════════════════════════════════════
    # 6. Compute aggregated stats for HTML report
    # ═══════════════════════════════════════════════════════════════
    mutated_seqs = 0
    mutated_reads = 0
    for r in successful:
        rc = r.get("readCount", 1)
        stats = r["stats"]
        if stats["mismatches"] > 0 or stats["gaps_in_query"] > 0 or stats["gaps_in_ref"] > 0:
            mutated_seqs += 1
            mutated_reads += rc
    unmutated_seqs = len(successful) - mutated_seqs
    unmutated_reads = total_reads_success - mutated_reads
    total_reads_all = sum(r["readCount"] for r in results)

    # ═══════════════════════════════════════════════════════════════
    # 7. read_to_allele.tsv (when original_read_names are available)
    # ═══════════════════════════════════════════════════════════════
    if read_to_allele_path and any("original_read_names" in r for r in successful):
        parent = os.path.dirname(read_to_allele_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        with open(read_to_allele_path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(["read_name", "allele_label"])
            for r in successful:
                label = _allele_label(r.get("mutations", []))
                for read_name in r.get("original_read_names", []):
                    writer.writerow([read_name, label])
        log.info("Read-to-allele mapping saved to %s", read_to_allele_path)

    return {
        "total_sequences": len(results),
        "total_reads_all": total_reads_all,
        "successful_alignments": len(successful),
        "failed_alignments": len(failed),
        "total_reads_successful": total_reads_success,
        "mutated_sequences": mutated_seqs,
        "unmutated_sequences": unmutated_seqs,
        "mutated_reads": mutated_reads,
        "unmutated_reads": unmutated_reads,
        "editing_efficiency_pct": round(mutated_seqs / len(successful) * 100, 2) if successful else 0.0,
        "editing_efficiency_reads_pct": round(mutated_reads / total_reads_success * 100, 2) if total_reads_success else 0.0,
        "del_length_reads": dict(del_len_reads),
        "ins_length_reads": dict(ins_len_reads),
        "per_target": per_target_data,
    }


def _save_event_level_table(successful: List[Dict], output_dir: str,
                            cutsites: Optional[list] = None,
                            target_region_left: int = 13,
                            target_region_right: int = 7) -> None:
    """Write event-level statistics: one row per unique mutation event.

    Each mutation event is uniquely identified by (type, start_pos, length).
    For events spanning multiple targets, records the start/end target names
    and the number of targets crossed.

    Args:
        successful: List of successful alignment result dicts.
        output_dir: Output directory path.
        cutsites: CutsiteRegion list for target interval computation.
        target_region_left: Bases left of cutsite start for target interval.
        target_region_right: Bases right of cutsite end for target interval.
    """
    from collections import defaultdict

    import re
    # Build target intervals: cutsite ± flanking regions
    target_intervals = []
    if cutsites:
        for cs in cutsites:
            if cs.name.startswith("Target"):
                target_intervals.append((cs.name, cs.start - target_region_left, cs.end + target_region_right))

    # Helper for numeric comparison of target names (e.g. "Target10" → 10)
    _target_num = lambda name: int(re.search(r'\d+', name).group()) if re.search(r'\d+', name) else 0

    # Aggregate events by signature: (type, start_pos, length)
    event_groups = defaultdict(lambda: {"sequences": 0, "reads": 0,
                                        "affected": set(), "min_target": None, "max_target": None})

    for r in successful:
        rc = r.get("readCount", 1)
        for m in r.get("mutations", []):
            typ = m.get("type", "")
            rp = m.get("ref_pos", -1)
            rl = m.get("length", 1)
            key = (typ, rp, rl)
            eg = event_groups[key]
            eg["sequences"] += 1
            eg["reads"] += rc

            # Event span (end position)
            if typ in ("deletion", "indel"):
                event_end = rp + rl - 1
            else:
                event_end = rp

            # Check which targets this event overlaps
            for tname, ts, te in target_intervals:
                if rp <= te and event_end >= ts:
                    eg["affected"].add(tname)
                    # Use numeric comparison for target name ordering
                    tnum = _target_num(tname)
                    if eg["min_target"] is None or tnum < _target_num(eg["min_target"]):
                        eg["min_target"] = tname
                    if eg["max_target"] is None or tnum > _target_num(eg["max_target"]):
                        eg["max_target"] = tname

    if not event_groups:
        return

    path_ev = os.path.join(output_dir, "event_level_details.tsv")
    fieldnames = ["Type", "Start_Pos", "End_Pos", "Length",
                  "Affected_Targets", "N_Targets",
                  "Target_Range", "Sequences", "Reads"]
    with open(path_ev, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        # Sort by reads descending, then position ascending
        sorted_events = sorted(event_groups.items(),
                               key=lambda x: (-x[1]["reads"], x[0][1]))
        for key, eg in sorted_events:
            typ, rp, rl = key
            if typ in ("deletion", "indel"):
                end_pos = rp + rl - 1
            else:
                end_pos = rp
            affected_str = ";".join(sorted(eg["affected"])) if eg["affected"] else ""
            n_targets = len(eg["affected"])
            target_range = f"{eg['min_target']}-{eg['max_target']}" if eg["min_target"] and eg["max_target"] and eg["min_target"] != eg["max_target"] else (affected_str if n_targets == 1 else "")
            writer.writerow({
                "Type": typ, "Start_Pos": rp, "End_Pos": end_pos,
                "Length": rl,
                "Affected_Targets": affected_str, "N_Targets": n_targets,
                "Target_Range": target_range,
                "Sequences": eg["sequences"], "Reads": eg["reads"],
            })

    log.info("Event-level detail table saved to %s", path_ev)
