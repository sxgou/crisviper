"""crisviper/summary.py — 总结性统计表格生成模块

产生 4 个 TSV 表格，保存到输出目录：
  - allele_frequency.tsv    按突变指纹聚合的 Allele 频率表
  - per_target_editing.tsv  每个 Target 的编辑类型统计
  - filter_reason.tsv       序列丢弃原因统计
  - indel_length.tsv        Indel 长度分布（按实际长度分组）
"""

import csv
import os
from typing import List, Dict, Optional
from collections import Counter

from crisviper.logging_config import get_logger

log = get_logger(__name__)


def save_summary_tables(results: List[Dict], output_dir: str,
                         ref_seq: str = "", cutsites: Optional[list] = None) -> None:
    """生成4个总结性 TSV 表格到 output_dir。

    参数:
        results: 比对结果列表（AlignmentResult.to_dict() 格式）
        output_dir: 输出目录（自动创建）
        ref_seq: 参考序列（用于 cutsite 检测）
        cutsites: CutsiteRegion 列表（可选，缺失时 per_target 表为空）
    """
    os.makedirs(output_dir, exist_ok=True)
    successful = [r for r in results if "error" not in r]
    failed = [r for r in results if "error" in r]
    total_reads_success = sum(r.get("readCount", 1) for r in successful)

    # ── 辅助: 生成 Allele 标签 ──
    def _split_indel(m: Dict) -> List[str]:
        """将 INDEL 事件解析为独立的 del/ins 描述列表。"""
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
                ins_len = 0
                while i < len(ref_base) and ref_base[i] == '-' and query_base[i] != '-':
                    ins_len += 1
                    i += 1
                parts.append(f"ins:{ins_start}+{ins_len}bp")
            else:
                if rb != '-':
                    ref_pos += 1
                i += 1
        return parts if parts else [f"indel:{ref_pos}-{m.get('length', 1)}bp"]

    def _allele_label(mutations: List[Dict]) -> str:
        if not mutations:
            return "wt"
        parts = []
        for m in sorted(mutations, key=lambda x: (x.get("ref_pos", 0), x.get("type", ""), x.get("length", 0))):
            typ = m.get("type", "")
            pos = m.get("ref_pos", -1)
            if typ == "deletion":
                parts.append(f"del:{pos}-{m.get('length', 1)}bp")
            elif typ == "insertion":
                parts.append(f"ins:{pos}+{m.get('length', 1)}bp")
            elif typ == "indel":
                parts.extend(_split_indel(m))
            elif typ == "substitution":
                parts.append(f"sub:{pos}{m.get('ref_base', '?')}>{m.get('query_base', '?')}")
            else:
                parts.append(f"{typ}:{pos}")
        return ";".join(parts) if parts else "wt"

    def _mutation_type_label(mutations: List[Dict]) -> str:
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
    log.info("Allele频率表已保存至 %s", path_af)

    # ═══════════════════════════════════════════════════════════════
    # 2. per_target_editing.tsv
    # ═══════════════════════════════════════════════════════════════
    def _mutation_overlaps_target(m: Dict, target_start: int, target_end: int) -> bool:
        """判断突变是否与Target的20bp区域重叠（全长参考坐标）

        deletion/complex 使用跨度重叠检查（[pos, pos+len-1] 与 [ts, te] 相交）
        substitution/insertion 使用点位置检查（pos in [ts, te]）
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

    path_pt = os.path.join(output_dir, "per_target_editing.tsv")
    with open(path_pt, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Target", "Total", "Edited", "Rate_Pct",
                                                "Del", "Ins", "Sub", "Avg_Mut_Length"],
                                delimiter='\t')
        writer.writeheader()
        if cutsites:
            for cs in cutsites:
                if not cs.name.startswith("Target"):
                    continue
                # Full target 20bp region: conserved 13bp + cutsite 7bp
                target_start = cs.start - 13
                target_end = cs.end
                covering = 0
                edited = 0
                del_c = ins_c = sub_c = 0
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
                        if mt == "deletion":
                            del_c += rc
                            total_len_weighted += ml * rc
                            total_len_weight += rc
                        elif mt == "insertion":
                            ins_c += rc
                            total_len_weighted += ml * rc
                            total_len_weight += rc
                        elif mt == "indel":
                            # indel = 删除 + 插入复合事件，同时计入两者
                            del_c += rc
                            ins_c += rc
                            total_len_weighted += ml * rc
                            total_len_weight += rc
                        elif mt == "substitution":
                            sub_c += rc
                    if has_mut:
                        edited += rc
                rate = f"{edited/covering*100:.1f}" if covering else "NA"
                avg_len = f"{total_len_weighted/total_len_weight:.1f}" if total_len_weight > 0 else "NA"
                writer.writerow({
                    "Target": cs.name, "Total": covering,
                    "Edited": edited, "Rate_Pct": rate,
                    "Del": del_c, "Ins": ins_c, "Sub": sub_c,
                    "Avg_Mut_Length": avg_len,
                })
        else:
            writer.writerow({"Target": "No cutsite data available",
                            "Total": "", "Edited": "", "Rate_Pct": "",
                            "Del": "", "Ins": "", "Sub": "", "Avg_Mut_Length": ""})
    log.info("Per-target编辑表已保存至 %s", path_pt)

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
    log.info("过滤原因统计表已保存至 %s", path_fr)

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
    log.info("Deletion长度分布表已保存至 %s", path_dl)

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
    log.info("Insertion长度分布表已保存至 %s", path_il)

    # ═══════════════════════════════════════════════════════════════
    # 5. event_level_details.tsv
    # ═══════════════════════════════════════════════════════════════
    _save_event_level_table(successful, output_dir, cutsites)


def _save_event_level_table(successful: List[Dict], output_dir: str,
                            cutsites: Optional[list] = None) -> None:
    """事件级统计表：每个突变事件一行，记录位置、覆盖target、reads数。

    一个 mutation event 由 (type, start_pos, length) 唯一标识。
    对于横跨多个 target 的事件，记录起止 target 名称和跨越数量。
    """
    from collections import defaultdict

    # Build target intervals (20bp region matching summary logic)
    target_intervals = []
    if cutsites:
        for cs in cutsites:
            if cs.name.startswith("Target"):
                target_intervals.append((cs.name, cs.start - 13, cs.end))

    # Aggregate events by signature
    # key: (type, start_pos, length)
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
                    if eg["min_target"] is None or tname < eg["min_target"]:
                        eg["min_target"] = tname
                    if eg["max_target"] is None or tname > eg["max_target"]:
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

        # Sort: by reads desc, then position
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

    log.info("事件级统计表已保存至 %s", path_ev)
