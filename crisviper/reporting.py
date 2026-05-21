"""crisviper/reporting.py — 报告生成和结果保存模块"""

import json
import csv
import os
import sys
from typing import List, Dict, Optional
from collections import Counter

from crisviper.logging_config import get_logger
from crisviper.plotting import generate_charts

log = get_logger(__name__)


def save_alignment_results(results: List[Dict], output_path: str, fmt: str = "json") -> None:
    """
    保存比对结果

    参数:
        results: 比对结果列表
        output_path: 输出文件路径
        fmt: 输出格式 (json, tsv, all)
    """
    if fmt == "json":
        path = _ensure_extension(output_path, ".json")
        with open(path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        log.info("JSON结果已保存至 %s", path)

    elif fmt == "tsv":
        path = _ensure_extension(output_path, ".tsv")
        _save_tsv_results(results, path)

    elif fmt == "all":
        # 同时输出JSON和TSV
        json_path = _ensure_extension(output_path, ".json")
        tsv_path = _ensure_extension(output_path, ".tsv")

        with open(json_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        log.info("JSON结果已保存至 %s", json_path)

        _save_tsv_results(results, tsv_path)

    else:
        log.error("不支持的输出格式: %s", fmt)
        sys.exit(1)


def _ensure_extension(path: str, ext: str) -> str:
    """确保文件路径有指定的扩展名"""
    if not path.endswith(ext):
        return path + ext
    return path


def _save_tsv_results(results: List[Dict], output_path: str) -> None:
    """将比对结果保存为TSV格式，包含对齐后的序列"""
    simplified = []
    for result in results:
        if "error" in result:
            simplified.append({
                "readName": result["readName"],
                "cellBC": result["cellBC"],
                "UMI": result["UMI"],
                "readCount": result["readCount"],
                "score": "NA",
                "matches": "NA",
                "mismatches": "NA",
                "gaps_in_query": "NA",
                "similarity": "NA",
                "aligned_ref": "NA",
                "aligned_query": "NA",
                "error": result["error"]
            })
        else:
            simplified.append({
                "readName": result["readName"],
                "cellBC": result["cellBC"],
                "UMI": result["UMI"],
                "readCount": result["readCount"],
                "score": result["score"],
                "matches": result["stats"]["matches"],
                "mismatches": result["stats"]["mismatches"],
                "gaps_in_query": result["stats"]["gaps_in_query"],
                "similarity": result["stats"]["similarity"],
                "aligned_ref": result.get("aligned_ref", ""),
                "aligned_query": result.get("aligned_query", ""),
                "error": ""
            })

    fieldnames = ["readName", "cellBC", "UMI", "readCount", "score",
                  "matches", "mismatches", "gaps_in_query", "similarity",
                  "aligned_ref", "aligned_query", "error"]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(simplified)

    log.info("TSV结果已保存至 %s", output_path)


# ═══════════════════════════════════════════════════════════════
# 分析报告生成
# ═══════════════════════════════════════════════════════════════

def generate_report(results: List[Dict], output_path: str, fmt: str = "json",
                     ref_length: int = 0, ref_seq: str = "",
                     cutsites: list = None,
                     allele_window_start: int = 0,
                     allele_window_end: int = None,
                     allele_top_n: int = 50,
                     version: str = "2.1.0") -> None:
    """
    生成突变分析报告

    参数:
        results: 比对结果列表
        output_path: 输出文件路径
        format: 报告格式 (json, html)
        ref_length: 参考序列长度（用于图表）
        ref_seq: 参考序列（用于allele热图）
        cutsites: cutsite区域列表（用于allele热图标注）
        allele_window_start: Allele热图显示起始位置
        allele_window_end: Allele热图显示结束位置（含）
        allele_top_n: Allele热图展示的top allele数（默认50）
        version: 工具版本号
    """
    total_sequences = len(results)
    successful = [r for r in results if "error" not in r]
    failed = [r for r in results if "error" in r]

    total_successful = len(successful)
    total_failed = len(failed)
    total_reads = sum(r.get("readCount", 1) for r in successful)
    total_reads_all = sum(r.get("readCount", 1) for r in results)

    # 统计突变
    mutated_seqs = []
    unmutated_seqs = []
    for r in successful:
        stats = r["stats"]
        has_mismatch = stats["mismatches"] > 0
        has_deletion = stats["gaps_in_query"] > 0
        has_insertion = stats["gaps_in_ref"] > 0
        if has_mismatch or has_deletion or has_insertion:
            mutated_seqs.append(r)
        else:
            unmutated_seqs.append(r)

    total_mutated = len(mutated_seqs)
    total_unmutated = len(unmutated_seqs)
    mutated_reads = sum(r.get("readCount", 1) for r in mutated_seqs)
    efficiency = total_mutated / total_successful * 100 if total_successful > 0 else 0.0

    # 突变类型统计（序列数 + Reads数）
    only_insertion = {"sequences": 0, "reads": 0}
    only_deletion = {"sequences": 0, "reads": 0}
    only_substitution = {"sequences": 0, "reads": 0}
    insertion_and_deletion = {"sequences": 0, "reads": 0}
    insertion_and_substitution = {"sequences": 0, "reads": 0}
    deletion_and_substitution = {"sequences": 0, "reads": 0}
    all_three = {"sequences": 0, "reads": 0}

    ins_lengths = []     # 插入长度列表（每条序列一次）
    del_lengths = []     # 删除长度列表（每条序列一次）
    ins_length_reads = Counter()   # 插入长度 → Reads数
    del_length_reads = Counter()   # 删除长度 → Reads数
    total_mismatches = 0

    for r in mutated_seqs:
        rc = r.get("readCount", 1)
        stats = r["stats"]
        has_ins = stats["gaps_in_ref"] > 0
        has_del = stats["gaps_in_query"] > 0
        has_sub = stats["mismatches"] > 0

        if has_ins and not has_del and not has_sub:
            only_insertion["sequences"] += 1
            only_insertion["reads"] += rc
        elif has_del and not has_ins and not has_sub:
            only_deletion["sequences"] += 1
            only_deletion["reads"] += rc
        elif has_sub and not has_ins and not has_del:
            only_substitution["sequences"] += 1
            only_substitution["reads"] += rc
        elif has_ins and has_del and not has_sub:
            insertion_and_deletion["sequences"] += 1
            insertion_and_deletion["reads"] += rc
        elif has_ins and has_sub and not has_del:
            insertion_and_substitution["sequences"] += 1
            insertion_and_substitution["reads"] += rc
        elif has_del and has_sub and not has_ins:
            deletion_and_substitution["sequences"] += 1
            deletion_and_substitution["reads"] += rc
        elif has_ins and has_del and has_sub:
            all_three["sequences"] += 1
            all_three["reads"] += rc

        # 收集插入长度（序列级 + Reads加权）
        for block in stats.get("gap_blocks_ref", []):
            ins_lengths.append(block)
            ins_length_reads[block] += rc

        # 收集删除长度（序列级 + Reads加权）
        for block in stats.get("gap_blocks_query", []):
            del_lengths.append(block)
            del_length_reads[block] += rc

        total_mismatches += stats["mismatches"]

    avg_insertion_len = sum(ins_lengths) / len(ins_lengths) if ins_lengths else 0.0
    avg_deletion_len = sum(del_lengths) / len(del_lengths) if del_lengths else 0.0
    max_insertion_len = max(ins_lengths) if ins_lengths else 0
    max_deletion_len = max(del_lengths) if del_lengths else 0

    # Reads加权的平均/最大长度
    total_ins_reads = sum(ins_length_reads.values())
    total_del_reads = sum(del_length_reads.values())
    avg_ins_len_reads = (sum(k * v for k, v in ins_length_reads.items()) / total_ins_reads
                         if total_ins_reads else 0.0)
    avg_del_len_reads = (sum(k * v for k, v in del_length_reads.items()) / total_del_reads
                         if total_del_reads else 0.0)

    # 构建报告字典
    report = {
        "tool": "CARLIN序列分析工具",
        "version": version,
        "summary": {
            "total_sequences": total_sequences,
            "total_reads_all": total_reads_all,
            "successful_alignments": total_successful,
            "failed_alignments": total_failed,
            "total_reads_successful": total_reads,
            "mutated_sequences": total_mutated,
            "unmutated_sequences": total_unmutated,
            "mutated_reads": mutated_reads,
            "editing_efficiency_pct": round(efficiency, 2)
        },
        "mutation_types": {
            "only_insertion": only_insertion,
            "only_deletion": only_deletion,
            "only_substitution": only_substitution,
            "insertion_and_deletion": insertion_and_deletion,
            "insertion_and_substitution": insertion_and_substitution,
            "deletion_and_substitution": deletion_and_substitution,
            "insertion_deletion_substitution": all_three
        },
        "mutation_stats": {
            "total_point_mutations": total_mismatches,
            "total_insertion_events": len(ins_lengths),
            "total_deletion_events": len(del_lengths),
            "avg_insertion_length": round(avg_insertion_len, 2),
            "avg_deletion_length": round(avg_deletion_len, 2),
            "max_insertion_length": max_insertion_len,
            "max_deletion_length": max_deletion_len,
            "all_insertion_lengths": ins_lengths,
            "all_deletion_lengths": del_lengths,
            "all_insertion_lengths_reads": dict(ins_length_reads),
            "all_deletion_lengths_reads": dict(del_length_reads),
        },
        "mutated_sequences_detail": [
            {
                "readName": r["readName"],
                "readCount": r.get("readCount", 1),
                "mismatches": r["stats"]["mismatches"],
                "gaps_in_ref": r["stats"]["gaps_in_ref"],
                "gaps_in_query": r["stats"]["gaps_in_query"],
                "similarity": r["stats"]["similarity"],
                "mutations": r.get("mutations", []),
            }
            for r in sorted(mutated_seqs, key=lambda x: -x.get("readCount", 1))
        ]
    }

    # Per-target mutation stats
    report["per_target"] = {}
    if cutsites:
        for cs in cutsites:
            if not cs.name.startswith("Target"):
                continue
            covering = 0
            mutated = 0
            del_c = ins_c = sub_c = 0
            for r in successful:
                ar = r.get("aligned_ref", "")
                if not ar:
                    continue
                ref_bases = sum(1 for c in ar if c != '-')
                if ref_bases <= cs.start:
                    continue
                covering += 1
                has_mut = False
                for m in r.get("mutations", []):
                    mp = m.get("ref_pos", -1)
                    if cs.start - 3 <= mp <= cs.end + 3:
                        has_mut = True
                        mt = m.get("type", "")
                        if mt == "deletion": del_c += 1
                        elif mt == "insertion": ins_c += 1
                        elif mt == "substitution": sub_c += 1
                if has_mut:
                    mutated += 1
            rate = str(round(mutated/covering*100, 1)) if covering else "N/A"
            report["per_target"][cs.name] = {
                "total": covering, "mutated": mutated, "rate": rate,
                "del": del_c, "ins": ins_c, "sub": sub_c,
            }

    if fmt == "json":
        path = _ensure_extension(output_path, ".json")
        with open(path, 'w') as f:
            json.dump(report, f, indent=2, default=str)
        log.info("JSON分析报告已保存至 %s", path)

    elif fmt == "html":
        path = _ensure_extension(output_path, ".html")
        charts = generate_charts(results, report, ref_length,
                                  ref_seq=ref_seq, cutsites=cutsites,
                                  allele_window_start=allele_window_start,
                                  allele_window_end=allele_window_end,
                                  allele_top_n=allele_top_n)
        _save_report_html(report, path, charts,
                            ref_seq=ref_seq, cutsites=cutsites)
        log.info("HTML分析报告已保存至 %s", path)

    else:
        log.error("不支持的报告格式: %s", fmt)
        sys.exit(1)


def _save_report_html(report: dict, output_path: str, charts: dict = None,
                       ref_seq: str = "", cutsites: list = None) -> None:
    """将报告保存为自包含的HTML文件（包含丰富的可视化元素）"""
    s = report["summary"]
    mt = report["mutation_types"]
    ms = report["mutation_stats"]
    detail = report["mutated_sequences_detail"]

    # 突变类型表格行
    type_rows = ""
    type_labels = [
        ("only_substitution", "仅点突变（替换）"),
        ("only_deletion", "仅删除"),
        ("only_insertion", "仅插入"),
        ("insertion_and_deletion", "插入 + 删除"),
        ("insertion_and_substitution", "插入 + 点突变"),
        ("deletion_and_substitution", "删除 + 点突变"),
        ("insertion_deletion_substitution", "插入 + 删除 + 点突变"),
    ]
    for key, label in type_labels:
        val = mt.get(key, {"sequences": 0, "reads": 0})
        if isinstance(val, dict):
            seq_count = val.get("sequences", 0)
            reads_count = val.get("reads", 0)
        else:
            seq_count = reads_count = val
        seq_pct = seq_count / s["mutated_sequences"] * 100 if s["mutated_sequences"] > 0 else 0
        reads_pct = reads_count / s["mutated_reads"] * 100 if s["mutated_reads"] > 0 else 0
        type_rows += (
            f"<tr>"
            f"<td>{label}</td>"
            f"<td>{seq_count}</td><td>{seq_pct:.1f}%</td>"
            f"<td>{reads_count}</td><td>{reads_pct:.1f}%</td>"
            f"</tr>\n"
        )

    # 详细序列表格
    MUT_CLASSES = {"deletion": "mut-del", "insertion": "mut-ins", "substitution": "mut-sub"}

    def _mut_summary(muts):
        parts = []
        for m in muts:
            typ = m.get("type", "")
            pos = m.get("ref_pos", -1) + 1
            if typ == "deletion":
                label = f"Δ{pos}-{m.get('length',1)}bp"
                cls = "mut-tag mut-tag-del"
            elif typ == "insertion":
                label = f"ins@{pos}+{m.get('query_base','')}"
                cls = "mut-tag mut-tag-ins"
            elif typ == "substitution":
                label = f"{m.get('ref_base','?')}{pos}{m.get('query_base','?')}"
                cls = "mut-tag mut-tag-sub"
            else:
                label = str(typ)
                cls = "mut-tag"
            parts.append(f'<span class="{cls}">{label}</span>')
        return " ".join(parts)

    detail_rows = ""
    for i, d in enumerate(detail[:100]):
        muts = d.get("mutations", [])
        mut_summary = _mut_summary(muts)
        mut_html = f'<span class="mut-summary">{mut_summary}</span>' if mut_summary else '<span style="color:#999;">-</span>'
        detail_rows += (
            f"<tr>"
            f"<td>{i+1}</td>"
            f"<td>{d['readName']}</td>"
            f"<td>{d['readCount']}</td>"
            f"<td>{d['mismatches']}</td>"
            f"<td>{d['gaps_in_ref']}</td>"
            f"<td>{d['gaps_in_query']}</td>"
            f"<td>{d['similarity']}</td>"
            f"<td>{mut_html}</td>"
            f"</tr>\n"
        )
    if len(detail) > 100:
        detail_rows += f"<tr><td colspan='8' style='text-align:center;color:#888;'>... 还有 {len(detail)-100} 条突变序列未显示</td></tr>\n"

    # 图表嵌入
    def _img_html(key):
        b64 = (charts or {}).get(key, '')
        if not b64:
            return '<p style="color:#999;font-style:italic;">图表不可用</p>'
        klass = 'chart-img' if key != 'allele_heatmap' else 'chart-img allele-img'
        return f'<div class="chart-scroll"><img class="{klass}" src="data:image/png;base64,{b64}" onclick="openModal(this)" /></div>'

    has_charts = bool(charts)
    chart_map = {k: _img_html(k) for k in [
        'indel_length', 'allele_heatmap',
    ]}


    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>CARLIN Analysis Report</title>
<style>
:root{{--pri:#2563eb;--sec:#059669;--red:#dc2626;--pur:#7c3aed;--amb:#f59e0b;
  --bg:#f8fafc;--card:#fff;--txt:#1e293b;--txt2:#64748b;--bdr:#e2e8f0;}}
*{{margin:0;padding:0;box-sizing:border-box;}}
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;
  background:var(--bg);color:var(--txt);line-height:1.6;}}
.sidebar{{position:fixed;top:0;left:0;width:200px;height:100vh;
  background:linear-gradient(180deg,#1e293b,#0f172a);color:#fff;overflow-y:auto;z-index:100;}}
.sidebar h1{{font-size:14px;padding:16px 12px 8px;border-bottom:1px solid rgba(255,255,255,.1);}}
.sidebar a{{display:block;padding:8px 12px;color:#cbd5e1;text-decoration:none;font-size:12px;
  border-left:3px solid transparent;}}
.sidebar a:hover{{background:rgba(255,255,255,.06);color:#fff;border-left-color:var(--pri);}}
.main{{margin-left:200px;padding:20px 24px;max-width:1400px;}}
h2{{font-size:18px;font-weight:700;margin:28px 0 12px;padding-bottom:4px;border-bottom:2px solid var(--pri);}}
.cards{{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:8px;margin:12px 0;}}
.card{{background:var(--card);border-radius:6px;padding:12px;
  box-shadow:0 1px 3px rgba(0,0,0,.06);border:1px solid var(--bdr);}}
.card .lbl{{font-size:9px;text-transform:uppercase;letter-spacing:.5px;color:var(--txt2);}}
.card .val{{font-size:20px;font-weight:700;margin-top:2px;}}
.crd-hl .val{{color:var(--pri);}}
table{{width:100%;border-collapse:collapse;background:var(--card);border-radius:6px;overflow:hidden;
  box-shadow:0 1px 3px rgba(0,0,0,.06);margin:12px 0;font-size:13px;}}
th{{background:var(--pri);color:#fff;padding:8px 10px;text-align:left;font-weight:600;}}
td{{padding:6px 10px;border-bottom:1px solid #f1f5f9;}}
tr:hover td{{background:#f8fafc;}}
.mut-summary{{font-size:11px;color:var(--txt);white-space:normal;word-break:break-all;}}
.mut-tag{{display:inline-block;padding:0 4px;border-radius:3px;font-size:10px;font-weight:600;margin:1px;}}
.mut-tag-del{{background:#fee2e2;color:#dc2626;}}
.mut-tag-ins{{background:#dbeafe;color:#2563eb;}}
.mut-tag-sub{{background:#fef3c7;color:#d97706;}}
.section{{background:var(--card);border-radius:6px;padding:16px;
  box-shadow:0 1px 3px rgba(0,0,0,.06);margin:12px 0;}}
.stat-grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(180px,1fr));gap:10px;}}
.stat-item{{padding:6px 10px;background:#f8f9fc;border-radius:4px;}}
.stat-item .sl{{font-size:11px;color:var(--txt2);}}
.stat-item .sv{{font-size:16px;font-weight:600;}}

/* Ruler */

/* Chart images */
.chart-box{{background:var(--card);border-radius:6px;padding:12px;
  box-shadow:0 1px 3px rgba(0,0,0,.06);margin:12px 0;text-align:center;}}
.chart-box .chart-img{{max-width:100%;height:auto;border-radius:4px;cursor:zoom-in;}}
.chart-row{{display:flex;flex-wrap:wrap;gap:12px;margin:12px 0;}}
.chart-row .chart-box{{flex:1;min-width:280px;}}
.chart-scroll{{overflow-x:auto;max-width:100%;}}
.chart-scroll .allele-img{{max-width:none;}}

footer{{margin-top:30px;padding-top:10px;border-top:1px solid var(--bdr);
  color:#aaa;font-size:11px;text-align:center;}}
@media(max-width:768px){{.sidebar{{display:none;}}.main{{margin-left:0;padding:10px;}}}}

/* Modal */
.modal-overlay{{display:none;position:fixed;top:0;left:0;width:100%;height:100%;
  background:rgba(0,0,0,0.85);z-index:9999;cursor:zoom-out;}}
.modal-overlay.active{{display:flex;align-items:center;justify-content:center;}}
.modal-content{{position:relative;max-width:95vw;max-height:95vh;overflow:auto;}}
.modal-content img{{display:block;max-width:none;border-radius:4px;box-shadow:0 4px 30px rgba(0,0,0,0.3);}}
</style>
</head>
<body>
<nav class="sidebar">
  <h1>CrisViper</h1>
  <a href="#summary">Summary</a>
  <a href="#indel-length">Indel Length</a>
  <a href="#target-breakdown">Targets</a>
  <a href="#heatmap">Allele Heatmap</a>
  <a href="#mutation-types">Mutation Types</a>
  <a href="#stats">Statistics</a>
  <a href="#details">Details</a>
</nav>
<div class="main">

<h2 id="summary">Summary</h2>
<div class="cards">
  <div class="card"><div class="lbl">Total Sequences</div><div class="val">{s["total_sequences"]}</div></div>
  <div class="card"><div class="lbl">Total Reads</div><div class="val">{s["total_reads_all"]:,}</div></div>
  <div class="card"><div class="lbl">Aligned</div><div class="val">{s["successful_alignments"]}</div></div>
  <div class="card crd-hl"><div class="lbl">Editing Efficiency</div><div class="val">{s["editing_efficiency_pct"]}%</div></div>
  <div class="card"><div class="lbl">Mutated Sequences</div><div class="val">{s["mutated_sequences"]}</div></div>
  <div class="card"><div class="lbl">Mutated Reads</div><div class="val">{s["mutated_reads"]:,}</div></div>
  <div class="card"><div class="lbl">Unmutated</div><div class="val">{s["unmutated_sequences"]}</div></div>
  <div class="card"><div class="lbl">Failed</div><div class="val">{s["failed_alignments"]}</div></div>
</div>

<h2 id="indel-length">Indel Length Distribution</h2>
{has_charts and f'''<div class="chart-box">{chart_map["indel_length"]}</div>''' or ''}

<h2 id="mutation-types">Mutation Type Distribution</h2>
<table>
<tr><th>Mutation Type</th><th>Sequences</th><th>Seq %</th><th>Reads</th><th>Reads %</th></tr>
{type_rows}
</table>

<h2 id="stats">Mutation Statistics</h2>
<div class="section"><div class="stat-grid">
  <div class="stat-item"><div class="sl">Point Mutations</div><div class="sv">{ms["total_point_mutations"]}</div></div>
  <div class="stat-item"><div class="sl">Insertion Events</div><div class="sv">{ms["total_insertion_events"]}</div></div>
  <div class="stat-item"><div class="sl">Deletion Events</div><div class="sv">{ms["total_deletion_events"]}</div></div>
  <div class="stat-item"><div class="sl">Avg Insertion Length</div><div class="sv">{ms["avg_insertion_length"]} bp</div></div>
  <div class="stat-item"><div class="sl">Avg Deletion Length</div><div class="sv">{ms["avg_deletion_length"]} bp</div></div>
  <div class="stat-item"><div class="sl">Max Insertion</div><div class="sv">{ms["max_insertion_length"]} bp</div></div>
  <div class="stat-item"><div class="sl">Max Deletion</div><div class="sv">{ms["max_deletion_length"]} bp</div></div>
</div></div>

<h2 id="target-breakdown">Target-specific Editing Efficiency</h2>
<div class="section"><div class="stat-grid">
  {report.get("per_target") and "".join(
    f'<div class=stat-item>'
    f'<div class=sl>{csname}</div>'
    f'<div class=sv>{d["rate"]}%</div>'
    f'<div style=font-size:10px;color:#94a3b8;>mutated {d["mutated"]}/{d["total"]} sequences</div></div>\n'
    for csname, d in report.get("per_target", {}).items()
  ) or '<p style=color:#94a3b8;>No cutsite data available.</p>'}
</div></div>

{has_charts and f'''<h2 id="heatmap">Allele Heatmap</h2>
<div class="chart-box">{chart_map["allele_heatmap"]}</div>''' or ''}

<h2 id="details">Mutated Sequence Details (Top 100)</h2>
<table id="detailTable">
<thead><tr><th onclick="sortTable(0)">#</th><th onclick="sortTable(1)">Name</th><th onclick="sortTable(2)">Reads</th><th onclick="sortTable(3)">Mismatches</th><th onclick="sortTable(4)">Insertions</th><th onclick="sortTable(5)">Deletions</th><th onclick="sortTable(6)">Similarity</th><th>Mutations</th></tr></thead>
<tbody>
{detail_rows}
</tbody>
</table>

<footer>Generated by CrisViper v{report["version"]}</footer>
</div>

<div id="modalOverlay" class="modal-overlay" onclick="document.getElementById('modalOverlay').classList.remove('active')">
  <div class="modal-content"><img id="modalImg" src="" /></div>
</div>
<script>
function openModal(img){{
  document.getElementById('modalImg').src=img.src;
  document.getElementById('modalOverlay').classList.add('active');
}}
function sortTable(col){{
  var t=document.getElementById('detailTable'),b=t.querySelector('tbody'),r=Array.from(b.querySelectorAll('tr'));
  var lr=r[r.length-1]; if(lr&&lr.cells.length===1) r.pop();
  var asc=t.getAttribute('data-sort')==col?t.getAttribute('data-asc')!=='true':true;
  r.sort(function(a,b){{
    var va=a.cells[col].textContent.trim(),vb=b.cells[col].textContent.trim();
    var na=parseFloat(va),nb=parseFloat(vb);
    return !isNaN(na)&&!isNaN(nb)?(asc?na-nb:nb-na):(asc?va.localeCompare(vb):vb.localeCompare(va));
  }});
  r.forEach(function(x){{b.appendChild(x);}});
  if(lr&&lr.cells.length===1) b.appendChild(lr);
  t.setAttribute('data-sort',col); t.setAttribute('data-asc',asc);
  var ths=t.querySelectorAll('th');ths.forEach(function(th,i){{if(i<7) th.textContent=th.textContent.replace(/[\u25B2\u25BC]/g,'');}});
  if(col<7) ths[col].textContent+=asc?' \u25B2':' \u25BC';
}}
</script>
</body>
</html>"""
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)


def save_text_report(
    results: List[Dict],
    output_path: str,
    ref_seq: str = "",
    version: str = "2.1.0",
    called_alleles: Optional[List] = None,
    cutsites: Optional[List] = None,
) -> str:
    """Generate MATLAB-compatible text report files.

    Produces in output_path/:
      - Results.txt: Summary statistics, read breakdown, mutation stats.
      - AlleleAnnotations.txt: HGVS annotations for each mutated allele.
      - Warnings.txt: Quality warnings (off-target amp, filtering, etc.).

    MATLAB equivalent: reports/generate_text_output.m + reports/generate_warnings.m

    Args:
        results: List of alignment result dicts.
        output_path: Output directory path.
        ref_seq: Reference sequence (for template matching).
        version: Tool version string.
        called_alleles: Optional list of called allele dicts.

    Returns:
        The output directory path.
    """
    from crisviper.mutation import extract_mutations, annotate_mutations

    os.makedirs(output_path, exist_ok=True)

    total_seqs = len(results)
    successful = [r for r in results if "error" not in r]
    failed = [r for r in results if "error" in r]
    total_successful = len(successful)
    total_failed = len(failed)
    total_reads = sum(r.get("readCount", 1) for r in successful)
    total_reads_all = sum(r.get("readCount", 1) for r in results)

    # Mutation classification
    mut_types = Counter()
    mut_types_reads = Counter()
    for r in successful:
        s = r["stats"]
        rc = r.get("readCount", 1)
        key = (
            ("I" if s["gaps_in_ref"] > 0 else "") +
            ("D" if s["gaps_in_query"] > 0 else "") +
            ("S" if s["mismatches"] > 0 else "")
        ) or "WT"
        mut_types[key] += 1
        mut_types_reads[key] += rc

    mutated = [r for r in successful if r["stats"]["mismatches"] > 0
               or r["stats"]["gaps_in_query"] > 0 or r["stats"]["gaps_in_ref"] > 0]
    edited_seqs = len(mutated)
    edited_reads = sum(r.get("readCount", 1) for r in mutated)

    # ── AlleleAnnotations.txt ──
    _write_allele_annotations(output_path, mutated, successful, ref_seq, cutsites)

    # ── Results.txt ──
    _write_results_txt(output_path, total_seqs, total_successful, total_failed,
                       total_reads, total_reads_all, edited_seqs, edited_reads,
                       mut_types, mut_types_reads, version, called_alleles)

    # ── Warnings.txt ──
    _write_warnings_txt(output_path, total_seqs, total_successful, total_failed,
                        total_reads, total_reads_all, edited_seqs)

    return output_path


def _write_allele_annotations(
    output_path: str,
    mutated: List[Dict],
    all_successful: List[Dict],
    ref_seq: str = "",
    cutsites: Optional[List] = None,
) -> None:
    """Write AlleleAnnotations.txt — HGVS annotations for each read."""
    from crisviper.mutation import extract_mutations, annotate_mutations

    path = os.path.join(output_path, "AlleleAnnotations.txt")
    with open(path, 'w') as f:
        for r in all_successful:
            aligned_ref = r.get("aligned_ref", "")
            aligned_query = r.get("aligned_query", "")
            if aligned_ref and aligned_query and ref_seq:
                events = extract_mutations(aligned_ref, aligned_query, cutsites=cutsites)
                ann = annotate_mutations(events, full=False)
                f.write(ann + "\n")
            else:
                f.write("[]\n")
    log.info("AlleleAnnotations saved: %s", path)


def _write_results_txt(
    output_path: str,
    total_seqs: int,
    total_successful: int,
    total_failed: int,
    total_reads: int,
    total_reads_all: int,
    edited_seqs: int,
    edited_reads: int,
    mut_types: Dict,
    mut_types_reads: Dict,
    version: str,
    called_alleles: Optional[List] = None,
) -> None:
    """Write Results.txt — MATLAB-compatible results summary."""
    path = os.path.join(output_path, "Results.txt")
    with open(path, 'w') as f:
        f.write("RESULTS\n\n")
        f.write(f"{'Tool Version:':<45} {version}\n")
        f.write(f"{'Total Sequences:':<45} {total_seqs}\n")
        f.write(f"{'Successful Alignments:':<45} {total_successful}\n")
        f.write(f"{'Failed Alignments:':<45} {total_failed}\n")
        f.write(f"{'Total Reads (all):':<45} {total_reads_all:,}\n")
        f.write(f"{'Total Reads (aligned):':<45} {total_reads:,}\n")

        efficiency = edited_seqs / total_successful * 100 if total_successful > 0 else 0.0
        f.write(f"{'Editing Efficiency:':<45} {efficiency:.1f}%\n")

        f.write(f"\n{'READ BREAKDOWN':<45} {'Reads':>10}\n\n")
        f.write(f"{'  in_fastq:':<45} {total_reads_all:>10}\n")
        f.write(f"{'  valid_alignment:':<45} {total_reads:>10}\n")

        pct_aligned = total_reads / total_reads_all * 100 if total_reads_all > 0 else 0
        f.write(f"{'  aligned_pct:':<45} {pct_aligned:>9.0f}%\n")
        f.write(f"{'  mutated:':<45} {edited_reads:>10}\n")
        f.write(f"{'  unmutated:':<45} {total_reads - edited_reads:>10}\n")

        f.write(f"\n{'MUTATION BREAKDOWN':<45} {'Seqs':>10} {'Reads':>10}\n\n")
        labels = {
            "I": "Insertion only", "D": "Deletion only", "S": "Substitution only",
            "ID": "Ins+Del", "IS": "Ins+Sub", "DS": "Del+Sub",
            "IDS": "All three", "WT": "Unmutated",
        }
        for key, label in labels.items():
            n = mut_types.get(key, 0)
            nr = mut_types_reads.get(key, 0)
            if n > 0 or nr > 0:
                f.write(f"  {label:<42} {n:>10} {nr:>10}\n")

        f.write(f"\n{'MUTATION STATS':<45}\n\n")
        ins_events = sum(v for k, v in mut_types.items() if "I" in k)
        del_events = sum(v for k, v in mut_types.items() if "D" in k)
        sub_events = mut_types.get("S", 0) + mut_types.get("IS", 0) + mut_types.get("DS", 0) + mut_types.get("IDS", 0)
        f.write(f"{'  Edited sequences:':<45} {edited_seqs}\n")
        f.write(f"{'  Edited reads:':<45} {edited_reads:,}\n")
        f.write(f"{'  Events with insertion:':<45} {ins_events}\n")
        f.write(f"{'  Events with deletion:':<45} {del_events}\n")

        if called_alleles:
            f.write(f"\n{'ALLELES':<45}\n\n")
            f.write(f"{'  Total alleles:':<45} {len(called_alleles)}\n")
            singletons = sum(1 for a in called_alleles if a.get("freq", 0) == 1)
            f.write(f"{'  Singletons:':<45} {singletons}\n")

    log.info("Results.txt saved: %s", path)


def _write_warnings_txt(
    output_path: str,
    total_seqs: int,
    total_successful: int,
    total_failed: int,
    total_reads: int,
    total_reads_all: int,
    edited_seqs: int = 0,
) -> None:
    """Write Warnings.txt — quality warnings for the run."""
    path = os.path.join(output_path, "Warnings.txt")
    with open(path, 'w') as f:
        f.write("WARNINGS\n\n")

        # Off-target amplification
        ota_pct = (1 - total_reads / total_reads_all) * 100 if total_reads_all > 0 else 0
        f.write("OFF-TARGET AMPLIFICATION\n")
        if ota_pct > 10:
            f.write(f"\nSignificant off-target amplification detected: "
                    f"{ota_pct:.0f}% of reads are not CARLIN.\n")
        else:
            f.write(f"\nInsignificant off-target amplification detected. "
                    f"Only {ota_pct:.0f}% of reads are off-target.\n")

        # Filtering warnings
        f.write("\nFILTERING\n")
        issues = False
        if total_successful / total_seqs * 100 < 75 if total_seqs > 0 else False:
            f.write(f"\nOnly {total_successful / total_seqs * 100:.0f}% of sequences "
                    f"aligned successfully. Consider checking primer anchoring.\n")
            issues = True
        if not issues:
            f.write("\nNo issues detected at filtering step.\n")

        # Analysis warnings
        f.write("\nANALYSIS\n")
        analysis_issues = False
        if edited_seqs < 500:
            f.write(f"\nNumber of edited sequences is low ({edited_seqs}).\n")
            analysis_issues = True
        if not analysis_issues:
            f.write("\nNo issues detected during analysis.\n")

        # Results warnings
        f.write("\nRESULTS\n")
        results_issues = False
        pct_edited = edited_seqs / total_successful * 100 if total_successful > 0 else 0
        if pct_edited < 30:
            f.write(f"\nLow editing efficiency detected. "
                    f"Only {pct_edited:.0f}% of sequences reported an edited allele.\n")
            results_issues = True
        if not results_issues:
            f.write("\nNo issues detected in results.\n")

    log.info("Warnings.txt saved: %s", path)
