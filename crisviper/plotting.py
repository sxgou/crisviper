"""Charting and visualization functions for alignment results."""

import io
import base64
import sys
from typing import List, Dict, Tuple, Optional
import numpy as np
from crisviper.logging_config import get_logger
from crisviper.mutation import _build_ref_pos_map_full

log = get_logger(__name__)

# Lazy-load matplotlib (slow import — defer until actual charting)
_HAS_MPL = False
plt = None
ticker = None

def _ensure_mpl():
    """Ensure matplotlib is loaded; returns True if available."""
    global _HAS_MPL, plt, ticker
    if _HAS_MPL:
        return True
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as _plt
        import matplotlib.ticker as _ticker
        plt = _plt
        ticker = _ticker
        _HAS_MPL = True
        return True
    except ImportError:
        log.warning("未安装matplotlib，图表功能将不可用。请安装: pip install matplotlib")
        return False


def _img_to_b64(fig) -> str:
    """将matplotlib图形转换为base64 HTML可嵌入字符串"""
    if fig is None:
        return ""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=120, bbox_inches='tight', facecolor='white')
    buf.seek(0)
    img_b64 = base64.b64encode(buf.read()).decode('utf-8')
    if plt is not None:
        plt.close(fig)
    return img_b64


def _gen_reads_distribution(results: List[Dict]) -> str:
    """Reads数分布直方图（仅统计成功比对的结果，对数bins + 截断尾部离群值）"""
    read_counts = [r.get("readCount", 1) for r in results if "error" not in r]
    if not read_counts:
        return ""
    fig, ax = plt.subplots(figsize=(6, 3))
    # log-spaced bins覆盖主要数据范围（covers 99.95%的数据）
    log_min, log_max = 0.5, 1000
    bins = np.logspace(np.log10(log_min), np.log10(log_max), 40)
    ax.hist(read_counts, bins=bins, color='#4a6cf7', edgecolor='white', alpha=0.8)
    ax.set_xscale('log')
    ax.set_xlim(log_min, log_max)
    ax.set_xlabel('Reads per Unique Sequence')
    ax.set_ylabel('Sequence Count')
    ax.set_title('Reads Count Distribution')
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_alignment_rate(results: List[Dict]) -> str:
    """比对成功率饼图"""
    successful = sum(1 for r in results if "error" not in r)
    failed = len(results) - successful
    if failed == 0 and successful == 0:
        return ""
    values = [successful, failed] if failed > 0 else [successful]
    labels = ['Successful', 'Failed'] if failed > 0 else ['Successful']
    colors = ['#4a6cf7', '#ff6b6b'] if failed > 0 else ['#4a6cf7']
    fig, ax = plt.subplots(figsize=(4, 3))
    wedges, texts, autotexts = ax.pie(
        values, labels=labels, autopct='%1.1f%%',
        colors=colors, startangle=90, textprops={'fontsize': 10}
    )
    for t in autotexts:
        t.set_fontweight('bold')
    ax.set_title('Alignment Rate')
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_length_distribution(results: List[Dict], ref_length: int = 0,
                              weight_by_reads: bool = False) -> str:
    """查询序列长度分布直方图（以reads或sequence计数）

    参数:
        weight_by_reads: True时以reads数量加权统计
    """
    lengths = []
    weights = []
    for r in results:
        if "error" in r:
            continue
        rc = r.get("readCount", 1)
        aq = r.get("aligned_query", "")
        if aq:
            l = len(aq) - aq.count('-')
            lengths.append(l)
            if weight_by_reads:
                weights.append(rc)
        else:
            stats = r.get("stats")
            if not stats:
                continue
            alen = stats.get("alignment_length", 0)
            if alen:
                l = int(alen)
                lengths.append(l)
                if weight_by_reads:
                    weights.append(rc)
    if not lengths:
        return ""

    n_total = sum(weights) if weight_by_reads else len(lengths)
    fig, ax = plt.subplots(figsize=(6, 3))
    bins = min(30, max(lengths) - min(lengths) + 1) if max(lengths) > min(lengths) else 1
    hist_kw = dict(bins=bins, color='#20c997', edgecolor='white', alpha=0.8)
    if weight_by_reads:
        hist_kw['weights'] = weights
    ax.hist(lengths, **hist_kw)
    if ref_length > 0:
        ax.axvline(ref_length, color='#dc3545', linestyle='--', linewidth=1.5, label=f'Ref ({ref_length}bp)')
        ax.legend(fontsize=9)
    ylbl = 'Reads' if weight_by_reads else 'Sequences'
    ax.set_xlabel('Sequence Length (bp)')
    ax.set_ylabel(f'{ylbl} (count)')
    ax.set_title(f'Sequence Length Distribution ({ylbl.lower()})')
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_mutation_overview(report_data: dict) -> Tuple[str, str]:
    """突变总览：序列角度和Reads角度的饼图"""
    s = report_data["summary"]
    mutated_seq = s.get("mutated_sequences", 0)
    unmutated_seq = s.get("unmutated_sequences", 0)
    mutated_reads = s.get("mutated_reads", 0)
    total_reads = s.get("total_reads_successful", 0)
    unmutated_reads = total_reads - mutated_reads

    def make_pie(values, labels, title):
        # 过滤掉零值
        nz = [(v, l) for v, l in zip(values, labels) if v > 0]
        if len(nz) <= 1:
            return ""
        nz_values = [v for v, l in nz]
        nz_labels = [l for v, l in nz]
        fig, ax = plt.subplots(figsize=(4, 3))
        colors_list = ['#ff6b6b', '#adb5bd', '#20c997', '#4a6cf7']
        wedges, texts, autotexts = ax.pie(
            nz_values, labels=nz_labels, autopct='%1.1f%%',
            colors=colors_list[:len(nz_values)], startangle=90,
            textprops={'fontsize': 10}
        )
        for t in autotexts:
            t.set_fontweight('bold')
        ax.set_title(title)
        fig.tight_layout()
        return _img_to_b64(fig)

    seq_pie = make_pie(
        [mutated_seq, unmutated_seq],
        [f'Mutated ({mutated_seq})', f'Unmutated ({unmutated_seq})'],
        'Mutation Rate (by Sequence)'
    )
    reads_pie = make_pie(
        [mutated_reads, unmutated_reads],
        [f'Mutated ({mutated_reads})', f'Unmutated ({unmutated_reads})'],
        'Mutation Rate (by Reads)'
    )
    return seq_pie, reads_pie


def _gen_mutation_type_chart(report_data: dict) -> str:
    """突变类型分布柱状图（序列数和Reads数并排）"""
    mt = report_data["mutation_types"]
    labels = ['Only\nSubstitution', 'Only\nDeletion', 'Only\nInsertion',
              'Ins+Del', 'Ins+Sub', 'Del+Sub', 'All Three']
    keys = ['only_substitution', 'only_deletion', 'only_insertion',
            'insertion_and_deletion', 'insertion_and_substitution',
            'deletion_and_substitution', 'insertion_deletion_substitution']

    seq_counts = []; reads_counts = []
    for k in keys:
        v = mt.get(k, {"sequences": 0, "reads": 0})
        if isinstance(v, dict):
            seq_counts.append(v.get("sequences", 0))
            reads_counts.append(v.get("reads", 0))
        else:
            seq_counts.append(v); reads_counts.append(v)

    if all(s == 0 for s in seq_counts):
        return ""

    fig, ax = plt.subplots(figsize=(7, 3.5))
    x = range(len(labels))
    w = 0.35
    ax.bar([i - w/2 for i in x], seq_counts, w, label='Sequences', color='#4a6cf7', alpha=0.85)
    ax.bar([i + w/2 for i in x], reads_counts, w, label='Reads', color='#ff6b6b', alpha=0.85)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('Count')
    ax.set_title('Mutation Type Distribution')
    ax.legend(fontsize=9)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_indel_length_chart(report_data: dict) -> str:
    """蝴蝶图：Insertion向上，Deletion向下的diverging bar chart (Reads加权, 百分比Y轴)"""
    ms = report_data.get("mutation_stats", {})
    ins_counts = ms.get("all_insertion_lengths_reads", {})
    del_counts = ms.get("all_deletion_lengths_reads", {})
    if not ins_counts and not del_counts:
        from collections import Counter
        ins_raw = ms.get("all_insertion_lengths", [])
        del_raw = ms.get("all_deletion_lengths", [])
        if not ins_raw and not del_raw:
            return ""
        ins_counts = Counter(ins_raw)
        del_counts = Counter(del_raw)
    else:
        ins_counts = {int(k): v for k, v in ins_counts.items()}
        del_counts = {int(k): v for k, v in del_counts.items()}

    all_keys = list(ins_counts.keys()) + list(del_counts.keys())
    max_len = min(max(all_keys) if all_keys else 1, 50)
    all_lengths = list(range(1, max_len + 1))
    ins_raw_vals = [ins_counts.get(l, 0) for l in all_lengths]
    del_raw_vals = [del_counts.get(l, 0) for l in all_lengths]

    total_ins = sum(ins_counts.values())
    total_del = sum(del_counts.values())
    # 转为百分比
    ins_pct = [v / total_ins * 100 if total_ins else 0 for v in ins_raw_vals]
    del_pct = [-v / total_del * 100 if total_del else 0 for v in del_raw_vals]

    max_pct = max(max(ins_pct) if ins_pct else 0, max(abs(x) for x in del_pct) if del_pct else 0)
    y_lim = max(10, ((max_pct // 10) + 1) * 10)  # 以10为单位向上取整

    fig, ax = plt.subplots(figsize=(8, 3.5))
    ax.bar(all_lengths, ins_pct, color='#20c997', width=0.8,
           label=f'Insertion ({total_ins} reads)')
    ax.bar(all_lengths, del_pct, color='#ff6b6b', width=0.8,
           label=f'Deletion ({total_del} reads)')
    ax.axhline(0, color='#333', linewidth=0.5)
    ax.set_xlabel('Indel Length (bp)')
    ax.set_ylabel('Reads (%)')
    ax.set_title('Indel Length Distribution (Butterfly Chart)')
    ax.set_xticks(all_lengths)
    ax.set_ylim(-y_lim, y_lim)
    ax.set_yticks(range(-int(y_lim), int(y_lim) + 1, 10))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{int(abs(x))}%'))
    ax.legend(fontsize=9, loc='upper right')
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_indel_position_chart(results: List[Dict], ref_length: int) -> str:
    """基于参考序列位置的Indel分布柱状图"""
    if not results or ref_length == 0:
        return ""

    ins_positions = []; del_positions = []
    for r in results:
        if "error" in r:
            continue
        ar = r.get("aligned_ref", ""); aq = r.get("aligned_query", "")
        if not ar or not aq:
            continue
        ref_pos = 0
        for a, b in zip(ar, aq):
            if a == '-' and b != '-':
                ins_positions.append(ref_pos)
            elif a != '-' and b == '-':
                del_positions.append(ref_pos)
            if a != '-':
                ref_pos += 1

    if not ins_positions and not del_positions:
        return ""

    fig, ax = plt.subplots(figsize=(10, 3.5))
    bins = min(100, ref_length)
    if del_positions:
        ax.hist(del_positions, bins=bins, alpha=0.6, color='#ff6b6b',
                label=f'Deletion ({len(del_positions)})', edgecolor=None)
    if ins_positions:
        ax.hist(ins_positions, bins=bins, alpha=0.6, color='#20c997',
                label=f'Insertion ({len(ins_positions)})', edgecolor=None)
    ax.set_xlabel('Reference Position (bp)'); ax.set_ylabel('Indel Count')
    ax.set_title('Position-based Indel Distribution')
    ax.legend(fontsize=9)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_mutation_composition_chart(report_data: dict) -> str:
    """点突变/插入/删除总体比例图"""
    ms = report_data.get("mutation_stats", {})
    values = [ms.get("total_point_mutations", 0),
              ms.get("total_insertion_events", 0),
              ms.get("total_deletion_events", 0)]
    labels = ['Substitution', 'Insertion', 'Deletion']
    colors = ['#ffd43b', '#20c997', '#ff6b6b']
    # 过滤掉零值
    non_zero = [(v, l, c) for v, l, c in zip(values, labels, colors) if v > 0]
    if len(non_zero) <= 1:
        return ""
    fig, ax = plt.subplots(figsize=(5, 3))
    nz_values = [v for v, l, c in non_zero]
    nz_labels = [l for v, l, c in non_zero]
    nz_colors = [c for v, l, c in non_zero]
    wedges, texts, autotexts = ax.pie(
        nz_values, labels=nz_labels, autopct='%1.1f%%',
        colors=nz_colors, startangle=90, textprops={'fontsize': 10}
    )
    for t in autotexts:
        t.set_fontweight('bold')
    centre = plt.Circle((0, 0), 0.5, fc='white')
    ax.add_artist(centre)
    ax.set_title('Mutation Composition')
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_allele_heatmap(results: List[Dict], ref_seq: str,
                         cutsites: list = None, top_n: int = 50,
                         window_start: int = 0,
                         window_end: int = None,
                         primer5_len: int = 23,
                         primer3_len: int = 33) -> str:
    """
    Allele热图：每行一个allele，每列一个碱基位置。

    展示top_n个最频繁allele的碱基级比对，支持自定义显示区域。
    window_start/window_end = None 时显示全长参考序列。

    参数:
        results: 比对结果列表
        ref_seq: 参考序列
        cutsites: 用于标注cutsite区域
        top_n: 展示的top allele数
        window_start: 显示起始位置（0-indexed，默认0）
        window_end: 显示结束位置（含，默认ref_seq全长）

    返回: base64 PNG
    """
    if not results or not ref_seq:
        return ""

    ref_len = len(ref_seq)
    if window_end is None:
        window_end = ref_len - 1
    window_end = min(window_end, ref_len - 1)
    window_start = max(0, window_start)
    window_len = window_end - window_start + 1
    if window_len <= 0:
        return ""

    # 1. 聚合allele（按内部区域去重，消除引物差异导致的假分裂）
    allele_counts = {}
    allele_data = {}
    for r in results:
        if "error" in r:
            continue
        aq = r.get("aligned_query", "")
        ar = r.get("aligned_ref", "")
        if not aq or not ar:
            continue
        rc = r.get("readCount", 1)

        # 提取内部区域（去掉引物区）
        pos_map, total_ref = _build_ref_pos_map_full(ar)
        internal_start = next((i for i, p in enumerate(pos_map) if p == primer5_len), primer5_len)
        internal_end = next((i for i, p in enumerate(pos_map) if p == total_ref - primer3_len),
                            len(ar) - primer3_len)
        internal_end = max(internal_start, min(internal_end, len(ar)))
        internal_key = aq[internal_start:internal_end]

        allele_counts[internal_key] = allele_counts.get(internal_key, 0) + rc
        if internal_key not in allele_data:
            # 保存第一份完整对齐序列用于展示
            allele_data[internal_key] = {
                "aligned_ref": ar,
                "aligned_query": aq,
                "internal_start": internal_start,
                "internal_end": internal_end,
                "readCount": rc
            }
    if not allele_counts:
        return ""

    # 2. 排序取top_n
    sorted_alleles = sorted(allele_counts.items(), key=lambda x: -x[1])
    top_alleles = sorted_alleles[:top_n]
    total_reads = sum(allele_counts.values())

    # 3. Insertion感知的显示列构建
    window_ref_bases = ref_seq[window_start:window_end + 1]
    n_cols = window_len  # 参考序列列数

    def _build_display_with_insertions(aligned_ref, aligned_query):
        """构建包含insertion的显示列（ref不开gap，query顺移展示插入碱基）。

        遍历比对列，参考碱基占一列，插入碱基（ref为gap但query有碱基时）
        额外占一列，后续碱基顺移。返回：
            seq:      每列对应的query碱基列表（长度 >= window_len）
            ref:      每列对应的ref碱基（''表示insertion列）
            ins_blocks: [(start, end), ...] 插入块的列范围
            width:    总显示列数
        """
        seq = []
        ref = []
        ins_blocks = []
        ref_pos = 0
        col = 0
        in_ins = False
        ins_start = 0
        for rb, qb in zip(aligned_ref, aligned_query):
            if rb != '-':   # ref有碱基 → 正常列
                if in_ins:
                    ins_blocks.append((ins_start, col - 1))
                    in_ins = False
                if window_start <= ref_pos <= window_end:
                    seq.append(qb)
                    ref.append(rb)
                    col += 1
                ref_pos += 1
                if ref_pos > window_end:
                    break
            elif qb != '-':  # ref为gap、query有碱基 → insertion列
                if window_start <= ref_pos <= window_end:
                    if not in_ins:
                        ins_start = col
                        in_ins = True
                    seq.append(qb)
                    ref.append('')
                    col += 1
        if in_ins:
            ins_blocks.append((ins_start, col - 1))
        # 填充至至少window_len个ref列
        ref_cnt = sum(1 for r in ref if r != '')
        while ref_cnt < window_len:
            seq.append('-')
            ref.append('')
            col += 1
            ref_cnt += 1
        return seq, ref, ins_blocks, col

    # 预计算所有allele行的显示列，截断至n_cols
    row_data = []  # [(seq, ref, ins_blocks, readCount), ...]
    for internal_key, total_rc in top_alleles:
        ad = allele_data.get(internal_key, {})
        ar = ad.get("aligned_ref", "")
        aq = ad.get("aligned_query", "")
        if ar and aq:
            sq, rf, blk, w = _build_display_with_insertions(ar, aq)
        else:
            sq = list(window_ref_bases)
            rf = list(window_ref_bases)
            blk = []
            w = len(sq)
        # 截断至n_cols，超出reference长度的部分丢弃
        if len(sq) > n_cols:
            # 修正insertion块坐标：丢弃超出n_cols的部分
            clipped_blk = []
            for start, end in blk:
                if start >= n_cols:
                    continue
                clipped_blk.append((start, min(end, n_cols - 1)))
            sq = sq[:n_cols]
            rf = rf[:n_cols]
            blk = clipped_blk
        # 填充至n_cols（短于参考序列的补gap）
        while len(sq) < n_cols:
            sq.append('-')
            rf.append('')
        row_data.append((sq, rf, blk, total_rc))

    grid_cols = n_cols  # 固定为参考序列列数

    # 4. 碱基配色
    base_colors = {
        'A': '#cbe9cb', 'C': '#fee5ce', 'G': '#ffffd6',
        'T': '#e5deed', 'U': '#e5deed', 'N': '#e9e9e9',
        'a': '#cbe9cb', 'c': '#fee5ce', 'g': '#ffffd6',
        't': '#e5deed', 'u': '#e5deed', 'n': '#e9e9e9',
        '-': '#f0f0f0',
    }

    n_rows = len(top_alleles)

    # 5. 自适应单元格尺寸（基于原始参考列数）
    if n_cols > 200:
        cell_size = 10
        font_size = 5.5
    elif n_cols > 100:
        cell_size = 14
        font_size = 6.5
    elif n_cols > 60:
        cell_size = 18
        font_size = 7
    else:
        cell_size = 22
        font_size = 9

    label_w = 180  # 右侧标注宽度(px)
    fig_w_px = grid_cols * cell_size + label_w
    fig_h_px = (n_rows + 5.0) * cell_size

    dpi = 100
    fig_w = max(8, fig_w_px / dpi)
    fig_h = max(3, fig_h_px / dpi)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, grid_cols)
    ax.set_ylim(-2.5, n_rows + 2.5)
    ax.set_aspect('equal')
    ax.axis('off')

    # 6. 参考序列行（grid_cols列，超出n_cols部分留白）
    ref_y = n_rows + 1.5
    for ci in range(grid_cols):
        if ci < n_cols:
            base = window_ref_bases[ci] if ci < len(window_ref_bases) else ' '
            color = base_colors.get(base, '#e9e9e9')
            ax.add_patch(plt.Rectangle((ci, ref_y - 0.5), 1, 1,
                                        facecolor=color, edgecolor='white', linewidth=0.3))
            ax.text(ci + 0.5, ref_y, base, ha='center', va='center',
                    fontsize=font_size, fontweight='normal')
        else:
            ax.add_patch(plt.Rectangle((ci, ref_y - 0.5), 1, 1,
                                        facecolor='white', edgecolor='white', linewidth=0))

    # 标注cutsite区域（灰色半透明方框，基于参考坐标n_cols）
    if cutsites:
        for cs in cutsites:
            if cs.start > window_end or cs.end < window_start:
                continue
            cs_c0 = max(cs.start - window_start, 0)
            cs_c1 = min(cs.end - window_start, n_cols - 1)
            if cs_c0 <= cs_c1:
                ax.add_patch(plt.Rectangle(
                    (cs_c0, ref_y - 0.5), cs_c1 - cs_c0 + 1, 1,
                    facecolor='#888888', alpha=0.3, edgecolor=None, zorder=0
                ))

    # 8. sgRNA区域示意（半高灰色方块，紧贴reference下方）
    if cutsites:
        ref_pos_to_col = {rp: rp - window_start for rp in range(window_start, window_end + 1)}
        grna_y = ref_y - 0.75
        for cs in cutsites:
            if not cs.name.startswith("Target"):
                continue
            tgt_start = cs.start - 13
            tgt_end = cs.end
            if tgt_start > window_end or tgt_end < window_start:
                continue
            col0 = ref_pos_to_col.get(max(tgt_start, window_start))
            col1 = ref_pos_to_col.get(min(tgt_end, window_end))
            if col0 is None or col1 is None:
                continue
            ax.add_patch(plt.Rectangle(
                (col0, grna_y - 0.3), col1 - col0 + 1, 0.6,
                facecolor='#bbbbbb', edgecolor='#999999', linewidth=0.3, zorder=0
            ))
            ax.text(col0 - 0.1, grna_y, cs.name.replace('Target', 'T'),
                    ha='right', va='center', fontsize=max(6, font_size - 0.5),
                    color='#666', fontweight='bold')

    # 7. Allele行（使用扩展显示列，insertion用红色框标出）
    ins_box_color = '#e94560'
    for ai, (sq, rf, blk, total_rc) in enumerate(row_data):
        y = n_rows - 1 - ai

        # 隔行背景色
        if ai % 2 == 1:
            ax.add_patch(plt.Rectangle((0, y - 0.5), grid_cols, 1,
                                        facecolor='#f8f8f8', edgecolor=None, zorder=0))

        for ci in range(grid_cols):
            base = sq[ci] if ci < len(sq) else '-'
            ref_base = rf[ci] if ci < len(rf) else ''
            is_del = (base == '-')
            is_ins_col = (ref_base == '')
            is_sub = (not is_del and not is_ins_col and base != ref_base)

            if is_ins_col:
                color = '#ffcccc'
            elif is_del:
                color = '#f0f0f0'
            else:
                color = base_colors.get(base, '#e9e9e9')

            ax.add_patch(plt.Rectangle((ci, y - 0.5), 1, 1,
                                        facecolor=color, edgecolor='white', linewidth=0.2))
            if is_ins_col:
                ax.text(ci + 0.5, y, base, ha='center', va='center',
                        fontsize=font_size, color='red', fontweight='bold')
            elif is_del:
                ax.text(ci + 0.5, y, '-', ha='center', va='center',
                        fontsize=font_size, color='#666', fontweight='bold')
            else:
                wgt = 'bold' if is_sub else 'normal'
                ax.text(ci + 0.5, y, base, ha='center', va='center',
                        fontsize=font_size, fontweight=wgt)

        # 红色框标记整体insertion块
        for start, end in blk:
            ax.add_patch(plt.Rectangle(
                (start, y - 0.5), end - start + 1, 1,
                facecolor='none', edgecolor=ins_box_color, linewidth=1.5, zorder=5
            ))

        # 右侧label
        pct = total_rc / total_reads * 100
        ax.text(grid_cols + 0.3, y, f"{pct:.1f}% ({total_rc})",
                ha='left', va='center', fontsize=max(8, font_size + 1))

    # 图例
    lgy = -1.2
    items = [('A', '#cbe9cb'), ('C', '#fee5ce'), ('G', '#ffffd6'),
             ('T', '#e5deed'), ('-', '#f0f0f0'), ('ins', '#ffcccc')]
    lx = 0
    leg_font = max(9, font_size + 2)
    for label, color in items:
        ax.add_patch(plt.Rectangle((lx, lgy - 0.35), 0.7, 0.7,
                                    facecolor=color, edgecolor='#666', linewidth=0.6))
        ax.text(lx + 0.9, lgy, label, ha='left', va='center', fontsize=leg_font)
        lx += 2.2

    ax.set_title(f'Top {top_n} Alleles — ref {window_start + 1}–{window_end + 1}bp ({n_cols}bp)',
                 fontsize=max(10, font_size + 3), pad=4)
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_per_target_stacked_bar(results: List[Dict], cutsites: list) -> str:
    """Per-target stacked bar: each bar split by mutation type (reads-weighted)."""
    if not results or not cutsites:
        return ""

    target_names = [cs.name for cs in cutsites if cs.name.startswith("Target")]
    n_targets = len(target_names)
    if n_targets == 0:
        return ""

    categories = ["Unedited", "Substitution", "Insertion", "Deletion", "Complex", "Mixed"]
    n_cats = len(categories)
    data = [[0] * n_cats for _ in range(n_targets)]

    successful = [r for r in results if "error" not in r]

    for ti, cs in enumerate(cutsites):
        if not cs.name.startswith("Target"):
            continue
        cs_lo, cs_hi = cs.start - 3, cs.end + 3
        for r in successful:
            ar = r.get("aligned_ref", "")
            if not ar:
                continue
            rc = r.get("readCount", 1)

            # --- Coverage check ---
            ref_bases = sum(1 for c in ar if c != '-')
            covers = ref_bases > cs.start
            # Also check if any mutation spans into the target
            for m in r.get("mutations", []):
                mp = m.get("ref_pos", -1)
                ml = m.get("length", 1)
                mt = m.get("type", "")
                if mt in ("deletion", "complex"):
                    if mp <= cs_hi and (mp + ml - 1) >= cs_lo:
                        covers = True
                        break
                elif cs_lo <= mp <= cs_hi:
                    covers = True
                    break
            if not covers:
                continue

            # --- Classify mutation types within target window ---
            has_del = has_ins = has_sub = has_complex = False
            for m in r.get("mutations", []):
                mp = m.get("ref_pos", -1)
                mt = m.get("type", "")
                if mt == "complex":
                    if mp <= cs_hi and (mp + m.get("length", 1) - 1) >= cs_lo:
                        has_complex = True
                elif cs_lo <= mp <= cs_hi:
                    if mt == "deletion":
                        has_del = True
                    elif mt == "insertion":
                        has_ins = True
                    elif mt == "substitution":
                        has_sub = True

            # Priority: Complex > Mixed > Deletion > Insertion > Substitution > Unedited
            if has_complex:
                data[ti][4] += rc
            elif has_del and has_ins:
                data[ti][5] += rc  # Mixed (both indels)
            elif has_del:
                data[ti][3] += rc
            elif has_ins:
                data[ti][2] += rc
            elif has_sub:
                data[ti][1] += rc
            else:
                data[ti][0] += rc

    n_reads_per_target = [sum(row) for row in data]
    pct_data = []
    for row, total_n in zip(data, n_reads_per_target):
        pct_data.append([v / total_n * 100 if total_n else 0 for v in row])

    fig, ax = plt.subplots(figsize=(7, 3.5))
    x = range(n_targets)
    bar_w = 0.75
    colors = ['#e6e6e6', '#ffd43b', '#4a6cf7', '#ff6b6b', '#e94560', '#b366b3']
    bottom = [0.0] * n_targets

    for ci, (color, label) in enumerate(zip(colors, categories)):
        vals = [pct_data[ti][ci] for ti in range(n_targets)]
        bars = ax.bar(x, vals, bar_w, bottom=bottom, color=color,
                      label=label, edgecolor='white', linewidth=0.3)
        for ti in range(n_targets):
            if vals[ti] > 4:
                ax.text(ti, bottom[ti] + vals[ti] / 2,
                        f'{vals[ti]:.0f}%', ha='center', va='center',
                        fontsize=6.5, color='white' if ci != 1 else '#333',
                        fontweight='bold')
        bottom = [b + v for b, v in zip(bottom, vals)]

    for ti in range(n_targets):
        ax.text(ti, -3.5, f'n={n_reads_per_target[ti]}',
                ha='center', va='top', fontsize=7, color='#888')

    ax.set_xticks(list(x))
    ax.set_xticklabels([cs.name.replace('Target', 'T') for cs in cutsites
                        if cs.name.startswith('Target')], fontsize=8)
    ax.set_ylabel('Reads (%)')
    ax.set_title('Per-Target Editing Type Distribution (weighted by reads)')
    ax.set_ylim(-8, 105)
    ax.legend(fontsize=7, loc='upper right')
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_mutation_segment_plot(results: List[Dict], ref_seq: str,
                                cutsites: list = None, top_n: int = 50) -> str:
    """Mutation segment plot: per-allele horizontal view of indels and substitutions."""
    if not results or not ref_seq:
        return ""

    ref_len = len(ref_seq)

    # 1. Aggregate alleles
    allele_data = {}
    for r in results:
        if "error" in r:
            continue
        aq = r.get("aligned_query", "")
        ar = r.get("aligned_ref", "")
        rc = r.get("readCount", 1)
        if not aq or not ar:
            continue
        if ar not in allele_data:
            allele_data[ar] = {"aq": aq, "ar": ar, "rc": 0}
        allele_data[ar]["rc"] += rc
    sorted_alleles = sorted(allele_data.values(), key=lambda x: -x["rc"])
    top = sorted_alleles[:top_n]

    # 2. Extract segments per allele
    rows = []
    for ad in top:
        ar, aq, rc = ad["ar"], ad["aq"], ad["rc"]
        del_segs, ins_segs, sub_positions = [], [], []
        ref_pos = 0
        i = 0
        while i < len(ar):
            if ar[i] != '-' and aq[i] == '-':            # deletion
                start = ref_pos
                while i < len(ar) and ar[i] != '-' and aq[i] == '-':
                    i += 1; ref_pos += 1
                del_segs.append((start, ref_pos - 1))
            elif ar[i] == '-' and aq[i] != '-':           # insertion
                ins_pos = ref_pos
                j = i
                while j < len(ar) and ar[j] == '-' and aq[j] != '-':
                    j += 1
                ins_segs.append((ins_pos, ins_pos))
                i = j
            else:
                if ar[i] != '-' and aq[i] != '-' and ar[i] != aq[i]:  # substitution
                    sub_positions.append(ref_pos)
                if ar[i] != '-':
                    ref_pos += 1
                i += 1
        rows.append((rc, del_segs, ins_segs, sub_positions))

    if not rows:
        return ""

    n_rows = len(rows)
    row_h = 0.55  # half the original height
    fig_h = max(3, n_rows * row_h + 1.5)
    fig, ax = plt.subplots(figsize=(12, fig_h))
    ax.set_xlim(-2, ref_len + 22)
    ax.set_ylim(-1.8, n_rows + 0.8)
    ax.set_xlabel('Reference Position (bp)')
    ax.set_ylabel('Top Alleles')
    ax.set_title(f'Mutation Segment Plot (Top {n_rows} Alleles by Read Count)')
    ax.yaxis.set_ticks([])

    for ri, (rc, del_segs, ins_segs, sub_positions) in enumerate(rows):
        y = n_rows - 1 - ri

        # Alternating row background
        if ri % 2 == 1:
            ax.axhspan(y - row_h / 2, y + row_h / 2, color='#f8f8f8', zorder=0)

        # Gray gRNA targets
        for cs in cutsites or []:
            if not cs.name.startswith("Target"):
                continue
            rect = plt.Rectangle((cs.start, y - row_h / 2 + 0.05),
                                 cs.end - cs.start + 1, row_h - 0.1,
                                 facecolor='#e0e0e0', alpha=0.4, edgecolor=None, zorder=1)
            ax.add_patch(rect)

        # Draw in order: insertion → deletion (on top) → substitution (on top)
        # Insertion markers (blue diamonds, lowest zorder)
        for pos, _ in ins_segs:
            ax.plot(pos, y, marker='D', color='#4a6cf7', markersize=5,
                    linewidth=0, zorder=2)

        # Deletion segments (red, drawn on top of insertions)
        for s, e in del_segs:
            ax.plot([s, e], [y, y], color='#ff6b6b', linewidth=4,
                    solid_capstyle='round', zorder=3)

        # Substitution point mutations (black asterisks, topmost)
        for pos in sub_positions:
            ax.plot(pos, y, marker='*', color='black', markersize=8,
                    linewidth=0, zorder=4)

        # Label
        ax.text(ref_len + 2, y, f'{ri+1}. ({rc})', ha='left', va='center',
                fontsize=7, color='#666')

    # Reference line
    ax.axhline(y=-0.5, color='#333', linewidth=0.5)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='#ff6b6b', linewidth=3, label='Deletion'),
        Line2D([0], [0], marker='D', color='#4a6cf7', label='Insertion',
               markersize=6, linewidth=0),
        Line2D([0], [0], marker='*', color='black', label='Substitution',
               markersize=8, linewidth=0),
        plt.Rectangle((0, 0), 1, 1, facecolor='#e0e0e0', alpha=0.4, label='gRNA Target'),
    ]
    ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(0, -0.9),
              ncol=4, fontsize=8, framealpha=0.8)

    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_cross_target_chord(results: List[Dict], cutsites: list) -> str:
    """Chord diagram: circular links showing deletions spanning multiple targets."""
    if not results or not cutsites:
        return ""

    target_names = [cs.name for cs in cutsites if cs.name.startswith("Target")]
    n_targets = len(target_names)
    if n_targets < 2:
        return ""

    from matplotlib.path import Path
    import matplotlib.patches as mpatches

    target_ranges = [(cs.start, cs.end) for cs in cutsites if cs.name.startswith("Target")]

    # 1. Collect cross-target events
    cross_counts = {}
    successful = [r for r in results if "error" not in r]
    for r in successful:
        ar = r.get("aligned_ref", "")
        rc = r.get("readCount", 1)
        if not ar:
            continue
        for m in r.get("mutations", []):
            mt = m.get("type", "")
            mp = m.get("ref_pos", -1)
            ml = m.get("length", 1)
            if mt not in ("deletion", "insertion"):
                continue
            event_end = mp + ml - 1 if mt == "deletion" else mp
            affected = []
            for ti, (ts, te) in enumerate(target_ranges):
                if max(mp, ts) <= min(event_end, te):
                    affected.append(ti)
            if len(affected) >= 2:
                for a in range(len(affected)):
                    for b in range(a + 1, len(affected)):
                        t1, t2 = sorted((affected[a], affected[b]))
                        key = (t1, t2, mt)
                        cross_counts[key] = cross_counts.get(key, 0) + rc

    if not cross_counts:
        return ""

    # 2. Cartesian layout
    angles = np.linspace(0, 2 * np.pi, n_targets, endpoint=False)
    gap = 0.02  # gap between target sectors (ratio of circumference)
    inner_r = 0.25
    outer_r = 1.0
    arc_span = (2 * np.pi / n_targets) * (1 - gap)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect('equal')
    ax.axis('off')

    # Draw target arcs
    for ti in range(n_targets):
        a0 = angles[ti] - arc_span / 2
        a1 = angles[ti] + arc_span / 2
        theta = np.linspace(a0, a1, 40)

        # Outer arc
        ax.plot(np.cos(theta) * outer_r, np.sin(theta) * outer_r,
                color='#333', linewidth=2, zorder=3)
        # Fill sector
        fill_theta = np.linspace(a0, a1, 20)
        for ft in fill_theta:
            ax.plot([np.cos(ft) * inner_r, np.cos(ft) * outer_r],
                    [np.sin(ft) * inner_r, np.sin(ft) * outer_r],
                    color='#e0e0e0', linewidth=0.8, zorder=1)

        # Target number label outside
        lr = outer_r + 0.1
        ax.text(np.cos(angles[ti]) * lr, np.sin(angles[ti]) * lr,
                str(ti + 1), ha='center', va='center',
                fontsize=9, fontweight='bold', zorder=4)

    del_count = sum(v for (_, _, t), v in cross_counts.items() if t == "deletion")
    ins_count = sum(v for (_, _, t), v in cross_counts.items() if t == "insertion")

    if cross_counts:
        max_count = max(cross_counts.values())
        min_width = 0.01
        max_width = 0.06

        merged = {}
        for (t1, t2, etype), count in cross_counts.items():
            key = (t1, t2)
            if key not in merged:
                merged[key] = {"deletion": 0, "insertion": 0}
            merged[key][etype] += count

        for (t1, t2), vals in merged.items():
            a1, a2 = angles[t1], angles[t2]
            # Skip adjacent targets
            if abs(t1 - t2) == 1 or (t1 == 0 and t2 == n_targets - 1):
                continue

            # Draw chord as bezier curves
            ctrl_factor = 0.3
            for etype, count in [("deletion", vals["deletion"]), ("insertion", vals["insertion"])]:
                if count == 0:
                    continue
                width = min_width + (max_width - min_width) * (count / max_count)
                color = '#ff6b6b' if etype == 'deletion' else '#4a6cf7'
                alpha = min(0.7, 0.15 + 0.55 * (count / max_count))

                # Two ends of the chord: at the inner edge of each sector
                r_start = inner_r
                r_end = inner_r

                p1 = np.array([np.cos(a1) * r_start, np.sin(a1) * r_start])
                p2 = np.array([np.cos(a2) * r_end, np.sin(a2) * r_end])

                # Control points: pull toward center then outward
                mid = (p1 + p2) / 2
                center_dist = np.linalg.norm(mid)
                if center_dist > 0:
                    inward = -mid / center_dist * 0.15
                else:
                    inward = np.array([0, 0])
                cp1 = p1 * ctrl_factor + mid * (1 - ctrl_factor) + inward
                cp2 = p2 * ctrl_factor + mid * (1 - ctrl_factor) + inward

                # For a ribbon, offset perpendicular to chord direction
                chord_dir = p2 - p1
                perp = np.array([-chord_dir[1], chord_dir[0]])
                perp_norm = np.linalg.norm(perp)
                if perp_norm > 0:
                    perp /= perp_norm

                offset = perp * width / 2

                # Ribbon vertices: outer curve and inner curve
                verts = [
                    p1 + offset, cp1 + offset, cp2 + offset, p2 + offset,
                    p2 - offset, cp2 - offset, cp1 - offset, p1 - offset,
                ]
                codes = [
                    Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                    Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                ]
                path = Path(verts, codes)
                patch = mpatches.PathPatch(path, facecolor=color, alpha=alpha,
                                           edgecolor='none', zorder=2)
                ax.add_patch(patch)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#ff6b6b', alpha=0.6, label=f'Deletion (n={del_count})'),
        mpatches.Patch(facecolor='#4a6cf7', alpha=0.6, label=f'Insertion (n={ins_count})'),
    ]
    ax.legend(handles=legend_elements, loc='lower left', framealpha=0.8, fontsize=9)

    ax.set_title('Cross-Target Deletion/Insertion Chord Diagram\n'
                 f'Deletions: {del_count}, Insertions: {ins_count}',
                 fontsize=11, pad=15)

    fig.tight_layout()
    return _img_to_b64(fig)


def generate_charts(results: List[Dict], report_data: dict = None,
                     ref_length: int = 0, ref_seq: str = "",
                     cutsites: list = None,
                     allele_window_start: int = 0,
                     allele_window_end: int = None,
                     allele_top_n: int = 50,
                     primer5_len: int = 23,
                     primer3_len: int = 33) -> Dict[str, str]:
    """
    生成所有图表，返回 {图表名: base64图片} 字典。

    参数:
        results: 比对结果列表
        report_data: 报告数据字典
        ref_length: 参考序列长度
        ref_seq: 参考序列字符串（用于allele热图）
        cutsites: CutsiteRegion列表（用于allele热图标注）
        allele_window_start: Allele热图显示起始位置（默认0）
        allele_window_end: Allele热图显示结束位置（含，默认全长）
        allele_top_n: Allele热图展示的top allele数（默认50）
    """
    charts = {}
    if not _ensure_mpl():
        return charts
    try:
        charts['reads_length_dist'] = _gen_length_distribution(results, ref_length, weight_by_reads=True)
        charts['align_rate'] = _gen_alignment_rate(results)
        charts['length_dist'] = _gen_length_distribution(results, ref_length)
        if report_data:
            seq_pie, reads_pie = _gen_mutation_overview(report_data)
            charts['mutation_seq_pie'] = seq_pie
            charts['mutation_reads_pie'] = reads_pie
            charts['mutation_type'] = _gen_mutation_type_chart(report_data)
            charts['indel_length'] = _gen_indel_length_chart(report_data)
            charts['indel_position'] = _gen_indel_position_chart(results, ref_length)
            charts['mutation_composition'] = _gen_mutation_composition_chart(report_data)
            # New charts from notebook integration
            if cutsites:
                charts['per_target_bar'] = _gen_per_target_stacked_bar(results, cutsites)
                charts['cross_target_chord'] = _gen_cross_target_chord(results, cutsites)
            if ref_seq:
                charts['segment_plot'] = _gen_mutation_segment_plot(
                    results, ref_seq, cutsites=cutsites, top_n=allele_top_n)
        # Allele热图（默认显示全长，窗口可调）
        if ref_seq and ref_length > 0:
            ws = allele_window_start
            we = allele_window_end if allele_window_end is not None else ref_length - 1
            charts['allele_heatmap'] = _gen_allele_heatmap(
                results, ref_seq, cutsites=cutsites, top_n=allele_top_n,
                window_start=ws, window_end=we,
                primer5_len=primer5_len, primer3_len=primer3_len,
            )
    except Exception as e:
        log.warning("图表生成失败 - %s", e)
    return charts
