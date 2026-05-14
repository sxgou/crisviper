"""Charting and visualization functions for alignment results."""

import io
import base64
import sys
from typing import List, Dict, Tuple, Optional
import numpy as np
from ltlib.logging_config import get_logger
from ltlib.mutation import _build_ref_pos_map_full

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
    """Reads数分布直方图"""
    read_counts = [r.get("readCount", 1) for r in results]
    if not read_counts:
        return ""
    fig, ax = plt.subplots(figsize=(6, 3))
    max_rc = max(read_counts)
    bins = min(50, max_rc) if max_rc > 1 else 1
    ax.hist(read_counts, bins=bins, color='#4a6cf7', edgecolor='white', alpha=0.8)
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


def _gen_length_distribution(results: List[Dict], ref_length: int = 0) -> str:
    """查询序列长度分布直方图"""
    lengths = []
    for r in results:
        if "error" in r:
            continue
        seq = r.get("seq", "")
        if seq:
            lengths.append(len(seq))
        else:
            stats = r.get("stats")
            if not stats:
                continue
            alen = stats.get("alignment_length", 0)
            if alen:
                lengths.append(int(alen))
    if not lengths:
        return ""
    fig, ax = plt.subplots(figsize=(6, 3))
    bins = min(30, max(lengths) - min(lengths) + 1) if max(lengths) > min(lengths) else 1
    ax.hist(lengths, bins=bins, color='#20c997', edgecolor='white', alpha=0.8)
    if ref_length > 0:
        ax.axvline(ref_length, color='#dc3545', linestyle='--', linewidth=1.5, label=f'Ref ({ref_length}bp)')
        ax.legend(fontsize=9)
    ax.set_xlabel('Sequence Length (bp)')
    ax.set_ylabel('Count')
    ax.set_title('Sequence Length Distribution')
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
    """Insertion和Deletion长度分布直方图并排"""
    ms = report_data.get("mutation_stats", {})
    ins_len = ms.get("avg_insertion_length", 0)
    del_len = ms.get("avg_deletion_length", 0)
    max_ins = ms.get("max_insertion_length", 0)
    max_del = ms.get("max_deletion_length", 0)

    ins_lengths = ms.get("all_insertion_lengths", [])
    del_lengths = ms.get("all_deletion_lengths", [])
    if not ins_lengths and not del_lengths:
        return ""

    def _do_hist(ax, lengths, color, title, avg, maxv):
        if lengths:
            max_len = max(lengths)
            bins = min(max_len, 50)
            ax.hist(lengths, bins=bins, color=color, edgecolor='white', alpha=0.8)
            ax.set_xlabel('Length (bp)'); ax.set_ylabel('Count')
            ax.set_title(f'{title} (avg={avg:.1f}bp, max={maxv}bp)')
            ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        else:
            ax.text(0.5, 0.5, f'No {title.lower()}', ha='center', va='center',
                    transform=ax.transAxes)
            ax.set_title(title)

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
    _do_hist(axes[0], ins_lengths, '#20c997', 'Insertion', ins_len, max_ins)
    _do_hist(axes[1], del_lengths, '#ff6b6b', 'Deletion', del_len, max_del)
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

    # 3. 从对齐序列中提取窗口碱基
    # 使用原始参考序列（无gap）确定列数，避免矫正管线引入的gap干扰坐标映射
    n_cols = window_len
    window_ref_bases = ref_seq[window_start:window_end + 1]

    def _extract_window_bases(aligned_ref, aligned_query):
        """通过对齐列遍历（跳过ref gap列）提取窗口各位置对应的碱基。"""
        result = []
        ref_pos = 0
        for col, c in enumerate(aligned_ref):
            if c != '-':
                if window_start <= ref_pos <= window_end:
                    base = aligned_query[col] if col < len(aligned_query) else '-'
                    result.append(base)
                ref_pos += 1
                if ref_pos > window_end:
                    break
        while len(result) < n_cols:
            result.append('-')
        return ''.join(result[:n_cols])

    # 4. 碱基配色
    base_colors = {
        'A': '#cbe9cb', 'C': '#fee5ce', 'G': '#ffffd6',
        'T': '#e5deed', 'U': '#e5deed', 'N': '#e9e9e9',
        'a': '#cbe9cb', 'c': '#fee5ce', 'g': '#ffffd6',
        't': '#e5deed', 'u': '#e5deed', 'n': '#e9e9e9',
        '-': '#f0f0f0',
    }

    n_rows = len(top_alleles)

    # 5. 自适应单元格尺寸
    if window_len > 200:
        cell_size = 10
        font_size = 5.5
    elif window_len > 100:
        cell_size = 14
        font_size = 6.5
    elif window_len > 60:
        cell_size = 18
        font_size = 7
    else:
        cell_size = 22
        font_size = 9

    label_w = 180  # 右侧标注宽度(px)
    fig_w_px = n_cols * cell_size + label_w
    fig_h_px = (n_rows + 5.0) * cell_size

    dpi = 100
    fig_w = max(8, fig_w_px / dpi)
    fig_h = max(3, fig_h_px / dpi)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, n_cols)
    ax.set_ylim(-2.5, n_rows + 2.5)
    ax.set_aspect('equal')
    ax.axis('off')

    # 6. 参考序列行
    ref_y = n_rows + 1.5
    for ci in range(n_cols):
        base = window_ref_bases[ci] if ci < len(window_ref_bases) else ' '
        x = ci + 0.5
        color = base_colors.get(base, '#e9e9e9')
        ax.add_patch(plt.Rectangle((ci, ref_y - 0.5), 1, 1,
                                    facecolor=color, edgecolor='white', linewidth=0.3))
        ax.text(x, ref_y, base, ha='center', va='center',
                fontsize=font_size, fontweight='normal')

    # 标注cutsite区域（灰色半透明方框）
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
        # 建立 ref_pos → display column 映射（基于无gap的参考序列坐标）
        ref_pos_to_col = {rp: rp - window_start for rp in range(window_start, window_end + 1)}

        grna_y = ref_y - 0.75  # 半高行中心
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
            # 绘制灰色方块
            ax.add_patch(plt.Rectangle(
                (col0, grna_y - 0.3), col1 - col0 + 1, 0.6,
                facecolor='#bbbbbb', edgecolor='#999999', linewidth=0.3, zorder=0
            ))
            # 在前方标注Target名称
            ax.text(col0 - 0.1, grna_y, cs.name.replace('Target', 'T'),
                    ha='right', va='center', fontsize=max(6, font_size - 0.5),
                    color='#666', fontweight='bold')

    # 7. Allele行
    for ai, (internal_key, total_rc) in enumerate(top_alleles):
        y = n_rows - 1 - ai
        ad = allele_data.get(internal_key, {})
        full_aligned_ref = ad.get("aligned_ref", "")
        full_aligned_query = ad.get("aligned_query", "")

        # 显示全长对齐序列（含引物），按窗口裁剪
        win_bases = _extract_window_bases(full_aligned_ref, full_aligned_query) if full_aligned_ref and full_aligned_query else ""
        win_ref_b = window_ref_bases

        while len(win_bases) < n_cols:
            win_bases += '-'
        while len(win_ref_b) < n_cols:
            win_ref_b += '-'
        win_bases = win_bases[:n_cols]
        win_ref_b = win_ref_b[:n_cols]

        for ci in range(n_cols):
            base = win_bases[ci]
            ref_base = win_ref_b[ci]
            is_del = (base == '-')
            is_ins = (base != '-' and ref_base == '-')
            is_sub = (not is_del and not is_ins and base != ref_base)

            if is_ins:
                color = '#ffcccc'
            elif is_del:
                color = '#f0f0f0'
            else:
                color = base_colors.get(base, '#e9e9e9')

            ax.add_patch(plt.Rectangle((ci, y - 0.5), 1, 1,
                                        facecolor=color, edgecolor='white', linewidth=0.2))

            if is_ins:
                ax.text(ci + 0.5, y, base, ha='center', va='center',
                        fontsize=font_size, color='red', fontweight='bold')
                ax.add_patch(plt.Rectangle(
                    (ci, y - 0.5), 1, 1,
                    facecolor='none', edgecolor='red', linewidth=1.2, zorder=5
                ))
            elif is_del:
                ax.text(ci + 0.5, y, '-', ha='center', va='center',
                        fontsize=font_size, color='#666', fontweight='bold')
            else:
                w = 'bold' if is_sub else 'normal'
                ax.text(ci + 0.5, y, base, ha='center', va='center',
                        fontsize=font_size, fontweight=w)

        # 右侧label
        pct = total_rc / total_reads * 100
        ax.text(n_cols + 0.3, y, f"{pct:.1f}% ({total_rc})",
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

    ax.set_title(f'Top {top_n} Alleles — ref {window_start + 1}–{window_end + 1}bp ({window_len}bp)',
                 fontsize=max(10, font_size + 3), pad=4)
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
        charts['reads_dist'] = _gen_reads_distribution(results)
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
