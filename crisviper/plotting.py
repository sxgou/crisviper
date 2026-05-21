"""Charting and visualization functions for alignment results."""

import io
import base64
from typing import List, Dict, Optional
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
        if report_data:
            charts['indel_length'] = _gen_indel_length_chart(report_data)
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
