"""Charting and visualization functions for alignment results."""

import io
import base64
from typing import List, Dict
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
        log.warning("matplotlib not installed — charting disabled. Install with: pip install matplotlib")
        return False


def _img_to_b64(fig, **savefig_kw) -> str:
    """Convert a matplotlib figure to a base64-encoded PNG string for HTML embedding.

    Args:
        fig: Matplotlib figure to convert.
        **savefig_kw: Extra kwargs forwarded to fig.savefig() (e.g., bbox_inches=None).
    """
    if fig is None:
        return ""
    buf = io.BytesIO()
    kwargs = dict(format='png', dpi=120, facecolor='white')
    kwargs.update(savefig_kw)
    fig.savefig(buf, **kwargs)
    buf.seek(0)
    img_b64 = base64.b64encode(buf.read()).decode('utf-8')
    if plt is not None:
        plt.close(fig)
    return img_b64


def _gen_indel_length_chart(ins_length_reads: dict = None,
                            del_length_reads: dict = None,
                            report_data: dict = None) -> str:
    """Butterfly chart: insertion lengths upward, deletion lengths downward (read-weighted, percentage Y-axis).

    Args:
        ins_length_reads: Dict of {insertion_length: read_count} (from summary data).
        del_length_reads: Dict of {deletion_length: read_count} (from summary data).
        report_data: Fallback report dict (reads from mutation_stats if above not provided).
    """
    ins_counts = ins_length_reads or {}
    del_counts = del_length_reads or {}
    if not ins_counts and not del_counts and report_data:
        ms = report_data.get("mutation_stats", {})
        ins_counts = ms.get("all_insertion_lengths_reads", {})
        del_counts = ms.get("all_deletion_lengths_reads", {})
        if not ins_counts and not del_counts:
            from collections import Counter
            ins_raw = ms.get("all_insertion_lengths", [])
            del_raw = ms.get("all_deletion_lengths", [])
            if ins_raw:
                ins_counts = Counter(ins_raw)
            if del_raw:
                del_counts = Counter(del_raw)
    if ins_counts:
        ins_counts = {int(k): v for k, v in ins_counts.items()}
    if del_counts:
        del_counts = {int(k): v for k, v in del_counts.items()}

    # Full range from 1 to max observed length (no 50bp cap)
    all_keys = list(ins_counts.keys()) + list(del_counts.keys())
    if not all_keys:
        return ""
    max_len = max(all_keys)
    all_lengths = list(range(1, max_len + 1))
    ins_raw_vals = [ins_counts.get(l, 0) for l in all_lengths]
    del_raw_vals = [del_counts.get(l, 0) for l in all_lengths]

    total_ins = sum(ins_counts.values())
    total_del = sum(del_counts.values())
    ins_pct = [v / total_ins * 100 if total_ins else 0 for v in ins_raw_vals]
    del_pct = [-v / total_del * 100 if total_del else 0 for v in del_raw_vals]

    max_pct = max(max(ins_pct) if ins_pct else 0, max(abs(x) for x in del_pct) if del_pct else 0)
    y_lim = max(10, ((max_pct // 10) + 1) * 10)

    # Fixed aspect ratio matching the pie chart (in chart-scroll for overflow)
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.bar(all_lengths, ins_pct, color='#20c997', width=0.8,
           label=f'Insertion ({total_ins} reads)')
    ax.bar(all_lengths, del_pct, color='#ff6b6b', width=0.8,
           label=f'Deletion ({total_del} reads)')
    ax.axhline(0, color='#333', linewidth=0.5)
    ax.set_xlabel('Indel Length (bp)')
    ax.set_ylabel('Reads (%)')
    ax.set_title('Indel Length Distribution')
    if max_len <= 50:
        ax.set_xticks(all_lengths)
    elif max_len <= 100:
        ax.set_xticks(range(0, max_len + 1, 10))
    elif max_len <= 200:
        ax.set_xticks(range(0, max_len + 1, 20))
    else:
        ax.set_xticks(range(0, max_len + 1, 50))
    ax.set_xlim(-1, max_len + 0.5)
    ax.set_ylim(-y_lim, y_lim)
    ax.set_yticks(range(-int(y_lim), int(y_lim) + 1, 10))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{int(abs(x))}%'))
    ax.legend(fontsize=9, loc='upper right')
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_editing_pie_charts(report_data: dict) -> str:
    """Dual pie charts: editing efficiency by sequence count and by read count."""
    s = report_data.get("summary", {})
    mutated_seqs = s.get("mutated_sequences", 0)
    unmutated_seqs = s.get("unmutated_sequences", 0)
    mutated_reads = s.get("mutated_reads", 0)
    total_reads = s.get("total_reads_successful", 0)
    unmutated_reads = total_reads - mutated_reads

    if mutated_seqs + unmutated_seqs == 0 or total_reads == 0:
        return ""

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 3))
    colors = ['#ff6b6b', '#e2e8f0']
    explode = (0.03, 0)

    # Pie 1: by sequence count
    ax1.pie([mutated_seqs, unmutated_seqs], labels=['Mutated', 'Unmutated'],
            autopct='%1.1f%%', colors=colors, explode=explode,
            startangle=90, textprops={'fontsize': 8})
    ax1.set_title(f'By Sequences\n({mutated_seqs:,} / {mutated_seqs + unmutated_seqs:,})',
                  fontsize=9)

    # Pie 2: by read count
    ax2.pie([mutated_reads, unmutated_reads], labels=['Mutated', 'Unmutated'],
            autopct='%1.1f%%', colors=colors, explode=explode,
            startangle=90, textprops={'fontsize': 8})
    ax2.set_title(f'By Reads\n({mutated_reads:,} / {total_reads:,})',
                  fontsize=9)

    fig.suptitle('Editing Efficiency', fontsize=11, y=0.96)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    return _img_to_b64(fig)


def _gen_allele_heatmap(results: List[Dict], ref_seq: str,
                         cutsites: list = None, top_n: int = 50,
                         window_start: int = 0,
                         window_end: int = None,
                         primer5_len: int = 23,
                         primer3_len: int = 33) -> str:
    """
    Allele heatmap: one row per allele, one column per base position.

    Displays the top N most frequent alleles as a base-level alignment grid,
    with color-coded bases, deletion/insertion markers, and cutsite annotations.
    Supports a configurable display window over the reference region.

    Args:
        results: List of alignment result dicts.
        ref_seq: Full reference sequence.
        cutsites: Cutsite regions for annotations on the heatmap.
        top_n: Number of top alleles to display.
        window_start: Start position for display window (0-indexed, default 0).
        window_end: End position for display window (inclusive, default full length).

    Returns:
        Base64-encoded PNG image string, or empty string on failure.
    """
    if not results or not ref_seq:
        return ""

    ref_len = len(ref_seq)
    if window_end is None:
        window_end = ref_len - 1
    window_end = min(window_end, ref_len - 1)
    window_start = 0 if window_start is None else max(0, window_start)
    window_len = window_end - window_start + 1
    if window_len <= 0:
        return ""

    # 1. Aggregate alleles (deduplicate by internal region to avoid
    #    primer-variation-induced allele fragmentation)
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

        # Extract internal region (trim primers)
        pos_map, total_ref = _build_ref_pos_map_full(ar)
        internal_start = next((i for i, p in enumerate(pos_map) if p == primer5_len), primer5_len)
        internal_end = next((i for i, p in enumerate(pos_map) if p == total_ref - primer3_len),
                            len(ar) - primer3_len)
        internal_end = max(internal_start, min(internal_end, len(ar)))
        internal_key = aq[internal_start:internal_end]

        allele_counts[internal_key] = allele_counts.get(internal_key, 0) + rc
        if internal_key not in allele_data:
            # Save first full alignment data for display
            allele_data[internal_key] = {
                "aligned_ref": ar,
                "aligned_query": aq,
                "internal_start": internal_start,
                "internal_end": internal_end,
                "readCount": rc
            }
    if not allele_counts:
        return ""

    # 2. Sort and take top N
    sorted_alleles = sorted(allele_counts.items(), key=lambda x: -x[1])
    top_alleles = sorted_alleles[:top_n]
    total_reads = sum(allele_counts.values())

    # 3. Insertion-aware display column construction
    window_ref_bases = ref_seq[window_start:window_end + 1]
    n_cols = window_len  # Number of reference columns in display window

    def _build_display_with_insertions(aligned_ref, aligned_query):
        """Build display columns that accommodate insertions.

        Reference bases each occupy one column. Insertion columns (reference
        gap, query base) create extra columns beyond the reference width.
        Subsequent bases shift right to make room.

        Returns:
            Tuple of (seq, ref, ins_blocks, width):
            - seq: query base per column (length >= window_len)
            - ref: reference base per column (empty string for insertion columns)
            - ins_blocks: list of [(start_col, end_col)] for insertion blocks
            - width: total number of display columns
        """
        seq = []
        ref = []
        ins_blocks = []
        ref_pos = 0
        col = 0
        in_ins = False
        ins_start = 0
        for rb, qb in zip(aligned_ref, aligned_query):
            if rb != '-':   # Reference has base → normal column
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
            elif qb != '-':  # Ref gap, query base → insertion column
                if window_start <= ref_pos <= window_end:
                    if not in_ins:
                        ins_start = col
                        in_ins = True
                    seq.append(qb)
                    ref.append('')
                    col += 1
        if in_ins:
            ins_blocks.append((ins_start, col - 1))
        # Pad to at least window_len reference-mapped columns
        ref_cnt = sum(1 for r in ref if r != '')
        while ref_cnt < window_len:
            seq.append('-')
            ref.append('')
            col += 1
            ref_cnt += 1
        return seq, ref, ins_blocks, col

    # Pre-compute display columns for all allele rows, truncate to n_cols
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
        # Clip to n_cols, discarding columns beyond reference width
        if len(sq) > n_cols:
            # Fix insertion block coordinates: discard blocks past n_cols
            clipped_blk = []
            for start, end in blk:
                if start >= n_cols:
                    continue
                clipped_blk.append((start, min(end, n_cols - 1)))
            sq = sq[:n_cols]
            rf = rf[:n_cols]
            blk = clipped_blk
        # Pad to n_cols (alleles shorter than reference get gaps)
        while len(sq) < n_cols:
            sq.append('-')
            rf.append('')
        row_data.append((sq, rf, blk, total_rc))

    grid_cols = n_cols  # Fixed to reference column count

    # 4. Base color map (hex → RGB for imshow)
    import numpy as np

    def _hex_to_rgb(h):
        h = h.lstrip('#')
        return (int(h[0:2], 16) / 255.0, int(h[2:4], 16) / 255.0, int(h[4:6], 16) / 255.0)

    base_colors_rgb = {
        'A': _hex_to_rgb('#cbe9cb'), 'C': _hex_to_rgb('#fee5ce'),
        'G': _hex_to_rgb('#ffffd6'), 'T': _hex_to_rgb('#e5deed'),
        'N': _hex_to_rgb('#e9e9e9'), '-': _hex_to_rgb('#f0f0f0'),
    }
    ins_bg_rgb = _hex_to_rgb('#ffcccc')   # insertion column background
    del_bg_rgb = _hex_to_rgb('#f0f0f0')   # deletion cell background
    alt_bg_rgb = _hex_to_rgb('#f8f8f8')   # alternating row tint

    n_rows = len(top_alleles)

    # 5. Adaptive cell sizing based on reference column count
    if n_cols > 300:
        cell_size = 14
        font_size = 6.5
    elif n_cols > 200:
        cell_size = 18
        font_size = 7
    elif n_cols > 100:
        cell_size = 22
        font_size = 8
    elif n_cols > 60:
        cell_size = 26
        font_size = 9
    else:
        cell_size = 30
        font_size = 10

    label_w = 180  # Right-side annotation width (px)
    fig_w_px = grid_cols * cell_size + label_w

    # Calculate fig_h to match set_aspect('equal') so no letterboxing occurs.
    # With subplots_adjust(top=0.94): ppu_X = fig_w_px / x_range
    #                                 ppu_Y = fig_h_px * 0.94 / y_range
    # Set ppu_X = ppu_Y → fig_h_px = fig_w_px * y_range / (0.94 * x_range)
    x_range = grid_cols + 15  # xlim(-5, grid_cols + 10)
    y_range = n_rows + 5       # ylim(-2.5, n_rows + 2.5)
    fig_h_px = fig_w_px * y_range / (0.94 * x_range)

    dpi = 100
    fig_w = max(8, fig_w_px / dpi)
    fig_h = max(3, fig_h_px / dpi)

    # Build RGB grid for imshow (single render call instead of thousands of patches+texts)
    grid = np.zeros((n_rows, n_cols, 3), dtype=np.float32)
    for ai, (sq, rf, blk, total_rc) in enumerate(row_data):
        for ci in range(n_cols):
            base = sq[ci] if ci < len(sq) else '-'
            ref_base = rf[ci] if ci < len(rf) else ''
            if ref_base == '':       # insertion column
                color = ins_bg_rgb
            elif base == '-':        # deletion
                color = del_bg_rgb
            else:
                color = base_colors_rgb.get(base.upper(), (0.91, 0.91, 0.91))
            if ai % 2 == 1:
                color = tuple(c * 0.92 + 0.08 for c in color)  # blend with alt row bg
            grid[ai, ci] = color

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=0.94)
    ax.set_xlim(-5, grid_cols + 10)
    ax.set_ylim(-2.5, n_rows + 2.5)
    ax.set_aspect('equal')
    ax.axis('off')

    # 6. Render the allele grid as a single image
    ax.imshow(grid, extent=[0, n_cols, 0, n_rows], aspect='equal',
              interpolation='nearest', origin='upper', zorder=1)

    # 7. Reference sequence row (individual patches+text, ~O(n_cols) ≈ fast)
    ref_y = n_rows + 1.5
    for ci in range(n_cols):
        base = window_ref_bases[ci] if ci < len(window_ref_bases) else ' '
        color = '#e9e9e9'
        for k, v in base_colors_rgb.items():
            if k == base or k == base.upper():
                color = '#{:02x}{:02x}{:02x}'.format(
                    int(v[0]*255), int(v[1]*255), int(v[2]*255))
                break
        ax.add_patch(plt.Rectangle((ci, ref_y - 0.5), 1, 1,
                                    facecolor=color, edgecolor='white', linewidth=0.3))
        ax.text(ci + 0.5, ref_y, base, ha='center', va='center',
                fontsize=font_size, fontweight='normal')

    # Cutsite annotations on reference row
    if cutsites:
        for cs in cutsites:
            if cs.start > window_end or cs.end < window_start:
                continue
            cs_c0 = max(cs.start - window_start, 0)
            cs_c1 = min(cs.end - window_start, n_cols - 1)
            if cs_c0 <= cs_c1:
                ax.add_patch(plt.Rectangle(
                    (cs_c0, ref_y - 0.5), cs_c1 - cs_c0 + 1, 1,
                    facecolor='#888888', alpha=0.3, edgecolor=None, zorder=2
                ))

    # sgRNA region indicators (half-height grey bars just below reference row)
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
                facecolor='#bbbbbb', edgecolor='#999999', linewidth=0.3, zorder=2
            ))
            ax.text(col0 - 0.1, grna_y, cs.name.replace('Target', 'T'),
                    ha='right', va='center', fontsize=max(6, font_size - 0.5),
                    color='#666', fontweight='bold')

    # 8. Overlay: insertion red boxes, right-side labels, optional cell text
    for ai, (sq, rf, blk, total_rc) in enumerate(row_data):
        y = n_rows - 0.5 - ai

        # Red box marking insertion blocks
        for start, end in blk:
            ax.add_patch(plt.Rectangle(
                (start, y - 0.5), end - start + 1, 1,
                facecolor='none', edgecolor='#e94560', linewidth=1.5, zorder=5
            ))

        # Cell text: always show bases, colored by mutation type
        for ci in range(n_cols):
            base = sq[ci] if ci < len(sq) else '-'
            ref_base = rf[ci] if ci < len(rf) else ''
            is_del = (base == '-')
            is_ins_col = (ref_base == '')
            is_sub = (not is_del and not is_ins_col and base != ref_base)
            if is_ins_col:
                ax.text(ci + 0.5, y, base, ha='center', va='center',
                        fontsize=font_size, color='red', fontweight='bold')
            elif is_del:
                ax.text(ci + 0.5, y, '-', ha='center', va='center',
                        fontsize=font_size, color='#666', fontweight='bold')
            elif is_sub:
                ax.text(ci + 0.5, y, base, ha='center', va='center',
                        fontsize=font_size, fontweight='bold')
            else:
                ax.text(ci + 0.5, y, base, ha='center', va='center',
                        fontsize=font_size, fontweight='normal')

        # Right-side label with percentage and read count
        pct = total_rc / total_reads * 100
        ax.text(grid_cols + 0.3, y, f"{pct:.1f}% ({total_rc})",
                ha='left', va='center', fontsize=max(8, font_size + 1),
                clip_on=False)

    # 9. Legend (A, C, G, T, -, ins)
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
    return _img_to_b64(fig, bbox_inches=None)


def generate_charts(results: List[Dict], report_data: dict = None,
                     ref_length: int = 0, ref_seq: str = "",
                     cutsites: list = None,
                     allele_window_start: int = 0,
                     allele_window_end: int = None,
                     allele_top_n: int = 50,
                     primer5_len: int = 23,
                     primer3_len: int = 33,
                     summary_data: dict = None) -> Dict[str, str]:
    """
    Generate all chart images and return as a {chart_name: base64} dictionary.

    Currently produces three charts:
      - editing_pie: dual pie charts of editing efficiency (by sequence and by read)
      - indel_length: butterfly chart of insertion/deletion length distribution
      - allele_heatmap: base-level heatmap of top N alleles

    Args:
        results: List of alignment result dicts.
        report_data: Report data dict containing mutation_stats.
        ref_length: Reference sequence length.
        ref_seq: Reference sequence string (required for allele heatmap).
        cutsites: CutsiteRegion list for heatmap annotations.
        allele_window_start: Start position for heatmap window (default 0).
        allele_window_end: End position for heatmap window (inclusive).
        allele_top_n: Number of top alleles to show in heatmap (default 50).
        summary_data: Optional summary_data dict from save_summary_tables.
                      Provides ins_length_reads/del_length_reads for the indel chart.
    """
    charts = {}
    if not _ensure_mpl():
        return charts
    try:
        if report_data:
            ins_r = summary_data.get("ins_length_reads") if summary_data else None
            del_r = summary_data.get("del_length_reads") if summary_data else None
            charts['indel_length'] = _gen_indel_length_chart(
                ins_length_reads=ins_r, del_length_reads=del_r,
                report_data=report_data)
            charts['editing_pie'] = _gen_editing_pie_charts(report_data)
        if ref_seq and ref_length > 0:
            ws = allele_window_start if allele_window_start is not None else 0
            we = allele_window_end if allele_window_end is not None else ref_length - 1
            charts['allele_heatmap'] = _gen_allele_heatmap(
                results, ref_seq, cutsites=cutsites, top_n=allele_top_n,
                window_start=ws, window_end=we,
                primer5_len=primer5_len, primer3_len=primer3_len,
            )
    except Exception as e:
        log.warning("Chart generation failed - %s", e)
    return charts
