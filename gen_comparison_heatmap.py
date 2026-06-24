"""Generate single unified comparison heatmap: old (left) vs new (right) side by side."""
import csv, sys, os, io, base64, multiprocessing as mp
try:
    mp.set_start_method('fork')
except RuntimeError:
    pass
from concurrent.futures import ProcessPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from crisviper.lineage import lineage_tracer_align, get_amplicon_structure

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ── Module-level globals for multiprocessing ──
raw = open("examples/reference.fa").read().strip()
lines = raw.split("\n")
_REF_SEQ = "".join(lines[1:]).replace(" ", "").replace("\n", "")
_CUTSITES = get_amplicon_structure(_REF_SEQ)

def find_isolated(aq, ar, max_len=4):
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
            if l <= max_len:
                left_gap = s > 0 and aq[s-1] == '-'
                right_gap = i < len(aq) and aq[i] == '-'
                if left_gap and right_gap and any(ar[j] != '-' for j in range(s, i)):
                    segs.append((s, l))
        else:
            i += 1
    return segs

def align_pair(args):
    name, seq, rc = args
    try:
        s0, aq0, ar0, st0 = lineage_tracer_align(
            _REF_SEQ, seq, _CUTSITES,
            base_gap_open=-2.0, base_gap_extend=-0.1,
            match_score=2.0, mismatch_penalty=-3.0,
            base_gap_exit=0.0, min_scale=1.0, max_scale=6.0)
        old_segs = find_isolated(aq0, ar0)
        s1, aq1, ar1, st1 = lineage_tracer_align(
            _REF_SEQ, seq, _CUTSITES,
            base_gap_open=-2.0, base_gap_extend=-0.1,
            match_score=2.0, mismatch_penalty=-3.0,
            base_gap_exit=-3.5, min_scale=1.0, max_scale=6.0)
        new_segs = find_isolated(aq1, ar1)
        return (name, rc, aq0, ar0, aq1, ar1, len(old_segs), len(new_segs))
    except:
        return (name, rc, "", "", "", "", -1, -1)

def extract_window(aq, ar, ref_seq, win_start, win_end):
    """Extract query bases at each reference position within window."""
    n = win_end - win_start + 1
    result = []
    ref_pos = 0
    for col, c in enumerate(ar):
        if c != '-':
            if win_start <= ref_pos <= win_end:
                base = aq[col] if col < len(aq) else '-'
                result.append(base)
            ref_pos += 1
            if ref_pos > win_end:
                break
    while len(result) < n:
        result.append('-')
    return ''.join(result[:n])

def main():
    # Scan TSV
    print("Reading queries...")
    queries = []
    with open("examples/output/all_queries.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            queries.append((row["readName"], row["seq"], int(row.get("readCount", 1))))
    print(f"Total: {len(queries)}")
    sample = queries[:2000]

    # Align
    print(f"Aligning {len(sample)} sequences...")
    fixed = []
    with ProcessPoolExecutor(max_workers=12) as pool:
        futures = [pool.submit(align_pair, q) for q in sample]
        for i, f in enumerate(as_completed(futures)):
            if (i+1) % 200 == 0:
                print(f"  {i+1}/{len(sample)}")
            name, rc, aq0, ar0, aq1, ar1, n_old, n_new = f.result()
            if n_old > 0 and n_new == 0:
                fixed.append({
                    "readName": name, "readCount": rc,
                    "aligned_query_old": aq0, "aligned_ref_old": ar0,
                    "aligned_query_new": aq1, "aligned_ref_new": ar1,
                })

    fixed.sort(key=lambda x: -x["readCount"])
    top20 = fixed[:20]
    print(f"\nFound {len(fixed)} fixed. Top 20:")
    for r in top20:
        print(f"  {r['readName']}: reads={r['readCount']}")

    if not top20:
        print("No fixed sequences!")
        return

    # ── Draw a SINGLE figure with old on left, new on right ──
    ref_len = len(_REF_SEQ)
    win_start, win_end = 0, ref_len - 1
    n_cols = win_end - win_start + 1
    window_ref = _REF_SEQ[win_start:win_end + 1]
    n_rows = len(top20)

    base_colors = {
        'A': '#cbe9cb', 'C': '#fee5ce', 'G': '#ffffd6',
        'T': '#e5deed', 'U': '#e5deed', 'N': '#e9e9e9',
        '-': '#f0f0f0',
    }
    gap_color = '#f0f0f0'

    # Extract display data for all rows
    old_display = []
    new_display = []
    old_highlights = []

    for rec in top20:
        old_display.append(extract_window(rec["aligned_query_old"], rec["aligned_ref_old"],
                                           _REF_SEQ, win_start, win_end))
        new_display.append(extract_window(rec["aligned_query_new"], rec["aligned_ref_new"],
                                           _REF_SEQ, win_start, win_end))
        segs = find_isolated(rec["aligned_query_old"], rec["aligned_ref_old"])
        old_highlights.append({p for p, _ in segs if win_start <= p <= win_end})

    # Layout: 2 side-by-side panels in one figure
    cell_size = 11
    font_size = 5
    label_w = 160

    # Each half: n_cols columns + label
    half_w_px = n_cols * cell_size + label_w
    total_w_px = half_w_px * 2 + cell_size * 2  # gap between panels
    row_h_px = cell_size * (n_rows + 4)

    dpi = 100
    fig_w = total_w_px / dpi
    fig_h = row_h_px / dpi

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')
    plt.axis('off')

    # Left panel
    ax_l = fig.add_axes([0, 0, 0.48, 1])
    ax_l.set_xlim(0, n_cols)
    ax_l.set_ylim(-3, n_rows + 3)
    ax_l.set_aspect('equal')
    ax_l.axis('off')
    ax_l.set_title("Without geb_profile (old)", fontsize=10, fontweight='bold', pad=8)

    # Right panel
    ax_r = fig.add_axes([0.52, 0, 0.48, 1])
    ax_r.set_xlim(0, n_cols)
    ax_r.set_ylim(-3, n_rows + 3)
    ax_r.set_aspect('equal')
    ax_r.axis('off')
    ax_r.set_title("With geb_profile (new, base=-3.5)", fontsize=10, fontweight='bold', pad=8)

    def draw_panel(ax, displays, highlights, show_hl=False):
        """Draw heatmap panel. All rows use SAME reference row."""
        ref_y = n_rows + 1.5

        # Reference row
        for ci in range(n_cols):
            base = window_ref[ci] if ci < len(window_ref) else 'N'
            in_cs = any(cs.start <= (win_start + ci) <= cs.end for cs in _CUTSITES)
            color = '#e0e0e0' if in_cs else base_colors.get(base, '#e9e9e9')
            ax.add_patch(Rectangle((ci, ref_y - 0.5), 1, 1,
                        facecolor=color, edgecolor='white', linewidth=0.3))
            ax.text(ci + 0.5, ref_y, base, ha='center', va='center',
                    fontsize=font_size, fontweight='normal')

        # Cutsite shading
        for cs in _CUTSITES:
            if cs.start > win_end or cs.end < win_start:
                continue
            c0 = max(cs.start - win_start, 0)
            c1 = min(cs.end - win_start, n_cols - 1)
            if c0 <= c1:
                ax.add_patch(Rectangle(
                    (c0, ref_y - 0.5), c1 - c0 + 1, 1,
                    facecolor='#888888', alpha=0.3, edgecolor=None, zorder=0))
                mid = (c0 + c1) / 2
                ax.text(mid, ref_y + 0.7, cs.name, ha='center', va='bottom',
                        fontsize=max(4, font_size-1), color='#666', fontstyle='italic')

        # Allele rows (use pre-sorted, fixed order)
        for ri in range(n_rows):
            y = ri + 0.5
            bases = displays[ri]
            hl = highlights[ri]

            # Alternating row background
            if ri % 2 == 1:
                ax.add_patch(Rectangle((0, y - 0.5), n_cols, 1,
                            facecolor='#f8f8f8', edgecolor=None, zorder=0))

            for ci in range(n_cols):
                base = bases[ci] if ci < len(bases) else '-'
                bg = base_colors.get(base, gap_color) if base != '-' else gap_color
                is_hl = show_hl and (win_start + ci) in hl
                ec = '#e94560' if is_hl else 'white'
                lw = 2.0 if is_hl else 0.3
                ax.add_patch(Rectangle((ci, y - 0.5), 1, 1,
                            facecolor=bg, edgecolor=ec, linewidth=lw))
                if base != '-':
                    ax.text(ci + 0.5, y, base, ha='center', va='center',
                            fontsize=font_size)

            # Right label
            name = top20[ri]['readName'][:14]
            rc = top20[ri]['readCount']
            ax.text(n_cols + 0.5, y, f"{ri+1}. {name} (n={rc})",
                    ha='left', va='center', fontsize=font_size)

    draw_panel(ax_l, old_display, old_highlights, show_hl=True)
    draw_panel(ax_r, new_display, [set() for _ in range(n_rows)], show_hl=False)

    # Vertical divider
    mid_x = 0.51
    fig.lines.append(plt.Line2D([mid_x, mid_x], [0.02, 0.98],
                                transform=fig.transFigure, color='#ccc',
                                linewidth=1, linestyle='-'))

    # ── Save ──
    fig.savefig("examples/output/comparison_heatmap.png", dpi=dpi,
                bbox_inches='tight', facecolor='white')
    print(f"Saved: examples/output/comparison_heatmap.png")

    # Also embed in HTML
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight', facecolor='white')
    buf.seek(0)
    img_b64 = base64.b64encode(buf.read()).decode()

    html = f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8">
<title>Geb-profile comparison — Allele Heatmap</title>
<style>
*{{margin:0;padding:0;box-sizing:border-box;}}
body{{background:#f8fafc;color:#1e293b;padding:20px;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;}}
h1{{font-size:18px;margin-bottom:10px;padding-bottom:6px;border-bottom:2px solid #2563eb;}}
.info{{background:#fff;padding:10px 14px;border-radius:6px;margin:0 0 16px;font-size:12px;line-height:1.5;border:1px solid #e2e8f0;}}
.img-wrap{{background:#fff;border-radius:8px;padding:12px;border:1px solid #e2e8f0;overflow-x:auto;text-align:center;}}
.img-wrap img{{max-width:none;}}
footer{{margin-top:16px;color:#94a3b8;font-size:11px;}}
</style>
</head><body>
<h1>Geb-profile comparison: Before vs After — Allele Heatmap</h1>
<div class="info">
<strong>Left:</strong> alignment without geb_profile (red border = isolated match fragment).<br>
<strong>Right:</strong> alignment with geb_profile (base=-3.5). Top {n_rows} alleles by readCount from {len(fixed)} corrected sequences (sampled {len(sample)}).
</div>
<div class="img-wrap">
<img src="data:image/png;base64,{img_b64}" alt="Old vs New comparison heatmap">
</div>
<footer>Reference: {_REF_SEQ[:50]}... | Both panels share the same reference coordinate system and row ordering.</footer>
</body></html>"""

    with open("examples/output/comparison_heatmap.html", "w") as f:
        f.write(html)
    print(f"Saved: examples/output/comparison_heatmap.html")
    print(f"Done!")

if __name__ == '__main__':
    main()
