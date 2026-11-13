"""Render one assay set's F/Probe/R binding sites on a single reference genome as an SVG
genome track (ruler + optional per-base sequence + colored arrows), mirroring the
'참조유전체비교' slides. Dependency-free (no matplotlib), consistent with svg_charts.py."""

ROLE_COLOR = {"F": "#4C78A8", "P": "#F2A65A", "R": "#54A24B"}
ROLE_LABEL = {"F": "Forward", "P": "Probe", "R": "Reverse"}


def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def genome_track_svg(set_id, ref_seq, hits, width=900, pad_bp=25, max_letter_window=420):
    """hits: list of dicts with keys role, oligo_id, found(bool), start, end, strand, mismatches.
    Coordinates are 0-based half-open, as reported by mfeprimer search."""
    found_hits = [h for h in hits if h["found"]]
    ref_len = len(ref_seq)

    if found_hits:
        win_start = max(0, min(h["start"] for h in found_hits) - pad_bp)
        win_end = min(ref_len, max(h["end"] for h in found_hits) + pad_bp)
    else:
        win_start, win_end = 0, min(ref_len, 200)

    win_len = max(1, win_end - win_start)
    left_margin, right_margin, top_margin = 20, 20, 30
    plot_w = width - left_margin - right_margin
    row_h = 46
    n_rows = len(hits)
    show_letters = win_len <= max_letter_window
    letters_h = 18 if show_letters else 0
    ruler_h = 22
    height = top_margin + ruler_h + letters_h + n_rows * row_h + 20

    def x_of(pos):
        return left_margin + (pos - win_start) / win_len * plot_w

    parts = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
             f'font-family="ui-monospace,monospace" font-size="11">']
    parts.append(f'<text x="{left_margin}" y="14" font-family="ui-sans-serif,system-ui,sans-serif" '
                 f'font-size="12" font-weight="600" fill="currentColor">{_esc(set_id)} '
                 f'&mdash; {_esc(hits[0]["chrom_label"]) if hits else ""}</text>')

    ruler_y = top_margin + ruler_h
    step = max(10, round(win_len / 8 / 10) * 10)
    tick = win_start - (win_start % step) + step if win_start % step else win_start
    while tick < win_end:
        tx = x_of(tick)
        parts.append(f'<line x1="{tx:.1f}" y1="{ruler_y-6}" x2="{tx:.1f}" y2="{ruler_y}" '
                     f'stroke="currentColor" stroke-opacity="0.4"/>')
        parts.append(f'<text x="{tx:.1f}" y="{ruler_y-8}" text-anchor="middle" fill="currentColor" '
                     f'opacity="0.6" font-size="9">{tick}</text>')
        tick += step
    parts.append(f'<line x1="{left_margin}" y1="{ruler_y}" x2="{left_margin+plot_w}" y2="{ruler_y}" '
                 f'stroke="currentColor" stroke-opacity="0.3"/>')

    if show_letters:
        window_seq = ref_seq[win_start:win_end]
        ly = ruler_y + 14
        char_w = plot_w / win_len
        # draw as a single text run with letter-spacing approximated by explicit x per char
        xs = " ".join(f'{(left_margin + i*char_w):.1f}' for i in range(len(window_seq)))
        parts.append(f'<text x="{xs}" y="{ly}" fill="currentColor" opacity="0.55" font-size="9">'
                     f'{_esc(window_seq)}</text>')

    row_top = ruler_y + letters_h + 10
    for i, h in enumerate(hits):
        ry = row_top + i * row_h
        color = ROLE_COLOR.get(h["role"], "#999")
        label = f'{h["oligo_id"]} ({ROLE_LABEL.get(h["role"], h["role"])})'
        parts.append(f'<text x="{left_margin}" y="{ry+10:.1f}" fill="currentColor" '
                     f'font-family="ui-sans-serif,system-ui,sans-serif" font-size="10" '
                     f'font-weight="600">{_esc(label)}</text>')
        if not h["found"]:
            bx0, bx1 = x_of(win_start + win_len * 0.3), x_of(win_start + win_len * 0.7)
            by = ry + 16
            text = "PROBE X" if h["role"] == "P" else "NO MATCH"
            parts.append(f'<rect x="{bx0:.1f}" y="{by:.1f}" width="{bx1-bx0:.1f}" height="20" '
                         f'fill="#E4575622" stroke="#E45756" stroke-dasharray="4,3" rx="4"/>')
            parts.append(f'<text x="{(bx0+bx1)/2:.1f}" y="{by+14:.1f}" text-anchor="middle" '
                         f'fill="#E45756" font-family="ui-sans-serif,system-ui,sans-serif" '
                         f'font-size="10" font-weight="700">{text}</text>')
            continue

        x0, x1 = x_of(max(win_start, h["start"])), x_of(min(win_end, h["end"]))
        y = ry + 16
        bar_h = 12
        if h["strand"] == "-":
            head = x0
            parts.append(f'<polygon points="{x0:.1f},{y-6:.1f} {x0+10:.1f},{y:.1f} {x0:.1f},{y+bar_h+6:.1f}" fill="{color}"/>')
            parts.append(f'<rect x="{x0+8:.1f}" y="{y:.1f}" width="{max(2,x1-x0-8):.1f}" height="{bar_h}" fill="{color}"/>')
        else:
            parts.append(f'<rect x="{x0:.1f}" y="{y:.1f}" width="{max(2,x1-x0-8):.1f}" height="{bar_h}" fill="{color}"/>')
            parts.append(f'<polygon points="{x1:.1f},{y-6:.1f} {x1-10:.1f},{y:.1f} {x1:.1f},{y+bar_h+6:.1f}" fill="{color}"/>')
        mm = h["mismatches"]
        mm_color = "#54A24B" if mm == 0 else ("#F2A65A" if mm <= 2 else "#E45756")
        parts.append(f'<text x="{(x0+x1)/2:.1f}" y="{y+bar_h+18:.1f}" text-anchor="middle" '
                     f'fill="{mm_color}" font-family="ui-sans-serif,system-ui,sans-serif" '
                     f'font-size="10" font-weight="600">mismatches: {mm}</text>')

    parts.append('</svg>')
    return "".join(parts)
