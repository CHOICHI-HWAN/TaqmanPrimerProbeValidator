"""Tiny dependency-free SVG chart helpers for the dashboard (no matplotlib needed)."""


def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def grouped_hbar_chart(categories, series, colors, max_value=100.0, unit="%",
                        row_height=54, bar_gap=4, left_margin=190, chart_width=640,
                        title=None):
    """categories: list[str] (rows, top to bottom)
    series: dict[str, list[float]] each list aligned with `categories`
    colors: dict[str, str] series name -> hex color
    """
    series_names = list(series.keys())
    n_series = len(series_names)
    bar_h = (row_height - bar_gap * (n_series + 1)) / n_series
    top_margin = 46 if title else 30
    plot_h = row_height * len(categories)
    height = top_margin + plot_h + 30
    width = left_margin + chart_width + 70

    parts = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
             f'font-family="ui-sans-serif,system-ui,sans-serif" font-size="12">']
    parts.append(f'<rect x="0" y="0" width="{width}" height="{height}" fill="none"/>')

    if title:
        parts.append(f'<text x="{width/2}" y="20" text-anchor="middle" font-size="14" '
                      f'font-weight="600" fill="currentColor">{_esc(title)}</text>')

    # legend
    lx = left_margin
    ly = top_margin - 14
    for name in series_names:
        parts.append(f'<rect x="{lx}" y="{ly-9}" width="10" height="10" fill="{colors[name]}" rx="2"/>')
        parts.append(f'<text x="{lx+14}" y="{ly}" fill="currentColor" opacity="0.8">{_esc(name)}</text>')
        lx += 16 + 8 * len(name) + 22

    # gridlines + x-axis ticks
    for frac in (0, 0.25, 0.5, 0.75, 1.0):
        gx = left_margin + chart_width * frac
        parts.append(f'<line x1="{gx}" y1="{top_margin}" x2="{gx}" y2="{top_margin+plot_h}" '
                      f'stroke="currentColor" stroke-opacity="0.12"/>')
        label = f"{int(max_value*frac)}{unit}"
        parts.append(f'<text x="{gx}" y="{top_margin+plot_h+18}" text-anchor="middle" '
                      f'fill="currentColor" opacity="0.6" font-size="10">{label}</text>')

    for i, cat in enumerate(categories):
        row_y = top_margin + i * row_height
        parts.append(f'<text x="{left_margin-10}" y="{row_y+row_height/2+4}" text-anchor="end" '
                      f'fill="currentColor" font-size="11">{_esc(cat)}</text>')
        by = row_y + bar_gap
        for name in series_names:
            val = series[name][i]
            bw = max(0.0, min(1.0, val / max_value)) * chart_width
            parts.append(f'<rect x="{left_margin}" y="{by:.1f}" width="{bw:.1f}" height="{bar_h:.1f}" '
                          f'fill="{colors[name]}" rx="2"/>')
            parts.append(f'<text x="{left_margin+bw+6:.1f}" y="{by+bar_h/2+3.5:.1f}" '
                          f'fill="currentColor" opacity="0.75" font-size="10">{val:.1f}{unit}</text>')
            by += bar_h + bar_gap
        if i < len(categories) - 1:
            parts.append(f'<line x1="{left_margin}" y1="{row_y+row_height}" x2="{left_margin+chart_width}" '
                          f'y2="{row_y+row_height}" stroke="currentColor" stroke-opacity="0.06"/>')

    parts.append('</svg>')
    return "".join(parts)


def simple_bar_chart(categories, values, colors=None, max_value=None, unit="", width=640, height=240,
                      title=None):
    """Simple vertical bar chart for a handful of categories (e.g. dataset sizes)."""
    max_value = max_value or (max(values) * 1.15 if values else 1)
    n = len(categories)
    left_margin, bottom_margin, top_margin = 50, 60, 40 if title else 20
    plot_w = width - left_margin - 20
    plot_h = height - top_margin - bottom_margin
    bar_w = plot_w / n * 0.6
    gap = plot_w / n

    parts = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
             f'font-family="ui-sans-serif,system-ui,sans-serif" font-size="12">']
    if title:
        parts.append(f'<text x="{width/2}" y="20" text-anchor="middle" font-size="14" '
                      f'font-weight="600" fill="currentColor">{_esc(title)}</text>')
    for frac in (0, 0.25, 0.5, 0.75, 1.0):
        gy = top_margin + plot_h * (1 - frac)
        parts.append(f'<line x1="{left_margin}" y1="{gy:.1f}" x2="{left_margin+plot_w}" y2="{gy:.1f}" '
                      f'stroke="currentColor" stroke-opacity="0.12"/>')
        parts.append(f'<text x="{left_margin-8}" y="{gy+3:.1f}" text-anchor="end" fill="currentColor" '
                      f'opacity="0.6" font-size="10">{int(max_value*frac)}{unit}</text>')
    for i, (cat, val) in enumerate(zip(categories, values)):
        bh = plot_h * (val / max_value) if max_value else 0
        bx = left_margin + gap * i + (gap - bar_w) / 2
        by = top_margin + plot_h - bh
        color = colors[i] if colors else "#4C78A8"
        parts.append(f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bar_w:.1f}" height="{bh:.1f}" '
                      f'fill="{color}" rx="3"/>')
        parts.append(f'<text x="{bx+bar_w/2:.1f}" y="{by-6:.1f}" text-anchor="middle" '
                      f'fill="currentColor" font-size="11">{val}{unit}</text>')
        parts.append(f'<text x="{bx+bar_w/2:.1f}" y="{top_margin+plot_h+16:.1f}" text-anchor="middle" '
                      f'fill="currentColor" opacity="0.8" font-size="11">{_esc(cat)}</text>')
    parts.append('</svg>')
    return "".join(parts)


def heatmap_table(row_labels, col_labels, matrix, color="76,120,168"):
    """matrix[row][col] -> int; returns an HTML <table> with intensity-shaded cells."""
    flat = [v for row in matrix for v in row if v]
    vmax = max(flat) if flat else 1
    html = ['<table class="heatmap"><thead><tr><th></th>']
    for c in col_labels:
        html.append(f'<th>{_esc(c)}</th>')
    html.append('</tr></thead><tbody>')
    for r, row in zip(row_labels, matrix):
        html.append(f'<tr><th>{_esc(r)}</th>')
        for v in row:
            alpha = 0.08 + 0.82 * (v / vmax) if vmax else 0.08
            bg = f"rgba({color},{alpha:.2f})"
            txt = v if v else ""
            html.append(f'<td style="background:{bg}">{txt}</td>')
        html.append('</tr>')
    html.append('</tbody></table>')
    return "".join(html)
