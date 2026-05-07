#!/usr/bin/env python3
import csv
import math
import os
import sys
from pathlib import Path

WIDTH = 1000
HEIGHT = 620
MARGIN_LEFT = 90
MARGIN_RIGHT = 30
MARGIN_TOP = 50
MARGIN_BOTTOM = 80
PLOT_W = WIDTH - MARGIN_LEFT - MARGIN_RIGHT
PLOT_H = HEIGHT - MARGIN_TOP - MARGIN_BOTTOM
COLORS = {"java": "#1f77b4", "rust": "#d62728"}
LABELS = {"java": "Java", "rust": "Rust auto"}


def load_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def nice_ticks(max_value, count=5):
    if max_value <= 0:
        return [0, 1]
    raw = max_value / count
    power = 10 ** math.floor(math.log10(raw))
    for mult in (1, 2, 2.5, 5, 10):
        step = mult * power
        if raw <= step:
            break
    top = math.ceil(max_value / step) * step
    ticks = [i * step for i in range(int(round(top / step)) + 1)]
    return ticks


def x_positions(xs):
    xmin, xmax = min(xs), max(xs)
    if xmin == xmax:
        return {xmin: MARGIN_LEFT + PLOT_W / 2}
    log_min = math.log10(xmin)
    log_max = math.log10(xmax)
    out = {}
    for x in xs:
        t = (math.log10(x) - log_min) / (log_max - log_min)
        out[x] = MARGIN_LEFT + t * PLOT_W
    return out


def y_pos(y, ymax):
    if ymax <= 0:
        return MARGIN_TOP + PLOT_H
    return MARGIN_TOP + PLOT_H - (y / ymax) * PLOT_H


def polyline(points, color):
    coords = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    return f'<polyline fill="none" stroke="{color}" stroke-width="3" points="{coords}" />'


def circles(points, color):
    return "\n".join(
        f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4.5" fill="{color}" />' for x, y in points
    )


def render_chart(series, y_key, y_label, title, out_path):
    xs = sorted({int(row['reads']) for rows in series.values() for row in rows})
    xmap = x_positions(xs)
    ymax = max(float(row[y_key]) for rows in series.values() for row in rows)
    yticks = nice_ticks(ymax)
    ymax = yticks[-1]

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
        '<rect width="100%" height="100%" fill="white" />',
        f'<text x="{WIDTH/2:.1f}" y="28" text-anchor="middle" font-size="24" font-family="sans-serif">{title}</text>',
        f'<text x="{WIDTH/2:.1f}" y="{HEIGHT-18}" text-anchor="middle" font-size="16" font-family="sans-serif">Read pairs processed</text>',
        f'<text x="22" y="{HEIGHT/2:.1f}" text-anchor="middle" font-size="16" font-family="sans-serif" transform="rotate(-90 22 {HEIGHT/2:.1f})">{y_label}</text>',
        f'<line x1="{MARGIN_LEFT}" y1="{MARGIN_TOP+PLOT_H}" x2="{MARGIN_LEFT+PLOT_W}" y2="{MARGIN_TOP+PLOT_H}" stroke="#111" stroke-width="1.5" />',
        f'<line x1="{MARGIN_LEFT}" y1="{MARGIN_TOP}" x2="{MARGIN_LEFT}" y2="{MARGIN_TOP+PLOT_H}" stroke="#111" stroke-width="1.5" />',
    ]

    for tick in yticks:
        y = y_pos(tick, ymax)
        parts.append(f'<line x1="{MARGIN_LEFT}" y1="{y:.2f}" x2="{MARGIN_LEFT+PLOT_W}" y2="{y:.2f}" stroke="#ddd" stroke-width="1" />')
        label = f"{tick:.2f}" if ymax <= 10 else f"{tick:.1f}" if ymax <= 100 else f"{tick:.0f}"
        parts.append(f'<text x="{MARGIN_LEFT-10}" y="{y+5:.2f}" text-anchor="end" font-size="13" font-family="sans-serif">{label}</text>')

    for x in xs:
        xp = xmap[x]
        parts.append(f'<line x1="{xp:.2f}" y1="{MARGIN_TOP}" x2="{xp:.2f}" y2="{MARGIN_TOP+PLOT_H}" stroke="#f0f0f0" stroke-width="1" />')
        parts.append(f'<text x="{xp:.2f}" y="{MARGIN_TOP+PLOT_H+22}" text-anchor="middle" font-size="13" font-family="sans-serif">{x:,}</text>')

    legend_x = WIDTH - 210
    legend_y = 52
    parts.append(f'<rect x="{legend_x}" y="{legend_y}" width="170" height="58" fill="white" stroke="#ccc" />')
    for i, tool in enumerate(("java", "rust")):
        yy = legend_y + 20 + i * 22
        parts.append(f'<line x1="{legend_x+12}" y1="{yy}" x2="{legend_x+42}" y2="{yy}" stroke="{COLORS[tool]}" stroke-width="3" />')
        parts.append(f'<circle cx="{legend_x+27}" cy="{yy}" r="4.5" fill="{COLORS[tool]}" />')
        parts.append(f'<text x="{legend_x+52}" y="{yy+5}" font-size="14" font-family="sans-serif">{LABELS[tool]}</text>')

    for tool, rows in series.items():
        pts = []
        for row in sorted(rows, key=lambda r: int(r['reads'])):
            pts.append((xmap[int(row['reads'])], y_pos(float(row[y_key]), ymax)))
        parts.append(polyline(pts, COLORS[tool]))
        parts.append(circles(pts, COLORS[tool]))

    parts.append('</svg>')
    Path(out_path).write_text("\n".join(parts))


def write_summary(rows, out_path):
    xs = sorted({int(r['reads']) for r in rows})
    by_reads = {x: {} for x in xs}
    for row in rows:
        by_reads[int(row['reads'])][row['tool']] = row
    with open(out_path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow([
            'reads', 'java_seconds', 'rust_seconds', 'java_rss_kb', 'rust_rss_kb',
            'rust_speed_vs_java', 'rust_memory_vs_java'
        ])
        for x in xs:
            j = by_reads[x]['java']
            r = by_reads[x]['rust']
            writer.writerow([
                x,
                j['elapsed_seconds'],
                r['elapsed_seconds'],
                j['max_rss_kb'],
                r['max_rss_kb'],
                r['speed_vs_java'],
                r['memory_vs_java'],
            ])


def main():
    if len(sys.argv) != 3:
        raise SystemExit('Usage: plot_scaling_java_rust_auto_svg.py <scaling.tsv> <outdir>')
    scaling_tsv = sys.argv[1]
    outdir = Path(sys.argv[2])
    rows = load_rows(scaling_tsv)
    series = {'java': [], 'rust': []}
    for row in rows:
        series[row['tool']].append(row)
    render_chart(series, 'elapsed_seconds', 'Elapsed seconds', 'Java vs Rust auto: speed vs read count', outdir / 'speed_vs_reads.svg')
    render_chart(series, 'max_rss_gib', 'Peak RSS (GiB)', 'Java vs Rust auto: memory vs read count', outdir / 'memory_vs_reads.svg')
    write_summary(rows, outdir / 'summary.tsv')


if __name__ == '__main__':
    main()
