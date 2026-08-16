#!/usr/bin/env python3
"""Render the h-refinement convergence study as a log-log SVG.

Reads the p,h,error CSVs written by the solver and produces a two panel plot,
one panel per spatial dimension, with a dashed reference line of the
theoretical slope O(h^(p+1)) behind each series.

Standard library only, so it runs anywhere the solver builds.

    python3 scripts/plot_convergence.py

Regenerate the inputs by setting `dimension` in main/main.cpp and copying
main/hp_error.csv over the corresponding file in docs/.
"""

import csv
import math
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PANELS = [
    ("One dimension", os.path.join(ROOT, "docs", "convergence_1d.csv")),
    ("Two dimensions", os.path.join(ROOT, "docs", "convergence_2d.csv")),
]
OUT = os.path.join(ROOT, "docs", "convergence.svg")

# vivid enough to read on both light and dark backgrounds
SERIES = ["#4c9be8", "#e8743b", "#19a979", "#bf5af2", "#e8c33b"]
# github uses this muted grey in both themes
CHROME = "#7d8590"

W, H = 900, 460
PAD_L, PAD_R, PAD_T, PAD_B = 66, 20, 46, 78
GAP = 46
PANEL_W = (W - GAP) // 2


def ticks_125(lo, hi):
    """Tick positions at 1, 2 and 5 times each power of ten, as log10 values."""
    out = []
    for decade in range(math.floor(lo) - 1, math.ceil(hi) + 1):
        for mult in (1.0, 2.0, 5.0):
            lv = decade + math.log10(mult)
            if lo <= lv <= hi:
                out.append((lv, mult, decade))
    return out


def fmt_tick(mult, decade):
    if mult == 1.0:
        return f'10<tspan font-size="9" dy="-5">{decade}</tspan>'
    return f'{mult:g}&#183;10<tspan font-size="9" dy="-5">{decade}</tspan>'


def load(path):
    series = defaultdict(list)
    with open(path) as fh:
        for row in csv.DictReader(fh):
            series[int(row["p"])].append((float(row["h"]), float(row["error"])))
    for p in series:
        series[p].sort()
    return series


def fitted_rate(points):
    xs = [math.log(h) for h, _ in points]
    ys = [math.log(e) for _, e in points]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    return num / den


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def panel(ox, title, series, out):
    xs, ys = [], []
    for pts in series.values():
        xs += [math.log10(h) for h, _ in pts]
        ys += [math.log10(e) for _, e in pts]
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    # pad so markers and reference lines are not clipped
    px, py = 0.08 * (x1 - x0 or 1), 0.10 * (y1 - y0 or 1)
    x0, x1 = x0 - px, x1 + px
    y0, y1 = y0 - py, y1 + py

    pw = PANEL_W - PAD_L - PAD_R
    ph = H - PAD_T - PAD_B

    def sx(lx):
        return ox + PAD_L + (lx - x0) / (x1 - x0) * pw

    def sy(ly):
        return PAD_T + (y1 - ly) / (y1 - y0) * ph

    out.append(
        f'<text x="{ox + PAD_L}" y="{PAD_T - 22}" fill="{CHROME}" '
        f'font-size="15" font-weight="600">{esc(title)}</text>'
    )

    # x ticks at 1-2-5 per decade: h often spans well under a decade
    for lv, mult, decade in ticks_125(x0, x1):
        out.append(
            f'<line x1="{sx(lv):.1f}" y1="{PAD_T}" x2="{sx(lv):.1f}" '
            f'y2="{PAD_T + ph}" stroke="{CHROME}" stroke-opacity="0.18"/>'
        )
        out.append(
            f'<text x="{sx(lv):.1f}" y="{PAD_T + ph + 20}" fill="{CHROME}" '
            f'font-size="11.5" text-anchor="middle">{fmt_tick(mult, decade)}</text>'
        )
    # y spans many decades, so decade lines are the right density there
    for ly in range(math.floor(y0), math.ceil(y1) + 1):
        if y0 <= ly <= y1:
            out.append(
                f'<line x1="{ox + PAD_L}" y1="{sy(ly):.1f}" x2="{ox + PAD_L + pw}" '
                f'y2="{sy(ly):.1f}" stroke="{CHROME}" stroke-opacity="0.18"/>'
            )
            out.append(
                f'<text x="{ox + PAD_L - 9}" y="{sy(ly) + 4:.1f}" fill="{CHROME}" '
                f'font-size="11.5" text-anchor="end">10<tspan font-size="9" dy="-5">{ly}</tspan></text>'
            )

    # axis frame
    out.append(
        f'<rect x="{ox + PAD_L}" y="{PAD_T}" width="{pw}" height="{ph}" fill="none" '
        f'stroke="{CHROME}" stroke-opacity="0.45"/>'
    )
    out.append(
        f'<text x="{ox + PAD_L + pw / 2:.1f}" y="{PAD_T + ph + 42}" fill="{CHROME}" '
        f'font-size="13" text-anchor="middle">mesh size h</text>'
    )
    lab_x = ox + 18
    out.append(
        f'<text x="{lab_x}" y="{PAD_T + ph / 2:.1f}" fill="{CHROME}" font-size="13" '
        f'text-anchor="middle" transform="rotate(-90 {lab_x} {PAD_T + ph / 2:.1f})">'
        f'L2 error</text>'
    )

    legend = []
    for idx, p in enumerate(sorted(series)):
        pts = series[p]
        colour = SERIES[idx % len(SERIES)]
        rate = fitted_rate(pts)

        # theoretical slope p+1, anchored at the finest mesh
        hf, ef = pts[0]
        ax, ay = math.log10(hf), math.log10(ef)
        ref = [(x0, ay + (p + 1) * (x0 - ax)), (x1, ay + (p + 1) * (x1 - ax))]
        ref = [(x, y) for x, y in ref if y0 - 1 <= y <= y1 + 1]
        if len(ref) == 2:
            out.append(
                f'<line x1="{sx(ref[0][0]):.1f}" y1="{sy(ref[0][1]):.1f}" '
                f'x2="{sx(ref[1][0]):.1f}" y2="{sy(ref[1][1]):.1f}" stroke="{colour}" '
                f'stroke-opacity="0.38" stroke-width="1.2" stroke-dasharray="5 4"/>'
            )

        pl = " ".join(f"{sx(math.log10(h)):.1f},{sy(math.log10(e)):.1f}" for h, e in pts)
        out.append(
            f'<polyline points="{pl}" fill="none" stroke="{colour}" stroke-width="2.1" '
            f'stroke-linejoin="round"/>'
        )
        for h, e in pts:
            out.append(
                f'<circle cx="{sx(math.log10(h)):.1f}" cy="{sy(math.log10(e)):.1f}" '
                f'r="3.3" fill="{colour}"/>'
            )

        legend.append((p, colour, rate))

    # drawn last so the labels sit above every curve
    for idx, (p, colour, rate) in enumerate(legend):
        ly = PAD_T + ph - 13 - (len(legend) - 1 - idx) * 18
        lx = ox + PAD_L + pw - 176
        out.append(
            f'<line x1="{lx}" y1="{ly - 4}" x2="{lx + 20}" y2="{ly - 4}" '
            f'stroke="{colour}" stroke-width="2.1"/>'
        )
        out.append(
            f'<text x="{lx + 27}" y="{ly}" fill="{CHROME}" font-size="12">'
            f'p={p}  slope {rate:.2f} / {p + 1}</text>'
        )


def main():
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" '
        f'height="{H}" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">',
        f'<title>L2 error under h-refinement</title>',
    ]
    for i, (title, path) in enumerate(PANELS):
        panel(i * (PANEL_W + GAP), title, load(path), out)
    out.append(
        f'<text x="{W / 2}" y="{H - 14}" fill="{CHROME}" font-size="12" '
        f'text-anchor="middle" fill-opacity="0.9">Solid lines are measured; dashed lines '
        f'show the theoretical rate O(h^(p+1)). Legend gives measured slope / theory.</text>'
    )
    out.append("</svg>")
    with open(OUT, "w") as fh:
        fh.write("\n".join(out) + "\n")
    print(f"wrote {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
