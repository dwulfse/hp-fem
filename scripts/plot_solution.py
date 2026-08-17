#!/usr/bin/env python3
"""Render a 2D solution field as an SVG.

Reads field.csv, which the solver writes with --field: a list of mesh vertices
carrying the solution value, followed by the triangles joining them.

Each triangle is subdivided and its pieces are flat filled by the linearly
interpolated value, which approximates smooth shading closely enough at this
scale while keeping the output a plain list of polygons.

Standard library only.

    cd main && ./FEM -d 2 -p 2 -m L.1 --problem const --field
    python3 scripts/plot_solution.py main/field.csv docs/solution.svg "title"
"""

import os
import sys

CHROME = "#7d8590"
W, H = 720, 620
PAD = 54
BAR_W, BAR_H = 16, 300

# viridis, sampled at eight points and interpolated between
VIRIDIS = [
    (68, 1, 84), (72, 40, 120), (62, 74, 137), (49, 104, 142),
    (38, 130, 142), (31, 158, 137), (53, 183, 121), (109, 205, 89),
    (180, 222, 44), (253, 231, 37),
]


def colour(t):
    """Map t in [0, 1] to a viridis colour."""
    t = min(max(t, 0.0), 1.0)
    x = t * (len(VIRIDIS) - 1)
    i = min(int(x), len(VIRIDIS) - 2)
    f = x - i
    a, b = VIRIDIS[i], VIRIDIS[i + 1]
    return "#%02x%02x%02x" % tuple(round(a[k] + f * (b[k] - a[k])) for k in range(3))


def load(path):
    verts, tris = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(",")
            if parts[0] == "v":
                verts.append(tuple(float(v) for v in parts[1:4]))
            elif parts[0] == "t":
                tris.append(tuple(int(v) for v in parts[1:4]))
    return verts, tris


def subdivide(tri, levels):
    """Split a triangle of (x, y, u) corners into 4^levels smaller ones."""
    out = [tri]
    for _ in range(levels):
        nxt = []
        for (a, b, c) in out:
            ab = tuple((a[k] + b[k]) / 2 for k in range(3))
            bc = tuple((b[k] + c[k]) / 2 for k in range(3))
            ca = tuple((c[k] + a[k]) / 2 for k in range(3))
            nxt += [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]
        out = nxt
    return out


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "main/field.csv"
    dst = sys.argv[2] if len(sys.argv) > 2 else "docs/solution.svg"
    title = sys.argv[3] if len(sys.argv) > 3 else "Finite element solution"
    note = sys.argv[4] if len(sys.argv) > 4 else ""
    levels = int(sys.argv[5]) if len(sys.argv) > 5 else 0

    verts, tris = load(src)
    if not verts or not tris:
        sys.exit(f"no mesh found in {src}")

    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]
    us = [v[2] for v in verts]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    u0, u1 = min(us), max(us)
    span = (u1 - u0) or 1.0

    # equal aspect, leaving room for the colour bar on the right
    plot_w = W - 2 * PAD - BAR_W - 46
    plot_h = H - 2 * PAD
    scale = min(plot_w / (x1 - x0), plot_h / (y1 - y0))
    ox = PAD + (plot_w - scale * (x1 - x0)) / 2
    oy = PAD + (plot_h - scale * (y1 - y0)) / 2

    def sx(x):
        return ox + (x - x0) * scale

    def sy(y):
        return oy + (y1 - y) * scale          # svg y runs downwards

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" '
        f'height="{H}" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">',
        f'<title>{title}</title>',
        f'<text x="{PAD}" y="{PAD - 20}" fill="{CHROME}" font-size="15" '
        f'font-weight="600">{title}</text>',
    ]

    # Triangles are bucketed by colour and emitted as one path per bucket. Each
    # triangle then costs only its own coordinates, rather than repeating a fill
    # attribute and element tags, which cuts the file by roughly four times.
    # Coordinates are rounded to whole pixels; adjacent triangles share vertex
    # values exactly, so they round identically and no seams open up.
    BUCKETS = 128
    paths = {}

    for (i, j, k) in tris:
        for piece in subdivide((verts[i], verts[j], verts[k]), levels):
            mean = sum(p[2] for p in piece) / 3.0
            bucket = min(int((mean - u0) / span * BUCKETS), BUCKETS - 1)
            a, b, c = ((round(sx(p[0])), round(sy(p[1]))) for p in piece)
            paths.setdefault(bucket, []).append(
                f"M{a[0]},{a[1]}L{b[0]},{b[1]}L{c[0]},{c[1]}Z"
            )

    # stroking each path in its own fill colour closes the hairline gaps that
    # antialiasing otherwise leaves between abutting triangles; because the
    # stroke sits on the path rather than on every triangle it is nearly free
    for bucket in sorted(paths):
        fill = colour((bucket + 0.5) / BUCKETS)
        out.append(
            f'<path d="{"".join(paths[bucket])}" fill="{fill}" stroke="{fill}" '
            f'stroke-width="1" stroke-linejoin="round"/>'
        )

    # colour bar
    bx = W - PAD - BAR_W
    by = PAD + (plot_h - BAR_H) / 2
    out.append('<defs><linearGradient id="bar" x1="0" y1="1" x2="0" y2="0">')
    for s in range(11):
        out.append(f'<stop offset="{s / 10:.1f}" stop-color="{colour(s / 10)}"/>')
    out.append("</linearGradient></defs>")
    out.append(
        f'<rect x="{bx}" y="{by:.1f}" width="{BAR_W}" height="{BAR_H}" fill="url(#bar)" '
        f'stroke="{CHROME}" stroke-opacity="0.45"/>'
    )
    for s in range(5):
        v = u0 + span * s / 4
        ty = by + BAR_H * (1 - s / 4)
        out.append(
            f'<text x="{bx - 7}" y="{ty + 4:.1f}" fill="{CHROME}" font-size="11.5" '
            f'text-anchor="end">{v:.3g}</text>'
        )
    out.append(
        f'<text x="{bx + BAR_W / 2}" y="{by - 12:.1f}" fill="{CHROME}" font-size="12" '
        f'text-anchor="middle">u</text>'
    )

    out.append(
        f'<text x="{PAD}" y="{H - 18}" fill="{CHROME}" font-size="11.5" '
        f'fill-opacity="0.9">{note}</text>'
    )
    out.append("</svg>")

    os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
    with open(dst, "w") as fh:
        fh.write("\n".join(out) + "\n")
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()
