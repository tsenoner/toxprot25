#!/usr/bin/env python3
"""ToxProt Definition Comparison — Venn→phylum→order flow figure.

Single-panel figure flowing from a Venn diagram (left) through phyla bars
(centre) to individual order dots (right): general to specific.

Usage:
    python -m src.analysis.combine_definition_panels
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse, PathPatch, Rectangle
from matplotlib.path import Path as MplPath

from .colors import DEFINITION_COLORS

DEFINITION_ORDER = ["both", "venom_tissue", "kw_toxin"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _crosstab(df: pd.DataFrame) -> pd.DataFrame:
    """Phylum x Definition cross-tabulation sorted by total descending."""
    clean = df.dropna(subset=["Phylum", "ToxProt definition"])
    ct = pd.crosstab(clean["Phylum"], clean["ToxProt definition"])
    for d in DEFINITION_ORDER:
        if d not in ct.columns:
            ct[d] = 0
    ct = ct[DEFINITION_ORDER]
    ct["total"] = ct.sum(axis=1)
    return ct.sort_values("total", ascending=False)


def _order_data(df: pd.DataFrame) -> pd.DataFrame:
    """Per-order entry count, phylum, and exclusivity classification.

    Returns a DataFrame with columns: Order, Phylum, count, exclusivity.
    Exclusivity is 'kw_toxin' (keyword-exclusive) or 'both' (shared).
    """
    from .analyze_definitions import get_criterion_subset

    clean = df.dropna(subset=["Order", "Phylum", "ToxProt definition"])

    # Entry count per order
    order_counts = clean.groupby("Order").size().rename("count")

    # Phylum for each order (mode — each order belongs to one phylum)
    order_phylum = clean.groupby("Order")["Phylum"].first()

    # Classify exclusivity: orders in venom tissue set vs keyword-only
    venom_orders = set(get_criterion_subset(clean, "venom_tissue")["Order"].dropna().unique())
    keyword_orders = set(get_criterion_subset(clean, "kw_toxin")["Order"].dropna().unique())
    shared_orders = venom_orders & keyword_orders

    orders = pd.DataFrame({"count": order_counts, "Phylum": order_phylum})
    orders["exclusivity"] = orders.index.map(lambda o: "both" if o in shared_orders else "kw_toxin")
    orders = orders.sort_values("count", ascending=False).reset_index()
    orders.rename(columns={"index": "Order"}, inplace=True)
    return orders


# ---------------------------------------------------------------------------
# Main plot
# ---------------------------------------------------------------------------


def plot_grouped_phylum_venn_flow(df: pd.DataFrame, ax: plt.Axes) -> None:
    """Draw Venn -> phylum -> order flow diagram on *ax*."""
    from .analyze_definitions import get_taxa_by_exclusivity

    ct = _crosstab(df)

    # --- classify phyla by exclusivity ------------------------------------
    phylum_excl = get_taxa_by_exclusivity(df, "Phylum")
    kw_exclusive_phyla = phylum_excl["Keyword"]  # set of phylum names
    shared_phyla = phylum_excl["Shared"]

    # Split phyla from crosstab into two groups, each sorted by total desc
    kw_group = [p for p in ct.index if p in kw_exclusive_phyla]
    shared_group = [p for p in ct.index if p in shared_phyla]
    # Include venom-exclusive phyla in shared group if any exist
    venom_exclusive = phylum_excl.get("Venom", set())
    venom_group = [p for p in ct.index if p in venom_exclusive]
    shared_group = shared_group + venom_group

    # Combined phyla list: kw-exclusive first (bottom y), shared after (top y)
    # Since y increases upward in matplotlib, items later in the list are higher.
    phyla = kw_group[::-1] + shared_group[::-1]

    # --- layout constants -------------------------------------------------
    box_w = 0.7
    x_venn = -2.5
    x_left = 3.0
    x_orders = 5.5
    gap_within = 0.4
    gap_between = 2.5
    min_bar = 0.5

    total_entries = ct.loc[phyla, "total"].sum()
    scale = 28.0 / total_entries

    # --- left: phylum positions (kw-exclusive at bottom, shared at top) ---
    n_kw_phyla = len(kw_group)
    left_pos = {}
    y = 0.0
    for i, phylum in enumerate(phyla):
        # Insert group gap between kw-exclusive and shared groups
        if i == n_kw_phyla and n_kw_phyla > 0:
            y += gap_between - gap_within
        cnt = ct.loc[phylum, "total"]
        h = max(cnt * scale, min_bar)
        left_pos[phylum] = {"bottom": y, "top": y + h, "count": cnt}
        y += h + gap_within
    total_left_h = y - gap_within

    # --- order strip plot (left column) -----------------------------------
    orders = _order_data(df)
    n_orders = len(orders)
    max_count = orders["count"].max()
    log_max = np.log10(max_count)

    # Dot size proportional to log(count)
    s_min, s_max = 60, 600
    orders["dot_s"] = s_min + (s_max - s_min) * np.log10(orders["count"]) / log_max

    # Estimate display→data conversion for separate x/y radii.
    # Scatter circles are round in display space but elliptical in data space.
    fig_w_in, fig_h_in = 16, 9
    ax_w_in = fig_w_in * 0.90
    ax_h_in = fig_h_in * 0.88
    y_range_est = total_left_h + 8.0
    x_range_pre = (x_orders + 3.0) - (x_venn - 6)
    pts_per_x = 72 * ax_w_in / x_range_pre
    pts_per_y = 72 * ax_h_in / y_range_est

    r_pts = np.sqrt(orders["dot_s"].values / np.pi)
    radii_x = r_pts / pts_per_x
    radii_y = r_pts / pts_per_y

    orders["dot_x"] = 0.0
    orders["dot_y"] = 0.0
    orders["arrival_y"] = 0.0
    pad = 0.10
    gap_pts = 3  # visual gap in points (uniform on screen)
    spacing_x = gap_pts / pts_per_x
    spacing_y = gap_pts / pts_per_y

    dot_x_max_global = x_orders  # track rightmost dot edge for xlim

    for phylum in left_pos:
        pos = left_pos[phylum]
        mask = orders["Phylum"] == phylum
        if not mask.any():
            continue
        ph_orders = orders.loc[mask].sort_values("count", ascending=False)
        idxs = ph_orders.index.tolist()
        n = len(idxs)
        dot_rx = [radii_x[i] for i in idxs]
        dot_ry = [radii_y[i] for i in idxs]

        y_top = pos["top"] - pad
        y_bot = pos["bottom"] + pad
        if y_top - y_bot < 0.1:
            y_top, y_bot = pos["top"], pos["bottom"]
        available = y_top - y_bot

        if n == 1:
            orders.at[idxs[0], "dot_x"] = x_orders
            orders.at[idxs[0], "dot_y"] = (y_top + y_bot) / 2
            orders.at[idxs[0], "arrival_y"] = (y_top + y_bot) / 2
            dot_x_max_global = max(dot_x_max_global, x_orders + dot_rx[0])
            continue

        # Start narrow (single-column), widen only when vertical space runs out
        strip_w = max(dot_rx) + 0.01
        max_strip_w = 6.0

        while True:
            rows = []
            cur_row = []
            row_x_cursor = 0.0
            row_max_ry = 0.0
            y_cursor = y_top - dot_ry[0]

            for k in range(n):
                rx = dot_rx[k]
                ry = dot_ry[k]
                if not cur_row:
                    cur_row.append((k, 0.0, y_cursor))
                    row_x_cursor = rx
                    row_max_ry = ry
                else:
                    new_x = row_x_cursor + spacing_x + rx
                    if new_x + rx <= strip_w:
                        cur_row.append((k, new_x, y_cursor))
                        row_x_cursor = new_x + rx
                        row_max_ry = max(row_max_ry, ry)
                    else:
                        rows.append((cur_row, row_max_ry))
                        y_cursor -= row_max_ry + ry + spacing_y
                        cur_row = [(k, 0.0, y_cursor)]
                        row_x_cursor = rx
                        row_max_ry = ry

            if cur_row:
                rows.append((cur_row, row_max_ry))

            # Check vertical fit
            last_row_dots, last_row_max_ry = rows[-1]
            lowest_y = last_row_dots[-1][2] - last_row_max_ry
            total_height = y_top - lowest_y

            if total_height <= available or strip_w >= max_strip_w:
                # Centre each row horizontally around x_orders
                packed = []
                for row_dots, _ in rows:
                    if len(row_dots) == 1:
                        k, _, cy = row_dots[0]
                        packed.append((k, x_orders, cy))
                    else:
                        k_last = row_dots[-1][0]
                        row_right = row_dots[-1][1] + dot_rx[k_last]
                        k_first = row_dots[0][0]
                        row_left_edge = row_dots[0][1] - dot_rx[k_first]
                        row_span = row_right - row_left_edge
                        offset = x_orders - row_span / 2 - row_left_edge
                        for k, local_x, cy in row_dots:
                            packed.append((k, local_x + offset, cy))

                # Vertically centre the block within the band
                actual_top = max(cy + dot_ry[k] for k, _, cy in packed)
                actual_bot = min(cy - dot_ry[k] for k, _, cy in packed)
                block_h = actual_top - actual_bot
                band_mid = (y_top + y_bot) / 2
                y_shift = band_mid - (actual_bot + block_h / 2)
                packed = [(k, cx, cy + y_shift) for k, cx, cy in packed]

                for k, cx, _ in packed:
                    dot_x_max_global = max(dot_x_max_global, cx + dot_rx[k])
                break

            strip_w = min(strip_w * 1.5, max_strip_w)

        # Stagger y within multi-dot rows so each dot gets its own level
        row_groups = []
        cur_group = [packed[0]]
        for item in packed[1:]:
            if abs(item[2] - cur_group[0][2]) < 1e-9:
                cur_group.append(item)
            else:
                row_groups.append(cur_group)
                cur_group = [item]
        row_groups.append(cur_group)

        row_centers = [g[0][2] for g in row_groups]
        staggered = []
        for ri, group in enumerate(row_groups):
            if len(group) == 1:
                staggered.extend(group)
                continue
            n_in_row = len(group)
            # Upper boundary: midpoint to row above, or band top
            if ri > 0:
                y_hi = (row_centers[ri] + row_centers[ri - 1]) / 2
            else:
                y_hi = y_top
            # Lower boundary: midpoint to row below, or band bottom
            if ri < len(row_groups) - 1:
                y_lo = (row_centers[ri] + row_centers[ri + 1]) / 2
            else:
                y_lo = y_bot
            # Place dot centers from top to bottom within bounds
            rys_in = [dot_ry[k] for k, _, _ in group]
            y_first = y_hi - rys_in[0]
            y_last = y_lo + rys_in[-1]
            if y_first < y_last:
                y_first = y_last = row_centers[ri]
            ys = np.linspace(y_first, y_last, n_in_row)
            for i, (k, cx, _) in enumerate(group):
                staggered.append((k, cx, ys[i]))
        packed = staggered

        # Assign packed x-positions (y will be overridden below)
        for k, cx, _cy in packed:
            idx = idxs[k]
            orders.at[idx, "dot_x"] = cx

    # Equal y-spacing: within each phylum bar, largest dot on top
    for phylum in left_pos:
        pos = left_pos[phylum]
        mask = orders["Phylum"] == phylum
        if not mask.any():
            continue
        ph_idx = orders.loc[mask].sort_values("count", ascending=False).index.tolist()
        n = len(ph_idx)
        bar_pad = (pos["top"] - pos["bottom"]) * 0.10
        y_top = pos["top"] - bar_pad
        y_bot = pos["bottom"] + bar_pad
        if y_top - y_bot < 0.1:
            y_top, y_bot = pos["top"], pos["bottom"]

        ys = np.linspace(y_top, y_bot, n) if n > 1 else [(y_top + y_bot) / 2]

        for k, idx in enumerate(ph_idx):
            orders.at[idx, "dot_y"] = ys[k]
            orders.at[idx, "arrival_y"] = ys[k]

    # X-axis repulsion: shift overlapping dots to the left
    dot_x = orders["dot_x"].values.copy()
    dot_y = orders["dot_y"].values.copy()

    for _iteration in range(200):
        moved = False
        for i in range(n_orders):
            for j in range(i + 1, n_orders):
                dx_pts = (dot_x[i] - dot_x[j]) * pts_per_x
                dy_pts = (dot_y[i] - dot_y[j]) * pts_per_y
                dist_pts = np.hypot(dx_pts, dy_pts)
                min_dist_pts = r_pts[i] + r_pts[j]  # just touch, no extra gap

                if dist_pts < min_dist_pts:
                    # Need more x separation (move the smaller dot left)
                    needed_dx_pts = np.sqrt(max(min_dist_pts**2 - dy_pts**2, 0))
                    current_dx = abs(dx_pts)
                    if current_dx < needed_dx_pts:
                        push = (needed_dx_pts - current_dx) / pts_per_x
                        # Move the smaller dot further right
                        if r_pts[i] <= r_pts[j]:
                            dot_x[i] += push
                        else:
                            dot_x[j] += push
                        moved = True
        if not moved:
            break

    orders["dot_x"] = dot_x
    dot_x_max_global = max(dot_x_max_global, (dot_x + radii_x).max())

    # Draw phylum→order bezier flows
    for _, row in orders.iterrows():
        if row["Phylum"] not in left_pos:
            continue
        x0 = x_left + box_w / 2  # right edge of phylum bar
        y0 = row["arrival_y"]
        x1 = row["dot_x"]  # order dot (on the right)
        y1 = row["dot_y"]
        lw = 0.5 + 0.7 * np.log10(row["count"]) / log_max
        color = DEFINITION_COLORS[row["exclusivity"]]

        ctrl_dx = (x1 - x0) * 0.45
        verts = [
            (x0, y0),
            (x0 + ctrl_dx, y0),
            (x1 - ctrl_dx, y1),
            (x1, y1),
        ]
        codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
        ax.add_patch(
            PathPatch(
                MplPath(verts, codes),
                fc="none",
                ec=color,
                lw=lw,
                alpha=0.3,
            )
        )

    # Draw order dots (scatter, size proportional to log entry count)
    for excl_key in ["kw_toxin", "both"]:
        subset = orders[orders["exclusivity"] == excl_key]
        ax.scatter(
            subset["dot_x"],
            subset["dot_y"],
            s=subset["dot_s"],
            c=DEFINITION_COLORS[excl_key],
            edgecolors="black",
            linewidths=0.4,
            alpha=0.85,
            zorder=5,
        )

    # Count labels on dots
    for _, row in orders.iterrows():
        ax.text(
            row["dot_x"],
            row["dot_y"],
            f"{int(row['count']):,}",
            ha="center",
            va="center",
            fontsize=5 + 2 * np.log10(row["count"]) / log_max,
            fontweight="bold",
            color="white",
            zorder=6,
        )

    # Column header — color-coded counts
    n_kw_orders = int((orders["exclusivity"] == "kw_toxin").sum())
    n_shared_orders = int((orders["exclusivity"] == "both").sum())
    header_y = total_left_h + 0.8
    ax.text(
        x_orders,
        header_y + 2.0,
        f"Orders (n={n_orders})",
        ha="center",
        va="bottom",
        fontsize=11,
        fontweight="bold",
        color="black",
    )
    ax.text(
        x_orders,
        header_y + 0.8,
        f"{n_shared_orders} both",
        ha="center",
        va="bottom",
        fontsize=9,
        fontweight="bold",
        color=DEFINITION_COLORS["both"],
    )
    ax.text(
        x_orders,
        header_y,
        f"{n_kw_orders} keyword",
        ha="center",
        va="bottom",
        fontsize=9,
        fontweight="bold",
        color=DEFINITION_COLORS["kw_toxin"],
    )

    # --- right: Venn geometry ---------------------------------------------
    n_kw = int(ct.loc[phyla, "kw_toxin"].sum())
    n_both = int(ct.loc[phyla, "both"].sum())
    n_vt = int(ct.loc[phyla, "venom_tissue"].sum())

    h_kw = max(n_kw * scale, min_bar)
    h_both = max(n_both * scale, min_bar)
    h_vt = max(n_vt * scale, min_bar)

    total_right_h = h_kw + h_both + h_vt + 2 * gap_within
    right_start = (total_left_h - total_right_h) / 2

    zone_kw = {"bottom": right_start, "top": right_start + h_kw, "count": n_kw}
    zone_both = {
        "bottom": zone_kw["top"] + gap_within,
        "top": zone_kw["top"] + gap_within + h_both,
        "count": n_both,
    }
    zone_vt = {
        "bottom": zone_both["top"] + gap_within,
        "top": zone_both["top"] + gap_within + h_vt,
        "count": n_vt,
    }
    right_zones = {
        "both": zone_both,
        "venom_tissue": zone_vt,
        "kw_toxin": zone_kw,
    }

    # --- Ellipse geometry -------------------------------------------------
    kw_lo = zone_kw["bottom"]
    kw_hi = zone_both["top"]
    vt_lo = zone_both["bottom"]
    vt_hi = zone_vt["top"]

    ry_vt = (vt_hi - vt_lo) / 2
    cy_vt = (vt_hi + vt_lo) / 2
    ry_kw = (kw_hi - kw_lo) / 2
    cy_kw = (kw_hi + kw_lo) / 2

    # Aspect correction for true circles
    fig = ax.get_figure()
    fig_w, fig_h = fig.get_size_inches()
    bbox = ax.get_position()
    ax_w = fig_w * bbox.width
    ax_h = fig_h * bbox.height
    margin_lbl = 2.5
    dot_y_min = orders["dot_y"].min() - radii_y.max()
    dot_y_max = orders["dot_y"].max() + radii_y.max()
    y_lo = min(zone_kw["bottom"], 0, dot_y_min) - margin_lbl
    y_hi = max(zone_vt["top"], total_left_h, dot_y_max) + margin_lbl
    y_range = y_hi - y_lo
    x_range_est = (dot_x_max_global + 1.0) - (x_venn - 5)
    aspect_corr = (ax_h * x_range_est) / (ax_w * y_range)
    rx_vt = ry_vt * aspect_corr
    rx_kw = ry_kw * aspect_corr

    flow_colors = {
        "both": DEFINITION_COLORS["both"],
        "kw_toxin": DEFINITION_COLORS["kw_toxin"],
        "venom_tissue": DEFINITION_COLORS["venom_tissue"],
    }

    # --- draw Venn circles (6-patch technique) ----------------------------
    # Filled ellipses (no borders — borders drawn separately on top)
    ell_vt = Ellipse(
        (x_venn, cy_vt),
        width=rx_vt * 2,
        height=ry_vt * 2,
        facecolor=DEFINITION_COLORS["venom_tissue"],
        edgecolor="none",
        alpha=0.15,
    )
    ax.add_patch(ell_vt)

    ell_kw = Ellipse(
        (x_venn, cy_kw),
        width=rx_kw * 2,
        height=ry_kw * 2,
        facecolor=DEFINITION_COLORS["kw_toxin"],
        edgecolor="none",
        alpha=0.15,
    )
    ax.add_patch(ell_kw)

    # White mask to erase blue+orange blend in intersection
    mask = Ellipse(
        (x_venn, cy_vt),
        width=rx_vt * 2,
        height=ry_vt * 2,
        facecolor="white",
        edgecolor="none",
        alpha=1.0,
    )
    ax.add_patch(mask)
    mask.set_clip_path(
        Ellipse(
            (x_venn, cy_kw),
            width=rx_kw * 2,
            height=ry_kw * 2,
            transform=ax.transData,
        )
    )

    # Clean purple intersection fill
    inter = Ellipse(
        (x_venn, cy_vt),
        width=rx_vt * 2,
        height=ry_vt * 2,
        facecolor=DEFINITION_COLORS["both"],
        edgecolor="none",
        alpha=0.20,
    )
    ax.add_patch(inter)
    inter.set_clip_path(
        Ellipse(
            (x_venn, cy_kw),
            width=rx_kw * 2,
            height=ry_kw * 2,
            transform=ax.transData,
        )
    )

    # --- Borders (drawn after fills so they appear on top) ----------------
    # KW border: subtle orange outside intersection, subtle purple inside
    kw_border_full = Ellipse(
        (x_venn, cy_kw),
        width=rx_kw * 2,
        height=ry_kw * 2,
        facecolor="none",
        edgecolor=DEFINITION_COLORS["kw_toxin"],
        linewidth=1.5,
        alpha=0.3,
    )
    ax.add_patch(kw_border_full)

    # Erase orange in intersection, then redraw as purple
    kw_border_erase = Ellipse(
        (x_venn, cy_kw),
        width=rx_kw * 2,
        height=ry_kw * 2,
        facecolor="none",
        edgecolor="white",
        linewidth=1.5,
        alpha=1.0,
    )
    ax.add_patch(kw_border_erase)
    kw_border_erase.set_clip_path(
        Ellipse(
            (x_venn, cy_vt),
            width=rx_vt * 2,
            height=ry_vt * 2,
            transform=ax.transData,
        )
    )

    kw_border_inter = Ellipse(
        (x_venn, cy_kw),
        width=rx_kw * 2,
        height=ry_kw * 2,
        facecolor="none",
        edgecolor=DEFINITION_COLORS["both"],
        linewidth=1.5,
        alpha=0.3,
    )
    ax.add_patch(kw_border_inter)
    kw_border_inter.set_clip_path(
        Ellipse(
            (x_venn, cy_vt),
            width=rx_vt * 2,
            height=ry_vt * 2,
            transform=ax.transData,
        )
    )

    # Highlighted VT border: blue outside intersection, purple inside
    # (drawn last so it's on top of everything)
    vt_border_full = Ellipse(
        (x_venn, cy_vt),
        width=rx_vt * 2,
        height=ry_vt * 2,
        facecolor="none",
        edgecolor=DEFINITION_COLORS["venom_tissue"],
        linewidth=2.5,
        alpha=0.7,
    )
    ax.add_patch(vt_border_full)

    vt_border_inter = Ellipse(
        (x_venn, cy_vt),
        width=rx_vt * 2,
        height=ry_vt * 2,
        facecolor="none",
        edgecolor=DEFINITION_COLORS["both"],
        linewidth=2.5,
        alpha=0.7,
    )
    ax.add_patch(vt_border_inter)
    vt_border_inter.set_clip_path(
        Ellipse(
            (x_venn, cy_kw),
            width=rx_kw * 2,
            height=ry_kw * 2,
            transform=ax.transData,
        )
    )

    # --- Ellipse right-edge helpers ---------------------------------------
    def ellipse_right_x(y, cy, ry, r_x):
        t = (y - cy) / ry
        t = np.clip(t, -1.0, 1.0)
        return x_venn + r_x * np.sqrt(1.0 - t * t)

    def flow_target_x(defn, y_mid):
        if defn == "kw_toxin":
            return ellipse_right_x(y_mid, cy_kw, ry_kw, rx_kw)
        elif defn == "venom_tissue":
            return ellipse_right_x(y_mid, cy_vt, ry_vt, rx_vt)
        else:  # both — follow outermost individual circle outline
            return max(
                ellipse_right_x(y_mid, cy_kw, ry_kw, rx_kw),
                ellipse_right_x(y_mid, cy_vt, ry_vt, rx_vt),
            )

    # --- Bezier flow helper (Venn right edge → phylum bar left edge) ------
    n_arc = 30
    flow_gap = 0.05  # gap between flow start and ellipse boundary

    def draw_flow(x1, y1b, y1t, y2b, y2t, defn, alpha=0.35):
        color = flow_colors[defn]
        # Trace ellipse right-edge curvature offset by gap
        ys_up = np.linspace(y2b, y2t, n_arc)
        arc_up = [(flow_target_x(defn, y) + flow_gap, y) for y in ys_up]
        x2b, _ = arc_up[0]
        x2t, _ = arc_up[-1]
        x2_mid = (x2b + x2t) / 2
        bar_left = x1 - box_w / 2
        ctrl = (bar_left - x2_mid) * 0.45

        # Bottom bezier: Venn right edge → bar left edge
        verts = [
            (x2b, y2b),
            (x2b + ctrl, y2b),
            (bar_left - ctrl, y1b),
            (bar_left, y1b),
        ]
        codes = [
            MplPath.MOVETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
        ]
        # Line up bar left edge
        verts.append((bar_left, y1t))
        codes.append(MplPath.LINETO)
        # Top bezier: bar left edge → Venn right edge
        verts += [
            (bar_left - ctrl, y1t),
            (x2t + ctrl, y2t),
            (x2t, y2t),
        ]
        codes += [MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
        # Arc back down along Venn ellipse
        for pt in reversed(arc_up[:-1]):
            verts.append(pt)
            codes.append(MplPath.LINETO)
        verts.append((x2b, y2b))
        codes.append(MplPath.CLOSEPOLY)

        ax.add_patch(PathPatch(MplPath(verts, codes), fc=color, ec="none", alpha=alpha))

    # --- draw flows -------------------------------------------------------
    stack_order = ["kw_toxin", "both", "venom_tissue"]

    left_off = {p: left_pos[p]["bottom"] for p in phyla}
    right_off = {d: right_zones[d]["bottom"] for d in stack_order}

    for phylum in phyla:
        pos = left_pos[phylum]
        bar_h = pos["top"] - pos["bottom"]
        tot = pos["count"]
        for defn in stack_order:
            cnt = ct.loc[phylum, defn]
            if cnt == 0:
                continue
            fh_left = bar_h * cnt / tot
            fh_right = cnt * scale

            y2b = right_off[defn]
            y2t = y2b + fh_right

            draw_flow(
                x_left,
                left_off[phylum],
                left_off[phylum] + fh_left,
                y2b,
                y2t,
                defn,
            )
            left_off[phylum] += fh_left
            right_off[defn] += fh_right

    # --- draw left phylum bars --------------------------------------------
    small_thr = 1.5
    for phylum in phyla:
        pos = left_pos[phylum]
        bar_h = pos["top"] - pos["bottom"]
        tot = pos["count"]

        active = [d for d in DEFINITION_ORDER if ct.loc[phylum, d] > 0]
        lbl_color = flow_colors[active[0]] if len(active) == 1 else "black"

        # stacked coloured segments
        seg_y = pos["bottom"]
        for defn in stack_order:
            sc = ct.loc[phylum, defn]
            if sc == 0:
                continue
            sh = bar_h * sc / tot
            ax.add_patch(
                Rectangle(
                    (x_left - box_w / 2, seg_y),
                    box_w,
                    sh,
                    fc=flow_colors[defn],
                    ec="none",
                )
            )
            seg_y += sh

        # border
        ax.add_patch(
            Rectangle(
                (x_left - box_w / 2, pos["bottom"]),
                box_w,
                bar_h,
                fc="none",
                ec="black",
                lw=1,
            )
        )

        # labels
        fs_name = 10 if bar_h > small_thr else 9
        fs_count = 9 if bar_h > small_thr else 7
        ax.text(
            x_left + box_w / 2 + 0.1,
            pos["bottom"] + bar_h / 2,
            phylum,
            ha="left",
            va="center",
            fontsize=fs_name,
            fontstyle="italic",
            color=lbl_color,
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=1),
            zorder=4,
        )
        ax.text(
            x_left,
            pos["bottom"] + bar_h / 2,
            f"{tot:,}",
            ha="center",
            va="center",
            fontsize=fs_count,
            fontweight="bold",
            color="white",
        )

    # --- group header annotations -----------------------------------------
    phylum_header_y = total_left_h + 0.8
    n_phyla_total = len(kw_group) + len(shared_group)
    ax.text(
        x_left,
        phylum_header_y + 2.0,
        f"Phyla (n={n_phyla_total})",
        ha="center",
        va="bottom",
        fontsize=11,
        fontweight="bold",
        color="black",
    )
    ax.text(
        x_left,
        phylum_header_y + 0.8,
        f"{len(shared_group)} both",
        ha="center",
        va="bottom",
        fontsize=9,
        fontweight="bold",
        color=DEFINITION_COLORS["both"],
    )
    ax.text(
        x_left,
        phylum_header_y,
        f"{len(kw_group)} keyword",
        ha="center",
        va="bottom",
        fontsize=9,
        fontweight="bold",
        color=DEFINITION_COLORS["kw_toxin"],
    )

    # --- Venn region labels -----------------------------------------------
    for defn, zone in right_zones.items():
        cy = (zone["bottom"] + zone["top"]) / 2
        ax.text(
            x_venn,
            cy,
            f"{zone['count']:,}",
            ha="center",
            va="center",
            fontsize=11,
            fontweight="bold",
            color=flow_colors[defn],
        )

    ax.text(
        x_venn,
        zone_vt["top"] + 2.6,
        "Venom Tissue",
        ha="center",
        va="bottom",
        fontsize=14,
        fontweight="bold",
        color=flow_colors["venom_tissue"],
    )
    ax.text(
        x_venn,
        zone_vt["top"] + 1.2,
        "cc_tissue_specificity:venom",
        ha="center",
        va="bottom",
        fontsize=12,
        fontstyle="italic",
        color=flow_colors["venom_tissue"],
    )
    ax.text(
        x_venn,
        zone_kw["bottom"] - 1.2,
        "Toxin Keyword",
        ha="center",
        va="top",
        fontsize=14,
        fontweight="bold",
        color=flow_colors["kw_toxin"],
    )
    ax.text(
        x_venn,
        zone_kw["bottom"] - 2.6,
        "keyword:KW-0800",
        ha="center",
        va="top",
        fontsize=12,
        fontstyle="italic",
        color=flow_colors["kw_toxin"],
    )

    # --- legend -----------------------------------------------------------
    legend_elements = [
        plt.Rectangle(
            (0, 0),
            1,
            1,
            fc=DEFINITION_COLORS["kw_toxin"],
            ec="black",
            label="Toxin Keyword",
        ),
        plt.Rectangle(
            (0, 0),
            1,
            1,
            fc=DEFINITION_COLORS["venom_tissue"],
            ec="black",
            label="Venom Tissue",
        ),
        plt.Rectangle(
            (0, 0),
            1,
            1,
            fc=DEFINITION_COLORS["both"],
            ec="black",
            label="Both",
        ),
    ]
    ax.legend(
        handles=legend_elements,
        loc="upper left",
        fontsize=12,
        framealpha=0.9,
    )

    # --- axis limits ------------------------------------------------------
    rx_max = max(rx_vt, rx_kw)
    ax.set_xlim(x_venn - rx_max - 1.0, dot_x_max_global + 0.3)
    ax.set_ylim(y_lo, y_hi)
    ax.axis("off")


# ---------------------------------------------------------------------------
# Figure creation
# ---------------------------------------------------------------------------


def create_combined_definition_figure(
    df: pd.DataFrame,
    output_dir: Path,
) -> Path:
    """Create Venn->phylum->order definition comparison figure."""
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(16, 9))
    fig.subplots_adjust(left=0.05, right=0.95, top=0.93, bottom=0.05)
    plot_grouped_phylum_venn_flow(df, ax)

    out = output_dir / "definition_comparison.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    return out


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data",
        type=Path,
        default=Path("data/processed/toxprot/toxprot_2025.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("figures/definitions"),
    )
    args = parser.parse_args()

    df = pd.read_csv(args.data)
    print(f"Loaded {len(df):,} entries")
    path = create_combined_definition_figure(df, args.output_dir)
    print(f"Saved {path}")
