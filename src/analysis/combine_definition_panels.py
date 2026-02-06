#!/usr/bin/env python3
"""Phylum → Definition Venn flow diagram.

Single unified visualization: alluvial flows from phyla on the left
into a vertically-oriented Venn diagram on the right.

  - Bottom circle  = Venom Tissue definition
  - Top circle     = Toxin Keyword definition
  - Overlap        = Both definitions

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
    """Phylum × Definition cross-tabulation sorted by total descending."""
    clean = df.dropna(subset=["Phylum", "ToxProt definition"])
    ct = pd.crosstab(clean["Phylum"], clean["ToxProt definition"])
    for d in DEFINITION_ORDER:
        if d not in ct.columns:
            ct[d] = 0
    ct = ct[DEFINITION_ORDER]
    ct["total"] = ct.sum(axis=1)
    return ct.sort_values("total", ascending=False)


# ---------------------------------------------------------------------------
# Main plot
# ---------------------------------------------------------------------------


def plot_phylum_venn_flow(df: pd.DataFrame, ax: plt.Axes) -> None:
    """Draw phylum → Venn flow diagram on *ax*."""
    ct = _crosstab(df)
    phyla = ct.index.tolist()

    # --- layout constants ---------------------------------------------------
    box_w = 0.7  # phylum bar width
    x_left = 0.0  # phylum column x
    x_venn = 6.0  # Venn centre x
    gap = 0.4  # vertical gap between phylum bars

    total_entries = ct["total"].sum()
    target_h = 30.0  # total visual height for bars
    scale = target_h / total_entries
    min_bar = 0.4

    # --- left: phylum positions ---------------------------------------------
    left_pos = {}
    y = 0.0
    for phylum in phyla:
        cnt = ct.loc[phylum, "total"]
        h = max(cnt * scale, min_bar)
        left_pos[phylum] = {"bottom": y, "top": y + h, "count": cnt}
        y += h + gap
    total_left_h = y - gap

    # --- right: Venn geometry -----------------------------------------------
    # Two overlapping circles stacked vertically.
    # Region heights (scaled) determine circle sizing / placement.
    n_kw = int(ct["kw_toxin"].sum())  # kw_toxin-only
    n_both = int(ct["both"].sum())  # both
    n_vt = int(ct["venom_tissue"].sum())  # venom_tissue-only

    h_kw = n_kw * scale
    h_both = n_both * scale
    h_vt = n_vt * scale

    # Stack the three zones vertically, centred on the left phyla stack.
    total_right_h = h_kw + h_both + h_vt + 2 * gap
    right_start = (total_left_h - total_right_h) / 2

    # Zone positions (bottom, top) — order: kw_toxin (top), both, venom_tissue
    zone_vt = {"bottom": right_start, "top": right_start + h_vt, "count": n_vt}
    zone_both = {
        "bottom": zone_vt["top"] + gap,
        "top": zone_vt["top"] + gap + h_both,
        "count": n_both,
    }
    zone_kw = {
        "bottom": zone_both["top"] + gap,
        "top": zone_both["top"] + gap + h_kw,
        "count": n_kw,
    }
    right_zones = {
        "both": zone_both,
        "venom_tissue": zone_vt,
        "kw_toxin": zone_kw,
    }

    # Ellipse geometry — each ellipse covers its exclusive zone + the overlap.
    # Bottom ellipse (venom tissue): spans zone_vt bottom → zone_both top
    # Top ellipse (kw_toxin):        spans zone_both bottom → zone_kw top
    vt_lo = zone_vt["bottom"]
    vt_hi = zone_both["top"]
    kw_lo = zone_both["bottom"]
    kw_hi = zone_kw["top"]

    ry_vt = (vt_hi - vt_lo) / 2  # vertical semi-axis
    cy_vt = (vt_hi + vt_lo) / 2
    ry_kw = (kw_hi - kw_lo) / 2
    cy_kw = (kw_hi + kw_lo) / 2

    # Compute rx so that ellipses render as circles on screen.
    # For a circle: rx_data / ry_data = (fig_h * x_range) / (fig_w * y_range)
    fig_w, fig_h = ax.get_figure().get_size_inches()
    margin_lbl = 2.5  # extra vertical room for labels above/below Venn
    y_lo = min(zone_vt["bottom"], 0) - margin_lbl
    y_hi = max(zone_kw["top"], total_left_h) + margin_lbl
    y_range = y_hi - y_lo
    # Estimate x_range (will finalise in axis limits)
    x_range_est = (x_venn + 5) - (x_left - 2.5)
    aspect_corr = (fig_h * x_range_est) / (fig_w * y_range)
    rx_vt = ry_vt * aspect_corr
    rx_kw = ry_kw * aspect_corr

    # --- "Both" color: purple as perceptual mix of orange + blue -------------
    both_color = "#7B5EA7"  # muted violet — reads as intersection of warm + cool
    flow_colors = {
        "both": both_color,
        "kw_toxin": DEFINITION_COLORS["kw_toxin"],
        "venom_tissue": DEFINITION_COLORS["venom_tissue"],
    }

    # --- draw Venn circles (exclusive regions + purple intersection) ---------
    # Draw both full circles
    ell_vt = Ellipse(
        (x_venn, cy_vt),
        width=rx_vt * 2,
        height=ry_vt * 2,
        facecolor=DEFINITION_COLORS["venom_tissue"],
        edgecolor=DEFINITION_COLORS["venom_tissue"],
        alpha=0.15,
        linewidth=2,
    )
    ax.add_patch(ell_vt)

    ell_kw = Ellipse(
        (x_venn, cy_kw),
        width=rx_kw * 2,
        height=ry_kw * 2,
        facecolor=DEFINITION_COLORS["kw_toxin"],
        edgecolor=DEFINITION_COLORS["kw_toxin"],
        alpha=0.15,
        linewidth=2,
    )
    ax.add_patch(ell_kw)

    # White mask to erase the blue+orange blend in the intersection
    mask = Ellipse(
        (x_venn, cy_vt),
        width=rx_vt * 2,
        height=ry_vt * 2,
        facecolor="white",
        edgecolor="none",
        alpha=1.0,
    )
    ax.add_patch(mask)
    clip_for_mask = Ellipse(
        (x_venn, cy_kw),
        width=rx_kw * 2,
        height=ry_kw * 2,
        transform=ax.transData,
    )
    mask.set_clip_path(clip_for_mask)

    # Clean purple intersection fill
    inter = Ellipse(
        (x_venn, cy_vt),
        width=rx_vt * 2,
        height=ry_vt * 2,
        facecolor=both_color,
        edgecolor="none",
        alpha=0.20,
    )
    ax.add_patch(inter)
    clip_for_inter = Ellipse(
        (x_venn, cy_kw),
        width=rx_kw * 2,
        height=ry_kw * 2,
        transform=ax.transData,
    )
    inter.set_clip_path(clip_for_inter)

    # Purple border around the intersection — softer than full, darker than fill
    border1 = Ellipse(
        (x_venn, cy_vt),
        width=rx_vt * 2,
        height=ry_vt * 2,
        facecolor="none",
        edgecolor=both_color,
        linewidth=2,
        alpha=0.4,
    )
    ax.add_patch(border1)
    clip_b1 = Ellipse(
        (x_venn, cy_kw),
        width=rx_kw * 2,
        height=ry_kw * 2,
        transform=ax.transData,
    )
    border1.set_clip_path(clip_b1)

    border2 = Ellipse(
        (x_venn, cy_kw),
        width=rx_kw * 2,
        height=ry_kw * 2,
        facecolor="none",
        edgecolor=both_color,
        linewidth=2,
        alpha=0.4,
    )
    ax.add_patch(border2)
    clip_b2 = Ellipse(
        (x_venn, cy_vt),
        width=rx_vt * 2,
        height=ry_vt * 2,
        transform=ax.transData,
    )
    border2.set_clip_path(clip_b2)

    # --- Ellipse left-edge helper -------------------------------------------
    def ellipse_left_x(y: float, cy: float, ry: float, r_x: float) -> float:
        """Return the x of the left edge of the ellipse at height *y*."""
        t = (y - cy) / ry
        t = np.clip(t, -1.0, 1.0)
        return x_venn - r_x * np.sqrt(1.0 - t * t)

    def flow_target_x(defn: str, y_mid: float) -> float:
        """Compute where a flow for *defn* should land at height *y_mid*.

        - kw_toxin:      left edge of the kw_toxin ellipse
        - venom_tissue:  left edge of the venom_tissue ellipse
        - both:          left edge of the intersection (the more-rightward
                         of the two ellipse left edges)
        """
        if defn == "kw_toxin":
            return ellipse_left_x(y_mid, cy_kw, ry_kw, rx_kw)
        elif defn == "venom_tissue":
            return ellipse_left_x(y_mid, cy_vt, ry_vt, rx_vt)
        else:  # both — intersection boundary
            return max(
                ellipse_left_x(y_mid, cy_kw, ry_kw, rx_kw),
                ellipse_left_x(y_mid, cy_vt, ry_vt, rx_vt),
            )

    # --- Bezier flow helper (polyline right edge on ellipse) -----------------
    n_arc = 30  # points to sample along the ellipse arc for the right edge

    def draw_flow(x1, y1b, y1t, y2b, y2t, defn, alpha=0.35):
        """Draw flow whose right edge is a polyline sampled on the ellipse."""
        color = flow_colors[defn]

        # Sample right edge along the ellipse from y2b → y2t
        ys_up = np.linspace(y2b, y2t, n_arc)
        arc_up = [(flow_target_x(defn, y), y) for y in ys_up]

        x2b, _ = arc_up[0]
        x2t, _ = arc_up[-1]
        x2_mid = (x2b + x2t) / 2
        ctrl = (x2_mid - x1) * 0.45

        verts = [
            # Bottom Bezier: left → right
            (x1 + box_w / 2, y1b),
            (x1 + box_w / 2 + ctrl, y1b),
            (x2b - ctrl, y2b),
            (x2b, y2b),
        ]
        codes = [
            MplPath.MOVETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
        ]

        # Right edge: polyline following the ellipse (skip first, already there)
        for pt in arc_up[1:]:
            verts.append(pt)
            codes.append(MplPath.LINETO)

        # Top Bezier: right → left
        verts += [
            (x2t - ctrl, y2t),
            (x1 + box_w / 2 + ctrl, y1t),
            (x1 + box_w / 2, y1t),
        ]
        codes += [
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
        ]

        # Close
        verts.append((x1 + box_w / 2, y1b))
        codes.append(MplPath.CLOSEPOLY)

        ax.add_patch(
            PathPatch(MplPath(verts, codes), fc=color, ec="none", alpha=alpha)
        )

    # --- draw flows ---------------------------------------------------------
    # Stacking order: venom_tissue (bottom), both (middle), kw_toxin (top)
    stack_order = ["venom_tissue", "both", "kw_toxin"]

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

    # --- draw left phylum bars ----------------------------------------------
    small_thr = 1.5
    for phylum in phyla:
        pos = left_pos[phylum]
        bar_h = pos["top"] - pos["bottom"]
        tot = pos["count"]

        active = [d for d in DEFINITION_ORDER if ct.loc[phylum, d] > 0]
        lbl_color = flow_colors[active[0]] if len(active) == 1 else "black"

        # stacked coloured segments (venom_tissue → both → kw_toxin)
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
        if bar_h > small_thr:
            ax.text(
                x_left - box_w / 2 - 0.1,
                pos["bottom"] + bar_h / 2,
                phylum,
                ha="right",
                va="center",
                fontsize=10,
                fontstyle="italic",
                color=lbl_color,
            )
            ax.text(
                x_left,
                pos["bottom"] + bar_h / 2,
                f"{tot:,}",
                ha="center",
                va="center",
                fontsize=9,
                fontweight="bold",
                color="white",
            )
        else:
            ax.text(
                x_left - box_w / 2 - 0.1,
                pos["bottom"] + bar_h / 2,
                f"{phylum}  ({tot:,})",
                ha="right",
                va="center",
                fontsize=9,
                fontstyle="italic",
                color=lbl_color,
            )

    # --- Venn region labels -------------------------------------------------
    # Count labels inside each zone
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

    # Definition labels above / below the Venn
    ax.text(
        x_venn,
        zone_kw["top"] + 1.0,
        "Toxin Keyword",
        ha="center",
        va="bottom",
        fontsize=14,
        fontweight="bold",
        color=flow_colors["kw_toxin"],
    )
    ax.text(
        x_venn,
        zone_vt["bottom"] - 1.0,
        "Venom Tissue",
        ha="center",
        va="top",
        fontsize=14,
        fontweight="bold",
        color=flow_colors["venom_tissue"],
    )

    # --- axis limits --------------------------------------------------------
    rx_max = max(rx_vt, rx_kw)
    ax.set_xlim(x_left - 2.5, x_venn + rx_max + 1.0)
    ax.set_ylim(y_lo, y_hi)
    ax.axis("off")


# ---------------------------------------------------------------------------
# Figure creation
# ---------------------------------------------------------------------------


def create_combined_definition_figure(
    df: pd.DataFrame,
    output_dir: Path,
) -> Path:
    """Create and save the phylum→Venn flow figure."""
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(14, 7))
    plot_phylum_venn_flow(df, ax)
    plt.tight_layout()

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
