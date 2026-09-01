#!/usr/bin/env python
"""Figure style for comparison/ figures — RE-EXPORT of charte_graphique (single source of truth).

Étape 6 (2026-08-10): this module no longer defines any colour or size of its own; it re-exports
charte_graphique/palette.py and applies charte_graphique/article.mplstyle. The historical local
definitions (C_BG '#1a4fb4', 16 pt panel letters) diverged from the charte and are gone.
Rules (charte §0): panel letters lowercase bold 12 pt; nothing < 7 pt; no line < 0.75 pt;
50 %-reduction test before validating any figure.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))     # project root
_CHARTE = os.path.join(_ROOT, "charte_graphique")
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from charte_graphique.palette import (   # noqa: F401  (re-exports)
    GROUP as ENTITY, LAYER, SEMANTIC, INK, seq_blue, seq_aqua, div_bwr,
    ORDER, RELIEF, GROUP_EDGE, rgba,
)

import matplotlib.pyplot as plt


def use_style():
    """Apply the article mplstyle (fonts, sizes, export settings). Call FIRST in every figure script."""
    plt.style.use(os.path.join(_CHARTE, "article.mplstyle"))


# --- the three biological groups of the comparison, mapped onto charte entities ---
C_FC = ENTITY["FC"]          # '#d62728' red   — Fañanas cell / displaced Bergmann (focus)
C_BG = ENTITY["BG"]          # '#2a78d6' blue  — canonical PCL Bergmann
C_OTHER = ENTITY["Astro_CB"]  # '#009e73' green — other cerebellar (parenchymal) astrocytes
GROUP = {"FC": C_FC, "BG": C_BG, "OtherCBastro": C_OTHER}
GROUP_LABEL = {"FC": "FC (Fañanas cell)", "BG": "BG (canonical Bergmann)",
               "OtherCBastro": "parenchymal astrocytes"}

# qualitative palette for gene/score series (never the three group colours)
GENE_PAL = ["#ff7f0e", "#9467bd", "#17becf", "#8c564b", "#bcbd22", "#e377c2"]

# canonical dataset order + display names (one form everywhere)
DS_ORDER = ["zhuang", "hao", "suyama"]
DSL = {"zhuang": "MERFISH", "hao": "Stereo-seq", "suyama": "Xenium"}

# figure widths (charte §1), inches
W1, W15, W2 = 3.50, 4.72, 7.20


def panel_letter(ax, letter, dx=-0.09, dy=1.06, fs=12):
    """Charte §2–3: lowercase, bold, 12 pt, outside the axes top-left."""
    ax.text(dx, dy, letter.lower(), transform=ax.transAxes, fontsize=fs,
            fontweight="bold", va="top", ha="right")


def qa_overlaps(fig, min_px=3.0):
    """Objective layout QA: report every pair of rendered text artists whose bounding boxes
    overlap by more than min_px in both directions, plus texts clipped by the figure edge.
    Returns a list of problem strings (empty = clean)."""
    fig.canvas.draw()
    items = []
    for ax in fig.get_axes():
        cands = list(ax.texts) + [ax.title, ax.xaxis.label, ax.yaxis.label]
        if ax.axison:  # ticks of axis('off') axes report visible=True but are never drawn
            xlo, xhi = sorted(ax.get_xlim()); ylo, yhi = sorted(ax.get_ylim())
            for t in ax.get_xticklabels():
                px = t.get_position()[0]
                if xlo - 1e-9 <= px <= xhi + 1e-9:   # out-of-range ticks are never drawn
                    cands.append(t)
            for t in ax.get_yticklabels():
                py = t.get_position()[1]
                if ylo - 1e-9 <= py <= yhi + 1e-9:
                    cands.append(t)
        leg = ax.get_legend()
        if leg is not None:
            cands += list(leg.get_texts())
        for t in cands:
            if t.get_visible() and t.get_text().strip():
                items.append(t)
    for leg in fig.legends:
        items += [t for t in leg.get_texts() if t.get_text().strip()]
    items += [t for t in fig.texts if t.get_text().strip()]
    boxes = []
    for t in items:
        try:
            boxes.append((t, t.get_window_extent()))
        except Exception:
            pass
    probs = []
    for i in range(len(boxes)):
        ti, bi = boxes[i]
        for j in range(i + 1, len(boxes)):
            tj, bj = boxes[j]
            ox = min(bi.x1, bj.x1) - max(bi.x0, bj.x0)
            oy = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
            if ox > min_px and oy > min_px:
                probs.append(f"OVERLAP: '{ti.get_text()[:40]}' <-> '{tj.get_text()[:40]}' ({ox:.0f}x{oy:.0f}px)")
    # charte §0 floors: no text < 7 pt; report text clipped by the figure edge (unless clip disabled on purpose)
    fw, fh = fig.get_window_extent().x1, fig.get_window_extent().y1
    for t, b in boxes:
        if t.get_fontsize() < 6.9:
            probs.append(f"FONT<7pt: '{t.get_text()[:40]}' ({t.get_fontsize():.1f} pt)")
        if (b.x0 < -2 or b.y0 < -2 or b.x1 > fw + 2 or b.y1 > fh + 2) and not t.get_clip_on():
            pass  # bbox_inches='tight' extends the canvas; outside-canvas text is kept by export
    # charte §0: no line thinner than 0.75 pt (data lines; spines excluded — mplstyle governs them)
    for ax in fig.get_axes():
        for ln in ax.get_lines():
            if ln.get_visible() and ln.get_linewidth() < 0.74 and ln.get_linestyle() != "None":
                probs.append(f"LINE<0.75pt: line in axes '{ax.get_title() or ax.get_ylabel()[:30]}' "
                             f"({ln.get_linewidth():.2f} pt)")
                break
    return probs


def save(fig, name, figdir=None, qa=True):
    """Uniform export: PDF (master) + SVG (editable) + PNG 600 dpi (charte §1) + layout QA."""
    figdir = figdir or os.path.join(_ROOT, "comparison", "figures")
    if qa:
        probs = qa_overlaps(fig)
        if probs:
            print(f"!! LAYOUT QA — {name}: {len(probs)} problem(s)")
            for p in probs:
                print("   ", p)
        else:
            print(f"LAYOUT QA — {name}: clean")
    for ext in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(figdir, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.06)
    print("saved", name, "->", figdir)
