#!/usr/bin/env python
"""Render comparison/MANUSCRIPT.md to the print PDF (review copy for O.G.).
Markdown subset: #/##/### headings, **bold**, *italic*, `code`, - bullets, | tables,
![img](path). Arial family registered for full Latin+Greek coverage; unicode
superscripts converted to <super> markup (reportlab built-in glyph limitation)."""
import re, sys, os
from PIL import Image as PILImage
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                Table, TableStyle, HRFlowable)

import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from project_paths import COMPARISON as _COMPARISON
ROOT = str(_COMPARISON)
SRC = ROOT + "/MANUSCRIPT.md"
OUT = ROOT + "/MANUSCRIPT_print.pdf"

F = "/System/Library/Fonts/Supplemental/"
pdfmetrics.registerFont(TTFont("Arial", F + "Arial.ttf"))
pdfmetrics.registerFont(TTFont("Arial-Bold", F + "Arial Bold.ttf"))
pdfmetrics.registerFont(TTFont("Arial-Italic", F + "Arial Italic.ttf"))
pdfmetrics.registerFont(TTFont("Arial-BoldItalic", F + "Arial Bold Italic.ttf"))
pdfmetrics.registerFontFamily("Arial", normal="Arial", bold="Arial-Bold",
                              italic="Arial-Italic", boldItalic="Arial-BoldItalic")
pdfmetrics.registerFont(TTFont("Courier10", F + "Courier New.ttf"))

BODY = ParagraphStyle("body", fontName="Arial", fontSize=9.3, leading=12.6,
                      spaceAfter=6, alignment=4)
H1 = ParagraphStyle("h1", parent=BODY, fontName="Arial-Bold", fontSize=13.5,
                    leading=16, spaceBefore=14, spaceAfter=8, alignment=0)
H2 = ParagraphStyle("h2", parent=BODY, fontName="Arial-Bold", fontSize=11,
                    leading=14, spaceBefore=10, spaceAfter=6, alignment=0)
TITLE = ParagraphStyle("title", parent=H1, fontSize=15.5, leading=19)
BULLET = ParagraphStyle("bullet", parent=BODY, leftIndent=14, bulletIndent=4)
CELL = ParagraphStyle("cell", parent=BODY, fontSize=8, leading=10, spaceAfter=0,
                      alignment=0)
CELLH = ParagraphStyle("cellh", parent=CELL, fontName="Arial-Bold")
CAPTION = ParagraphStyle("cap", parent=BODY, fontSize=8.6, leading=11.4)

SUP = {"⁰":"0","¹":"1","²":"2","³":"3","⁴":"4","⁵":"5","⁶":"6","⁷":"7","⁸":"8",
       "⁹":"9","⁻":"-","⁺":"+"}

def inline(s):
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    s = re.sub("[" + "".join(SUP) + "]+",
               lambda m: "<super>" + "".join(SUP[c] for c in m.group(0)) + "</super>", s)
    s = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", s, flags=re.S)
    # emphase : l\'astérisque ouvrante doit être suivie d\'un non-espace et la fermante
    # précédée d\'un non-espace — sinon une astérisque de significativité (« * padj < 0.05 »)
    # ouvre une italique qui contamine tout le reste du bloc (bug constaté sur les blurbs S6-S15)
    s = re.sub(r"(?<!\*)\*(?!\s)([^*\n]+?)(?<!\s)\*(?!\*)", r"<i>\1</i>", s)
    s = re.sub(r"`([^`]+)`", r'<font face="Courier10" size="8">\1</font>', s)
    return s

PAGE_W = A4[0] - 2 * 1.7 * cm

def image_flowable(path):
    im = PILImage.open(path)
    dpi = im.info.get("dpi", (600, 600))[0] or 600
    w_in, h_in = im.size[0] / dpi, im.size[1] / dpi
    w = min(w_in * 72, PAGE_W)
    h = h_in * 72 * (w / (w_in * 72))
    max_h = A4[1] - 2 * 1.7 * cm - 20
    if h > max_h:
        w *= max_h / h; h = max_h
    return Image(path, width=w, height=h)

story = []
lines = open(SRC).read().split("\n")
i, para = 0, []

def flush(style=BODY):
    global para
    if para:
        story.append(Paragraph(inline(" ".join(para)), style))
        para = []

while i < len(lines):
    l = lines[i]
    if l.startswith("# "):
        flush(); story.append(Paragraph(inline(l[2:]), TITLE))
    elif l.startswith("## "):
        flush(); story.append(HRFlowable(width="100%", thickness=0.6, color=colors.HexColor("#999999"), spaceBefore=8))
        story.append(Paragraph(inline(l[3:]), H1))
    elif l.startswith("### "):
        flush(); story.append(Paragraph(inline(l[4:]), H2))
    elif l.startswith("!["):
        flush()
        m = re.match(r"!\[[^\]]*\]\(([^)]+)\)", l)
        p = os.path.join(ROOT, m.group(1))
        story.append(Spacer(1, 6)); story.append(image_flowable(p)); story.append(Spacer(1, 4))
    elif l.startswith("|"):
        flush()
        rows = []
        while i < len(lines) and lines[i].startswith("|"):
            cells = [c.strip() for c in lines[i].strip("|").split("|")]
            if not all(re.fullmatch(r":?-{2,}:?", c or "-") for c in cells):
                rows.append(cells)
            i += 1
        i -= 1
        ncol = max(len(r) for r in rows)
        data = [[Paragraph(inline(c), CELLH if k == 0 else CELL)
                 for c in r + [""] * (ncol - len(r))] for k, r in enumerate(rows)]
        tab = Table(data, colWidths=[PAGE_W / ncol] * ncol, repeatRows=1)
        tab.setStyle(TableStyle([
            ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#bbbbbb")),
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#efefef")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
            ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]))
        story.append(Spacer(1, 4)); story.append(tab); story.append(Spacer(1, 6))
    elif l.startswith("- "):
        flush()
        story.append(Paragraph(inline(l[2:]), BULLET, bulletText="•"))
    elif not l.strip():
        flush()
    else:
        # legend paragraphs (bold-start Figure lines) in caption style
        if l.startswith("**Figure") and not para:
            para = [l]
            while i + 1 < len(lines) and lines[i + 1].strip() and not lines[i + 1].startswith(("#", "|", "-", "!")):
                i += 1; para.append(lines[i])
            flush(CAPTION)
        else:
            para.append(l)
    i += 1
flush()

doc = SimpleDocTemplate(OUT, pagesize=A4,
                        leftMargin=1.7 * cm, rightMargin=1.7 * cm,
                        topMargin=1.7 * cm, bottomMargin=1.7 * cm,
                        title="FC transcriptomic comparison — print review copy")
doc.build(story)
from pypdf import PdfReader
r = PdfReader(OUT)
print(f"OK: {OUT} — {len(r.pages)} pages, {os.path.getsize(OUT)/1e6:.1f} MB")
