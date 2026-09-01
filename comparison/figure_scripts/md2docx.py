#!/usr/bin/env python
"""Render comparison/MANUSCRIPT.md to a Word .docx (editable copy for co-authors).

Same Markdown subset as md2pdf.py: #/##/### headings, **bold**, *italic*, `code`,
- bullets, | tables |, ![img](path), --- rule. Built directly as OOXML (a .docx is a
ZIP of XML) because this machine has no pandoc / node-docx / python-docx / LibreOffice.

Real Heading1-3 styles are emitted so Word's navigation pane and any TOC work, and
bullets use a proper numbering definition rather than literal bullet characters.
"""
import os
import re
import struct
import sys
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from project_paths import COMPARISON

SRC = COMPARISON / "MANUSCRIPT.md"
OUT = COMPARISON / "MANUSCRIPT.docx"

EMU_PER_IN = 914400
CONTENT_IN = 6.27                      # A4 minus 1" margins
CONTENT_DXA = int(CONTENT_IN * 1440)
FONT = "Arial"


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def png_size(path):
    """Width/height from the PNG IHDR (no image library needed)."""
    with open(path, "rb") as fh:
        head = fh.read(24)
    assert head[:8] == b"\x89PNG\r\n\x1a\n", path
    return struct.unpack(">II", head[16:24])


# ---------- inline runs ----------
TOKEN = re.compile(r"(\*\*.+?\*\*|(?<!\*)\*[^*]+?\*(?!\*)|`[^`]+`)", re.S)


def runs(text, bold=False, size=None, color=None):
    out = []
    for part in TOKEN.split(text):
        if not part:
            continue
        b, i, mono = bold, False, False
        if part.startswith("**") and part.endswith("**"):
            part, b = part[2:-2], True
        elif part.startswith("*") and part.endswith("*"):
            part, i = part[1:-1], True
        elif part.startswith("`") and part.endswith("`"):
            part, mono = part[1:-1], True
        rpr = ""
        if b:
            rpr += "<w:b/>"
        if i:
            rpr += "<w:i/>"
        if mono:
            rpr += '<w:rFonts w:ascii="Courier New" w:hAnsi="Courier New"/>'
        if size:
            rpr += f'<w:sz w:val="{size}"/><w:szCs w:val="{size}"/>'
        if color:
            rpr += f'<w:color w:val="{color}"/>'
        rpr = f"<w:rPr>{rpr}</w:rPr>" if rpr else ""
        out.append(f'<w:r>{rpr}<w:t xml:space="preserve">{esc(part)}</w:t></w:r>')
    return "".join(out)


def para(text, style=None, size=None, jc=None, bullet=False, spacing=None):
    ppr = ""
    if style:
        ppr += f'<w:pStyle w:val="{style}"/>'
    if bullet:
        ppr += '<w:numPr><w:ilvl w:val="0"/><w:numId w:val="1"/></w:numPr>'
    if jc:
        ppr += f'<w:jc w:val="{jc}"/>'
    if spacing:
        ppr += f'<w:spacing w:before="{spacing[0]}" w:after="{spacing[1]}"/>'
    ppr = f"<w:pPr>{ppr}</w:pPr>" if ppr else ""
    return f"<w:p>{ppr}{runs(text, size=size)}</w:p>"


def image_para(rid, w_emu, h_emu, idx):
    return (
        '<w:p><w:pPr><w:jc w:val="center"/><w:spacing w:before="160" w:after="160"/></w:pPr>'
        f'<w:r><w:drawing><wp:inline distT="0" distB="0" distL="0" distR="0">'
        f'<wp:extent cx="{w_emu}" cy="{h_emu}"/>'
        f'<wp:docPr id="{idx}" name="Figure {idx}"/><wp:cNvGraphicFramePr/>'
        '<a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
        '<a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">'
        '<pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">'
        f'<pic:nvPicPr><pic:cNvPr id="{idx}" name="Figure {idx}"/><pic:cNvPicPr/></pic:nvPicPr>'
        f'<pic:blipFill><a:blip r:embed="{rid}"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill>'
        f'<pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{w_emu}" cy="{h_emu}"/></a:xfrm>'
        '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr>'
        "</pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>"
    )


def table(rows):
    ncol = max(len(r) for r in rows)
    widths = [CONTENT_DXA // ncol] * ncol
    widths[-1] = CONTENT_DXA - sum(widths[:-1])
    grid = "".join(f'<w:gridCol w:w="{w}"/>' for w in widths)
    body = []
    for ri, row in enumerate(rows):
        cells = []
        for ci in range(ncol):
            txt = row[ci] if ci < len(row) else ""
            shade = ('<w:shd w:val="clear" w:color="auto" w:fill="EDEDED"/>' if ri == 0 else "")
            cells.append(
                f'<w:tc><w:tcPr><w:tcW w:w="{widths[ci]}" w:type="dxa"/>{shade}</w:tcPr>'
                f'<w:p><w:pPr><w:spacing w:before="20" w:after="20"/></w:pPr>'
                f"{runs(txt, bold=(ri == 0), size=16)}</w:p></w:tc>")
        body.append(f"<w:tr>{''.join(cells)}</w:tr>")
    return (
        f'<w:tbl><w:tblPr><w:tblW w:w="{CONTENT_DXA}" w:type="dxa"/>'
        '<w:tblBorders>' + "".join(
            f'<w:{e} w:val="single" w:sz="4" w:space="0" w:color="AAAAAA"/>'
            for e in ("top", "left", "bottom", "right", "insideH", "insideV")) +
        "</w:tblBorders></w:tblPr>"
        f"<w:tblGrid>{grid}</w:tblGrid>{''.join(body)}</w:tbl><w:p/>")


# ---------- parse ----------
lines = SRC.read_text(encoding="utf-8").split("\n")
body, images = [], []
i = 0
while i < len(lines):
    line = lines[i].rstrip()
    s = line.strip()
    if not s:
        i += 1
        continue

    if s.startswith("!["):                                   # figure
        m = re.match(r"!\[([^\]]*)\]\(([^)]+)\)", s)
        path = COMPARISON / m.group(2)
        pw, ph = png_size(path)
        w = int(CONTENT_IN * EMU_PER_IN)
        h = int(w * ph / pw)
        rid = f"rIdImg{len(images) + 1}"
        images.append((rid, path))
        body.append(image_para(rid, w, h, len(images)))
    elif s.startswith("### "):
        body.append(para(s[4:], style="Heading3"))
    elif s.startswith("## "):
        body.append(para(s[3:], style="Heading2"))
    elif s.startswith("# "):
        body.append(para(s[2:], style="Heading1"))
    elif s.startswith("|"):                                   # table block
        rows = []
        while i < len(lines) and lines[i].strip().startswith("|"):
            cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
            if not all(set(c) <= set("-: ") and c for c in cells):   # skip separator
                rows.append(cells)
            i += 1
        body.append(table(rows))
        continue
    elif s.startswith("- "):
        body.append(para(s[2:], bullet=True))
    elif set(s) <= set("-") and len(s) > 2:
        body.append('<w:p><w:pPr><w:pBdr><w:bottom w:val="single" w:sz="6" w:space="1" '
                    'w:color="999999"/></w:pBdr></w:pPr></w:p>')
    else:
        cap = s.startswith("**Figure") or s.startswith("**Supplementary")
        body.append(para(s, size=18 if cap else None, jc="both" if not cap else None))
    i += 1

# ---------- package ----------
NS = ('xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
      'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
      'xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"')
document = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    f"<w:document {NS}><w:body>" + "".join(body) +
    '<w:sectPr><w:pgSz w:w="11906" w:h="16838"/>'
    '<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" '
    'w:header="708" w:footer="708" w:gutter="0"/></w:sectPr></w:body></w:document>')


def heading_style(sid, name, sz, before, after):
    return (f'<w:style w:type="paragraph" w:styleId="{sid}"><w:name w:val="{name}"/>'
            '<w:basedOn w:val="Normal"/><w:next w:val="Normal"/><w:qFormat/>'
            f'<w:pPr><w:keepNext/><w:spacing w:before="{before}" w:after="{after}"/>'
            f'<w:outlineLvl w:val="{sid[-1]}"/></w:pPr>'
            f'<w:rPr><w:b/><w:sz w:val="{sz}"/><w:szCs w:val="{sz}"/></w:rPr></w:style>')


styles = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    f"<w:styles {NS}>"
    '<w:docDefaults><w:rPrDefault><w:rPr>'
    f'<w:rFonts w:ascii="{FONT}" w:hAnsi="{FONT}" w:cs="{FONT}"/>'
    '<w:sz w:val="20"/><w:szCs w:val="20"/></w:rPr></w:rPrDefault>'
    '<w:pPrDefault><w:pPr><w:spacing w:after="120" w:line="276" w:lineRule="auto"/>'
    "</w:pPr></w:pPrDefault></w:docDefaults>"
    '<w:style w:type="paragraph" w:default="1" w:styleId="Normal">'
    '<w:name w:val="Normal"/><w:qFormat/></w:style>'
    + heading_style("Heading0", "heading 1", 32, 240, 160).replace("Heading0", "Heading1")
    + heading_style("Heading1", "heading 2", 26, 280, 140).replace(
        'w:styleId="Heading1"', 'w:styleId="Heading2"').replace('w:outlineLvl w:val="1"', 'w:outlineLvl w:val="1"')
    + heading_style("Heading2", "heading 3", 22, 240, 120).replace(
        'w:styleId="Heading2"', 'w:styleId="Heading3"')
    + "</w:styles>")

numbering = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    f"<w:numbering {NS}>"
    '<w:abstractNum w:abstractNumId="0"><w:lvl w:ilvl="0">'
    '<w:start w:val="1"/><w:numFmt w:val="bullet"/><w:lvlText w:val="•"/>'
    '<w:lvlJc w:val="left"/><w:pPr><w:ind w:left="360" w:hanging="240"/></w:pPr>'
    '<w:rPr><w:rFonts w:ascii="Symbol" w:hAnsi="Symbol" w:hint="default"/></w:rPr>'
    "</w:lvl></w:abstractNum>"
    '<w:num w:numId="1"><w:abstractNumId w:val="0"/></w:num></w:numbering>')

img_rels = "".join(
    f'<Relationship Id="{rid}" '
    'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" '
    f'Target="media/{os.path.basename(p)}"/>' for rid, p in images)
doc_rels = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
    '<Relationship Id="rIdStyles" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
    '<Relationship Id="rIdNum" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/numbering" Target="numbering.xml"/>'
    f"{img_rels}</Relationships>")

content_types = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
    '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
    '<Default Extension="xml" ContentType="application/xml"/>'
    '<Default Extension="png" ContentType="image/png"/>'
    '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
    '<Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>'
    '<Override PartName="/word/numbering.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.numbering+xml"/>'
    "</Types>")

root_rels = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>'
    "</Relationships>")

with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
    z.writestr("[Content_Types].xml", content_types)
    z.writestr("_rels/.rels", root_rels)
    z.writestr("word/document.xml", document)
    z.writestr("word/styles.xml", styles)
    z.writestr("word/numbering.xml", numbering)
    z.writestr("word/_rels/document.xml.rels", doc_rels)
    for _, p in images:
        z.write(p, f"word/media/{os.path.basename(p)}")

print(f"OK: {OUT} — {len(body)} blocs, {len(images)} figures, "
      f"{OUT.stat().st_size / 1e6:.1f} MB")
