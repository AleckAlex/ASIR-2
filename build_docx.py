"""Build Alexopoulos_ASIR_paper_v4.docx from manuscript.md (pandoc + table styling)."""
import subprocess
from docx import Document
from docx.shared import Pt
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

OUT = "Alexopoulos_ASIR_paper_v4.docx"
subprocess.run(["pandoc", "manuscript.md", "-o", OUT, "--reference-doc=ref_paper.docx"], check=True)

d = Document(OUT)


TBLPR = ["tblStyle", "tblpPr", "tblOverlap", "bidiVisual", "tblStyleRowBandSize", "tblStyleColBandSize",
         "tblW", "jc", "tblCellSpacing", "tblInd", "tblBorders", "shd", "tblLayout", "tblCellMar",
         "tblLook", "tblCaption", "tblDescription"]
TCPR = ["cnfStyle", "tcW", "gridSpan", "hMerge", "vMerge", "tcBorders", "shd", "noWrap", "tcMar",
        "textDirection", "tcFitText", "vAlign", "hideMark"]


def reorder(parent, order):
    rank = {n: i for i, n in enumerate(order)}
    kids = list(parent)
    for k in kids:
        parent.remove(k)
    for k in sorted(kids, key=lambda e: rank.get(e.tag.split("}")[1], len(order))):
        parent.append(k)


def add_borders(parent, tag, sides, sz):
    b = parent.find(qn(tag))
    if b is None:
        b = OxmlElement(tag); parent.append(b)
        reorder(parent, TBLPR if tag == "w:tblBorders" else TCPR)
    for s in sides:
        e = OxmlElement(f"w:{s}")
        for k, v in (("val", "single"), ("sz", str(sz)), ("space", "0"), ("color", "000000")):
            e.set(qn(f"w:{k}"), v)
        b.append(e)


for t in d.tables:                                   # three-rule journal tables
    add_borders(t._tbl.tblPr, "w:tblBorders", ["top", "bottom"], 10)
    for cell in t.rows[0].cells:
        add_borders(cell._tc.get_or_add_tcPr(), "w:tcBorders", ["bottom"], 6)
    for ri, row in enumerate(t.rows):
        for cell in row.cells:
            for p in cell.paragraphs:
                pf = p.paragraph_format
                pf.space_after = Pt(1); pf.space_before = Pt(1); pf.alignment = 0
                for r in p.runs:
                    r.font.size = Pt(9)
                    if ri == 0:
                        r.font.bold = True

for p in d.paragraphs:                               # keep table captions with tables
    if p.text.startswith(("Table ", "Table A")) and p.runs and p.runs[0].bold:
        p.paragraph_format.keep_with_next = True

d.save(OUT)

# --- schema tidy-ups for pandoc's default parts (settings order, page margins)
import re, shutil, zipfile, os
from lxml import etree
SETTINGS = ["writeProtection", "view", "zoom", "removePersonalInformation", "embedTrueTypeFonts",
            "embedSystemFonts", "saveSubsetFonts", "proofState", "stylePaneFormatFilter", "doNotTrackMoves",
            "defaultTabStop", "drawingGridHorizontalSpacing", "drawingGridVerticalSpacing",
            "displayHorizontalDrawingGridEvery", "displayVerticalDrawingGridEvery", "characterSpacingControl",
            "savePreviewPicture", "footnotePr", "endnotePr", "compat", "docVars", "rsids", "mathPr",
            "themeFontLang", "clrSchemeMapping", "shapeDefaults", "decimalSymbol", "listSeparator"]
tmp = OUT + ".tmp"
with zipfile.ZipFile(OUT) as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        if item.filename.endswith("/"):
            continue
        data = zin.read(item.filename)
        if item.filename == "word/settings.xml":
            root = etree.fromstring(data)
            reorder(root, SETTINGS)
            data = etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)
        elif item.filename == "word/document.xml":
            x = data.decode("utf-8")
            x = re.sub(r'<w:pgMar (?![^>]*w:header)([^>]*)/>',
                       r'<w:pgMar \1 w:header="708" w:footer="708" w:gutter="0"/>', x)
            def fix_dpr(m):                         # OMML order: begChr, sepChr, endChr, grow
                body = m.group(1)
                parts = re.findall(r'<m:(begChr|sepChr|endChr|grow|shp|ctrlPr)\b[^>]*/>', body)
                els = re.findall(r'(<m:(?:begChr|sepChr|endChr|grow|shp|ctrlPr)\b[^>]*/>)', body)
                rank = {"begChr": 0, "sepChr": 1, "endChr": 2, "grow": 3, "shp": 4, "ctrlPr": 5}
                return "<m:dPr>" + "".join(e for _, e in sorted(zip(parts, els), key=lambda p: rank[p[0]])) + "</m:dPr>"
            x = re.sub(r'<m:dPr>((?:<m:(?:begChr|sepChr|endChr|grow|shp|ctrlPr)\b[^>]*/>)*)</m:dPr>', fix_dpr, x)
            x = re.sub(r'<w:pPr>(<w:numPr>.*?</w:numPr>)(<w:pStyle [^>]*/>)', r'<w:pPr>\2\1', x)
            data = x.encode("utf-8")
        if item.filename == "word/styles.xml":
            root = etree.fromstring(data)
            ORD = {
                "style": ["name", "aliases", "basedOn", "next", "link", "autoRedefine", "uiPriority", "semiHidden",
                          "unhideWhenUsed", "qFormat", "locked", "personal", "personalCompose", "personalReply",
                          "rsid", "pPr", "rPr", "tblPr", "trPr", "tcPr", "tblStylePr"],
                "pPr": ["pStyle", "keepNext", "keepLines", "pageBreakBefore", "framePr", "widowControl", "numPr",
                        "suppressLineNumbers", "pBdr", "shd", "tabs", "suppressAutoHyphens", "kinsoku", "wordWrap",
                        "overflowPunct", "topLinePunct", "autoSpaceDE", "autoSpaceDN", "bidi", "adjustRightInd",
                        "snapToGrid", "spacing", "ind", "contextualSpacing", "mirrorIndents", "suppressOverlap", "jc",
                        "textDirection", "textAlignment", "textboxTightWrap", "outlineLvl", "divId", "cnfStyle",
                        "rPr", "sectPr", "pPrChange"],
                "rPr": ["rStyle", "rFonts", "b", "bCs", "i", "iCs", "caps", "smallCaps", "strike", "dstrike",
                        "outline", "shadow", "emboss", "imprint", "noProof", "snapToGrid", "vanish", "webHidden",
                        "color", "spacing", "w", "kern", "position", "sz", "szCs", "highlight", "u", "effect", "bdr",
                        "shd", "fitText", "vertAlign", "rtl", "cs", "em", "lang", "eastAsianLayout", "specVanish",
                        "oMath"],
                "tcPr": TCPR}
            for el in root.iter():
                ln = etree.QName(el).localname if isinstance(el.tag, str) else None
                if ln in ORD:
                    if ln == "rPr":
                        el.text = None
                        for c in el:
                            c.tail = None
                    reorder(el, ORD[ln])
            data = etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)
        if item.filename == "word/numbering.xml":
            data = re.sub(rb'(<w:nsid w:val=")([0-9A-Fa-f]{1,7})(")',
                          lambda m: m.group(1) + m.group(2).zfill(8) + m.group(3), data)
        zout.writestr(item, data)
shutil.move(tmp, OUT)
print("built", OUT, len(d.tables), "tables")
