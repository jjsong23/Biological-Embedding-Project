#!/usr/bin/env python3
"""Fix pandoc docx for review: Table 1 column widths and font size, US Letter page with 1-inch margins."""
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

from lxml import etree

src, dst = sys.argv[1:3]
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
M = "http://schemas.openxmlformats.org/officeDocument/2006/math"
w = lambda t: f"{{{W}}}{t}"
m = lambda t: f"{{{M}}}{t}"

# rPr child order (ECMA-376); sz/szCs must be inserted before anything listed after them
AFTER_SZ = ["highlight", "u", "effect", "bdr", "shd", "fitText", "vertAlign", "rtl", "cs", "em", "lang",
            "eastAsianLayout", "specVanish", "oMath"]
HALF_PTS = "17"


def set_size(rpr):
    for tag in ("sz", "szCs"):
        for old in rpr.findall(w(tag)):
            rpr.remove(old)
    later = [c for c in rpr if etree.QName(c).localname in AFTER_SZ]
    pos = list(rpr).index(later[0]) if later else len(rpr)
    for i, tag in enumerate(("sz", "szCs")):
        el = etree.Element(w(tag))
        el.set(w("val"), HALF_PTS)
        rpr.insert(pos + i, el)


tmp = Path(tempfile.mkdtemp())
with zipfile.ZipFile(src) as z:
    z.extractall(tmp)
doc_path = tmp / "word" / "document.xml"
tree = etree.parse(str(doc_path))
root = tree.getroot()

table1 = [t for t in root.iter(w("tbl"))
          if (c := t.find(f"{w('tblPr')}/{w('tblCaption')}")) is not None and c.get(w("val"), "").startswith("Table 1.")][0]
tblpr = table1.find(w("tblPr"))
tblw = tblpr.find(w("tblW"))
tblw.set(w("type"), "dxa")
tblw.set(w("w"), "9360")
layout = etree.Element(w("tblLayout"))
layout.set(w("type"), "fixed")
look = tblpr.find(w("tblLook"))
(look.addprevious(layout) if look is not None else tblpr.append(layout))
grid = table1.find(w("tblGrid"))
for c in list(grid):
    grid.remove(c)
for width in [1560] + [1300] * 6:
    gc = etree.SubElement(grid, w("gridCol"))
    gc.set(w("w"), str(width))

for r in table1.iter(w("r")):
    rpr = r.find(w("rPr"))
    if rpr is None:
        rpr = etree.Element(w("rPr"))
        r.insert(0, rpr)
    set_size(rpr)
body = root.find(w("body"))
sect = body.find(w("sectPr"))
if sect is None:
    sect = etree.SubElement(body, w("sectPr"))
for tag in ("pgSz", "pgMar"):
    for old in sect.findall(w(tag)):
        sect.remove(old)
pgsz = etree.Element(w("pgSz"))
pgsz.set(w("w"), "12240")
pgsz.set(w("h"), "15840")
pgmar = etree.Element(w("pgMar"))
for k, v in (("top", "1440"), ("right", "1440"), ("bottom", "1440"), ("left", "1440"),
             ("header", "720"), ("footer", "720"), ("gutter", "0")):
    pgmar.set(w(k), v)
# pgSz/pgMar follow header/footer references, footnotePr, endnotePr, and type
lead = [c for c in sect if etree.QName(c).localname in
        ("headerReference", "footerReference", "footnotePr", "endnotePr", "type")]
pos = list(sect).index(lead[-1]) + 1 if lead else 0
sect.insert(pos, pgsz)
sect.insert(pos + 1, pgmar)

tree.write(str(doc_path), xml_declaration=True, encoding="UTF-8", standalone=True)
with zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED) as z:
    for f in sorted(tmp.rglob("*")):
        if f.is_file():
            z.write(f, f.relative_to(tmp))
shutil.rmtree(tmp)
print("wrote", dst)
