#!/usr/bin/env python3
"""Build a pandoc reference.docx with Times New Roman 12 pt body text and black headings."""
import subprocess
import sys
import tempfile
import zipfile
from io import BytesIO

from lxml import etree

out = sys.argv[1]
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
w = lambda t: f"{{{W}}}{t}"
FONT, BODY_HALF_PTS = "Times New Roman", "24"

with tempfile.TemporaryDirectory() as d:
    path = f"{d}/default.docx"
    subprocess.run(["pandoc", "-o", path, "--print-default-data-file", "reference.docx"], check=True)
    default = open(path, "rb").read()
zin = zipfile.ZipFile(BytesIO(default))
styles = etree.fromstring(zin.read("word/styles.xml"))

for rf in styles.iter(w("rFonts")):
    for attr in list(rf.attrib):
        del rf.attrib[attr]
    for a in ("ascii", "hAnsi", "eastAsia", "cs"):
        rf.set(w(a), FONT)
for col in styles.iter(w("color")):
    col.set(w("val"), "000000")
    for a in ("themeColor", "themeShade", "themeTint"):
        col.attrib.pop(w(a), None)

defaults = styles.find(f"{w('docDefaults')}/{w('rPrDefault')}/{w('rPr')}")
for tag in ("sz", "szCs"):
    el = defaults.find(w(tag))
    if el is None:
        el = etree.SubElement(defaults, w(tag))
    el.set(w("val"), BODY_HALF_PTS)
if defaults.find(w("rFonts")) is None:
    rf = etree.Element(w("rFonts"))
    for a in ("ascii", "hAnsi", "eastAsia", "cs"):
        rf.set(w(a), FONT)
    defaults.insert(0, rf)

with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == "word/styles.xml":
            data = etree.tostring(styles, xml_declaration=True, encoding="UTF-8", standalone=True)
        zout.writestr(item, data)
print("wrote", out)
