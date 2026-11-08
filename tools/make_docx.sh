#!/usr/bin/env bash
# Build a Word copy of the manuscript for tracked changes and comments.
# Run from the repository root after a full pdflatex/bibtex build (needs manuscript.aux and manuscript.bbl).
# Requires pandoc and Python with lxml.
set -euo pipefail
out=${1:-manuscript_review.docx}
tmp=$(mktemp -d)
python3 tools/make_reference_docx.py "$tmp/reference.docx"
python3 tools/tex_to_docx_prep.py manuscript.tex manuscript.aux manuscript.bbl "$tmp/prep.tex"
pandoc "$tmp/prep.tex" -f latex -o "$tmp/raw.docx" --number-sections --resource-path=.:figures \
  --reference-doc="$tmp/reference.docx"
python3 tools/fix_docx.py "$tmp/raw.docx" "$out"
