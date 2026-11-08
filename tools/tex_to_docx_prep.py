#!/usr/bin/env python3
"""Prepare the bioRxiv manuscript.tex for pandoc -> docx: resolve refs/cites from .aux, number captions, flatten layout."""
import re
import sys

tex_path, aux_path, bbl_path, out_path = sys.argv[1:5]
s = open(tex_path).read()
aux = open(aux_path).read()
bbl = open(bbl_path).read()


def brace_group(text, start):
    """Return (content, end_index) for the balanced {...} group starting at text[start] == '{'."""
    assert text[start] == "{"
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start + 1:i], i + 1
    raise ValueError("unbalanced braces")


labels = {}
for m in re.finditer(r"\\newlabel\{([^}]*)\}\{\{([^}]*)\}", aux):
    labels[m.group(1)] = m.group(2)
cites = {m.group(1): m.group(2) for m in re.finditer(r"\\bibcite\{([^}]*)\}\{\{(\d+)\}", aux)}

# ── front matter (authors and affiliations read from the .tex) ──
title = re.search(r"\\title\{\\textbf\{(.*?)\}\}", s, re.S).group(1).replace("\n", " ")
abstract = s[s.index("\\begin{abstract}"):s.index("\\end{abstract}") + len("\\end{abstract}")]
authors = re.findall(r"\\author\[([^\]]*)\]\{([^}]*)\}", s)
names = [f"{n}\\textsuperscript{{{mark}}}" for mark, n in authors]
author_line = ", ".join(names[:-1]) + ", and " + names[-1] if len(names) > 1 else names[0]
affil_lines = []
for m in re.finditer(r"\\affil\[([^\]]*)\]", s):
    text, _ = brace_group(s, m.end())
    text = " ".join(text.split())
    mark = m.group(1).strip()
    affil_lines.append((f"\\textsuperscript{{{mark}}}" if mark else "") + text)
front = (f"\\title{{{title}}}\n\\begin{{document}}\n\\maketitle\n" + author_line + "\n\n"
         + "\n\n".join(affil_lines) + "\n\n"
         + abstract.replace("\\begin{abstract}", "\\section*{Abstract}").replace("\\end{abstract}", "") + "\n\n")

# main text runs to the bibliography; supplementary material (if any) follows it
SUPP_MARK = "\n%%SUPPLEMENT%%\n"
main = s[s.index("\\noindent\\textbf{Keywords:}"):s.index("\\bibliographystyle")]
supp = ""
if "\\section*{Supplementary" in s:
    supp = s[s.index("\\section*{Supplementary"):s.index("\\end{document}")]
    supp = re.sub(r"(?m)^\\(setcounter|renewcommand)\{.*\n", "", supp)
body = main + SUPP_MARK + supp
body = body.replace("\\resizebox{\\textwidth}{!}{%", "").replace("\\end{tabular}}", "\\end{tabular}")

# ── split side-by-side minipages into separate floats ──
def split_minipages(text, env):
    out, pos = [], 0
    pat = re.compile(r"\\begin\{" + env + r"\}\[!htbp\](.*?)\\end\{" + env + r"\}", re.S)
    for m in pat.finditer(text):
        out.append(text[pos:m.start()])
        inner = m.group(1)
        if "\\begin{minipage}" in inner:
            for mp in re.findall(r"\\begin\{minipage\}\[t\]\{[^}]*\}(.*?)\\end\{minipage\}", inner, re.S):
                out.append(f"\\begin{{{env}}}[!htbp]\n\\centering{mp}\\end{{{env}}}\n\n")
        else:
            out.append(m.group(0))
        pos = m.end()
    out.append(text[pos:])
    return "".join(out)


body = split_minipages(body, "figure")
body = split_minipages(body, "table")

# ── prefix captions with their numbers ──
def number_captions(text):
    out, pos = [], 0
    for m in re.finditer(r"\\begin\{(figure|table)\}.*?\\end\{\1\}", text, re.S):
        out.append(text[pos:m.start()])
        block = m.group(0)
        lab = re.search(r"\\label\{([^}]*)\}", block).group(1)
        kind = "Figure" if m.group(1) == "figure" else "Table"
        i = block.index("\\caption{")
        cap, end = brace_group(block, i + len("\\caption"))
        block = block[:i] + f"\\caption{{\\textbf{{{kind} {labels[lab]}.}} {cap}}}" + block[end:]
        out.append(block)
        pos = m.end()
    out.append(text[pos:])
    return "".join(out)


body = number_captions(body)


def table_math_to_text(text):
    """Inside tabular environments, turn simple inline math into plain text so Word cells hold ordinary runs."""
    def conv(mm):
        t = mm.group(1)
        t = t.replace("\\geq", "≥").replace("-", "−")
        t = re.sub(r"\^\{([^}]*)\}", r"\\textsuperscript{\1}", t)
        return t
    def fix_tab(m):
        return re.sub(r"\$([^$]*)\$", conv, m.group(0))
    return re.sub(r"\\begin\{tabular\}.*?\\end\{tabular\}", fix_tab, text, flags=re.S)


body = table_math_to_text(body)
body = re.sub(r"\$\^\{([^}]*)\}\s*(.*?)\$", lambda m: "\\textsuperscript{" + m.group(1) + "}" + (("$" + m.group(2) + "$") if m.group(2).strip() else ""), body)

# ── Table 1: flatten the two-row grouped header into one header row ──
hdr = re.search(r"\\toprule\n(.*?)\\midrule", body, re.S)
if "multicolumn" in hdr.group(1):
    flat = ("    Embedding & Atacformer, cell type (≥5K) & Atacformer, experiment (≥5K) [95\\% CI] "
            "& Atacformer, experiment (full) & RP, cell type (≥5K) & RP, experiment (≥5K) [95\\% CI] "
            "& RP, experiment (full) \\\\\n")
    body = body[:hdr.start(1)] + flat + body[hdr.end(1):]
body = body.replace("\\end{tabular}\n  \\\\[3pt]\n  {\\footnotesize", "\\end{tabular}\n\n  {\\footnotesize")

# ── cross-references and citations ──
body = re.sub(r"\\ref\{([^}]*)\}", lambda m: labels[m.group(1)], body)
body = re.sub(r"\\cite\{([^}]*)\}", lambda m: "[" + ", ".join(cites[k.strip()] for k in m.group(1).split(",")) + "]", body)
body = re.sub(r"\\label\{[^}]*\}", "", body)
body = body.replace("\\nolinkurl{", "\\texttt{")

# ── reference list from the .bbl, in citation order ──
items = re.split(r"\\bibitem(?:\[[^\]]*\])?\{[^}]*\}", bbl.split("\\end{thebibliography}")[0])[1:]
refs = ["\\section*{References}"]
for n, it in enumerate(items, 1):
    t = it.replace("\\newblock", " ").replace("\\penalty0", "").replace("~", " ")
    t = re.sub(r"\\doi\{([^}]*)\}", r"doi: \1", t)
    t = re.sub(r"\\natexlab\{([^}]*)\}", r"\1", t)
    t = re.sub(r"\s+", " ", t).strip()
    refs.append(f"[{n}] {t}\n")

main_out, supp_out = body.split(SUPP_MARK)
out = front + main_out + "\n" + "\n\n".join(refs) + "\n\n" + supp_out + "\n\\end{document}\n"
open(out_path, "w").write(out)
print("figures:", len(re.findall(r"\\begin\{figure\}", out)), "tables:", len(re.findall(r"\\begin\{table\}", out)),
      "refs:", len(items))
