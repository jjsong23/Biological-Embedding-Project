#!/usr/bin/env python3
"""Figure: promoter fraction vs peak count, and text-chromatin Spearman by peak-count stratum."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

S = Path("/home/jsong205/atac_out/strengthen")
OUT = Path("/home/jsong205/atac_submission/figures")
plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.labelsize": 9,
                     "axes.titlesize": 9, "axes.spines.top": False, "axes.spines.right": False})

MODELS = [("pubmedbert", "PubMedBERT", "#E63946", "o"), ("sapbert", "SapBERT", "#F4A261", "s"),
          ("medcpt", "MedCPT", "#457B9D", "^"), ("general", "General", "#9B5DE5", "D"),
          ("tfidf", "TF-IDF", "#555555", "v")]

pf = pd.read_csv(S / "I_promoter_fraction.csv").dropna()
edges = np.logspace(1, np.log10(pf.n_peaks.max() + 1), 25)
pf["b"] = pd.cut(pf.n_peaks, edges)
g = pf.groupby("b", observed=True).promoter_frac
mid = np.array([iv.mid for iv in g.median().index])
cnt = g.size().values
q = g.quantile([0.25, 0.5, 0.75]).unstack().values[cnt >= 30]
mid = mid[cnt >= 30]

st = pd.read_csv(S / "primary_ge5000/peak_count_strata.csv")
bins = list(dict.fromkeys(st.bin))
labels = ["<1K", "1–5K", "5–10K", "≥10K"]
x = np.arange(len(bins))

fig, ax = plt.subplots(1, 3, figsize=(7.16, 2.4), gridspec_kw={"width_ratios": [1.1, 1, 1]})
a = ax[0]
a.fill_between(mid, q[:, 0], q[:, 2], color="#457B9D", alpha=0.25, lw=0)
a.plot(mid, q[:, 1], color="#1D3557", lw=2)
a.axvline(5000, color="#888888", lw=1, ls="--")
a.text(5600, 0.93, "5,000 peaks", fontsize=7, color="#555555")
a.set_xscale("log")
a.set_ylim(0, 1)
a.set_xlabel("Peaks per experiment")
a.set_ylabel("Peaks ≤1 kb from TSS (frac.)")
a.set_title("(a) Promoter fraction", loc="left")

for k, (ref, title) in enumerate((("atacformer", "(b) Text vs. Atacformer"), ("rp", "(c) Text vs. RP"))):
    a = ax[k + 1]
    for key, name, col, mk in MODELS:
        y = [st[(st.bin == b) & (st.reference == ref) & (st.condition == key)].spearman_r.iloc[0] for b in bins]
        a.plot(x, y, color=col, lw=2, marker=mk, ms=6, mec="#222222", mew=0.6, label=name)
    a.axhline(0, color="#999999", lw=0.8)
    a.set_xticks(x, labels)
    a.set_ylim(-0.02, 0.45)
    a.set_xlabel("Peak-count stratum")
    a.set_title(title, loc="left")
    a.grid(axis="y", color="#e5e5e5", lw=0.6)
ax[1].set_ylabel("Spearman ρ")
ax[1].legend(frameon=False, fontsize=7, loc="upper right")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"fig_coverage.{ext}", dpi=300)
print("saved")
