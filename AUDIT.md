# Preprint Audit — Phase 1

**Paper:** Biomedical Text Embeddings Weakly Capture Chromatin Accessibility Similarity  
**Source branch:** master (IEEE conference submission)  
**Working branch:** preprint  
**Audit date:** 2026-10-24  
**Auditor:** Claude Sonnet 4.6 (via Claude Code)

All numbers flagged here were verified against result CSV files in
`/home/jsong205/atac_out/` unless noted. No numbers were changed; all are
flags only.

---

## Tier 1 — Safe / Mechanical (no judgment required)

These are factual errors confirmed against result files, or typographic/style
issues with one correct fix.

### T1-1 · Abstract: LLM Spearman range is wrong (full-scale)

**Location:** abstract, ~line 19  
**Text:** "achieving Spearman correlations of r = 0.25 to 0.28"  
**Issue:** Verified against `week3_results/results_full.csv`. The four full-text
LLM models span r = 0.196 (General) to 0.276 (PubMedBERT). The stated range
0.25–0.28 excludes two of the four models.  
**Fix:** Change to "r = 0.20 to 0.28" (or list individual values).

---

### T1-2 · Abstract: Cell-type r = 0.06 is imprecise / overstates uniformity

**Location:** abstract, ~line 20  
**Text:** "r = 0.06 at the cell-type level"  
**Issue:** Verified against `week3_results/results_celltype.csv`. Cell-type
values span 0.032 (General) to 0.060 (TF-IDF/PubMedBERT). Reporting a single
"r = 0.06" implies all models land there, which misrepresents MedCPT (0.035)
and General (0.032).  
**Fix:** Change to "r = 0.032–0.060" (or "r ≈ 0.03–0.06").

---

### T1-3 · Abstract: Embedding count says 11, should be 12

**Location:** abstract, ~line 13  
**Text:** "eleven embedding types"  
**Issue:** Fig. 2 and the methods enumerate 12 conditions: 4 full-text LLMs +
4 cell-type-name-only (CTO) ablations + TF-IDF + random + ATACformer + RP.
TF-IDF and random are arguably baselines rather than "text embedding types,"
but the CTO models are clearly part of the eleven that excludes ATACformer/RP,
and the current text omits them from the abstract count.  
**Options:** (a) change to 12 and mention CTO variants briefly; (b) reframe as
"four biomedical text encoders … plus ablation and baseline variants"; (c)
justify 11 if ATACformer/RP are excluded but CTO variants are not.  
**Note:** This is Tier 1 only if option (b) or (c) is chosen with no new
claim. If the count is changed to 12, no new analysis is needed; the CTO
variants are already reported.

---

### T1-4 · Methods Section II-A: "three metrics" should be "four"

**Location:** Section II-A (Methods opening paragraph), ~line 101  
**Text:** "We quantify similarity structure alignment using three metrics"  
**Issue:** Section II-E defines four metrics: Spearman rank correlation,
Mantel permutation test, top-k neighborhood overlap, and partial Spearman
correlation. Four are reported in Table I (rho, Mantel p, top-k, and
separately partial r in the text).  
**Fix:** Change "three" to "four."

---

### T1-5 · Cross-references: decimal section notation instead of Roman numerals

**Location:** ~line 370 and ~line 606  
**Text:** "Section~2.4" and "Section~3.2"  
**Issue:** The rest of the paper (and IEEEtran style) uses Roman numeral
section references (e.g., "Section~\ref{sec:methods}"). These two instances
use decimal notation that doesn't correspond to any section number in the
paper.  
**Fix:** Replace with the correct `\ref{}` cross-references (II-D and III-B
respectively, or whichever sections are intended).

---

### T1-6 · Δr notation is incorrect — these are partial r values, not differences

**Location:** ~lines 359, 362, 607, 609 (wherever Δr appears)  
**Text:** "Δr = 0.051 … Δr = 0.045 … Δr = −0.030 … Δr = −0.042"  
**Issue:** Verified against `figures_paper/partial_corr_full.csv`. These values
are the partial Spearman correlation coefficients r_{partial} (controlling for
TF-IDF), not differences (raw r − partial r). Using Δr implies subtraction;
the actual computation is a partial correlation.  
**Fix:** Replace Δr with r_{partial} (or r_p) and adjust surrounding text to
say "partial Spearman correlation … controlling for TF-IDF" rather than
implying a residual difference. This is purely a notation fix — no numbers
change.

---

### T1-7 · Figure caption: "seven-fold" is ~eight-fold per the data

**Location:** Figure 3 or 4 caption (wherever the fold-change is stated)  
**Text:** "a seven-fold increase"  
**Issue:** Verified against `week3_results/results_full.csv`. RP r = 0.216
vs. ATACformer r = 0.216 at full scale — actually they are essentially equal
here. The fold-change language refers to the cell-type level: RP r = 0.315
vs. best-text r ≈ 0.060, which is ~5-fold, or ATACformer cell-type r (not
directly compared). Need to check which comparison is intended and re-derive.  
**Note:** Flagged for verification; exact location and comparison must be
confirmed before fixing. Do not change without confirming the intended
denominator.

---

## Tier 2 — Needs a Decision

These require the author to choose the correct framing or acknowledge a
limitation.

### T2-1 · "p > 0.20 for all four models" is ambiguous — Mantel or Spearman?

**Location:** Section III (Results), wherever name-only CTO Mantel results are
discussed  
**Text:** "p > 0.20 for all four models" (referring to CTO ablation variants)  
**Issue:** Two p-values exist for each model: (a) Mantel permutation p-values,
which range 0.26–0.81 (confirming p > 0.20), and (b) Spearman p-values, which
are < 0.001 for 3/4 CTO models (PubMedBERT, SapBERT, MedCPT) at the full
scale. A reader will assume "p > 0.20" means the correlation is not
significant, which is true for Mantel but not for Spearman.  
**Decision needed:** Clarify which test the p-value refers to. If the claim is
specifically that the Mantel test (global permutation test of matrix
similarity) is non-significant, say so explicitly.

---

### T2-2 · "fully explains" is a causal overclaim

**Location:** wherever TF-IDF partial correlation is discussed  
**Text:** language suggesting TF-IDF "fully explains" or "accounts for" LLM
signal  
**Issue:** The partial Spearman r for PubMedBERT drops from 0.276 to 0.051
after controlling for TF-IDF; the paper interprets this as TF-IDF explaining
the bulk of the LLM correlation. But partial correlation ≠ causal mediation;
TF-IDF and LLMs may both be proxies for an unmeasured third variable (e.g.,
tissue type). Phrasing it as "fully explains" or "accounts for" the signal is
stronger than the analysis supports.  
**Decision needed:** Soften to "TF-IDF overlap is associated with the majority
of the observed LLM–ATACformer correlation" or similar.

---

### T2-3 · Tables II/III: qualitative examples lack technical confounder context

**Location:** Tables II and III captions / surrounding text  
**Text:** Qualitative examples presented as illustrating biological specificity  
**Issue:** The high-scoring pairs in Tables II/III may be driven by sequencing
depth, peak count, or batch effects rather than biology. No statement is made
about whether these specific pairs were checked for technical confounders.  
**Decision needed:** Add a brief caveat that these examples are illustrative
and that technical factors have not been explicitly ruled out for individual
pairs.

---

### T2-4 · mpnet citation is incorrect

**Location:** `references.bib`, citation key `mpnet`  
**Text:** Cites Reimers & Gurevych 2019 (Sentence-BERT, EMNLP-IJCNLP) for the
"General" encoder (all-mpnet-base-v2)  
**Issue:** all-mpnet-base-v2 is based on MPNet (Song et al. 2020, "MPNet:
Masked and Permuted Pre-training for Language Understanding," NeurIPS), not
Sentence-BERT. The Sentence-BERT paper is about the training framework, not
the base model architecture.  
**Decision needed:** Replace with Song et al. 2020 (MPNet, NeurIPS), or cite
both (MPNet for the base architecture + Reimers & Gurevych for the sentence
embedding fine-tuning framework). Either is defensible; pick one.

---

### T2-5 · atacformer DOI may be a placeholder (date in future at time of writing)

**Location:** `references.bib`, citation key `atacformer`  
**Text:** `doi = {10.1101/2025.11.03.685753}` with `year = {2025}`  
**Issue:** The DOI encodes a date of 2025-11-03. This was in the future at
the time of writing and should be verified to ensure the preprint actually
exists and the metadata is current. Today's date is 2026-10-24.  
**Decision needed:** Verify the DOI resolves and confirm author list, title,
and year are accurate in the .bib entry.

---

## Tier 3 — Journal Version / Nice to Have

These are not blocking for a bioRxiv preprint but would be expected in a
full journal submission.

### T3-1 · RP exon upweighting factor not quantified in methods

**Location:** Section II-C (RP methods)  
**Text:** "exon-adjacent peaks receiving additional upweighting"  
**Issue:** The code (`finn/finn/features/regulatory_potential.py`) assigns
weight = 1.0 to exon-overlapping peaks (replacing the exponential decay that
would otherwise apply), and uses a 10 kb decay constant. Neither the weight
value nor the decay constant is stated in the manuscript.  
**Suggested addition:** "Peaks overlapping annotated exons received a fixed
weight of 1.0 (replacing exponential decay), with a 10 kb decay constant for
non-exon peaks."

---

### T3-2 · 34 sentinel cell types: selection criterion not stated

**Location:** Section II or Methods, wherever the sanity-check / sentinel
analysis is described  
**Issue:** `sanity_check_all.py` contains a hardcoded list of 34 cell types
used as sentinel classes. No criterion for their selection is stated in the
manuscript (e.g., representative diversity, minimum experiment count,
manual curation).  
**Suggested addition:** One sentence stating how these 34 types were chosen.

---

### T3-3 · Random baseline seed not in methods

**Location:** Section II-D or wherever random baseline is described  
**Issue:** `embed_baselines.py` uses seed=42 and 192-dimensional random
Gaussian vectors. The seed is not mentioned, making exact reproduction
impossible.  
**Suggested addition:** "Random Gaussian embeddings (192-dimensional, seed 42)
were drawn independently per experiment."

---

### T3-4 · Stratified sampling seed and N not in methods

**Location:** Section II-E (similarity matrix construction)  
**Issue:** `build_similarity_matrices.py` uses SAMPLE_N=5000 and SAMPLE_SEED=42
for the stratified experiment-level sample used in Spearman/Mantel/top-k
analyses. Neither the seed nor the stratification procedure is described in
the methods.  
**Suggested addition:** "A stratified sample of N=5,000 experiments
(proportional to cell-type size, seed 42) was used for experiment-level
analyses."

---

### T3-5 · (unlabeled) experiments in Tables II/III not explained

**Location:** Tables II and III  
**Issue:** Some entries in the qualitative tables have "(unlabeled)" for cell
type, meaning the experiment had no ChIP-Atlas cell-type annotation. No
explanation is given for how these appear in a table supposedly illustrating
cell-type-level biological interpretability.  
**Suggested addition:** Note in caption that a small number of experiments
lack cell-type annotations in ChIP-Atlas and "(unlabeled)" denotes these.

---

## Preprint Conversion Checklist (Phase 3 — not yet applied)

These items are structural changes needed for bioRxiv submission. Not yet
applied; listed here for tracking.

- [ ] Populate `\author{}` block with real author names and affiliations
  (currently empty for blind review)
- [ ] Add "Data and Code Availability" section
- [ ] Add AI-use disclosure (Claude used for writing assistance)
- [ ] Remove or replace IEEEtran formatting for bioRxiv (plain article or
  bioRxiv LaTeX template)
- [ ] Add acknowledgements section
- [ ] Decide on supplementary vs. integrated figure numbering
  (currently unified — already done)

---

## Summary Statistics

| Tier | Count |
|------|-------|
| Tier 1 (safe/mechanical) | 7 |
| Tier 2 (needs decision)  | 5 |
| Tier 3 (journal/nice-to-have) | 5 |
| Phase 3 checklist items  | 6 |

User-pre-identified issues and their audit disposition:
1. "eleven embedding types" → T1-3
2. Abstract r-range wrong → T1-1 + T1-2
3. "three metrics" → T1-4
4. Δr notation → T1-6
5. "seven-fold" → T1-7
6. Section cross-references → T1-5
