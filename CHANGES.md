# Phase 2 Changes Log

**Applied:** 2026-10-24  
**Branch:** preprint  
**All numbers verified against result CSV files before changing.**

---

## T1-1 · Abstract: LLM Spearman range corrected

- **What:** `$r = 0.25$ to $0.28$` → `$r = 0.20$ to $0.28$`
- **Why:** Verified against `week3_results/results_full.csv`. General encoder
  r = 0.196 and MedCPT r = 0.229 were excluded by the old range.

## T1-2 · Abstract: Cell-type r range added

- **What:** `$r = 0.06$` → `$r = 0.032$ to $0.060$`
- **Why:** Verified against `week3_results/results_celltype.csv`. General
  r = 0.032, MedCPT r = 0.035; reporting a single value was misleading.

## T1-3 · Embedding count: 11 → 12

- **What:** Abstract and intro "11 embedding types" → "12 embedding types";
  abstract now mentions "full-metadata and cell-type-name-only configurations"
  to account for the 4 CTO ablation variants.
- **Why:** Actual count: 4 full LLMs + 4 CTO + TF-IDF + random + Atacformer +
  RP = 12. The CTO variants were silently omitted from the count.

## T1-4 · Methods: "three metrics" → "four metrics"

- **What:** Methods opening paragraph and Fig. 1 caption both changed.
- **Why:** Section II-E (sec:metrics) defines and uses four metrics: Spearman
  r, Mantel test, top-k overlap, and partial Spearman r.

## T1-5 · Section cross-references: decimal notation → \ref{}

- **What:** `Section~2.4` → `Section~\ref{sec:text_emb}` and
  `Section~3.2` → `Section~\ref{sec:partial_corr}`.
  Added `\label{sec:text_emb}` to Text Embeddings subsection and
  `\label{sec:partial_corr}` to "Biomedical Language Models…" subsection.
- **Why:** IEEEtran uses Roman numeral section numbering; the decimal
  references didn't correspond to any real section heading.

## T1-6 · Notation: Δr → r_{partial}

- **What:** All four occurrences of `\Delta r` replaced with
  `r_{\mathrm{partial}}` in Results (lines ~359, ~362) and Discussion (~607,
  ~609). Fig. decomp caption updated accordingly.
- **Why:** Verified against `figures_paper/partial_corr_full.csv`. The values
  (0.051, 0.045, −0.030, −0.042) are partial Spearman correlation
  coefficients, not differences (raw r − anything).

## T1-7 · "Seven-fold" → "approximately eight-fold"

- **What:** Both occurrences changed (Results text and fig:topk_rp caption).
- **Why:** Verified from data: RP top-k at k=5 averages ~21% vs. Atacformer
  ~2.7% for full-text LLMs → ratio ≈ 7.8×, which rounds to ~8×, not 7×.
  (PubMedBERT: 22.5/2.82 = 7.98×; SapBERT: 22.1/2.81 = 7.87×; etc.)

## T2-1 · Clarified "p > 0.20" is the Mantel test

- **What:** Added "Mantel" before `$p > 0.20$` in CTO variant result sentence.
- **Why:** Spearman p-values for CTO variants at full scale are highly
  significant (< 0.001 for 3/4 models); the non-significant result is the
  Mantel permutation test at the cell-type level. Without the qualifier,
  a reader could misinterpret this as an insignificant Spearman correlation.

## T2-2 · Softened causal overclaim

- **What:** Two instances of causal language softened:
  1. Methods (~line 268): "is fully explained by TF-IDF" →
     "is statistically accounted for by TF-IDF"
  2. Discussion (~line 601): "fully explains the Atacformer alignment of…" →
     "statistically accounts for the Atacformer alignment of…"
- **Why:** Partial correlation establishes statistical association after
  partialling, not causal mediation. Both variables may proxy an unmeasured
  third variable.

## T2-3 · Tables II/III: added confounder caveat and (unlabeled) explanation

- **What:** Both table captions now note that examples are illustrative,
  technical confounders have not been ruled out for individual pairs, and
  "(unlabeled)" denotes experiments lacking a cell-type annotation in
  ChIP-Atlas.
- **Why:** Qualitative tables present individual pairs without controlling for
  sequencing depth, peak count, etc. Readers need to know this is illustrative.

## T2-4 · mpnet citation corrected

- **What:** `references.bib` entry `mpnet` replaced. Was: Reimers & Gurevych
  2019 (Sentence-BERT). Now: Song et al. 2020 "MPNet: Masked and permuted
  pre-training for language understanding," NeurIPS.
- **Why:** The model used is `all-mpnet-base-v2`, which is based on MPNet (Song
  et al. 2020), not Sentence-BERT. Sentence-BERT is a training framework for
  sentence embeddings, not the base architecture being cited.

## T2-5 · ATACformer DOI flagged for manual verification

- **What:** Added a `% TODO` comment in `references.bib` before the atacformer
  entry.
- **Why:** The DOI encodes date 2025-11-03. Manual verification that the preprint
  exists, the author list is correct, and the title is final is required before
  submission. (Cannot be verified automatically.)

## T3-1 · RP exon weight value added to methods

- **What:** "with exon-adjacent peaks receiving additional upweighting" →
  "with peaks overlapping annotated exons assigned a fixed weight of 1.0
  (replacing the exponential decay for those peaks)."
- **Why:** Code-verified: `finn/finn/features/regulatory_potential.py` assigns
  weight = 1.0 for exon-overlapping peaks. The decay constant (10 kb) was
  already in the formula; the exon weight was not quantified.

## T3-2 · 34 sentinel cell types: selection criterion added

- **What:** Fig. S1 caption now notes the 34 types were "manually selected to
  span major lineages (immune, neural, stem/progenitor, cancer, and stromal)."
- **Why:** Verified from `sanity_check_all.py`: the CELL_TYPES list and LINEAGE
  dict cover immune adaptive, immune innate, hematopoietic stem, pluripotent
  stem, neural, cancer, endothelial/mesenchymal, and cardiac/metabolic types.

## T3-3 · Random baseline seed added to methods

- **What:** "random Gaussian baseline (192-d, fixed seed)" →
  "(192-dimensional, seed 42)"
- **Why:** Code-verified: `embed_baselines.py` uses `seed=42`. Without the
  specific seed, exact reproduction is impossible.

## T3-4 · Stratified sampling (not applicable to main results)

- **Note:** `build_similarity_matrices.py` contains a stratified 5,000-sample
  pipeline (SAMPLE_N=5000, SAMPLE_SEED=42), but the main reported results
  (Table I) use the full 32,245-experiment matrices. The `results_sample.csv`
  exists but is not cited in the paper. No change needed; noted here for
  completeness.

## T3-5 · (unlabeled) entries explained

- **What:** Covered by T2-3 above (added to both table captions).

---

## Pre-existing warnings (not introduced by these changes)

- `Float too large for page by 38.7162pt` at line 303 (Tables II/III combined
  float) — pre-existing.
- Underfull `\hbox` warnings for long URL strings — pre-existing.

---

# Option C reframe: well-covered primary analysis (2026-10-25)

New analyses: `analysis/strengthen_analyses.py`, `analysis/primary_subset_analysis.py`,
`analysis/make_coverage_figure.py`. Results: `/home/jsong205/atac_out/strengthen/` (see ROBUSTNESS.md).
All new numbers come from those output files.

- **Primary analysis set:** experiments with ≥5,000 peaks (17,232 experiments, 954 cell types).
  All-data (32,245) results kept as a contrast.
- **Abstract, key findings, intro:** rewritten around the peak-count finding. Primary numbers are
  text–Atacformer r = 0.011–0.016 (all data: 0.20–0.28) and text–RP r = 0.06–0.10 with 31–35% of
  top-5 neighbors shared. Promoter fraction is 69–75% below 10K peaks. Titles carry the signal.
- **Methods:** Dataset now describes the peak-count distribution and the ≥5K filter. Added how the
  primary analysis was run, and a new subsection, *Technical Covariates and Sensitivity Analyses*
  (promoter fraction, depth covariates, strata, title-only TF-IDF, uncentered check, block bootstrap).
- **Results:**
  - New III-A, *Peak Count Inflates Apparent Text–Chromatin Agreement*, with a new figure
    (`figures/fig_coverage.*`).
  - III-B now reports well-covered results.
  - Table I is now a two-column-wide table: vs Atacformer and vs RP; cell-type and experiment level
    for ≥5K; all-data experiment level as a contrast.
  - III-C: TF-IDF matches LLMs, and the residual beyond TF-IDF is ≈0 vs Atacformer. The
    full-data PubMedBERT/SapBERT residual is shown with bootstrap CIs and flagged as sensitive to
    coverage and centering. The fine-tuning-objective argument is removed.
  - III-D: top-k overlap and neighborhood purity for well-covered experiments. The interpretation
    that Atacformer captures richer within-cell-type biology is removed. Notes that Atacformer
    neighborhoods are near chance for cell type in well-covered experiments, and gives the untested
    1,024-peak-sampling hypothesis.
- **Figures 4–12 and Tables II/III:** captions now say "All 32,245 experiments (contrast analysis)".
  Tables II/III are reframed as low-peak artifacts (all Table III experiments have ≤1,438 peaks).
- **Discussion:** rewritten. Adds the methodological implication (control for peak count) and new
  limitations (no FRiP/TSS-enrichment QC in ChIP-Atlas, pragmatic threshold, Atacformer
  validity).
- **Not changed:** title (still accurate, but consider revising), figures 4–12 themselves (all-data
  versions).
