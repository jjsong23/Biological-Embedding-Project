# Robustness analyses (2026-10-24)

Code: `/home/jsong205/strengthen_analyses.py` (stages: extract, labels, purity, bootstrap, depth,
centering, subsets, ties). Outputs: `/home/jsong205/atac_out/strengthen/` (file prefixes A–H below).
All analyses reuse the saved 32,245 × 32,245 corpus-centered cosine matrices and the paper's
2M-random-pair design (seed 42). Recomputing the centered cosines from raw embeddings reproduces the
saved matrices to ≤ 5×10⁻⁶ (F_centering_recompute_check.csv), and the 2M-pair sample reproduces
Table I (e.g., PubMedBERT 0.275 vs 0.276; TF-IDF 0.300 vs 0.300; RP 0.215 vs 0.216).

Peak counts = lines in each experiment's ChIP-Atlas BED (q < 1e-10); median 6,287, IQR 780–19,836.

## Headline

1. **Both chromatin references are dominated by peak count.** Atacformer similarity correlates
   r = 0.66 with depth match (−|Δlog10 peaks|) and RP r = 0.59, which is larger than any text–chromatin
   correlation in the paper. Controlling for peak count (depth match + mean depth), Atacformer–RP
   agreement goes from 0.215 to −0.018.
2. **The full-scale text–Atacformer correlation comes from low-coverage experiments.** Restricted to
   pairs where both experiments have ≥ 10,000 peaks (345K pairs), every text model's r with Atacformer
   is ≈ 0.01 (it is 0.20–0.30 on all pairs). The ≥ 5,000-peak subset gives the same result, and so do
   uncentered similarities.
3. **Sample titles carry the text signal.** A TF-IDF on titles alone matches PubMedBERT on full
   metadata (r = 0.275 vs Atacformer). Within a cell type, titles beat every full-text model
   (0.24 vs ~0.21). Titles are SRA/GEO sample names ("GSM…: JH285 S48; Homo sapiens; ATAC-seq",
   "Illumina HiSeq 2500 sequencing of SAMD…"), so they encode provenance (study, platform, naming
   convention), not biology.
4. **Text–RP agreement is the only text–chromatin signal that survives** in high-coverage experiments
   (r ≈ 0.08–0.15, centered or uncentered).
5. **Neighborhood purity supports the paper's RP interpretation:** 72% of RP top-5 neighbors share the
   query's cell type, versus 21% for Atacformer (random 0.7%, label oracle 98%).

## A/B. Cell-type label controls (A_*.csv, B_stratified_spearman.csv)

Partial Spearman vs Atacformer (2M random pairs):

| Model | raw r | \| TF-IDF | \| same cell type + same class | \| labels + title template |
|---|---|---|---|---|
| PubMedBERT | 0.275 | 0.052 | 0.254 | 0.243 |
| SapBERT | 0.246 | 0.042 | 0.222 | 0.210 |
| MedCPT | 0.228 | −0.030 | 0.201 | 0.188 |
| General | 0.196 | −0.042 | 0.163 | 0.150 |
| TF-IDF | 0.300 | — | 0.281 | 0.271 |

Categorical labels explain little of the text–Atacformer correlation. Only 0.65% of pairs share a
cell type and 13.6% share a class. For RP, labels + title template take text–RP partial r from about
0.09 to about 0.03.

Stratified Spearman vs Atacformer:

| Model | same cell type (1M pairs) | same class, diff. cell type (258K) | diff. class (1.35M) |
|---|---|---|---|
| RP | 0.409 | 0.259 | 0.189 |
| PubMedBERT | 0.216 | 0.083 | 0.111 |
| SapBERT | 0.205 | 0.087 | 0.101 |
| MedCPT | 0.204 | 0.054 | 0.100 |
| General | 0.215 | 0.025 | 0.080 |
| TF-IDF | 0.203 | 0.053 | 0.148 |

Name-only models are constant within a cell type, since identical strings give identical embeddings.
Their within-cell-type values (−0.16 to +0.07) are float-rounding noise and should not be reported.

## C. Neighborhood label purity (C_neighborhood_purity.csv)

Fraction of top-k neighbors sharing the query's label (queries whose label has ≥ 2 experiments).

| | cell type k=5 | k=10 | class k=5 | k=10 |
|---|---|---|---|---|
| Random expectation | 0.007 | 0.007 | 0.151 | 0.151 |
| Label oracle | 0.980 | 0.953 | 1.000 | 1.000 |
| Atacformer | 0.208 | 0.191 | 0.379 | 0.363 |
| RP | 0.722 | 0.633 | 0.867 | 0.822 |
| Full-text LLMs | 0.956–0.975 | 0.908–0.946 | 0.991–0.998 | 0.981–0.996 |
| TF-IDF | 0.955 | 0.927 | 0.994 | 0.988 |

## D. Bootstrap CIs for partial correlations (D_partial_bootstrap_ci.csv)

Experiment-level block bootstrap: 1,000 resamples of experiments, 200K pairs each, Atacformer reference.

| Model | raw r [95% CI] | partial \| TF-IDF [95% CI] | partial \| cell type + class [95% CI] |
|---|---|---|---|
| PubMedBERT | 0.276 [0.268, 0.283] | 0.053 [0.048, 0.059] | 0.255 [0.248, 0.263] |
| SapBERT | 0.248 [0.241, 0.255] | 0.044 [0.039, 0.049] | 0.224 [0.217, 0.231] |
| MedCPT | 0.229 [0.222, 0.236] | −0.030 [−0.035, −0.025] | 0.202 [0.195, 0.209] |
| General | 0.196 [0.190, 0.202] | −0.042 [−0.047, −0.037] | 0.163 [0.157, 0.171] |
| TF-IDF | 0.300 [0.293, 0.308] | — | 0.282 [0.274, 0.290] |

On all pairs the sign pattern (PubMedBERT/SapBERT > 0, MedCPT/General < 0) holds under
experiment-level resampling, and the CIs could go straight into Section III-B. However, the pattern
vanishes among high-coverage experiments (G) and depends on centering (F).

## E. Peak-count confound (E_*.csv)

| | r vs depth match | r vs mean depth | partial r vs Atacformer \| depth |
|---|---|---|---|
| Atacformer | 0.659 | 0.539 | — |
| RP | 0.591 | −0.278 | −0.018 |
| PubMedBERT | 0.256 | 0.233 | 0.101 |
| SapBERT | 0.225 | 0.217 | 0.089 |
| MedCPT | 0.213 | 0.183 | 0.087 |
| General | 0.179 | 0.145 | 0.083 |
| TF-IDF | 0.279 | 0.271 | 0.103 |

- Within the same cell type, Atacformer–RP agreement goes from 0.409 to 0.065 after controlling for
  depth.
- Nearest neighbors are depth-matched. Median |Δlog10 peaks| to the top-5 neighbors is 0.20 for
  Atacformer and 0.10 for RP, versus 0.83 for random pairs.
- Within-cell-type similarity tracks peak-count spread. Across 406 cell types with ≥ 10 experiments,
  ρ(within-cell-type similarity, SD of log10 peaks) is −0.77 for Atacformer and −0.71 for RP. Much of
  the "within-cell-type heterogeneity" in Figs. 10–11 may therefore be technical.
- **Table III:** all 10 pairs are two low-peak experiments (≤ 1,438 peaks each, bottom ~30%; corpus
  median 6,287).
- **Table II:** the lead example (NCI-H1299 / CD4+ T cells) is 1,181 vs 1,272 peaks. The other Table II
  pairs are mixed (E_qual_table_peak_counts.csv).

## F. Corpus centering (F_*.csv)

| vs Atacformer | centered r | uncentered r | partial \| TF-IDF, centered | partial \| TF-IDF, uncentered |
|---|---|---|---|---|
| RP | 0.215 | 0.743 | — | — |
| PubMedBERT | 0.275 | 0.292 | 0.052 | 0.257 |
| SapBERT | 0.246 | 0.273 | 0.042 | 0.228 |
| MedCPT | 0.228 | 0.205 | −0.030 | 0.139 |
| General | 0.196 | 0.075 | −0.042 | −0.030 |
| TF-IDF | 0.300 | 0.153 | — | — |

"TF-IDF matches or outperforms all LLMs" holds only with centering. Uncentered, both the main
correlations and the Atacformer–RP agreement still collapse after depth control:

| vs Atacformer, uncentered | r | partial \| depth | r, both ≥ 10K peaks |
|---|---|---|---|
| RP | 0.743 | 0.227 | 0.005 |
| PubMedBERT | 0.292 | 0.033 | −0.003 |
| SapBERT | 0.273 | 0.011 | 0.003 |
| TF-IDF | 0.153 | 0.008 | 0.008 |

## G. High-coverage subsets and title-only TF-IDF (G_*.csv)

Spearman r (centered):

| | vs Atacformer, all | both ≥ 5K | both ≥ 10K | vs RP, all | both ≥ 10K |
|---|---|---|---|---|---|
| n pairs | 2,000,000 | 571,523 | 345,214 | 2,000,000 | 345,214 |
| RP / Atacformer | 0.215 | 0.106 | 0.060 | 0.215 | 0.060 |
| PubMedBERT | 0.275 | 0.009 | 0.007 | 0.090 | 0.080 |
| SapBERT | 0.246 | 0.012 | 0.009 | 0.072 | 0.103 |
| MedCPT | 0.228 | 0.014 | 0.012 | 0.088 | 0.136 |
| General | 0.196 | 0.012 | 0.008 | 0.087 | 0.138 |
| TF-IDF | 0.300 | 0.016 | 0.011 | 0.088 | 0.113 |

In the ≥ 10K subset, partial \| TF-IDF vs Atacformer is −0.001 to 0.006 for all four LLMs, so the
PubMedBERT/SapBERT residual disappears.

Title-only TF-IDF vs Atacformer: r = 0.275 on all pairs (0.092 after depth control) and 0.240 within
a cell type. Controlling for title similarity, within-cell-type text–Atacformer r falls from ~0.21 to
~0.05.

## H. Tie sensitivity of the paper's top-k overlap (H_topk_tie_sensitivity.csv)

Name-only embeddings are exactly constant within a cell type, and 16% of titles are duplicates, so
many similarities are tied. Re-running the paper's `argpartition` top-k reproduces the reported k=5
values exactly (e.g., PubMedBERT 0.0282 vs Atacformer, 0.2253 vs RP; name-only 0.0087–0.0093 and
0.0688–0.0700). Rounding similarities to 4 decimals and breaking ties randomly changes every value by
≤ 0.001. **The top-k results are robust to ties.**

## Caveats

- Peak count partly reflects biology across cell types, but within a cell type it is mostly
  technical (depth, library quality). The ≥ 5K and ≥ 10K thresholds are pragmatic, not standards.
- Partial correlations use linear residualization on ranks. Binary label covariates give only
  first-order control, so the stratified analysis is the cleaner test.
- Title template (digits → "#") is a loose provenance proxy. The two biggest templates are generic
  Illumina boilerplate.

## I. Promoter fraction vs peak count (I_promoter_fraction.csv)

Median fraction of peaks within 1 kb of a RefSeq TSS: <1K peaks 0.69, 1–2K 0.75, 2–5K 0.72,
5–10K 0.69, 10–20K 0.57, 20–50K 0.39, ≥50K 0.24.

## J. Primary analysis on ≥5,000-peak experiments (primary_ge5000/)

17,232 experiments; 954 cell types. Text vs Atacformer: cell-type r −0.074 to −0.023 (Mantel p > 0.98);
experiment r 0.011–0.016; top-5 overlap < 0.1%. Text vs RP: cell-type r 0.038–0.068 (p = 0.001);
experiment r 0.061–0.098; top-5 overlap 31–35%. Atacformer–RP: cell-type r 0.365, experiment r 0.107.
Partial | TF-IDF vs Atacformer: 0.000–0.006. Neighborhood purity (top-5, cell type): Atacformer 1.3%
(chance 0.8%), RP 76%. Atacformer r with depth match: 0.06 (RP 0.57). Peak-count strata are in
peak_count_strata.csv.
