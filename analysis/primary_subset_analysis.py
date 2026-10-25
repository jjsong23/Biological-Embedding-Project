#!/usr/bin/env python3
"""Re-run the paper's analyses on well-covered experiments (>= MIN_PEAKS peaks) plus a peak-count-stratified table."""

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from strengthen_analyses import (EMB, LLM, MODELS, OUT, SIM, corr, load_meta, partial,
                                 rank, sample_pairs, topk_neighbors)

MIN_PEAKS = 5000
N_BOOT = 1000
N_BOOT_PAIRS = 50_000
N_PERM = 1000
KS = (5, 10, 20, 50)
SEED = 42
DIRS = {
    "atacformer": EMB / "atacformer_embeddings", "rp": EMB / "rp_embeddings",
    **{m: EMB / "text_embeddings/full" / m for m in LLM},
    **{f"cto_{m}": EMB / "text_embeddings/cell_type_only" / m for m in LLM},
    "tfidf": EMB / "text_embeddings/tfidf", "random": EMB / "text_embeddings/random",
}
P = OUT / "primary_ge5000"


def well_covered(df):
    pk = pd.read_csv(OUT / "peak_counts.csv").set_index("experiment_id").loc[df.experiment_id]
    return pk.n_peaks.values >= MIN_PEAKS, pk.n_peaks.values


def upper(mat):
    return mat[np.triu_indices(mat.shape[0], 1)]


def mantel_p(ref_mat, oth_mat, rng, n_perm=N_PERM):
    """One-sided permutation p for Spearman between upper triangles (row/col labels of `oth` permuted)."""
    n = ref_mat.shape[0]
    iu = np.triu_indices(n, 1)
    rr = rank(ref_mat[iu])
    ro = np.zeros((n, n))
    ro[iu] = rank(oth_mat[iu])
    ro = ro + ro.T
    obs = corr(rr, ro[iu])
    cnt = 0
    for _ in range(n_perm):
        p = rng.permutation(n)
        cnt += corr(rr, ro[np.ix_(p, p)][iu]) >= obs
    return obs, (cnt + 1) / (n_perm + 1)


def boot_ci(mats, ref, n, rng, n_boot=N_BOOT, n_pairs=N_BOOT_PAIRS):
    out = {m: [] for m in mats if m != ref}
    for _ in range(n_boot):
        s = rng.integers(0, n, n)
        a, c = s[rng.integers(0, n, n_pairs)], s[rng.integers(0, n, n_pairs)]
        k = a != c
        a, c = a[k], c[k]
        rr = rank(mats[ref][a, c])
        for m in out:
            out[m].append(corr(rr, rank(mats[m][a, c])))
    return {m: (float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))) for m, v in out.items()}


def stage_celltype():
    df, ct, _, _ = load_meta()
    keep, _ = well_covered(df)
    sub = df[keep & df.cell_type.notna()]
    cts = sorted(sub.cell_type.unique())
    mats = {}
    for m, d in DIRS.items():
        t0 = time.time()
        paths = list(d.glob("*.npy"))[:10000]
        mu = np.stack([np.load(p) for p in paths]).astype(np.float32).mean(0)
        means = np.stack([np.stack([np.load(d / f"{e}.npy") for e in g.experiment_id]).astype(np.float32).mean(0)
                          for _, g in sub.groupby("cell_type", sort=True)])
        X = means - mu
        X /= np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
        mats[m] = (X @ X.T).astype(np.float32)
        print(f"[celltype] {m}: {time.time() - t0:.0f}s", flush=True)
    rng = np.random.default_rng(SEED)
    rows = []
    for ref in ("atacformer", "rp"):
        ci = boot_ci(mats, ref, len(cts), rng)
        for m in [x for x in mats if x != ref]:
            r, p = mantel_p(mats[ref], mats[m], rng)
            rows.append({"reference": ref, "condition": m, "n_celltypes": len(cts), "spearman_r": r,
                         "mantel_p": p, "ci_low": ci[m][0], "ci_high": ci[m][1]})
        print(f"[celltype] {ref} done", flush=True)
    pd.DataFrame(rows).to_csv(P / "celltype_results.csv", index=False)


def stage_full():
    df, ct, cls, tmpl = load_meta()
    keep, npk = well_covered(df)
    idx = np.where(keep)[0]
    n = len(idx)
    rng = np.random.default_rng(SEED)
    pi, pj = sample_pairs(n, 2_000_000, rng)
    perm = rng.permutation(n)
    mats, vals, nbrs = {}, {}, {}
    for m in MODELS:
        t0 = time.time()
        full = np.load(SIM / f"full_{m}.npy")
        mats[m] = full[np.ix_(idx, idx)]
        del full
        vals[m] = mats[m][pi, pj]
        nbrs[m] = topk_neighbors(mats[m], max(KS), perm)
        print(f"[full] {m}: {time.time() - t0:.0f}s", flush=True)
    R = {m: rank(v) for m, v in vals.items()}
    rows = []
    for ref in ("atacformer", "rp"):
        ci = boot_ci(mats, ref, n, rng)
        for m in [x for x in MODELS if x != ref]:
            row = {"reference": ref, "condition": m, "n_experiments": n, "spearman_r": corr(R[ref], R[m]),
                   "ci_low": ci[m][0], "ci_high": ci[m][1]}
            for k in KS:
                row[f"topk_{k}"] = float(np.mean([len(set(a[:k]) & set(b[:k])) / k
                                                  for a, b in zip(nbrs[ref], nbrs[m])]))
            if m in LLM:
                row["partial|tfidf"] = partial(R[ref], R[m], [R["tfidf"]])
                row["r_with_tfidf"] = corr(R[m], R["tfidf"])
                row["cross_topk5_with_tfidf"] = float(np.mean([len(set(a[:5]) & set(b[:5])) / 5
                                                               for a, b in zip(nbrs[m], nbrs["tfidf"])]))
            rows.append(row)
        print(f"[full] {ref} done", flush=True)
    pd.DataFrame(rows).to_csv(P / "full_results.csv", index=False)

    # partial-correlation bootstrap CIs (experiment-level), Atacformer reference
    pb = {m: [] for m in LLM}
    for _ in range(N_BOOT):
        s = rng.integers(0, n, n)
        a, c = s[rng.integers(0, n, 200_000)], s[rng.integers(0, n, 200_000)]
        k = a != c
        a, c = a[k], c[k]
        ra, rt = rank(mats["atacformer"][a, c]), rank(mats["tfidf"][a, c])
        for m in LLM:
            pb[m].append(partial(ra, rank(mats[m][a, c]), [rt]))
    pd.DataFrame([{"condition": m, "ci_low": np.quantile(v, 0.025), "ci_high": np.quantile(v, 0.975),
                   "frac_le_0": float(np.mean(np.array(v) <= 0))} for m, v in pb.items()]
                 ).to_csv(P / "partial_bootstrap_ci.csv", index=False)

    # within-cell-type self-similarity (Atacformer vs RP) and neighborhood purity
    lab = ct[idx]
    wrows = []
    for g in np.unique(lab[lab >= 0]):
        mem = np.where(lab == g)[0]
        if len(mem) < 2:
            continue
        iu = np.triu_indices(len(mem), 1)
        wrows.append({"cell_type": df.cell_type.values[idx[mem[0]]], "n": len(mem),
                      "atacformer": float(mats["atacformer"][np.ix_(mem, mem)][iu].mean()),
                      "rp": float(mats["rp"][np.ix_(mem, mem)][iu].mean())})
    w = pd.DataFrame(wrows)
    w.to_csv(P / "within_celltype.csv", index=False)
    sizes = np.bincount(lab[lab >= 0])
    n_c = np.where(lab >= 0, sizes[np.clip(lab, 0, None)], 0)
    valid = (lab >= 0) & (n_c >= 2)
    pur = {m: float((lab[nbrs[m][:, :5]] == lab[:, None])[valid].mean()) for m in MODELS}
    pur["random_expectation"] = float(((n_c[valid] - 1) / (n - 1)).mean())
    logp = np.log10(npk[idx].astype(float))
    dm = rank(-np.abs(logp[pi] - logp[pj]))
    summ = {"n_experiments": n, "n_celltypes_ge2": len(w),
            "within_median_atacformer": float(w.atacformer.median()), "within_median_rp": float(w.rp.median()),
            "within_frac_atac_lt_rp": float((w.atacformer < w.rp).mean()),
            "wilcoxon_p": float(wilcoxon(w.atacformer, w.rp).pvalue),
            "atac_r_vs_depth_match": corr(R["atacformer"], dm), "rp_r_vs_depth_match": corr(R["rp"], dm),
            **{f"purity5_{m}": v for m, v in pur.items()}}
    pd.Series(summ).to_csv(P / "summary.csv", header=["value"])
    print("[full] done")


def stage_strata():
    df, *_ = load_meta()
    _, npk = well_covered(df)
    v = np.load(OUT / "pair_values.npz")
    pi, pj = v["pi"], v["pj"]
    bins = [(0, 1000), (1000, 5000), (5000, 10000), (10000, 10**9)]
    b = np.full(len(npk), -1)
    for k, (lo, hi) in enumerate(bins):
        b[(npk >= lo) & (npk < hi)] = k
    rows = []
    for k, (lo, hi) in enumerate(bins):
        mk = (b[pi] == k) & (b[pj] == k)
        R = {m: rank(v[f"all_{m}"][mk]) for m in MODELS}
        for ref in ("atacformer", "rp"):
            for m in [x for x in MODELS if x != ref]:
                rows.append({"bin": f"{lo}-{hi if hi < 10**9 else 'inf'}", "n_pairs": int(mk.sum()),
                             "n_experiments": int((b == k).sum()), "reference": ref, "condition": m,
                             "spearman_r": corr(R[ref], R[m])})
    pd.DataFrame(rows).to_csv(P / "peak_count_strata.csv", index=False)
    print("[strata] done")


STAGES = {"celltype": stage_celltype, "full": stage_full, "strata": stage_strata}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stages", nargs="+", choices=list(STAGES))
    P.mkdir(parents=True, exist_ok=True)
    for s in ap.parse_args().stages:
        STAGES[s]()
