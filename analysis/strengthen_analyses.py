#!/usr/bin/env python3
"""Robustness analyses: label controls, neighborhood purity, partial-corr CIs, depth confound, centering."""

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

SIM = Path("/home/jsong205/atac_out/similarity_matrices")
OUT = Path("/home/jsong205/atac_out/strengthen")
META = Path("/home/jsong205/atac_out/atac_metadata_final.csv")
EMB = Path("/home/jsong205/atac_out")

MODELS = ["atacformer", "rp", "pubmedbert", "sapbert", "medcpt", "general",
          "cto_pubmedbert", "cto_sapbert", "cto_medcpt", "cto_general", "tfidf", "random"]
TEXT = [m for m in MODELS if m not in ("atacformer", "rp")]
LLM = ["pubmedbert", "sapbert", "medcpt", "general"]
EMB_DIRS = {
    "atacformer": EMB / "atacformer_embeddings", "rp": EMB / "rp_embeddings",
    **{m: EMB / "text_embeddings/full" / m for m in LLM},
    "tfidf": EMB / "text_embeddings/tfidf",
}
N_PAIRS = 2_000_000
N_WITHIN = 1_000_000
K_MAX = 10
SEED = 42


def load_meta():
    df = pd.read_csv(META)
    cls = df.cell_type_class.where(df.cell_type_class != "No description")
    tmpl = df.title.fillna("").str.lower().str.replace(r"\d+", "#", regex=True)
    ct_c, _ = pd.factorize(df.cell_type)      # NaN -> -1
    cls_c, _ = pd.factorize(cls)
    tmpl_c, _ = pd.factorize(tmpl)
    return df, ct_c, cls_c, tmpl_c


def same(lab, i, j):
    return ((lab[i] == lab[j]) & (lab[i] >= 0)).astype(np.float32)


def sample_pairs(n, n_pairs, rng):
    i = rng.integers(0, n, int(n_pairs * 1.01) + 10)
    j = rng.integers(0, n, len(i))
    keep = i != j
    i, j = i[keep][:n_pairs], j[keep][:n_pairs]
    return np.minimum(i, j), np.maximum(i, j)


def sample_within_pairs(lab, n_pairs, rng):
    """Uniform sample (with replacement) from all same-label pairs."""
    idx = np.arange(len(lab))
    groups = [idx[lab == g] for g in np.unique(lab[lab >= 0])]
    groups = [g for g in groups if len(g) >= 2]
    sizes = np.array([len(g) for g in groups], float)
    w = sizes * (sizes - 1)
    gsel = rng.choice(len(groups), n_pairs, p=w / w.sum())
    I = np.empty(n_pairs, np.int64)
    J = np.empty(n_pairs, np.int64)
    for gi in np.unique(gsel):
        sel = np.where(gsel == gi)[0]
        mem = groups[gi]
        a = rng.integers(0, len(mem), len(sel))
        b = rng.integers(0, len(mem) - 1, len(sel))
        b = b + (b >= a)
        I[sel], J[sel] = mem[a], mem[b]
    return np.minimum(I, J), np.maximum(I, J)


def topk_neighbors(mat, k, perm, chunk=1024):
    # Columns are visited in a random order so tied similarities are broken randomly.
    n = mat.shape[0]
    inv = np.argsort(perm)
    out = np.empty((n, k), np.int32)
    for s in range(0, n, chunk):
        r = np.arange(s, min(s + chunk, n))
        rows = mat[r][:, perm]
        rows[np.arange(len(r)), inv[r]] = -np.inf
        part = np.argpartition(rows, -k, axis=1)[:, -k:]
        pv = np.take_along_axis(rows, part, 1)
        order = np.argsort(-pv, axis=1, kind="stable")
        out[r] = perm[np.take_along_axis(part, order, 1)]
    return out


def rank(x):
    return rankdata(x).astype(np.float64)


def corr(a, b):
    a = a - a.mean()
    b = b - b.mean()
    return float((a @ b) / np.sqrt((a @ a) * (b @ b)))


def residualize(r, Z):
    if not Z:
        return r - r.mean()
    X = np.column_stack([np.ones(len(r))] + Z)
    beta = np.linalg.lstsq(X, r, rcond=None)[0]
    return r - X @ beta


def partial(rx, ry, Z):
    """Partial Spearman given pre-ranked x, y and list of pre-ranked covariates."""
    return corr(residualize(rx, Z), residualize(ry, Z))


# ── Stage 1: extract sampled pair values + top-k neighbors from each full matrix ──
def stage_extract():
    df, ct, cls, tmpl = load_meta()
    rng = np.random.default_rng(SEED)
    pi, pj = sample_pairs(len(df), N_PAIRS, rng)
    wi, wj = sample_within_pairs(ct, N_WITHIN, rng)
    perm = rng.permutation(len(df))
    vals = {"pi": pi, "pj": pj, "wi": wi, "wj": wj}
    for m in MODELS:
        t0 = time.time()
        mat = np.load(SIM / f"full_{m}.npy")
        vals[f"all_{m}"] = mat[pi, pj]
        vals[f"within_{m}"] = mat[wi, wj]
        np.save(OUT / f"nbrs_{m}.npy", topk_neighbors(mat, K_MAX, perm))
        n_nan = int(np.isnan(vals[f"all_{m}"]).sum())
        del mat
        print(f"[extract] {m}: {time.time() - t0:.0f}s (NaN in sample: {n_nan})", flush=True)
    np.savez(OUT / "pair_values.npz", **vals)


# ── Stage 2: label controls (baselines, stratified, partial) ──
def stage_labels():
    df, ct, cls, tmpl = load_meta()
    v = np.load(OUT / "pair_values.npz")
    pi, pj = v["pi"], v["pj"]
    R = {m: rank(v[f"all_{m}"]) for m in MODELS}
    ind = {"same_celltype": same(ct, pi, pj), "same_class": same(cls, pi, pj),
           "same_title_template": same(tmpl, pi, pj)}
    Ri = {k: rank(x) for k, x in ind.items()}
    print({k: float(x.mean()) for k, x in ind.items()})

    rows = []
    for ref in ("atacformer", "rp"):
        for m in [x for x in MODELS if x != ref] + list(ind):
            rm = R[m] if m in R else Ri[m]
            rows.append({"reference": ref, "condition": m, "spearman_r": corr(R[ref], rm)})
    base = pd.DataFrame(rows)
    base.to_csv(OUT / "A_label_baselines.csv", index=False)

    # how much of each text model's own similarity structure is cell-type labels
    lab_fit = []
    for m in TEXT:
        lab_fit.append({"condition": m,
                        "r_vs_same_celltype": corr(R[m], Ri["same_celltype"]),
                        "r_vs_same_class": corr(R[m], Ri["same_class"]),
                        "r_vs_same_title_template": corr(R[m], Ri["same_title_template"])})
    pd.DataFrame(lab_fit).to_csv(OUT / "A_text_vs_labels.csv", index=False)

    covsets = {
        "tfidf": [R["tfidf"]],
        "labels(ct+class)": [Ri["same_celltype"], Ri["same_class"]],
        "labels+template": [Ri["same_celltype"], Ri["same_class"], Ri["same_title_template"]],
        "tfidf+labels+template": [R["tfidf"], Ri["same_celltype"], Ri["same_class"],
                                  Ri["same_title_template"]],
    }
    prow = []
    for ref in ("atacformer", "rp"):
        for m in TEXT:
            row = {"reference": ref, "condition": m, "raw_r": corr(R[ref], R[m])}
            for name, Z in covsets.items():
                if name.startswith("tfidf") and m == "tfidf":
                    row[f"partial|{name}"] = np.nan
                    continue
                row[f"partial|{name}"] = partial(R[ref], R[m], Z)
            prow.append(row)
    pd.DataFrame(prow).to_csv(OUT / "A_partial_label_controls.csv", index=False)

    # stratified Spearman: same cell type / same class-different cell type / different class
    both_cls = (cls[pi] >= 0) & (cls[pj] >= 0)
    both_ct = (ct[pi] >= 0) & (ct[pj] >= 0)
    strata = {
        "same_celltype": None,
        "same_class_diff_celltype": both_ct & (ind["same_class"] == 1) & (ind["same_celltype"] == 0),
        "diff_class": both_cls & (ind["same_class"] == 0),
    }
    srows = []
    for sname, mask in strata.items():
        for ref in ("atacformer", "rp"):
            for m in [x for x in MODELS if x != ref]:
                if mask is None:
                    a, b, n = v[f"within_{ref}"], v[f"within_{m}"], len(v["wi"])
                else:
                    a, b, n = v[f"all_{ref}"][mask], v[f"all_{m}"][mask], int(mask.sum())
                srows.append({"stratum": sname, "n_pairs": n, "reference": ref, "condition": m,
                              "spearman_r": corr(rank(a), rank(b))})
    pd.DataFrame(srows).to_csv(OUT / "B_stratified_spearman.csv", index=False)
    print("[labels] done")


# ── Stage 3: neighborhood label purity ──
def stage_purity():
    df, ct, cls, tmpl = load_meta()
    n = len(df)
    rows = []
    for labname, lab in (("cell_type", ct), ("cell_type_class", cls)):
        sizes = np.bincount(lab[lab >= 0])
        n_c = np.where(lab >= 0, sizes[np.clip(lab, 0, None)], 0)
        valid = (lab >= 0) & (n_c >= 2)
        for k in (5, 10):
            rows.append({"label": labname, "k": k, "condition": "random_expectation",
                         "purity": float(((n_c[valid] - 1) / (n - 1)).mean())})
            rows.append({"label": labname, "k": k, "condition": "label_oracle",
                         "purity": float((np.minimum(k, n_c[valid] - 1) / k).mean())})
            for m in MODELS:
                nb = np.load(OUT / f"nbrs_{m}.npy")[:, :k]
                hit = lab[nb] == lab[:, None]
                rows.append({"label": labname, "k": k, "condition": m,
                             "purity": float(hit[valid].mean())})
    pd.DataFrame(rows).to_csv(OUT / "C_neighborhood_purity.csv", index=False)
    print("[purity] done")


# ── Stage 4: experiment-level block bootstrap CIs for partial correlations ──
def stage_bootstrap(n_boot=1000, n_pairs=200_000):
    df, ct, cls, tmpl = load_meta()
    n = len(df)
    need = ["atacformer", "tfidf"] + LLM
    mats = {}
    for m in need:
        mats[m] = np.load(SIM / f"full_{m}.npy")
        print(f"[boot] loaded {m}", flush=True)
    rng = np.random.default_rng(SEED + 1)
    rec = []
    t0 = time.time()
    for b in range(n_boot):
        samp = rng.integers(0, n, n)
        a = rng.integers(0, n, int(n_pairs * 1.01))
        c = rng.integers(0, n, len(a))
        I, J = samp[a], samp[c]
        keep = I != J
        I, J = I[keep][:n_pairs], J[keep][:n_pairs]
        R = {m: rank(mats[m][I, J]) for m in need}
        L = [rank(same(ct, I, J)), rank(same(cls, I, J))]
        for m in LLM:
            rec.append({"boot": b, "condition": m,
                        "raw_r": corr(R["atacformer"], R[m]),
                        "partial|tfidf": partial(R["atacformer"], R[m], [R["tfidf"]]),
                        "partial|labels(ct+class)": partial(R["atacformer"], R[m], L)})
        rec.append({"boot": b, "condition": "tfidf", "raw_r": corr(R["atacformer"], R["tfidf"]),
                    "partial|tfidf": np.nan,
                    "partial|labels(ct+class)": partial(R["atacformer"], R["tfidf"], L)})
        if b % 100 == 99:
            print(f"[boot] {b + 1}/{n_boot} ({time.time() - t0:.0f}s)", flush=True)
    bt = pd.DataFrame(rec)
    bt.to_csv(OUT / "D_bootstrap_draws.csv", index=False)
    summ = []
    for m, g in bt.groupby("condition", sort=False):
        for col in ("raw_r", "partial|tfidf", "partial|labels(ct+class)"):
            x = g[col].dropna()
            if len(x) == 0:
                continue
            summ.append({"condition": m, "statistic": col, "boot_mean": x.mean(),
                         "ci_low": x.quantile(0.025), "ci_high": x.quantile(0.975),
                         "frac_boot_le_0": float((x <= 0).mean())})
    pd.DataFrame(summ).to_csv(OUT / "D_partial_bootstrap_ci.csv", index=False)
    print("[boot] done")


# ── Stage 5: peak-count (sequencing depth / coverage) confound ──
def stage_depth():
    df, ct, cls, tmpl = load_meta()
    pk = pd.read_csv(OUT / "peak_counts.csv").set_index("experiment_id").loc[df.experiment_id]
    logp = np.log10(pk.n_peaks.values.astype(float))
    v = np.load(OUT / "pair_values.npz")
    out = []
    for sname, (i, j, pre) in {"all_pairs": (v["pi"], v["pj"], "all"),
                               "same_celltype_pairs": (v["wi"], v["wj"], "within")}.items():
        mismatch = rank(-np.abs(logp[i] - logp[j]))     # higher = more similar depth
        meandepth = rank((logp[i] + logp[j]) / 2)
        R = {m: rank(v[f"{pre}_{m}"]) for m in MODELS}
        for m in MODELS:
            row = {"pairs": sname, "condition": m,
                   "r_vs_depth_match": corr(R[m], mismatch),
                   "r_vs_mean_depth": corr(R[m], meandepth)}
            for ref in ("atacformer", "rp"):
                if m == ref:
                    continue
                row[f"raw_r_vs_{ref}"] = corr(R[ref], R[m])
                row[f"partial_r_vs_{ref}|depth"] = partial(R[ref], R[m], [mismatch, meandepth])
            out.append(row)
    pd.DataFrame(out).to_csv(OUT / "E_depth_confound.csv", index=False)

    # are each model's nearest neighbors depth-matched?
    rng = np.random.default_rng(SEED)
    rnd = np.abs(logp - logp[rng.permutation(len(logp))])
    nrows = [{"condition": "random_pairs", "median_abs_dlog10_peaks": float(np.median(rnd))}]
    for m in MODELS:
        nb = np.load(OUT / f"nbrs_{m}.npy")[:, :5]
        nrows.append({"condition": m,
                      "median_abs_dlog10_peaks": float(np.median(np.abs(logp[nb] - logp[:, None])))})
    pd.DataFrame(nrows).to_csv(OUT / "E_neighbor_depth_match.csv", index=False)

    # per-cell-type: does within-cell-type heterogeneity track depth heterogeneity?
    wv = pd.read_csv(EMB / "week3_results/within_celltype_variation.csv")
    sd = pd.Series(logp).groupby(df.cell_type.values).agg(["std", "size"])
    wv = wv.merge(sd, left_on=wv.columns[0], right_index=True, how="inner")
    wv.to_csv(OUT / "E_within_celltype_vs_depth_sd.csv", index=False)

    # peak counts for the pairs shown in manuscript Tables II/III (first 10 rows of each source CSV)
    npk = pk.n_peaks
    pct = npk.rank(pct=True)
    qrows = []
    for tab, fn in (("Table II (Atacformer)", "table_qual_high_atacformer_low_text.csv"),
                    ("Table III (RP)", "table_qual_high_rp_low_text.csv")):
        q = pd.read_csv(EMB / "figures_paper" / fn).head(10)
        for _, r in q.iterrows():
            e1, e2 = r.experiment_1, r.experiment_2
            qrows.append({"table": tab, "experiment_1": e1, "cell_type_1": r.cell_type_1,
                          "n_peaks_1": int(npk[e1]), "pctile_1": round(float(pct[e1]), 3),
                          "experiment_2": e2, "cell_type_2": r.cell_type_2,
                          "n_peaks_2": int(npk[e2]), "pctile_2": round(float(pct[e2]), 3)})
    pd.DataFrame(qrows).to_csv(OUT / "E_qual_table_peak_counts.csv", index=False)
    print("[depth] done")


# ── Stage 6: corpus-centering robustness ──
def load_emb(model, ids):
    d = EMB_DIRS[model]
    X = np.stack([np.load(d / f"{e}.npy") for e in ids]).astype(np.float32)
    paths = list(d.glob("*.npy"))[:10000]
    mu = np.stack([np.load(p) for p in paths]).astype(np.float32).mean(0)
    return X, mu


def pair_cos(X, i, j, chunk=5000):
    nrm = np.linalg.norm(X, axis=1)
    nrm[nrm == 0] = 1
    out = np.empty(len(i), np.float32)
    for s in range(0, len(i), chunk):
        a, b = X[i[s:s + chunk]], X[j[s:s + chunk]]
        out[s:s + chunk] = (a * b).sum(1) / (nrm[i[s:s + chunk]] * nrm[j[s:s + chunk]])
    return out


def stage_centering():
    df, *_ = load_meta()
    ids = df.experiment_id.tolist()
    v = np.load(OUT / "pair_values.npz")
    i, j = v["pi"], v["pj"]
    unc, chk = {}, []
    for m in EMB_DIRS:
        t0 = time.time()
        X, mu = load_emb(m, ids)
        unc[m] = pair_cos(X, i, j)
        cen = pair_cos(X - mu, i, j)
        chk.append({"condition": m,
                    "max_abs_diff_recomputed_vs_matrix": float(np.nanmax(np.abs(cen - v[f"all_{m}"])))})
        del X
        print(f"[center] {m}: {time.time() - t0:.0f}s", flush=True)
    R = {m: rank(x) for m, x in unc.items()}
    Rc = {m: rank(v[f"all_{m}"]) for m in EMB_DIRS}
    rows = []
    for ref in ("atacformer", "rp"):
        for m in [x for x in EMB_DIRS if x != ref]:
            row = {"reference": ref, "condition": m,
                   "r_centered": corr(Rc[ref], Rc[m]), "r_uncentered": corr(R[ref], R[m])}
            if m in LLM:
                row["partial|tfidf_centered"] = partial(Rc[ref], Rc[m], [Rc["tfidf"]])
                row["partial|tfidf_uncentered"] = partial(R[ref], R[m], [R["tfidf"]])
            rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "F_centering_robustness.csv", index=False)
    pd.DataFrame(chk).to_csv(OUT / "F_centering_recompute_check.csv", index=False)
    np.savez(OUT / "pair_values_uncentered.npz", **unc)

    pk = pd.read_csv(OUT / "peak_counts.csv").set_index("experiment_id").loc[df.experiment_id]
    npk = pk.n_peaks.values
    logp = np.log10(npk.astype(float))
    dep = [rank(-np.abs(logp[i] - logp[j])), rank((logp[i] + logp[j]) / 2)]
    hi = (npk[i] >= 10000) & (npk[j] >= 10000)
    Rh = {m: rank(x[hi]) for m, x in unc.items()}
    rows = []
    for ref in ("atacformer", "rp"):
        for m in [x for x in EMB_DIRS if x != ref]:
            rows.append({"reference": ref, "condition": m,
                         "uncentered_r_all": corr(R[ref], R[m]),
                         "uncentered_partial|depth": partial(R[ref], R[m], dep),
                         "uncentered_r_both>=10000": corr(Rh[ref], Rh[m])})
    pd.DataFrame(rows).to_csv(OUT / "F_uncentered_depth_check.csv", index=False)
    print("[center] done")


# ── Stage 7: high-coverage subsets + title-only TF-IDF ──
def stage_subsets():
    from sklearn.feature_extraction.text import TfidfVectorizer
    df, ct, cls, tmpl = load_meta()
    pk = pd.read_csv(OUT / "peak_counts.csv").set_index("experiment_id").loc[df.experiment_id]
    npk = pk.n_peaks.values
    logp = np.log10(npk.astype(float))
    v = np.load(OUT / "pair_values.npz")
    pi, pj = v["pi"], v["pj"]
    refs = ("atacformer", "rp")
    rows = []
    masks = {"all_pairs": np.ones(len(pi), bool),
             "both_>=5000_peaks": (npk[pi] >= 5000) & (npk[pj] >= 5000),
             "both_>=10000_peaks": (npk[pi] >= 10000) & (npk[pj] >= 10000),
             "depth_matched_|dlog10|<=0.1": np.abs(logp[pi] - logp[pj]) <= 0.1}
    for name, mk in masks.items():
        R = {m: rank(v[f"all_{m}"][mk]) for m in MODELS}
        for ref in refs:
            for m in [x for x in MODELS if x != ref]:
                row = {"subset": name, "n_pairs": int(mk.sum()), "reference": ref, "condition": m,
                       "spearman_r": corr(R[ref], R[m])}
                if m in LLM:
                    row["partial|tfidf"] = partial(R[ref], R[m], [R["tfidf"]])
                rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "G_high_coverage_subsets.csv", index=False)

    titles = df.title.fillna("").tolist()
    X = TfidfVectorizer(max_features=2048, sublinear_tf=True).fit_transform(titles).toarray().astype(np.float32)
    X -= X.mean(0)
    trow = []
    for pre, (i, j) in {"all": (pi, pj), "within": (v["wi"], v["wj"])}.items():
        title_sim = pair_cos(X, i, j)
        rt = rank(title_sim)
        R = {m: rank(v[f"{pre}_{m}"]) for m in ("atacformer", "rp", "tfidf") + tuple(LLM)}
        dep = [rank(-np.abs(logp[i] - logp[j])), rank((logp[i] + logp[j]) / 2)]
        for m in ("tfidf",) + tuple(LLM):
            trow.append({"pairs": pre, "condition": m, "r_with_title_tfidf": corr(R[m], rt)})
        for ref in refs:
            trow.append({"pairs": pre, "condition": f"title_tfidf_vs_{ref}", "spearman_r": corr(R[ref], rt),
                         "partial|depth": partial(R[ref], rt, dep)})
            for m in ("tfidf",) + tuple(LLM):
                trow.append({"pairs": pre, "condition": f"{m}_vs_{ref}", "spearman_r": corr(R[ref], R[m]),
                             "partial|title_tfidf": partial(R[ref], R[m], [rt])})
    pd.DataFrame(trow).to_csv(OUT / "G_title_only_tfidf.csv", index=False)
    print("[subsets] done")


# ── Stage 8: sensitivity of the paper's top-k overlap to tied similarities ──
def topk_sets(mat, k, perm=None, decimals=None, chunk=1024):
    n = mat.shape[0]
    p = np.arange(n) if perm is None else perm
    inv = np.argsort(p)
    out = np.empty((n, k), np.int32)
    for s in range(0, n, chunk):
        r = np.arange(s, min(s + chunk, n))
        rows = mat[r][:, p]
        if decimals is not None:
            rows = np.round(rows, decimals)
        rows[np.arange(len(r)), inv[r]] = -np.inf
        out[r] = p[np.argpartition(rows, -k, axis=1)[:, -k:]]
    return out


def overlap(a, b):
    k = a.shape[1]
    return float(np.mean([len(set(x) & set(y)) / k for x, y in zip(a, b)]))


def stage_ties(k=5):
    rng = np.random.default_rng(SEED + 2)
    ref_nb = {}
    for ref in ("atacformer", "rp"):
        mat = np.load(SIM / f"full_{ref}.npy")
        ref_nb[ref] = topk_sets(mat, k)
        del mat
    rows = []
    v = np.load(OUT / "pair_values.npz")
    for m in TEXT:
        if m == "random":
            continue
        t0 = time.time()
        mat = np.load(SIM / f"full_{m}.npy")
        paper = topk_sets(mat, k)
        robust = topk_sets(mat, k, perm=rng.permutation(mat.shape[0]), decimals=4)
        del mat
        w = v[f"within_{m}"]
        for ref in ("atacformer", "rp"):
            rows.append({"condition": m, "reference": ref, f"topk{k}_paper_method": overlap(ref_nb[ref], paper),
                         f"topk{k}_rounded_random_ties": overlap(ref_nb[ref], robust),
                         "within_celltype_sim_std": float(w.std())})
        print(f"[ties] {m}: {time.time() - t0:.0f}s", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "H_topk_tie_sensitivity.csv", index=False)
    print("[ties] done")


# ── Stage 9: promoter fraction of peaks vs peak count ──
REFSEQ = Path("/home/jsong205/.venv_atac/lib/python3.12/site-packages/lisa/genomes/hg38.refseq")
_TSS = None


def _load_tss():
    rs = pd.read_csv(REFSEQ, sep="\t", header=None, usecols=[1, 2, 3, 4])
    rs.columns = ["chrom", "strand", "start", "end"]
    rs["tss"] = np.where(rs.strand == "+", rs.start, rs.end)
    return {c: np.unique(g.tss.values) for c, g in rs.groupby("chrom")}


def _promoter_frac(path, window=1000):
    global _TSS
    if _TSS is None:
        _TSS = _load_tss()
    try:
        b = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], comment="#")
    except pd.errors.EmptyDataError:
        return np.nan
    if len(b) == 0:
        return np.nan
    mid = ((b[1].values + b[2].values) // 2)
    hit = np.zeros(len(b), bool)
    for c, idx in b.groupby(0).indices.items():
        t = _TSS.get(c)
        if t is None:
            continue
        m = mid[idx]
        k = np.searchsorted(t, m)
        d = np.minimum(np.abs(m - t[np.clip(k - 1, 0, len(t) - 1)]), np.abs(t[np.clip(k, 0, len(t) - 1)] - m))
        hit[idx] = d <= window
    return float(hit.mean())


def stage_promoter():
    from multiprocessing import Pool
    df, *_ = load_meta()
    with Pool(16) as p:
        fr = p.map(_promoter_frac, df.bed_path.tolist(), chunksize=64)
    pk = pd.read_csv(OUT / "peak_counts.csv").set_index("experiment_id").loc[df.experiment_id]
    out = pd.DataFrame({"experiment_id": df.experiment_id, "n_peaks": pk.n_peaks.values, "promoter_frac": fr})
    out.to_csv(OUT / "I_promoter_fraction.csv", index=False)
    bins = [0, 1000, 2000, 5000, 10000, 20000, 50000, 1e9]
    out["bin"] = pd.cut(out.n_peaks, bins, right=False)
    print(out.groupby("bin", observed=True).promoter_frac.describe()[["count", "25%", "50%", "75%"]].round(3))
    print("[promoter] done")


STAGES = {"promoter": stage_promoter, "extract": stage_extract, "labels": stage_labels, "purity": stage_purity,
          "bootstrap": stage_bootstrap, "depth": stage_depth, "centering": stage_centering,
          "subsets": stage_subsets, "ties": stage_ties}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stages", nargs="+", choices=list(STAGES))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    for s in args.stages:
        STAGES[s]()
