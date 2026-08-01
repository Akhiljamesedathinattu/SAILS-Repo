#!/usr/bin/env python3
"""
PIPELINE STEPS 8-10 : Hierarchical Clustering -> Dendrogram -> Cluster Assignment

Ward linkage on the z-scored filtered matrix, exactly as your pipeline specifies
(clustering on the expression matrix itself, not on principal components).

Two additions inside step 10 that your pipeline implies but does not name, and
which a reviewer will ask for:

  a stated rule for choosing k   silhouette and consensus stability across a
                                 range of cuts, so the number of clusters is a
                                 decision with a criterion behind it rather
                                 than an eyeball judgement
  a stability check              consensus clustering over resampled subsets.
                                 PAC = proportion of ambiguous clustering, the
                                 fraction of sample pairs that are neither
                                 reliably together nor reliably apart. Lower is
                                 better. Without this you cannot tell a real
                                 subgroup from an arbitrary cut of a continuum.

CHOOSING k — why the default rule is PAC, not silhouette.

  Silhouette rewards compact, well-separated, roughly spherical clusters. High-
  dimensional expression data almost never has that shape, so silhouette values
  are low across the board (on GSE13159, 0.06-0.20 for every k tested) and the
  maximum lands on whichever k merges the data into the fewest large blobs.
  On GSE13159 it picks k=3 for 2096 samples spanning 17 diagnostic categories,
  which collapses distinct leukaemia lineages together.

  PAC asks a different question: when the cohort is resampled, do the same
  samples keep landing together? That is stability, not compactness, and it is
  the property that matters for "is this a real subgroup". On GSE13159 PAC
  falls sharply to a minimum at k=5 and rises again at k=6.

  So --k-rule defaults to "pac". Both metrics are always computed and both are
  written to 05_clustering_summary.csv, and when they disagree the script says
  so loudly — the disagreement is a finding to report, not something to hide.
  Use --k-rule silhouette or --k N to override.

k-means on the same matrix is also run, purely as a second opinion.

Outputs
  results/05_k_selection.csv, 05_consensus_pac.csv
  results/05_cluster_assignments.csv
  results/05_consensus_matrix.csv, 05_consensus_matrix_anno.csv
  results/05_clustering_summary.csv
  results/dendro_{merge,height,order,labels}.csv
"""
import argparse

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import pdist
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score

from common import RES_DIR, SEED, load_matrix, ensure_dirs, log


def export_dendrogram(Z, labels, tag=""):
    """Convert a scipy linkage matrix into R hclust components."""
    n = len(labels)
    merge = np.zeros((n - 1, 2), dtype=np.int64)
    for i in range(n - 1):
        pair = [(-(int(x) + 1) if int(x) < n else int(x) - n + 1) for x in Z[i, :2]]
        pair.sort(key=lambda v: (v > 0, abs(v)))
        merge[i] = pair
    pd.DataFrame(merge, columns=["a", "b"]).to_csv(
        RES_DIR / f"dendro_merge{tag}.csv", index=False)
    pd.DataFrame({"height": Z[:, 2]}).to_csv(
        RES_DIR / f"dendro_height{tag}.csv", index=False)
    pd.DataFrame({"order": leaves_list(Z) + 1}).to_csv(
        RES_DIR / f"dendro_order{tag}.csv", index=False)
    pd.DataFrame({"label": labels}).to_csv(
        RES_DIR / f"dendro_labels{tag}.csv", index=False)


def consensus_matrix(M, k, reps, frac, seed=SEED):
    n = M.shape[0]
    co = np.zeros((n, n), dtype=np.float32)
    cnt = np.zeros((n, n), dtype=np.float32)
    rng = np.random.default_rng(seed + k)
    m = max(int(round(frac * n)), k + 1)
    for _ in range(reps):
        idx = rng.choice(n, m, replace=False)
        sub = M[idx]
        lab = fcluster(linkage(sub, method="ward"), t=k, criterion="maxclust")
        same = (lab[:, None] == lab[None, :]).astype(np.float32)
        co[np.ix_(idx, idx)] += same
        cnt[np.ix_(idx, idx)] += 1.0
    with np.errstate(invalid="ignore", divide="ignore"):
        C = np.where(cnt > 0, co / cnt, 0.0)
    np.fill_diagonal(C, 1.0)
    return C


def per_sample_consensus(C, labels):
    """Mean co-clustering of each sample with the OTHER members of its cluster.

    The diagonal of C is 1 by construction, so including self inflates the score
    — badly for small clusters, where self is a large share of the mean. A
    singleton cluster has no other members and scores NaN, which is the honest
    answer rather than a perfect 1.0.
    """
    out = np.full(len(labels), np.nan)
    for i in range(len(labels)):
        peers = np.where(labels == labels[i])[0]
        peers = peers[peers != i]
        if peers.size:
            out[i] = float(C[i, peers].mean())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmin", type=int, default=2)
    ap.add_argument("--kmax", type=int, default=20)
    ap.add_argument("--k", type=int, default=0,
                    help="explicit k. 0 = decide with --k-rule")
    ap.add_argument("--k-rule", choices=["pac", "silhouette"], default="pac",
                    help="criterion when --k is 0. Default pac (stability); "
                         "see the CHOOSING k note in this file.")
    ap.add_argument("--pac-kmin", type=int, default=3,
                    help="ignore k below this when minimising PAC. k=2 is often "
                         "trivially stable because splitting anything in two is "
                         "reproducible without being informative.")
    ap.add_argument("--consensus-kmax", type=int, default=12)
    ap.add_argument("--consensus-reps", type=int, default=20)
    ap.add_argument("--consensus-frac", type=float, default=0.8)
    ap.add_argument("--dendro-n", type=int, default=400)
    ap.add_argument("--heatmap-n", type=int, default=300)
    args = ap.parse_args()
    ensure_dirs()

    z = load_matrix("expr_zscore")
    samples = z.columns.to_numpy()
    M = np.ascontiguousarray(z.values.T, dtype=np.float64)   # samples x probes
    log(f"hierarchical clustering {M.shape[0]} samples on {M.shape[1]} probes")

    # STEP 8-9: linkage on the full cohort, then the dendrogram export
    D = pdist(M, metric="euclidean")
    Z = linkage(D, method="ward")
    log(f"linkage built (cophenetic distances: {D.size} pairs)")

    rng = np.random.default_rng(SEED)
    sub = np.sort(rng.choice(len(samples),
                             min(args.dendro_n, len(samples)), replace=False))
    export_dendrogram(linkage(M[sub], method="ward"), samples[sub])
    log(f"dendrogram exported for {len(sub)} samples (a 2000-leaf tree is unreadable)")

    # STEP 10a: silhouette scan
    rows = []
    for k in range(args.kmin, args.kmax + 1):
        lab = fcluster(Z, t=k, criterion="maxclust")
        if len(np.unique(lab)) < 2:
            continue
        sil = silhouette_score(M, lab, sample_size=min(2000, len(samples)),
                               random_state=SEED)
        rows.append({"k": k, "n_clusters_realised": len(np.unique(lab)),
                     "silhouette": sil})
        log(f"  k={k:>2}  silhouette={sil:.3f}")
    scan = pd.DataFrame(rows)
    scan.to_csv(RES_DIR / "05_k_selection.csv", index=False)

    # STEP 10b: consensus stability scan
    pac_rows = []
    for k in range(args.kmin, min(args.consensus_kmax, args.kmax) + 1):
        C = consensus_matrix(M, k, args.consensus_reps, args.consensus_frac)
        v = C[np.triu_indices_from(C, k=1)]
        pac = float(np.mean((v > 0.1) & (v < 0.9)))
        pac_rows.append({"k": k, "pac": pac})
        log(f"  consensus k={k:>2}  PAC={pac:.3f}")
    pac_df = pd.DataFrame(pac_rows)
    pac_df.to_csv(RES_DIR / "05_consensus_pac.csv", index=False)

    # STEP 10c: decide k
    k_sil = int(scan.loc[scan.silhouette.idxmax(), "k"])
    eligible = pac_df[pac_df.k >= args.pac_kmin]
    if eligible.empty:
        eligible = pac_df
    k_pac = int(eligible.loc[eligible.pac.idxmin(), "k"])

    if args.k:
        k_final, rule = args.k, "user specified"
    elif args.k_rule == "silhouette":
        k_final, rule = k_sil, "maximum mean silhouette"
    else:
        k_final, rule = k_pac, f"minimum PAC over k>={args.pac_kmin}"

    log(f"k = {k_final} ({rule})")
    log(f"  silhouette optimum k={k_sil} (sil={scan.silhouette.max():.3f})")
    log(f"  stability optimum  k={k_pac} (PAC={eligible.pac.min():.3f})")
    if k_sil != k_pac:
        log("  NOTE: the two criteria disagree. Silhouette favours few compact "
            "clusters; PAC favours reproducible ones. Report both and say which "
            "you used and why — this disagreement belongs in the write-up.")

    # STEP 10d: assign
    hclust_lab = fcluster(Z, t=k_final, criterion="maxclust")
    km_lab = KMeans(n_clusters=k_final, n_init=25,
                    random_state=SEED).fit_predict(M) + 1

    ari = adjusted_rand_score(hclust_lab, km_lab)
    log(f"hclust vs k-means agreement (ARI): {ari:.3f}"
        + ("  — low agreement, the partition is not method-robust"
           if ari < 0.5 else ""))

    C = consensus_matrix(M, k_final, args.consensus_reps, args.consensus_frac)
    order = np.argsort(hclust_lab, kind="mergesort")
    hn = min(args.heatmap_n, len(samples))
    pick = np.sort(order[np.linspace(0, len(order) - 1, hn).astype(int)])
    pd.DataFrame(C[np.ix_(pick, pick)], index=samples[pick],
                 columns=samples[pick]).round(3) \
        .to_csv(RES_DIR / "05_consensus_matrix.csv")
    pd.DataFrame({"sample": samples[pick], "cluster": hclust_lab[pick]}) \
        .to_csv(RES_DIR / "05_consensus_matrix_anno.csv", index=False)

    cons = per_sample_consensus(C, hclust_lab)
    out = pd.DataFrame({
        "sample": samples,
        "hclust_cluster": hclust_lab,
        "kmeans_cluster": km_lab,
        "cluster_consensus": cons,
    })
    out.to_csv(RES_DIR / "05_cluster_assignments.csv", index=False)

    sizes = out.hclust_cluster.value_counts().sort_index()
    log("cluster sizes: " + ", ".join(f"{i}:{n}" for i, n in sizes.items()))

    n_weak = int(np.nansum(cons < 0.6))
    log(f"per-sample consensus: mean {np.nanmean(cons):.3f}, "
        f"{n_weak} of {len(cons)} samples below 0.6")
    if n_weak > 0.25 * len(cons):
        log("  WARNING: over a quarter of samples cluster ambiguously. Treat "
            "these clusters as provisional.")

    pd.DataFrame({"k_final": [k_final], "k_rule": [rule],
                  "k_min_pac": [k_pac], "k_max_silhouette": [k_sil],
                  "criteria_agree": [k_sil == k_pac],
                  "silhouette": [scan.loc[scan.k == k_final, "silhouette"].iloc[0]],
                  "pac": [pac_df.loc[pac_df.k == k_final, "pac"].iloc[0]
                          if (pac_df.k == k_final).any() else np.nan],
                  "hclust_kmeans_ari": [ari],
                  "mean_consensus": [float(np.nanmean(cons))],
                  "n_below_0.6_consensus": [n_weak]}) \
        .to_csv(RES_DIR / "05_clustering_summary.csv", index=False)
    log("clusters are NOT compared to diagnosis yet — that is step 15")


if __name__ == "__main__":
    main()