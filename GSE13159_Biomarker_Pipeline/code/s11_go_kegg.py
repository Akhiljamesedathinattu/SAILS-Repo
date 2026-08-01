#!/usr/bin/env python3
"""
STEP 6 - GO and KEGG enrichment

Two complementary tests are run for every cluster (or diagnosis) and every
library:

  1. Over-representation (ORA), hypergeometric. Asks: among my significant
     genes, is this pathway present more often than chance? Needs a threshold,
     and the BACKGROUND is the filtered gene universe from step 3 - never the
     whole genome. Using a genome-wide background against a filtered foreground
     is the single most common way to manufacture enrichment.

  2. Rank-based competitive test, threshold-free. Asks: are this pathway's genes
     shifted towards the top of my ranked t-statistic list relative to all other
     genes? Implemented as a Mann-Whitney U on ranks with a normal
     approximation. Reports direction (up or down) and an AUC effect size.

Report both. Agreement between them is much stronger evidence than either alone,
and the rank test rescues pathways that are coherently shifted but never cross
an arbitrary fold-change cutoff.

Prerequisites
    python3 fetch_genesets.py          # puts GMT libraries in raw/genesets/
    python3 05_differential_expression.py

Outputs
  results/enrichment_ora_<library>.csv
  results/enrichment_rank_<library>.csv
  results/enrichment_top.csv           consolidated, for the R dotplot
  results/enrichment_summary.csv       one row per library and group
"""
import argparse
import re

import numpy as np
import pandas as pd
from scipy import stats

from common import RES_DIR, GENESET_DIR, load_matrix, bh_fdr, ensure_dirs, log

GO_ID = re.compile(r"\(?(GO:\d{7})\)?")
KEGG_ID = re.compile(r"\b(hsa\d{5})\b")


def read_gmt(path):
    """
    Parse a GMT / Enrichr library.

    Handles both layouts, because Enrichr omits the description column that the
    Broad GMT spec includes:
        term <TAB> description <TAB> gene <TAB> gene ...
        term <TAB> gene <TAB> gene ...
    A token is treated as a description (and dropped) if it contains whitespace
    or looks like a URL. Gene symbols never do. Trailing ",1.0" weights, which
    Enrichr sometimes appends, are stripped.
    """
    sets = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").rstrip("\t").split("\t")
            if len(parts) < 2:
                continue
            term, genes = parts[0].strip(), set()
            for tok in parts[1:]:
                tok = tok.strip()
                if not tok:
                    continue
                if " " in tok or tok.lower().startswith("http"):
                    continue
                g = tok.split(",")[0].strip().upper()
                if g:
                    genes.add(g)
            if term and len(genes) >= 2:
                sets[term] = genes
    return sets


def parse_term(term, library):
    """Split an Enrichr term name into a clean label plus an ontology id."""
    tid, label = "", term
    m = GO_ID.search(term)
    if m:
        tid = m.group(1)
        label = GO_ID.sub("", term).strip(" ()")
    else:
        m = KEGG_ID.search(term)
        if m:
            tid = m.group(1)
            label = KEGG_ID.sub("", term).strip(" ()")

    lib = library.lower()
    if "biological_process" in lib:
        domain = "GO:BP"
    elif "molecular_function" in lib:
        domain = "GO:MF"
    elif "cellular_component" in lib:
        domain = "GO:CC"
    elif "kegg" in lib:
        domain = "KEGG"
    elif "reactome" in lib:
        domain = "Reactome"
    elif "hallmark" in lib:
        domain = "Hallmark"
    else:
        domain = library
    return label.strip(), tid, domain


def ora(query, background, sets, min_size, max_size):
    """Hypergeometric over-representation against a filtered background."""
    bg = set(background)
    fg = set(query) & bg
    N, n = len(bg), len(fg)
    if n < 5:
        return pd.DataFrame()
    rows = []
    for term, genes in sets.items():
        gs = genes & bg
        K = len(gs)
        if not (min_size <= K <= max_size):
            continue
        hits = gs & fg
        x = len(hits)
        if x < 2:
            continue
        expected = n * K / N
        rows.append({
            "term": term, "n_overlap": x, "n_term_in_bg": K,
            "n_query": n, "n_background": N,
            "expected": expected,
            "fold_enrichment": x / expected if expected else np.nan,
            "p_value": stats.hypergeom.sf(x - 1, N, K, n),
            "genes": ";".join(sorted(hits)[:40]),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df["fdr"] = bh_fdr(df.p_value.to_numpy())
    return df


def rank_test(stat, genes, sets, min_size, max_size):
    """
    Threshold-free competitive test. Mann-Whitney U of set members against all
    other genes, on the ranks of the per-gene statistic, normal approximation.
    """
    genes = np.asarray([str(g).upper() for g in genes])
    ranks = stats.rankdata(np.asarray(stat, dtype=float))
    pos = {g: i for i, g in enumerate(genes)}
    n = len(genes)
    rows = []
    for term, members in sets.items():
        idx = np.fromiter((pos[g] for g in members if g in pos), dtype=np.int64)
        K = idx.size
        if not (min_size <= K <= max_size) or K >= n:
            continue
        R = float(ranks[idx].sum())
        U = R - K * (K + 1) / 2.0
        mu = K * (n - K) / 2.0
        sigma = np.sqrt(K * (n - K) * (n + 1) / 12.0)
        if sigma == 0:
            continue
        z = (U - mu) / sigma
        rows.append({
            "term": term, "n_term_in_data": K,
            "auc": U / (K * (n - K)),
            "z_score": z,
            "direction": "up" if z > 0 else "down",
            "p_value": 2.0 * stats.norm.sf(abs(z)),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df["fdr"] = bh_fdr(df.p_value.to_numpy())
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", default="kmeans_cluster",
                    help="matches the grouping used in step 05")
    ap.add_argument("--fdr-in", type=float, default=0.05,
                    help="FDR cutoff defining the ORA query set")
    ap.add_argument("--lfc-in", type=float, default=0.2154,
                    help="min |mean difference| for the ORA query set. NOTE: the "
                         "data is 0-1 scaled, not log2 — a cutoff of 1.0 selects "
                         "nothing. Must match the cutoff used in step 10.")
    ap.add_argument("--min-size", type=int, default=10)
    ap.add_argument("--max-size", type=int, default=500)
    ap.add_argument("--fdr-out", type=float, default=0.05,
                    help="FDR cutoff for reporting a term as enriched")
    ap.add_argument("--top", type=int, default=10,
                    help="terms per group carried into enrichment_top.csv")
    ap.add_argument("--skip-rank-test", action="store_true")
    args = ap.parse_args()
    ensure_dirs()

    GENESET_DIR.mkdir(parents=True, exist_ok=True)
    gmts = sorted(p for p in GENESET_DIR.iterdir()
                  if p.suffix.lower() in (".gmt", ".txt"))
    if not gmts:
        log(f"no gene set libraries in {GENESET_DIR}")
        log("run:  python3 fetch_genesets.py")
        return

    de_file = RES_DIR / "10_de_results.csv"
    if not de_file.exists():
        raise FileNotFoundError(f"{de_file} not found - run step 05 first")
    de = pd.read_csv(de_file)
    de["feature_u"] = de.gene.astype(str).str.upper()

    background = sorted({str(g).upper()
                         for g in load_matrix("expr_genes").index})
    log(f"background universe: {len(background)} filtered genes")
    if len(background) < 2000:
        log("warning: background looks small - if you are still at probe level, "
            "run optional_annotate_probes.R and repeat step 03")

    ora_all, rank_all, summary = [], [], []

    for gmt in gmts:
        sets = read_gmt(gmt)
        lib = gmt.stem
        sizes = np.array([len(v) for v in sets.values()])
        log(f"{lib}: {len(sets)} sets, median size {int(np.median(sizes))}")
        overlap = len(set().union(*sets.values()) & set(background)) if sets else 0
        log(f"  {overlap} library genes present in the background")
        if overlap < 500:
            log("  warning: poor symbol overlap - check gene identifier types")

        for g, sub in de.groupby("group"):
            query = sub.loc[(sub.fdr < args.fdr_in) &
                            (sub.log2FC >= args.lfc_in), "feature_u"].tolist()

            res = ora(query, background, sets, args.min_size, args.max_size)
            n_ora = 0
            if not res.empty:
                res.insert(0, "group", g)
                res.insert(0, "library", lib)
                parsed = res.term.map(lambda t: parse_term(t, lib))
                res["label"] = [p[0] for p in parsed]
                res["term_id"] = [p[1] for p in parsed]
                res["domain"] = [p[2] for p in parsed]
                ora_all.append(res)
                n_ora = int((res.fdr < args.fdr_out).sum())

            n_rank = 0
            if not args.skip_rank_test:
                rres = rank_test(sub.t_stat.to_numpy(), sub.feature_u.to_numpy(),
                                 sets, args.min_size, args.max_size)
                if not rres.empty:
                    rres.insert(0, "group", g)
                    rres.insert(0, "library", lib)
                    parsed = rres.term.map(lambda t: parse_term(t, lib))
                    rres["label"] = [p[0] for p in parsed]
                    rres["term_id"] = [p[1] for p in parsed]
                    rres["domain"] = [p[2] for p in parsed]
                    rank_all.append(rres)
                    n_rank = int((rres.fdr < args.fdr_out).sum())

            summary.append({"library": lib, "group": g, "n_query_genes": len(query),
                            "n_ora_significant": n_ora,
                            "n_rank_significant": n_rank})
            log(f"  group {str(g)[:26]:<28} query={len(query):>5}  "
                f"ORA={n_ora:>4}  rank={n_rank:>4}")

    if not ora_all and not rank_all:
        log("no enrichment results produced")
        return

    pd.DataFrame(summary).to_csv(RES_DIR / "enrichment_summary.csv", index=False)

    cols_keep = ["library", "domain", "group", "label", "term_id", "term"]
    top_frames = []

    if ora_all:
        full = pd.concat(ora_all, ignore_index=True)
        for lib, sub in full.groupby("library"):
            sub.sort_values(["group", "p_value"]).to_csv(
                RES_DIR / f"enrichment_ora_{lib}.csv", index=False)
        sig = full[full.fdr < args.fdr_out].copy()
        sig["test"] = "ora"
        sig["score"] = -np.log10(sig.p_value.clip(lower=1e-300))
        sig["effect"] = sig.fold_enrichment
        top_frames.append(
            sig.sort_values(["library", "group", "p_value"])
               .groupby(["library", "group"]).head(args.top)
               [cols_keep + ["test", "score", "effect", "fdr",
                             "n_overlap", "n_term_in_bg"]])

    if rank_all:
        full = pd.concat(rank_all, ignore_index=True)
        for lib, sub in full.groupby("library"):
            sub.sort_values(["group", "p_value"]).to_csv(
                RES_DIR / f"enrichment_rank_{lib}.csv", index=False)
        sig = full[(full.fdr < args.fdr_out) & (full.direction == "up")].copy()
        sig["test"] = "rank"
        sig["score"] = -np.log10(sig.p_value.clip(lower=1e-300))
        sig["effect"] = sig.auc
        sig["n_overlap"] = sig.n_term_in_data
        sig["n_term_in_bg"] = sig.n_term_in_data
        top_frames.append(
            sig.sort_values(["library", "group", "p_value"])
               .groupby(["library", "group"]).head(args.top)
               [cols_keep + ["test", "score", "effect", "fdr",
                             "n_overlap", "n_term_in_bg"]])

    if top_frames:
        top = pd.concat(top_frames, ignore_index=True)
        top.to_csv(RES_DIR / "enrichment_top.csv", index=False)
        log(f"consolidated -> {RES_DIR/'enrichment_top.csv'} ({len(top)} rows)")

        both = (top.groupby(["library", "group", "label"])["test"]
                   .nunique().reset_index())
        n_both = int((both.test > 1).sum())
        log(f"{n_both} terms are significant under BOTH tests - lead with these")


if __name__ == "__main__":
    main()
