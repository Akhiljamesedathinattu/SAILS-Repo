#!/usr/bin/env python3
"""
PIPELINE STEP 20 : Candidate Biomarker Identification

Filters the differential expression output into a defensible candidate list.
A good diagnostic biomarker is not simply the gene with the smallest p-value; it
needs to be:

  significant   FDR below the cutoff
  large         normalised mean difference above the cutoff
  specific      up in ONE group, not in many. A gene up in six of eight groups
                separates nothing. Specificity is scored as the effect size in
                its own group divided by the next-highest group's effect size.
  expressed     mean expression above the cohort median, so the assay would
                actually have signal to measure

Genes are ranked by a composite of significance, effect size and specificity.

EFFECT-SIZE SCALE

  The DE step reports `mean_diff`, a difference of group means on the [0,1]
  rescaled matrix. It is NOT a log2 fold change. Earlier versions of this
  script read a `log2FC` column that duplicated those values under a
  misleading name; that column has been removed upstream.

  The cutoff is no longer hardcoded here. It is read from
  results/10_effect_size_cutoff.csv, which the DE step writes as the single
  authoritative record, so the two steps cannot drift apart. Pass --lfc to
  override.

SPECIFICITY FLOOR — check this before reporting specificity values.

  Specificity divides a gene's own effect size by the next-best group's. When
  the next-best is near zero the ratio explodes, so the denominator is floored.
  The default floor of 0.1 was chosen for a log2 fold-change scale, where it is
  small relative to typical effects.

  On the [0,1] scale it is not small: the entire significant range is roughly
  0.2 to 0.6, so for many genes the floor IS the denominator and specificity
  degenerates into the effect size divided by a constant. The script now
  reports what fraction of candidates are floored. If that fraction is high,
  specificity is not measuring what its name claims and should either be
  rescaled (try --specificity-floor 0.02) or dropped from the composite score.

  The default is left at 0.1 so that existing results are reproduced exactly.
  Changing it changes the candidate ranking.

Outputs
  results/12_candidate_biomarkers.csv
  results/12_candidate_summary.csv
"""
import argparse

import numpy as np
import pandas as pd

from common import RES_DIR, load_matrix, ensure_dirs, log


def resolve_cutoff(explicit):
    """Prefer an explicit --lfc; otherwise read the DE step's recorded cutoff."""
    if explicit is not None:
        log(f"effect-size cutoff: |mean_diff| >= {explicit:g}  (given via --lfc)")
        return explicit
    path = RES_DIR / "10_effect_size_cutoff.csv"
    if not path.exists():
        raise SystemExit(
            f"{path} not found. Run the differential expression step first, or "
            "pass --lfc explicitly. Do not guess a cutoff: the effect-size "
            "scale is dataset-specific and a conventional value such as 1.0 "
            "silently returns nothing.")
    row = pd.read_csv(path).iloc[0]
    cutoff = float(row["cutoff"])
    log(f"effect-size cutoff: |mean_diff| >= {cutoff:.4f}  "
        f"(from 10_effect_size_cutoff.csv, rule: {row.get('rule', 'unknown')})")
    return cutoff


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--lfc", type=float, default=None,
                    help="absolute |mean_diff| cutoff. Omit to read the value "
                         "recorded by the DE step (recommended).")
    ap.add_argument("--min-specificity", type=float, default=1.5,
                    help="own-group mean_diff / next-best group mean_diff")
    ap.add_argument("--specificity-floor", type=float, default=0.1,
                    help="floor on the specificity denominator. Default 0.1 is "
                         "calibrated for a log2FC scale and may dominate on a "
                         "[0,1] scale — see the note in this file.")
    ap.add_argument("--top-per-group", type=int, default=15)
    args = ap.parse_args()
    ensure_dirs()

    de = pd.read_csv(RES_DIR / "10_de_results.csv")
    if "mean_diff" not in de.columns:
        raise SystemExit(
            "10_de_results.csv has no `mean_diff` column. If it has `log2FC`, "
            "it was written by an older version of the DE step — rerun that "
            "step rather than renaming the column here.")

    expr = load_matrix("expr_genes")
    gene_mean = expr.mean(axis=1)
    expr_thr = float(np.median(expr.values))
    cutoff = resolve_cutoff(args.lfc)

    sig = de[(de.fdr < args.fdr) & (de.mean_diff >= cutoff)].copy()
    log(f"{len(sig)} gene-group pairs pass FDR<{args.fdr} and "
        f"mean_diff>={cutoff:.4f}")
    if sig.empty:
        log("nothing passes — loosen the thresholds")
        return

    # specificity: own effect size vs the best that gene achieves elsewhere
    best_two = (de.sort_values("mean_diff", ascending=False)
                  .groupby("gene").head(2)[["gene", "group", "mean_diff"]])
    top1 = best_two.groupby("gene").nth(0).set_index("gene")
    top2 = best_two.groupby("gene").nth(1).set_index("gene")

    sig["next_best_diff"] = top2.reindex(sig.gene)["mean_diff"].to_numpy()
    sig["next_best_group"] = top2.reindex(sig.gene)["group"].to_numpy()
    sig["is_own_best"] = (top1.reindex(sig.gene)["group"].to_numpy()
                          == sig.group.to_numpy())

    raw_denom = np.abs(sig.next_best_diff.to_numpy())
    denom = np.maximum(raw_denom, args.specificity_floor)
    n_floored = int((raw_denom < args.specificity_floor).sum())
    sig["specificity"] = sig.mean_diff.to_numpy() / denom
    sig["specificity_floored"] = raw_denom < args.specificity_floor

    if n_floored:
        pct = 100.0 * n_floored / len(sig)
        msg = (f"specificity denominator floored at {args.specificity_floor:g} "
               f"for {n_floored}/{len(sig)} pairs ({pct:.0f}%)")
        if pct >= 25.0:
            log("WARNING: " + msg + ". For those genes specificity is just "
                "mean_diff / floor and does not measure specificity. Consider "
                "--specificity-floor 0.02, or drop specificity from the score.")
        else:
            log(msg)

    sig["mean_expression"] = gene_mean.reindex(sig.gene).to_numpy()
    sig["n_groups_up"] = (de[(de.fdr < args.fdr) & (de.mean_diff >= cutoff)]
                          .groupby("gene").size().reindex(sig.gene).to_numpy())

    cand = sig[(sig.is_own_best) &
               (sig.specificity >= args.min_specificity) &
               (sig.mean_expression > expr_thr)].copy()
    log(f"{len(cand)} pass specificity >= {args.min_specificity} and "
        f"expression > {expr_thr:.2f}")
    if cand.empty:
        log("no candidates survive — try --min-specificity 1.2")
        return

    def norm(s):
        s = np.asarray(s, dtype=float)
        rng = s.max() - s.min()
        return (s - s.min()) / rng if rng else np.ones_like(s)

    cand["score_significance"] = norm(-np.log10(cand.p_value.clip(lower=1e-300)))
    cand["score_effect"] = norm(cand.mean_diff)
    cand["score_specificity"] = norm(np.clip(cand.specificity, 0, 10))
    cand["candidate_score"] = (0.4 * cand.score_significance +
                               0.3 * cand.score_effect +
                               0.3 * cand.score_specificity)

    out = (cand.sort_values(["group", "candidate_score"], ascending=[True, False])
               .groupby("group").head(args.top_per_group)
               .sort_values("candidate_score", ascending=False)
               .reset_index(drop=True))
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    keep = ["rank", "gene", "group", "mean_diff", "fdr", "specificity",
            "specificity_floored", "next_best_group", "next_best_diff",
            "n_groups_up", "mean_expression", "candidate_score"]
    out[keep].to_csv(RES_DIR / "12_candidate_biomarkers.csv", index=False)

    summ = (out.groupby("group")
               .agg(n_candidates=("gene", "size"),
                    best_gene=("gene", "first"),
                    best_score=("candidate_score", "max"),
                    median_specificity=("specificity", "median"),
                    n_specificity_floored=("specificity_floored", "sum"))
               .reset_index())
    summ.to_csv(RES_DIR / "12_candidate_summary.csv", index=False)
    log("\n" + summ.to_string(index=False))
    log(f"top candidates -> {RES_DIR/'12_candidate_biomarkers.csv'}")
    log("these are hypotheses, not findings — steps 21-24 test whether a model "
        "actually relies on them")


if __name__ == "__main__":
    main()