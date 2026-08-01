#!/usr/bin/env python3
"""
Generate a small synthetic dataset in GEO series matrix format, plus a matching
probe annotation file and gene set library.

Run this first, against a scratch project root, to verify the whole pipeline
works before committing hours to the real 2 GB download:

    export SAILS_BASE=/tmp/sails_test
    python3 make_synthetic_data.py
    bash run_all.sh

The synthetic data has known structure: N classes, each with its own block of
up-regulated marker genes, plus a simulated laboratory effect. So you know what
the answers should look like — clusters should recover the classes, and the
marker gene sets should come out enriched.
"""
import argparse
import gzip

import numpy as np

from common import RAW_DIR, GENESET_DIR, ANNOT_FILE, log

CLASSES = ["ALL with t(12;21)", "AML with t(15;17)", "CLL",
           "CML", "non-leukemia bone marrow"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=180)
    ap.add_argument("--probes", type=int, default=3000)
    ap.add_argument("--markers", type=int, default=40,
                    help="marker genes per class")
    ap.add_argument("--labs", type=int, default=3)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    GENESET_DIR.mkdir(parents=True, exist_ok=True)

    n_s, n_p = args.samples, args.probes
    labels = np.array([CLASSES[i % len(CLASSES)] for i in range(n_s)])
    labs = np.array([f"Lab_{chr(65 + (i % args.labs))}" for i in range(n_s)])
    gsm = [f"GSM{9000000 + i}" for i in range(n_s)]

    # log2-scale baseline, then class marker blocks, then a lab offset
    X = rng.normal(7.0, 1.6, size=(n_p, n_s)).astype(np.float64)
    marker_map = {}
    cursor = 0
    for c in CLASSES:
        idx = np.arange(cursor, min(cursor + args.markers, n_p))
        cursor += args.markers
        if idx.size == 0:
            break
        m = labels == c
        X[np.ix_(idx, np.where(m)[0])] += rng.normal(3.0, 0.4, size=(idx.size, m.sum()))
        marker_map[c] = idx

    for i, lab in enumerate(sorted(set(labs))):
        X[:, labs == lab] += (i - 1) * 0.35

    X = np.clip(X, 0, None)
    linear = np.power(2.0, X) - 1.0          # back to a MAS5-like linear scale

    probes = [f"{1000 + i}_at" for i in range(n_p)]
    genes = [f"GENE{i:05d}" for i in range(n_p)]

    out = RAW_DIR / "GSE13159_series_matrix.txt.gz"
    with gzip.open(out, "wt", encoding="utf-8") as fh:
        fh.write("!Series_title\t\"Synthetic test series\"\n")
        fh.write("!Series_geo_accession\t\"GSE00000\"\n")
        fh.write("!Sample_title\t" + "\t".join(f'"sample_{i}"' for i in range(n_s)) + "\n")
        fh.write("!Sample_geo_accession\t" + "\t".join(f'"{g}"' for g in gsm) + "\n")
        fh.write("!Sample_platform_id\t" + "\t".join('"GPL570"' for _ in range(n_s)) + "\n")
        fh.write("!Sample_characteristics_ch1\t" +
                 "\t".join(f'"leukemia class: {c}"' for c in labels) + "\n")
        fh.write("!Sample_contact_institute\t" +
                 "\t".join(f'"{l}"' for l in labs) + "\n")
        fh.write("!series_matrix_table_begin\n")
        fh.write('"ID_REF"\t' + "\t".join(f'"{g}"' for g in gsm) + "\n")
        for i, p in enumerate(probes):
            fh.write(f'"{p}"\t' + "\t".join(f"{v:.3f}" for v in linear[i]) + "\n")
        fh.write("!series_matrix_table_end\n")
    log(f"series matrix -> {out}  ({n_p} probes x {n_s} samples)")

    with open(ANNOT_FILE, "w", encoding="utf-8") as fh:
        fh.write("probe,symbol\n")
        for p, g in zip(probes, genes):
            fh.write(f"{p},{g}\n")
    log(f"annotation -> {ANNOT_FILE}")

    gmt = GENESET_DIR / "SYNTHETIC_Pathways.txt"
    with open(gmt, "w", encoding="utf-8") as fh:
        for c, idx in marker_map.items():
            name = c.replace(" ", "_").replace("(", "").replace(")", "")
            fh.write(f"true markers {name} (GO:0000001)\t\t" +
                     "\t".join(genes[i] for i in idx) + "\n")
        for j in range(20):     # decoy sets that should NOT come out enriched
            pick = rng.choice(n_p, 30, replace=False)
            fh.write(f"decoy set {j} (GO:9{j:06d})\t\t" +
                     "\t".join(genes[i] for i in pick) + "\n")
    log(f"gene sets -> {gmt}")
    log(f"expected: {len(CLASSES)} classes, ~{args.markers} markers each")


if __name__ == "__main__":
    main()
