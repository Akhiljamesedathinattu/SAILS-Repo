#!/usr/bin/env python3
"""
PIPELINE STEPS 1-2 : Raw GSE13159 Dataset -> Data Understanding

Inspects the download without loading the full expression table into memory.
Reports platform, sample count, probe count, value scale, and every metadata
field present, so you know what you are working with before committing to it.

SCALE DETECTION

  An earlier version classified scale with a single test — vmax > 50 meant
  linear, anything else meant log2. That has no branch for an already-scaled
  matrix, so a [0,1] series matrix was silently reported as log2. GSE13159 as
  distributed by GEO is exactly that case, and the misreport propagated into
  the written methods.

  Scale is now classified into four cases with the evidence recorded alongside
  the label, so a downstream reader can check the call rather than trust it:

    scaled_0_1   max <= ~1, min >= ~0     differences are NOT fold changes
    log2         max <= ~30               differences ARE log2 fold changes
    linear       max > ~50                needs log2 transform before testing
    ambiguous    anything else            inspect manually

  The peek also samples probes from across the file rather than taking the
  first N rows. On GPL570 the leading rows are AFFX- control probes, whose
  range is not representative of the biological probes.

QUANTILE NORMALISATION

  If every sample shares a near-identical maximum, the arrays have been forced
  onto a common distribution. That is worth recording: it means cross-sample
  comparison is sound, and it distinguishes quantile normalisation from
  per-sample min-max scaling, which would not be.

Outputs
  results/01_dataset_overview.csv
  results/01_metadata_fields.csv
"""
import argparse
import gzip

import numpy as np
import pandas as pd

from common import SERIES_MATRIX, RES_DIR, ensure_dirs, log


def classify_scale(vmin, vmax):
    """Return (label, note). Ordered most-specific first."""
    if vmax <= 1.01 and vmin >= -0.01:
        return ("scaled_0_1",
                "values already rescaled to [0,1]; differences between group "
                "means are NOT log2 fold changes and cannot be converted to "
                "fold changes without the pre-scaling range")
    if vmax <= 30.0:
        return ("log2",
                "consistent with log2 expression; differences between group "
                "means are log2 fold changes")
    if vmax > 50.0:
        return ("linear",
                "consistent with linear intensities (MAS5-like); apply a log2 "
                "transform before differential testing")
    return ("ambiguous",
            f"max {vmax:.3g} falls between the log2 and linear ranges; "
            "inspect the distribution manually before proceeding")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(SERIES_MATRIX))
    ap.add_argument("--peek-rows", type=int, default=2000,
                    help="probes sampled across the table to characterise the "
                         "value scale (default 2000)")
    args = ap.parse_args()
    ensure_dirs()

    series, fields, n_header = {}, {}, 0
    with gzip.open(args.input, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if line.startswith("!series_matrix_table_begin"):
                n_header = i
                break
            parts = line.rstrip("\n").split("\t")
            key = parts[0].lstrip("!")
            vals = [p.strip('"') for p in parts[1:]]
            if key.startswith("Series_"):
                series.setdefault(key, []).append(" | ".join(vals)[:300])
            elif key.startswith("Sample_"):
                fields.setdefault(key, []).append(vals)

    n_samples = len(fields.get("Sample_geo_accession", [[]])[0])
    log(f"{n_samples} samples, {n_header} metadata lines")

    n_probes = 0
    with gzip.open(args.input, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i <= n_header + 1:
                continue
            if line.startswith("!"):
                break
            n_probes += 1

    # Sample probes spread across the table rather than the first N. The
    # leading rows of a GPL570 matrix are AFFX- control probes and are not
    # representative of the biological range.
    if n_probes > args.peek_rows:
        step = max(1, n_probes // args.peek_rows)
        skip = [r for r in range(n_header + 2, n_header + 2 + n_probes)
                if (r - n_header - 2) % step != 0]
    else:
        skip = []
    peek = pd.read_csv(args.input, sep="\t", header=0, index_col=0,
                       low_memory=False,
                       skiprows=list(range(n_header + 1)) + skip,
                       comment="!")
    peek = peek.apply(pd.to_numeric, errors="coerce")

    affx = peek.index.astype(str).str.startswith("AFFX")
    if affx.any():
        log(f"excluding {int(affx.sum())} AFFX control probes from scale check")
        peek = peek[~affx]

    vmin, vmax = float(np.nanmin(peek.values)), float(np.nanmax(peek.values))
    scale, scale_note = classify_scale(vmin, vmax)

    # Distinguish quantile normalisation from per-sample min-max scaling.
    col_max = np.nanmax(peek.values, axis=0)
    col_min = np.nanmin(peek.values, axis=0)
    max_spread = float(np.nanmax(col_max) - np.nanmin(col_max))
    quantile_like = bool(max_spread < 0.01 * max(abs(vmax), 1e-9))
    per_sample_minmax = bool(
        quantile_like
        and np.allclose(col_min, np.nanmin(col_min), atol=1e-6)
        and scale == "scaled_0_1")

    if quantile_like:
        log("every sample shares a near-identical maximum — arrays appear "
            "quantile normalised onto a common distribution")
    if scale == "scaled_0_1":
        log("value scale: rescaled to [0,1]. Effect sizes from this matrix are "
            "NOT fold changes. See the SCALE DETECTION note in this file.")

    platform = fields.get("Sample_platform_id", [[""]])[0]
    rows = [
        ("file", args.input),
        ("series_title", series.get("Series_title", [""])[0]),
        ("accession", series.get("Series_geo_accession", [""])[0]),
        ("platform", sorted(set(platform))[0] if platform else "unknown"),
        ("n_samples", n_samples),
        ("n_probes", n_probes),
        ("value_scale", scale),
        ("value_scale_note", scale_note),
        ("effect_sizes_are_fold_changes", scale == "log2"),
        ("quantile_normalised_likely", quantile_like),
        ("per_sample_minmax_likely", per_sample_minmax),
        ("value_min_sampled", round(vmin, 6)),
        ("value_max_sampled", round(vmax, 6)),
        ("n_probes_sampled_for_scale", int(peek.shape[0])),
        ("metadata_lines", n_header),
        ("probe_id_example", str(peek.index[0])),
    ]
    pd.DataFrame(rows, columns=["property", "value"]).to_csv(
        RES_DIR / "01_dataset_overview.csv", index=False)

    frows = []
    for key, occurrences in fields.items():
        for j, vals in enumerate(occurrences):
            uniq = pd.Series(vals).nunique()
            frows.append({
                "field": key.replace("Sample_", ""),
                "occurrence": j + 1,
                "n_unique_values": uniq,
                "constant": uniq == 1,
                "example": str(vals[0])[:120] if vals else "",
            })
    fdf = pd.DataFrame(frows).sort_values("n_unique_values", ascending=False)
    fdf.to_csv(RES_DIR / "01_metadata_fields.csv", index=False)

    log(f"platform={rows[3][1]}  probes={n_probes}  scale={scale}")
    log("metadata fields that vary across samples (candidate group labels):")
    for _, r in fdf[(~fdf.constant) & (fdf.n_unique_values < 100)].head(8).iterrows():
        log(f"    {r.field:<28} {r.n_unique_values:>4} values   e.g. {r.example[:60]}")
    log("STEP 2 complete — inspect results/01_metadata_fields.csv to pick the "
        "field that defines your patient groups in step 09")


if __name__ == "__main__":
    main()