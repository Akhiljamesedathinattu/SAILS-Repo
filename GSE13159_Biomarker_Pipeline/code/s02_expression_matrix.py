#!/usr/bin/env python3
"""
PIPELINE STEP 3 : Expression Matrix Loading

Reads the expression table as float32 and caches it. Metadata is parsed too but
is deliberately NOT used yet — clinical integration is step 15 in this pipeline.
It is only stored so step 09 does not have to re-read a 2 GB file.

Outputs
  work/expr_raw.{parquet,npz}
  work/sample_metadata_raw.csv
"""
import argparse

import numpy as np

from common import SERIES_MATRIX, WORK_DIR, ensure_dirs, save_matrix, log
from geo import scan_header, build_metadata, read_expression


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(SERIES_MATRIX))
    args = ap.parse_args()
    ensure_dirs()

    log(f"loading {args.input}")
    meta, n_skip = scan_header(args.input)
    expr = read_expression(args.input, n_skip)
    md = build_metadata(meta)
    log(f"expression matrix: {expr.shape[0]} probes x {expr.shape[1]} samples")

    shared = [s for s in expr.columns if s in set(md.index)]
    expr = expr[shared]
    md.loc[shared].to_csv(WORK_DIR / "sample_metadata_raw.csv")
    log(f"metadata parked for step 09 -> {WORK_DIR/'sample_metadata_raw.csv'}")
    log(f"cached -> {save_matrix(expr.astype(np.float32), 'expr_raw')}")


if __name__ == "__main__":
    main()
