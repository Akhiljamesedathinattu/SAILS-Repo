"""
common.py — shared paths, caching and small statistics helpers.

Project root is taken from $SAILS_BASE, defaulting to the SAILS repo path.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(os.environ.get(
    "SAILS_BASE", "/home/sails/SAILS-Repo/Gene_Expression_Clustering"))

RAW_DIR = BASE / "raw"
WORK_DIR = BASE / "work"
RES_DIR = BASE / "results"
FIG_DIR = BASE / "figures"
MODEL_DIR = BASE / "models"

SERIES_MATRIX = RAW_DIR / "GSE13159_series_matrix.txt.gz"
ANNOT_FILE = RAW_DIR / "GPL570_annot.csv"
GENESET_DIR = RAW_DIR / "genesets"

SEED = 42


def ensure_dirs():
    for d in (WORK_DIR, RES_DIR, FIG_DIR, MODEL_DIR):
        d.mkdir(parents=True, exist_ok=True)


def log(msg):
    print(f"[{Path(sys.argv[0]).stem}] {msg}", flush=True)


# ---------------------------------------------------------------- matrix I/O
def save_matrix(df, name):
    """Cache a features x samples matrix. Parquet when pyarrow is present."""
    ensure_dirs()
    try:
        import pyarrow  # noqa: F401
        p = WORK_DIR / f"{name}.parquet"
        df.to_parquet(p)
    except Exception:
        p = WORK_DIR / f"{name}.npz"
        np.savez_compressed(
            p,
            values=np.asarray(df.values, dtype=np.float32),
            index=df.index.to_numpy().astype(str),
            columns=df.columns.to_numpy().astype(str),
        )
    return p


def load_matrix(name):
    pq = WORK_DIR / f"{name}.parquet"
    if pq.exists():
        return pd.read_parquet(pq)
    npz = WORK_DIR / f"{name}.npz"
    if npz.exists():
        z = np.load(npz, allow_pickle=False)
        return pd.DataFrame(z["values"], index=z["index"], columns=z["columns"])
    raise FileNotFoundError(
        f"cached matrix '{name}' not found in {WORK_DIR} — run the earlier step")


def load_metadata():
    f = RES_DIR / "sample_metadata.csv"
    if not f.exists():
        raise FileNotFoundError("run 01_load.py first")
    return pd.read_csv(f, index_col=0)


# ------------------------------------------------------------------ statistics
def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, dtype=np.float64)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    adj = p[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out


def load_probe_annotation():
    """probe -> gene symbol map, or None when the annotation file is absent."""
    if not ANNOT_FILE.exists():
        return None
    ann = pd.read_csv(ANNOT_FILE)
    low = {c.lower(): c for c in ann.columns}
    pcol = low.get("probe") or low.get("probe_id") or ann.columns[0]
    scol = low.get("symbol") or low.get("gene_symbol") or ann.columns[1]
    s = ann[[pcol, scol]].dropna()
    s = s[s[scol].astype(str).str.len() > 0]
    return s.set_index(pcol)[scol].astype(str)
