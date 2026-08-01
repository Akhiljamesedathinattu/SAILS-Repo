"""
geo.py — GEO series matrix parsing, shared by 01_load.py and 09_validate_external.py.

A series matrix file is two documents in one: '!'-prefixed metadata lines, then
an expression table delimited by !series_matrix_table_begin / _end.
"""
import gzip
import re

import numpy as np
import pandas as pd


def scan_header(path):
    """Return (metadata dict, number of lines to skip before the table header)."""
    meta, n_skip = {}, None
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if line.startswith("!series_matrix_table_begin"):
                n_skip = i + 1
                break
            if line.startswith("!Sample_"):
                parts = line.rstrip("\n").split("\t")
                meta.setdefault(parts[0].lstrip("!"), []).append(
                    [p.strip('"') for p in parts[1:]])
    if n_skip is None:
        raise ValueError(f"{path}: !series_matrix_table_begin not found — "
                         "truncated or non-series-matrix file?")
    return meta, n_skip


def clean(name):
    return re.sub(r"[^0-9a-zA-Z]+", "_", str(name).strip().lower()).strip("_") or "field"


def build_metadata(meta):
    """Pivot !Sample_* lines into one row per GSM, splitting 'key: value' fields."""
    gsm = meta["Sample_geo_accession"][0]
    df = pd.DataFrame(index=pd.Index(gsm, name="gsm"))

    for key in ("Sample_title", "Sample_source_name_ch1", "Sample_platform_id",
                "Sample_scan_date", "Sample_contact_institute"):
        if key in meta and len(meta[key][0]) == len(gsm):
            df[clean(key.replace("Sample_", ""))] = meta[key][0]

    for i, vals in enumerate(meta.get("Sample_characteristics_ch1", [])):
        if len(vals) != len(gsm):
            continue
        keys = {v.split(":", 1)[0] for v in vals if ":" in v}
        col = clean(keys.pop()) if len(keys) == 1 else f"characteristic_{i + 1}"
        df[col] = [v.split(":", 1)[1].strip() if ":" in v else v.strip()
                   for v in vals]

    if "leukemia_class" not in df.columns:
        for cand in ("source_name_ch1", "characteristic_1", "title"):
            if cand in df.columns:
                df["leukemia_class"] = df[cand].astype(str)
                break
        else:
            df["leukemia_class"] = "unknown"
    df["leukemia_class"] = df["leukemia_class"].astype(str).str.strip()
    return df


def read_expression(path, n_skip):
    """Read the expression table as float32, dropping the table_end sentinel row."""
    head = pd.read_csv(path, sep="\t", skiprows=n_skip, nrows=0)
    dtypes = {c: np.float32 for c in head.columns[1:]}
    expr = pd.read_csv(path, sep="\t", skiprows=n_skip, index_col=0,
                       dtype=dtypes, na_values=["", "NA", "null", "NULL"])
    expr = expr[~expr.index.astype(str).str.startswith("!")]
    expr.index = expr.index.astype(str)
    expr.index.name = "probe_id"
    expr.columns = [str(c).strip('"') for c in expr.columns]
    return expr


def load_series_matrix(path):
    """Convenience wrapper: returns (expression DataFrame, metadata DataFrame)."""
    meta, n_skip = scan_header(path)
    return read_expression(path, n_skip), build_metadata(meta)
