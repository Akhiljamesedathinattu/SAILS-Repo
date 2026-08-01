#!/usr/bin/env python3
"""
Download GO and KEGG gene set libraries into raw/genesets/.

Enrichr serves each library as tab-delimited text:
    <term name>\t<optional description>\t<gene>\t<gene>\t...

Usage
    python3 fetch_genesets.py --list              # show available library names
    python3 fetch_genesets.py                     # download the default set
    python3 fetch_genesets.py --libraries KEGG_2021_Human Reactome_2022

If your machine has no outbound access, download the same files by hand from
https://maayanlab.cloud/Enrichr → Libraries, and drop them into raw/genesets/
with a .txt or .gmt extension. Nothing else in the pipeline cares how they got
there.
"""
import argparse
import json
import urllib.error
import urllib.request

from common import GENESET_DIR, log

BASE = "https://maayanlab.cloud/Enrichr"
LIB_URL = BASE + "/geneSetLibrary?mode=text&libraryName={}"
STATS_URL = BASE + "/datasetStatistics"

DEFAULT_LIBRARIES = [
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023",
    "KEGG_2021_Human",
]

EXTRA_SUGGESTIONS = [
    "MSigDB_Hallmark_2020",
    "Reactome_2022",
    "WikiPathways_2024_Human",
]


def fetch(url, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": "sails-project/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", errors="replace")


def list_libraries():
    try:
        stats = json.loads(fetch(STATS_URL))
    except Exception as e:
        log(f"could not reach Enrichr ({e}). Browse the Libraries tab manually.")
        return
    names = sorted(s["libraryName"] for s in stats.get("statistics", []))
    for n in names:
        if n.startswith(("GO_", "KEGG", "Reactome", "WikiPath", "MSigDB")):
            print(n)
    log(f"{len(names)} libraries available in total")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--libraries", nargs="+", default=DEFAULT_LIBRARIES)
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    if args.list:
        list_libraries()
        return

    GENESET_DIR.mkdir(parents=True, exist_ok=True)
    ok = 0
    for lib in args.libraries:
        out = GENESET_DIR / f"{lib}.txt"
        if out.exists() and out.stat().st_size > 1000:
            log(f"{lib}: already present, skipping")
            ok += 1
            continue
        try:
            text = fetch(LIB_URL.format(lib))
        except urllib.error.HTTPError as e:
            log(f"{lib}: HTTP {e.code} — name may have changed, try --list")
            continue
        except Exception as e:
            log(f"{lib}: {type(e).__name__}: {e}")
            continue
        if len(text) < 1000 or "\t" not in text:
            log(f"{lib}: response looks wrong ({len(text)} bytes), not saved")
            continue
        out.write_text(text, encoding="utf-8")
        n_sets = sum(1 for line in text.splitlines() if line.strip())
        log(f"{lib}: {n_sets} gene sets -> {out}")
        ok += 1

    log(f"{ok}/{len(args.libraries)} libraries in {GENESET_DIR}")
    if ok:
        log("next: python3 06_enrichment.py")
    else:
        log("nothing downloaded. Other libraries worth trying: "
            + ", ".join(EXTRA_SUGGESTIONS))


if __name__ == "__main__":
    main()
