#!/usr/bin/env python3

# =============================================================
# DOWNLOAD THE GENE SET FILES (GO AND KEGG)
#
# Run this ONCE, BEFORE simple_pipeline.py.
#
# WHY WE NEED IT
# Step 13 of the pipeline asks "what do my disease genes actually
# DO?". To answer that it needs the official lists that biologists
# have built up over decades, such as:
#
#   "T cell activation"  ->  CD3D, CD3E, LCK, ZAP70, ...
#   "DNA repair"         ->  BRCA1, BRCA2, ATM, TP53, ...
#
# This script downloads four of those libraries from Enrichr, a
# free public website that collects them.
#
# NO INTERNET? Go to https://maayanlab.cloud/Enrichr, click the
# Libraries tab, download the four names listed below by hand, and
# put them in the genesets folder. The pipeline does not care how
# the files got there.
# =============================================================

import os
import urllib.request


# Use the same project folder as simple_pipeline.py
BASE_FOLDER = "/home/sails/SAILS-Repo/Gene_Expression_Clustering"
GENESET_FOLDER = BASE_FOLDER + "/raw/genesets"

WEBSITE = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="

# The four libraries the original project used.
#   GO Biological Process = what a gene DOES
#   GO Molecular Function = what a gene IS (an enzyme, a receptor...)
#   GO Cellular Component = WHERE in the cell the gene works
#   KEGG                  = chains of genes working in sequence
WANTED_LIBRARIES = [
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023",
    "KEGG_2021_Human",
]


def download_one_library(library_name):
    # Download one library and save it as a text file.
    # Returns True if we ended up with a usable file.

    save_to = GENESET_FOLDER + "/" + library_name + ".txt"

    # Already have it? Do not download it again.
    if os.path.exists(save_to):
        size = os.path.getsize(save_to)
        if size > 1000:
            print(library_name, "- already here, skipping")
            return True

    web_address = WEBSITE + library_name

    print(library_name, "- downloading...")

    try:
        # Some websites refuse requests with no name attached, so we
        # politely say who we are.
        request = urllib.request.Request(
            web_address,
            headers={"User-Agent": "simple-pipeline/1.0"})

        connection = urllib.request.urlopen(request, timeout=120)
        raw_bytes = connection.read()
        connection.close()

        text = raw_bytes.decode("utf-8", errors="replace")

    except Exception as problem:
        # Print what actually went wrong, not a vague guess
        print("  FAILED:", type(problem).__name__, "-", problem)
        return False

    # ---- Check we got something sensible before saving ----
    # A wrong library name gives a tiny error page, not gene sets.
    # Saving that would break Step 13 in a confusing way later.
    if len(text) < 1000:
        print("  the reply was only", len(text), "characters long.")
        print("  That is too small to be real. Not saving it.")
        print("  The library name may have changed.")
        return False

    if "\t" not in text:
        print("  the reply has no tab characters in it, so it is not")
        print("  a gene set file. Not saving it.")
        return False

    # Save it
    output_file = open(save_to, "w", encoding="utf-8")
    output_file.write(text)
    output_file.close()

    # Count the lines, so the user can see it worked
    how_many_sets = 0
    for one_line in text.split("\n"):
        if one_line.strip() != "":
            how_many_sets = how_many_sets + 1

    print("  saved", how_many_sets, "gene sets")
    return True


def get_all_libraries():
    # Make the folder, then download each library in turn.

    os.makedirs(GENESET_FOLDER, exist_ok=True)
    print("Saving into:", GENESET_FOLDER)
    print("")

    how_many_worked = 0

    for one_library in WANTED_LIBRARIES:
        worked = download_one_library(one_library)
        if worked == True:
            how_many_worked = how_many_worked + 1

    print("")
    print("Got", how_many_worked, "of", len(WANTED_LIBRARIES), "libraries")

    if how_many_worked == 0:
        print("")
        print("Nothing downloaded. If you have no internet access, get")
        print("the files by hand from:")
        print("  https://maayanlab.cloud/Enrichr  (Libraries tab)")
        print("and put them in:", GENESET_FOLDER)
    else:
        print("Now you can run: python3 simple_pipeline.py")


def mymain():
    get_all_libraries()


if __name__ == "__main__":
    mymain()
