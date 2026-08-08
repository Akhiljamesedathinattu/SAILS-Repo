#!/usr/bin/env python3

# =============================================================
# SIMPLE LEUKAEMIA GENE EXPRESSION PIPELINE
#
# This program does the same science as the original 16-script
# project, but it is written in a very simple beginner style.
#
# Read it from top to bottom. Each function does ONE job.
# mymain() at the bottom calls them in order.
# =============================================================

# ---- Libraries we use ----
# gzip  : lets us read files that end in .gz (compressed files)
# os    : lets us make folders
# numpy : works with big grids of numbers (fast maths)
# pandas: works with tables of data (like Excel sheets)
# scipy : has ready-made statistics tools (t-test, clustering)
# sklearn: has ready-made machine learning tools

import gzip
import os

import numpy as np
import pandas as pd

from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import pdist, squareform

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.metrics import normalized_mutual_info_score
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.metrics import confusion_matrix, roc_curve
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.base import BaseEstimator
from sklearn.feature_selection._base import SelectorMixin
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


# ---- Settings (change these if you want) ----

# THIS IS THE ONLY LINE MOST PEOPLE NEED TO CHANGE.
#
# It is the full path to your project folder. We write out the WHOLE
# path, starting from "/", instead of a short name like "raw".
#
# Why? A short name like "raw" means "a folder called raw, next to
# wherever I happen to be standing right now". If you run the program
# from a different folder, Python looks in the wrong place and says
# "No such file or directory". A full path always points to the same
# place, no matter where you run the program from.

BASE_FOLDER = "/home/sails/SAILS-Repo/Gene_Expression_Clustering"

RAW_FOLDER = BASE_FOLDER + "/raw"           # where the downloaded data lives
RESULT_FOLDER = BASE_FOLDER + "/results"    # where our answers are saved
MODEL_FOLDER = BASE_FOLDER + "/models"      # where the trained models are saved

DATA_FILE = RAW_FOLDER + "/GSE13159_series_matrix.txt.gz"
ANNOTATION_FILE = RAW_FOLDER + "/GPL570_annot.csv"

SEED = 42                       # a fixed number so results repeat exactly
TOP_PROBES = 5000               # how many most-variable probes to keep
MIN_GROUP_SIZE = 20             # a disease group needs at least 20 patients
FDR_CUTOFF = 0.05               # how strict we are about p-values
SPECIFICITY_CUTOFF = 1.5        # how "one-group-only" a biomarker must be
SPECIFICITY_FLOOR = 0.1         # smallest allowed bottom of the specificity sum

# CANDIDATES_PER_DISEASE stops one big disease taking over the list.
# CLL has 448 patients, so its genes get the strongest statistics and
# would fill every top-10 place. Keeping only the best few from EACH
# disease means the final panel can actually tell the diseases apart.
CANDIDATES_PER_DISEASE = 15     # best N biomarkers kept from EACH disease

# ---- How do we decide the number of patient groups? ----
#
# "stability"  = which answer keeps coming back when we redo the
#                clustering on random parts of the data?  (gives k = 5)
# "tightness"  = which answer makes the neatest, roundest groups?
#                (gives k = 3)
#
# The original study used "stability". See the long explanation
# inside choose_number_of_groups() for why.

HOW_TO_CHOOSE_K = "stability"

FORCE_K = 0                     # 0 = decide automatically. Or put 5 here.
CONSENSUS_REPEATS = 20          # how many times we redo the clustering
CONSENSUS_MAX_K = 12            # biggest k we test for stability

# ---- Machine learning settings ----
#
# We train TWO models and compare them, like the original study:
#   random_forest        = hundreds of decision trees voting together
#   logistic_regression  = one equation weighing up all the genes
#
# The one whose gene ranking goes to Step 11 is decided below.

MODEL_NAMES = ["random_forest", "logistic_regression"]

# WHICH MODEL'S GENE RANKING FEEDS STEP 11?
#
# Do NOT just leave this on a name because someone else used it.
# On this dataset logistic regression scores clearly better on
# balanced accuracy (about 0.90 versus 0.77) while the two models tie
# on AUC (about 0.987 each). So "which is best" depends on the
# measure, and the two models lean on very different genes.
#
# Options:
#   "best"                = let the cross-validation scores decide
#   "random_forest"       = force it (what the original project used)
#   "logistic_regression" = force it
#
# "best" is the default because it is the only one of the three that
# is an argument rather than an assumption. Whatever it picks is
# printed and saved, so you can quote the reason in your write-up.

MAIN_MODEL = "best"

# Which measure decides, when MAIN_MODEL is "best"?
#   "cv_balanced_accuracy" = average over the 5 training splits.
#                            More trustworthy: it is an average of
#                            five numbers, not one lucky split.
#   "test_balanced_accuracy" = the single hidden-patient score
#   "test_auc"               = ranking ability rather than accuracy
BEST_MODEL_MEASURE = "cv_balanced_accuracy"

USE_SHAP = True                 # False = use the model's quick built-in score
SHAP_PATIENTS = 300             # how many test patients to explain
GENES_FOR_MODEL = 2000          # genes kept when HOW_TO_PICK_GENES = "global"

# How should the model choose which genes to use?
#   "per_disease" = best genes for EACH disease, then combine (recommended)
#   "global"      = best genes for telling all 17 apart at once
HOW_TO_PICK_GENES = "per_disease"
GENES_PER_DISEASE = 150         # genes kept per disease in "per_disease" mode

CV_FOLDS = 5                    # how many times we split the training data
PANEL_SIZE = 25                 # how many genes go in the biomarker panel

# ---- Gene network settings ----
NETWORK_GENES = 300             # most variable genes to put in the network
NETWORK_STRENGTH = 0.7          # how strongly two genes must move together
NETWORK_TEAMS = 8               # how many gene teams to look for

# ---- GO / KEGG enrichment settings ----
GENESET_FOLDER = RAW_FOLDER + "/genesets"
ENRICHMENT_MIN_SET = 10         # ignore gene sets smaller than this
ENRICHMENT_MAX_SET = 500        # ignore gene sets bigger than this
ENRICHMENT_TOP = 10             # best N results kept per disease


def make_folders():
    # Create the folders we will save things into.
    # exist_ok=True means "do not complain if it is already there".
    os.makedirs(RAW_FOLDER, exist_ok=True)
    os.makedirs(RESULT_FOLDER, exist_ok=True)
    os.makedirs(MODEL_FOLDER, exist_ok=True)
    os.makedirs(GENESET_FOLDER, exist_ok=True)
    print("Folders ready")


def check_the_data_file():
    # Look for the data file BEFORE we start, so that if it is missing
    # we can print a helpful message instead of a scary Python error.
    #
    # os.path.exists() answers a simple question: is there a file
    # at this exact path? It gives back True or False.

    print("Looking for the data file...")
    print("  " + DATA_FILE)

    if os.path.exists(DATA_FILE):
        # Show the size in gigabytes, so a half-finished download is obvious
        size_in_bytes = os.path.getsize(DATA_FILE)
        size_in_gb = size_in_bytes / 1024 / 1024 / 1024
        print("Found it. Size:", round(size_in_gb, 2), "GB")
        return

    # ---- The file is not there. Help the user find out why. ----
    print("")
    print("ERROR: I cannot find the data file.")
    print("")

    if os.path.exists(RAW_FOLDER) == False:
        print("The raw folder does not exist either:")
        print("  " + RAW_FOLDER)
        print("Check that BASE_FOLDER at the top of this file is correct.")
        raise SystemExit(1)

    # The folder exists, so list what IS inside it. Very often the file
    # is there but with a slightly different name.
    print("The raw folder exists. Here is what is inside it:")

    file_names = os.listdir(RAW_FOLDER)

    if len(file_names) == 0:
        print("  (the folder is empty)")
    else:
        for one_name in file_names:
            print("  " + one_name)

    print("")
    print("If you can see your data file in that list, copy its exact name")
    print("into the DATA_FILE line at the top of this program.")
    raise SystemExit(1)


# =============================================================
# STEP 1 - LOAD THE DATA
# =============================================================

def find_table_start():
    # The data file has two parts stuck together:
    #   1) patient information lines, each starting with "!"
    #   2) the big table of gene expression numbers
    #
    # A special line tells us where the table starts. We count
    # lines until we find it.

    line_number = 0
    found = False

    open_file = gzip.open(DATA_FILE, "rt", errors="replace")

    for line in open_file:
        if line.startswith("!series_matrix_table_begin"):
            found = True
            break
        line_number = line_number + 1

    open_file.close()

    if found == False:
        print("ERROR: the file does not look like a GEO series matrix.")
        print("It may be broken or only half downloaded.")
        raise SystemExit(1)

    # We must skip the information lines AND the marker line itself.
    lines_to_skip = line_number + 1
    print("The number table starts after", lines_to_skip, "lines")
    return lines_to_skip


def load_patient_information():
    # Read the "!" lines and pull out the disease name for each patient.
    #
    # Biology note:
    #   Sample = one patient's blood or bone marrow measurement
    #   leukemia class = the doctor's diagnosis for that patient
    #
    # CAREFUL: the file has SEVERAL "!Sample_characteristics_ch1" lines,
    # not just one. Each holds a different fact about the patient:
    #
    #   !Sample_characteristics_ch1   "tissue: bone marrow"      ...
    #   !Sample_characteristics_ch1   "leukemia class: CLL"      ...
    #
    # We want the DIAGNOSIS line. If we just grab the first one we get
    # the tissue instead, and the whole analysis compares bone marrow
    # against blood rather than one leukaemia type against another.
    # So we collect them all, then choose the right one on purpose.

    patient_ids = []

    # A dictionary: the label of each line -> the list of its values
    all_facts = {}

    open_file = gzip.open(DATA_FILE, "rt", errors="replace")

    for line in open_file:
        if line.startswith("!series_matrix_table_begin"):
            break

        # Split the line into pieces separated by tabs
        pieces = line.strip().split("\t")
        first_piece = pieces[0]

        # The patient ID line
        if first_piece == "!Sample_geo_accession":
            for one_piece in pieces[1:]:
                clean_text = one_piece.replace('"', "")
                patient_ids.append(clean_text)

        # One of the fact lines. Each value looks like "label: value".
        if first_piece == "!Sample_characteristics_ch1":
            this_label = ""
            these_values = []

            for one_piece in pieces[1:]:
                clean_text = one_piece.replace('"', "").strip()

                if ":" in clean_text:
                    two_parts = clean_text.split(":", 1)
                    this_label = two_parts[0].strip()
                    these_values.append(two_parts[1].strip())
                else:
                    these_values.append(clean_text)

            # If the line had no label at all, invent one so we can
            # still tell the lines apart
            if this_label == "":
                this_label = "unnamed_line_" + str(len(all_facts) + 1)

            all_facts[this_label] = these_values

    open_file.close()

    # ---- Show the user what we found ----
    print("Loaded information for", len(patient_ids), "patients")
    print("The file describes each patient with these facts:")

    for one_label in all_facts:
        how_many_different = len(set(all_facts[one_label]))
        print("  '" + one_label + "' -", how_many_different, "different values")

    # ---- Choose the diagnosis line ----
    chosen_label = ""

    # First choice: a label that mentions leukaemia
    for one_label in all_facts:
        if "leukemia" in one_label.lower():
            chosen_label = one_label
            break

    # Second choice: the line with the MOST different values.
    # The diagnosis has ~18 values; tissue has only 2. So the line
    # with the most variety is almost always the diagnosis.
    if chosen_label == "":
        print("No label mentions leukaemia. Using the most varied line instead.")
        most_so_far = -1
        for one_label in all_facts:
            how_many_different = len(set(all_facts[one_label]))
            if how_many_different > most_so_far:
                most_so_far = how_many_different
                chosen_label = one_label

    if chosen_label == "":
        print("ERROR: no patient facts found in the file at all.")
        raise SystemExit(1)

    print("Using '" + chosen_label + "' as the diagnosis")

    disease_names = all_facts[chosen_label]

    # ---- Safety check ----
    # If the two lists are different lengths, something is wrong with
    # the file and lining them up would silently mislabel patients.
    if len(disease_names) != len(patient_ids):
        print("ERROR: found", len(patient_ids), "patients but",
              len(disease_names), "diagnoses. They must match.")
        raise SystemExit(1)

    # Put the two lists side by side in a small table
    patients = pd.DataFrame()
    patients["sample"] = patient_ids
    patients["disease"] = disease_names

    # ALSO keep every other fact as its own column. We need them later
    # to check for batch effects: if a technical detail like which
    # tissue was sampled explains more of the data than the diagnosis
    # does, that is a problem we must report.
    for one_label in all_facts:
        if len(all_facts[one_label]) != len(patient_ids):
            continue
        column_name = one_label.replace(" ", "_")
        if column_name == "disease":
            continue
        patients[column_name] = all_facts[one_label]

    patients = patients.set_index("sample")

    # Show a few real diagnosis names so the user can check them
    print("Example diagnoses found:")
    counts = patients["disease"].value_counts()
    for i in range(min(5, len(counts))):
        print("  ", counts.index[i], "-", counts.iloc[i], "patients")

    return patients


def load_gene_data(lines_to_skip):
    # Read the big table of numbers.
    #
    # Biology note:
    #   Probe = a tiny spot on the chip that measures one gene
    #   Gene expression = how active a gene is (a bigger number
    #                     means the gene is more active)
    #
    # The table has one row per probe and one column per patient.

    print("Reading the big number table... this takes a few minutes")

    data = pd.read_csv(DATA_FILE,
                       sep="\t",
                       skiprows=lines_to_skip,
                       index_col=0,
                       low_memory=False)

    # The very last row is a marker line, not data. Remove any row
    # whose name starts with "!".
    good_rows = []
    for probe_name in data.index:
        if str(probe_name).startswith("!"):
            good_rows.append(False)
        else:
            good_rows.append(True)
    data = data[good_rows]

    # Make sure every value is a number, and use float32 to save memory
    data = data.astype("float32")

    # Tidy up the column names (remove quote marks)
    clean_names = []
    for name in data.columns:
        clean_names.append(str(name).replace('"', ""))
    data.columns = clean_names

    print("Table size:", data.shape[0], "probes x", data.shape[1], "patients")
    return data


# =============================================================
# STEP 2 - CLEAN THE DATA
# =============================================================

def fill_missing_values(data):
    # Sometimes a measurement failed and the value is missing (NaN).
    # We replace a missing value with the average of that probe
    # across all patients. This is fine because very few are missing.

    number_missing = int(data.isna().sum().sum())
    total_values = data.shape[0] * data.shape[1]
    percent_missing = 100.0 * number_missing / total_values

    print("Missing values:", number_missing,
          "(", round(percent_missing, 4), "% )")

    if number_missing > 0:
        # Average of each row (each probe)
        row_averages = data.mean(axis=1)

        # Fill the gaps row by row
        for probe_name in data.index:
            one_row = data.loc[probe_name]
            if one_row.isna().any():
                data.loc[probe_name] = one_row.fillna(row_averages[probe_name])

        print("Missing values filled with the probe average")

    # Save a small report
    report = pd.DataFrame()
    report["item"] = ["missing_values", "percent_missing"]
    report["value"] = [number_missing, round(percent_missing, 6)]
    report.to_csv(RESULT_FOLDER + "/step2_missing_values.csv", index=False)

    return data


def remove_control_probes(data):
    # Some probes are not real genes. They are control spots used to
    # check that the chip worked. Their names start with "AFFX".
    # We remove them because they carry no biological meaning.

    keep_rows = []
    number_removed = 0

    for probe_name in data.index:
        if str(probe_name).upper().startswith("AFFX"):
            keep_rows.append(False)
            number_removed = number_removed + 1
        else:
            keep_rows.append(True)

    data = data[keep_rows]

    print("Removed", number_removed, "control probes")
    print("Probes left:", data.shape[0])
    return data


def normalize_data(data):
    # Different chips can be brighter or dimmer overall, even for the
    # same patient biology. Quantile normalisation forces every
    # patient to have the SAME shape of value distribution, so we are
    # comparing biology and not chip brightness.
    #
    # The idea in three steps:
    #   1. Sort each patient's values, then average those sorted lists.
    #   2. That average list is our "target shape".
    #   3. For each patient, replace their smallest value with the
    #      smallest target value, second smallest with second, and so on.

    print("Normalising the chips...")

    values = data.values.astype("float32")
    number_of_probes = values.shape[0]
    number_of_patients = values.shape[1]

    # Step 1: add up all the sorted lists
    total = np.zeros(number_of_probes)

    for patient in range(number_of_patients):
        one_patient = values[:, patient]
        sorted_values = np.sort(one_patient)
        total = total + sorted_values

    # Step 2: the target shape is the average sorted list
    target_shape = total / number_of_patients

    # Step 3: give each patient the target shape, keeping their order
    for patient in range(number_of_patients):
        one_patient = values[:, patient]

        # argsort tells us the position of the smallest, 2nd smallest...
        order = np.argsort(one_patient)

        # ranks[position] = "this value is the Nth smallest"
        ranks = np.empty(number_of_probes, dtype=int)
        ranks[order] = np.arange(number_of_probes)

        # Look up the target value for each rank
        values[:, patient] = target_shape[ranks]

    # Put the numbers back into a table with the same names
    data = pd.DataFrame(values, index=data.index, columns=data.columns)
    print("Normalising done")
    return data


def check_quality(data):
    # A good chip should look similar to all the other chips.
    # We measure this with correlation (a number from -1 to 1;
    # closer to 1 means more similar).
    #
    # Any patient far below average is FLAGGED but NOT deleted,
    # because a rare disease type can look different for real
    # biological reasons.

    print("Checking chip quality...")

    # Using all 54,000 probes would be slow, so we take a random
    # sample of 4000 probes. random_state keeps it repeatable.
    random_maker = np.random.default_rng(SEED)
    how_many = min(4000, data.shape[0])
    chosen_rows = random_maker.choice(data.shape[0], how_many, replace=False)
    small_data = data.values[chosen_rows]

    # Correlation between every pair of patients
    correlation_grid = np.corrcoef(small_data.T.astype("float64"))

    # A patient always correlates 1.0 with itself. That is not
    # informative, so we blank out the diagonal.
    np.fill_diagonal(correlation_grid, np.nan)

    # Average similarity of each patient to everybody else
    average_similarity = np.nanmean(correlation_grid, axis=1)

    overall_average = average_similarity.mean()
    spread = average_similarity.std()

    print("Average similarity between chips:", round(overall_average, 3))
    print("Lowest similarity:", round(average_similarity.min(), 3))

    # How many standard deviations below average is each patient?
    z_score = (average_similarity - overall_average) / spread

    quality = pd.DataFrame()
    quality["sample"] = data.columns
    quality["mean_correlation"] = average_similarity
    quality["z_score"] = z_score

    # Mark as unusual if more than 5 standard deviations below average
    is_outlier = []
    for one_z in z_score:
        if one_z < -5.0:
            is_outlier.append(True)
        else:
            is_outlier.append(False)
    quality["outlier"] = is_outlier

    quality.to_csv(RESULT_FOLDER + "/step2_quality.csv", index=False)

    print("Flagged", sum(is_outlier), "unusual chips (flagged, NOT removed)")
    return quality


# =============================================================
# STEP 3 - PICK THE USEFUL PROBES
# =============================================================

def filter_probes(data):
    # We have ~54,000 probes but most are boring: either switched off
    # in everybody, or the same in everybody. We keep two kinds:
    #
    #   1. Probes that are switched ON in at least some patients.
    #   2. Of those, the TOP_PROBES most variable ones, because a
    #      probe that never changes cannot tell diseases apart.
    #
    # We measure variability with MAD (median absolute deviation).
    # MAD is like standard deviation but it is not fooled by one
    # single strange value.

    print("Filtering probes...")

    values = data.values
    starting_number = data.shape[0]

    # ---- Filter 1: keep probes with real signal ----
    # For each probe, find the value that 95% of patients are below.
    high_value = np.percentile(values, 95, axis=1)

    # Our cut-off is the middle value of the whole table
    cutoff = float(np.median(values))

    keep_rows = []
    for one_value in high_value:
        if one_value > cutoff:
            keep_rows.append(True)
        else:
            keep_rows.append(False)

    data = data[keep_rows]
    print("After signal filter:", data.shape[0], "of", starting_number, "kept")

    # ---- Filter 2: keep the most variable probes ----
    values = data.values

    # MAD, step by step
    middle_of_each_probe = np.median(values, axis=1, keepdims=True)
    distance_from_middle = np.abs(values - middle_of_each_probe)
    mad = np.median(distance_from_middle, axis=1)

    # Sort MAD from big to small and take the top ones
    order_big_to_small = np.argsort(mad)[::-1]
    how_many = min(TOP_PROBES, data.shape[0])
    chosen = order_big_to_small[0:how_many]
    chosen = np.sort(chosen)          # keep the original row order

    data = data.iloc[chosen]
    print("After variability filter:", data.shape[0], "probes kept")

    return data


def zscore_data(data):
    # IMPORTANT: this must happen AFTER the variability filter.
    # Z-scoring makes every probe have the same variability, so if we
    # did it first the variability filter would become meaningless.
    #
    # Z-score means: for each probe, subtract its average and divide
    # by its spread. Then every probe is on the same scale and no
    # single loud probe dominates the clustering.

    print("Z-scoring the probes...")

    probe_average = data.mean(axis=1)
    probe_spread = data.std(axis=1)

    # If a probe never changes its spread is 0, and dividing by 0
    # is not allowed. Replace those zeros with 1.
    probe_spread = probe_spread.replace(0, 1)

    centred = data.sub(probe_average, axis=0)
    scaled = centred.div(probe_spread, axis=0)

    print("Check: average should be near 0 ->", round(float(scaled.values.mean()), 6))
    print("Check: spread should be near 1  ->", round(float(scaled.values.std()), 3))

    return scaled


# =============================================================
# STEP 4 - GROUP SIMILAR PATIENTS (CLUSTERING)
# =============================================================

def measure_tightness(patient_grid, tree):
    # Score each possible number of groups by TIGHTNESS.
    #
    # The silhouette score asks: is each patient closer to its own
    # group than to the nearest other group? Higher is better.
    #
    # This works well for data that forms neat round balls. Gene
    # expression data does not, so read choose_number_of_groups()
    # before trusting the winner here.

    print("Measuring tightness (silhouette) for each k...")

    k_list = []
    score_list = []

    for k in range(2, 21):
        # Cut the tree so we get exactly k groups
        labels = fcluster(tree, t=k, criterion="maxclust")

        # Score this answer. sample_size makes it faster.
        score = silhouette_score(patient_grid, labels,
                                 sample_size=min(2000, len(labels)),
                                 random_state=SEED)

        k_list.append(k)
        score_list.append(score)
        print("  k =", k, " tightness =", round(score, 3))

    scores = pd.DataFrame()
    scores["k"] = k_list
    scores["silhouette"] = score_list
    scores.to_csv(RESULT_FOLDER + "/step4_tightness.csv", index=False)

    return scores


def measure_stability(patient_grid):
    # Score each possible number of groups by STABILITY.
    #
    # THE IDEA
    # If a group of patients is real, then throwing away a random
    # fifth of the data should not change it much. The same patients
    # should keep landing together.
    #
    # So we do this many times:
    #   1. Pick a random 80% of the patients.
    #   2. Cluster just those, into k groups.
    #   3. Write down which pairs of patients ended up together.
    #
    # Afterwards, for every pair of patients we can work out:
    #   "out of all the times BOTH of them were picked,
    #    what fraction of the time did they land together?"
    #
    # A good answer gives us fractions near 1 (always together) or
    # near 0 (never together). Fractions in the middle, like 0.5,
    # mean "we cannot decide about this pair" - that is instability.
    #
    # PAC = Proportion of Ambiguous Clustering
    #     = the fraction of pairs stuck in the middle (0.1 to 0.9).
    # LOWER PAC IS BETTER.

    print("Measuring stability (PAC) for each k...")
    print("This is the slow part - about 5 to 10 minutes. Please wait.")

    number_of_patients = patient_grid.shape[0]
    how_many_to_take = int(0.8 * number_of_patients)

    k_list = []
    pac_list = []

    # We only want each pair once, not twice. triu_indices gives us
    # the top-right half of the square, skipping the diagonal.
    top_half = np.triu_indices(number_of_patients, k=1)

    for k in range(2, CONSENSUS_MAX_K + 1):

        # IMPORTANT: each k gets its OWN fresh set of random samples.
        #
        # It is tempting to draw 20 samples once and reuse them for
        # every k, because building the trees is the slow part. Do not
        # do that. With only 20 repeats, one lucky draw can make one
        # particular k look better than it really is, and because
        # every k shares that same lucky draw the mistake does not
        # average out. We tried it and it picked the wrong k.
        #
        # SEED + k means "a different but repeatable set for each k".
        random_maker = np.random.default_rng(SEED + k)

        # How many times each pair landed in the same group
        times_together = np.zeros((number_of_patients, number_of_patients),
                                 dtype="float32")

        # How many times each pair was even picked together. This is
        # the fair denominator: a pair that was rarely picked should
        # not look unstable just because of that.
        times_picked = np.zeros((number_of_patients, number_of_patients),
                               dtype="float32")

        for repeat in range(CONSENSUS_REPEATS):
            # Pick a random 80%. replace=False means no patient twice.
            chosen = random_maker.choice(number_of_patients,
                                         how_many_to_take, replace=False)
            small_grid = patient_grid[chosen]

            # Build a tree for just these patients and cut it into k
            tree = linkage(small_grid, method="ward")
            labels = fcluster(tree, t=k, criterion="maxclust")

            # Compare every patient's label with every other patient's
            # label, all at once. The result is a True/False square.
            same_group = (labels[:, None] == labels[None, :])
            same_group = same_group.astype("float32")

            # np.ix_ picks out the little square of the big grid that
            # belongs to the chosen patients, so we can add to all
            # those pairs at once instead of looping over millions.
            times_together[np.ix_(chosen, chosen)] += same_group
            times_picked[np.ix_(chosen, chosen)] += 1.0

        # ---- Turn the counts into a PAC score ----
        # Fraction of the time each pair landed together.
        # np.maximum avoids dividing by zero for pairs never picked.
        fraction = times_together / np.maximum(times_picked, 1.0)

        pair_values = fraction[top_half]

        # Count the undecided pairs
        is_ambiguous = (pair_values > 0.1) & (pair_values < 0.9)
        pac = float(np.mean(is_ambiguous))

        k_list.append(k)
        pac_list.append(pac)
        print("  k =", k, " PAC =", round(pac, 4), "(lower is better)")

    scores = pd.DataFrame()
    scores["k"] = k_list
    scores["pac"] = pac_list
    scores.to_csv(RESULT_FOLDER + "/step4_stability.csv", index=False)

    return scores


def measure_cluster_confidence(patient_grid, group_labels, best_k):
    # HOW SURE ARE WE ABOUT EACH INDIVIDUAL PATIENT?
    #
    # A cluster can look fine overall while some patients inside it sit
    # right on the border. We measure that per patient: redo the
    # clustering on random 80% samples, then ask how often this patient
    # landed with the OTHER members of its final cluster.
    #
    # Near 1.0 = this patient always belongs here. Confident.
    # Near 0.5 = it could easily have gone elsewhere. Uncertain.
    #
    # IMPORTANT DETAIL: we exclude the patient itself from the average.
    # A patient always lands with itself, so including it would push
    # every score up - badly for small clusters, where "itself" is a
    # large share of the group.

    print("Measuring how confident we are about each patient...")

    number_of_patients = patient_grid.shape[0]
    how_many_to_take = int(0.8 * number_of_patients)
    random_maker = np.random.default_rng(SEED + best_k)

    times_together = np.zeros((number_of_patients, number_of_patients),
                             dtype="float32")
    times_picked = np.zeros((number_of_patients, number_of_patients),
                            dtype="float32")

    for repeat in range(CONSENSUS_REPEATS):
        chosen = random_maker.choice(number_of_patients,
                                     how_many_to_take, replace=False)
        small_grid = patient_grid[chosen]
        tree = linkage(small_grid, method="ward")
        labels = fcluster(tree, t=best_k, criterion="maxclust")

        same_group = (labels[:, None] == labels[None, :]).astype("float32")
        times_together[np.ix_(chosen, chosen)] += same_group
        times_picked[np.ix_(chosen, chosen)] += 1.0

    fraction = times_together / np.maximum(times_picked, 1.0)

    confidence = np.zeros(number_of_patients)

    for i in range(number_of_patients):
        # Who else is in my final cluster?
        my_cluster = group_labels[i]
        peer_positions = []

        for j in range(number_of_patients):
            if j != i and group_labels[j] == my_cluster:
                peer_positions.append(j)

        if len(peer_positions) == 0:
            confidence[i] = np.nan          # alone in its cluster
        else:
            confidence[i] = fraction[i, peer_positions].mean()

    average_confidence = float(np.nanmean(confidence))

    how_many_unsure = 0
    for one_value in confidence:
        if one_value < 0.6:
            how_many_unsure = how_many_unsure + 1

    print("Average confidence:", round(average_confidence, 3))
    print(how_many_unsure, "of", number_of_patients,
          "patients are below 0.6 (uncertain)")

    return confidence


def export_tree_for_r(tree, patient_names):
    # Save the family tree in the form R's plotting code expects.
    #
    # WHY THIS IS FIDDLY. Python (scipy) and R (hclust) both store a
    # tree as a list of "these two things joined here", but they
    # number things differently:
    #
    #   Python: patients are 0,1,2...  joins are n, n+1, n+2...
    #   R:      patients are NEGATIVE (-1,-2,-3...), joins are POSITIVE
    #
    # If we get this wrong R draws a tree that looks fine but is
    # wrong, so we convert carefully.
    #
    # We also only save 400 patients. A tree with 2096 branches is an
    # unreadable black smear on paper. The saved picture is a sample;
    # the real cluster numbers still come from all the patients.

    number_of_patients = len(patient_names)

    join_a = []
    join_b = []

    for i in range(number_of_patients - 1):
        this_pair = []

        for side in range(2):
            item = int(tree[i, side])

            if item < number_of_patients:
                # It is a patient. R wants a negative 1-based number.
                this_pair.append(-(item + 1))
            else:
                # It is an earlier join. R wants a positive number.
                this_pair.append(item - number_of_patients + 1)

        # R lists patients before joins on each row
        first = this_pair[0]
        second = this_pair[1]

        if first > 0 and second < 0:
            first, second = second, first
        elif first > 0 and second > 0 and first > second:
            first, second = second, first
        elif first < 0 and second < 0 and abs(first) > abs(second):
            first, second = second, first

        join_a.append(first)
        join_b.append(second)

    merge_table = pd.DataFrame()
    merge_table["a"] = join_a
    merge_table["b"] = join_b
    merge_table.to_csv(RESULT_FOLDER + "/step4_tree_merge.csv", index=False)

    height_table = pd.DataFrame()
    height_table["height"] = tree[:, 2]
    height_table.to_csv(RESULT_FOLDER + "/step4_tree_height.csv", index=False)

    order_table = pd.DataFrame()
    order_table["order"] = leaves_list(tree) + 1
    order_table.to_csv(RESULT_FOLDER + "/step4_tree_order.csv", index=False)

    label_table = pd.DataFrame()
    label_table["label"] = patient_names
    label_table.to_csv(RESULT_FOLDER + "/step4_tree_labels.csv", index=False)

    print("Tree saved for R (", number_of_patients, "patients )")


def choose_number_of_groups(patient_grid, tree):
    # WHY THERE ARE TWO ANSWERS, AND WHY WE PICK THE SECOND ONE
    #
    # Tightness (silhouette) says k = 3.
    # Stability (PAC) says k = 5.
    #
    # They disagree because they measure different things.
    #
    # Tightness wants groups shaped like neat round balls, all far
    # apart. Gene expression data is almost never shaped like that,
    # so every tightness score comes out low (0.06 to 0.20 here) and
    # the "winner" is just whichever k mashes the data into the
    # fewest big lumps. With 17 real leukaemia types in this data,
    # 3 lumps is clearly too few.
    #
    # Stability asks a more useful question: if I redo this on
    # different patients, do I get the same groups again? That is
    # what we actually mean by "is this a real group?"
    #
    # We ignore k = 2 when looking for the best stability. Splitting
    # anything into two halves is very repeatable, but it tells us
    # almost nothing.
    #
    # The disagreement is not a problem to hide - it is a finding.
    # Report both numbers and say which one you used and why.

    tightness_scores = measure_tightness(patient_grid, tree)

    # ---- Find the best k by tightness ----
    best_tightness_k = 0
    best_tightness = -999.0

    for row in tightness_scores.itertuples(index=False):
        if row.silhouette > best_tightness:
            best_tightness = row.silhouette
            best_tightness_k = row.k

    # ---- If we only want tightness, stop here ----
    if HOW_TO_CHOOSE_K == "tightness":
        print("")
        print("Chosen k =", best_tightness_k, "(by tightness)")
        return best_tightness_k

    # ---- Otherwise also measure stability ----
    stability_scores = measure_stability(patient_grid)

    best_stability_k = 0
    best_pac = 999.0

    for row in stability_scores.itertuples(index=False):
        if row.k < 3:
            continue                    # skip the boring k = 2 answer
        if row.pac < best_pac:
            best_pac = row.pac
            best_stability_k = row.k

    # ---- Report both, then choose ----
    print("")
    print("Two ways of choosing k, two answers:")
    print("  tightness  says k =", best_tightness_k,
          "(score", round(best_tightness, 3), ")")
    print("  stability  says k =", best_stability_k,
          "(PAC", round(best_pac, 3), ")")

    if best_tightness_k != best_stability_k:
        print("They disagree. This belongs in your write-up.")
        print("We use stability, because tightness is unreliable")
        print("on gene expression data. See the comments above.")

    print("Chosen k =", best_stability_k, "(by stability)")

    # Save both answers side by side for the report
    summary = pd.DataFrame()
    summary["measure"] = ["best_k_tightness", "best_k_stability", "k_used"]
    summary["value"] = [best_tightness_k, best_stability_k, best_stability_k]
    summary.to_csv(RESULT_FOLDER + "/step4_k_choice.csv", index=False)

    return best_stability_k


def do_clustering(data):
    # Hierarchical clustering builds a family tree of patients:
    # the two most similar patients join first, then the next two,
    # and so on until everybody is in one big tree.
    #
    # Then we "cut" the tree to get a chosen number of groups.

    print("Clustering the patients...")

    # Our table is probes x patients, but clustering wants
    # patients x probes, so we flip it with .T (transpose)
    patient_grid = data.values.T.astype("float64")

    print("Measuring the distance between every pair of patients...")
    distances = pdist(patient_grid, metric="euclidean")

    print("Building the tree (Ward method)...")
    tree = linkage(distances, method="ward")

    # ---- Decide how many groups to cut the tree into ----
    if FORCE_K > 0:
        best_k = FORCE_K
        print("Using k =", best_k, "because FORCE_K was set at the top")
    else:
        best_k = choose_number_of_groups(patient_grid, tree)

    # ---- Save a small version of the tree so R can draw it ----
    # 400 patients only: a tree with 2096 branches is unreadable.
    random_maker = np.random.default_rng(SEED)
    how_many_to_draw = min(400, patient_grid.shape[0])
    chosen = np.sort(random_maker.choice(patient_grid.shape[0],
                                        how_many_to_draw, replace=False))
    small_tree = linkage(patient_grid[chosen], method="ward")
    export_tree_for_r(small_tree, data.columns.values[chosen])

    # Cut the tree into groups
    tree_groups = fcluster(tree, t=best_k, criterion="maxclust")

    # Do it a second way (k-means) as a double check.
    # If two different methods agree, we trust the answer more.
    print("Running k-means as a second opinion...")
    kmeans_model = KMeans(n_clusters=best_k, n_init=25, random_state=SEED)
    kmeans_groups = kmeans_model.fit_predict(patient_grid) + 1

    # ARI compares two groupings. 1.0 = identical, 0 = random.
    agreement = adjusted_rand_score(tree_groups, kmeans_groups)
    print("Agreement between the two methods (ARI):", round(agreement, 3))

    confidence = measure_cluster_confidence(patient_grid, tree_groups, best_k)

    clusters = pd.DataFrame()
    clusters["sample"] = data.columns
    clusters["tree_cluster"] = tree_groups
    clusters["kmeans_cluster"] = kmeans_groups
    clusters["confidence"] = confidence
    clusters.to_csv(RESULT_FOLDER + "/step4_clusters.csv", index=False)

    # Print how big each group is
    for group_number in range(1, best_k + 1):
        count = 0
        for one_label in tree_groups:
            if one_label == group_number:
                count = count + 1
        print("  cluster", group_number, "has", count, "patients")

    # ---- Save one summary row, so the final report can quote it ----
    summary = pd.DataFrame()
    summary["item"] = ["k_used", "tree_vs_kmeans_ari", "average_confidence"]
    summary["value"] = [best_k, round(agreement, 4),
                        round(float(np.nanmean(confidence)), 4)]
    summary.to_csv(RESULT_FOLDER + "/step4_summary.csv", index=False)

    return clusters


# =============================================================
# STEP 5 - PCA (SEE THE DATA IN 2D)
# =============================================================

def do_pca(data, clusters):
    # We have 5000 numbers per patient. Nobody can picture 5000
    # dimensions. PCA squashes them down to a few new numbers
    # (PC1, PC2, PC3...) that keep as much of the variation as
    # possible, so we can draw a normal 2D scatter plot.

    print("Running PCA...")

    patient_grid = data.values.T.astype("float64")

    pca_model = PCA(n_components=10, random_state=SEED)
    new_coordinates = pca_model.fit_transform(patient_grid)

    # How much of the story does each new number tell?
    variance_share = pca_model.explained_variance_ratio_
    first_three = 100.0 * (variance_share[0] + variance_share[1] + variance_share[2])
    print("PC1 + PC2 + PC3 explain", round(first_three, 1), "% of the variation")

    # Save the coordinates so R can plot them
    pca_table = pd.DataFrame()
    pca_table["sample"] = data.columns
    for i in range(10):
        column_name = "PC" + str(i + 1)
        pca_table[column_name] = new_coordinates[:, i]

    # Attach the cluster number of each patient
    pca_table = pca_table.merge(clusters, on="sample", how="left")
    pca_table.to_csv(RESULT_FOLDER + "/step5_pca.csv", index=False)

    variance_table = pd.DataFrame()
    variance_table["pc"] = range(1, 11)
    variance_table["variance_explained"] = variance_share
    variance_table.to_csv(RESULT_FOLDER + "/step5_pca_variance.csv", index=False)

    # ---- Which probes DRIVE each new direction? ----
    # PC1 is a mixture of all 5000 probes, but a few contribute most.
    # Those are called the loadings, and they tell us what biological
    # story PC1 is telling.

    probe_names = data.index.values
    loading_rows = []

    for pc_number in range(3):
        weights = pca_model.components_[pc_number]

        # Biggest weights, ignoring plus or minus
        size_of_weight = np.abs(weights)
        order_big_to_small = np.argsort(size_of_weight)[::-1]
        top_50 = order_big_to_small[0:50]

        for position in top_50:
            loading_rows.append({"pc": pc_number + 1,
                                 "probe": probe_names[position],
                                 "loading": weights[position]})

    pd.DataFrame(loading_rows).to_csv(
        RESULT_FOLDER + "/step5_pca_loadings.csv", index=False)
    print("Top probes driving PC1-PC3 saved")

    return pca_table


# =============================================================
# STEP 6 - TURN PROBES INTO GENES
# =============================================================

def map_probes_to_genes(data):
    # Several probes can measure the SAME gene. If we leave them all
    # in, one gene appears many times in our results.
    #
    # Rule: for each gene, keep the probe with the highest average
    # value (the strongest, most reliable measurement).
    #
    # We also need gene names because pathway databases (GO, KEGG)
    # only know gene names, not probe codes like "205548_s_at".

    print("Turning probes into genes...")

    if os.path.exists(ANNOTATION_FILE) == False:
        print("ERROR:", ANNOTATION_FILE, "not found.")
        print("Make it first by running: Rscript simple_make_annotation.R")
        raise SystemExit(1)

    # The annotation file says which gene each probe measures
    annotation = pd.read_csv(ANNOTATION_FILE)
    probe_to_gene = {}
    for row in annotation.itertuples(index=False):
        probe_to_gene[str(row.probe)] = str(row.symbol)

    # Go through our probes and keep only the ones we have a name for
    average_value = data.mean(axis=1)

    best_probe_for_gene = {}     # gene name -> probe name
    best_average_for_gene = {}   # gene name -> that probe's average

    for probe_name in data.index:
        probe_name = str(probe_name)

        if probe_name not in probe_to_gene:
            continue

        gene_name = probe_to_gene[probe_name]
        this_average = average_value[probe_name]

        # First time we see this gene, or a better probe than before?
        if gene_name not in best_average_for_gene:
            best_probe_for_gene[gene_name] = probe_name
            best_average_for_gene[gene_name] = this_average
        elif this_average > best_average_for_gene[gene_name]:
            best_probe_for_gene[gene_name] = probe_name
            best_average_for_gene[gene_name] = this_average

    print("Found gene names for", len(best_probe_for_gene), "genes")

    # Build the new gene-level table
    chosen_probes = []
    chosen_genes = []
    for gene_name in best_probe_for_gene:
        chosen_genes.append(gene_name)
        chosen_probes.append(best_probe_for_gene[gene_name])

    gene_data = data.loc[chosen_probes]
    gene_data.index = chosen_genes
    gene_data = gene_data.sort_index()

    # Remove genes that are switched off everywhere
    values = gene_data.values
    high_value = np.percentile(values, 95, axis=1)
    cutoff = float(np.median(values))

    keep_rows = []
    for one_value in high_value:
        if one_value > cutoff:
            keep_rows.append(True)
        else:
            keep_rows.append(False)

    gene_data = gene_data[keep_rows]

    print("Final gene table:", gene_data.shape[0], "genes x",
          gene_data.shape[1], "patients")
    print("These genes are the 'background list' for pathway analysis")

    # Save the list of genes we are working with
    gene_list = pd.DataFrame()
    gene_list["gene"] = gene_data.index
    gene_list.to_csv(RESULT_FOLDER + "/step6_gene_list.csv", index=False)

    return gene_data


# =============================================================
# STEP 7 - MAKE THE PATIENT GROUPS (THE MOST IMPORTANT STEP)
# =============================================================

def make_groups(gene_data, patients, clusters):
    # THIS IS THE MOST IMPORTANT STEP. Read the reason carefully.
    #
    # We could group patients in two ways:
    #   (a) by the DOCTOR'S DIAGNOSIS  <- we use this
    #   (b) by the clusters we found ourselves in Step 4
    #
    # If we used (b), we would be testing our own answer against
    # itself. Every p-value and accuracy score would look amazing
    # and mean nothing. This is called circular reasoning.
    #
    # Using (a) means the labels came from outside our analysis,
    # so our results are honest.

    print("Making patient groups from the DOCTOR'S DIAGNOSIS")

    # Line up the patient information with our gene table
    patients = patients.reindex(gene_data.columns)
    clusters = clusters.set_index("sample").reindex(gene_data.columns)

    disease = patients["disease"].astype(str)

    # Count how many patients have each disease
    disease_counts = disease.value_counts()
    print("Found", len(disease_counts), "different diagnoses")

    # A group needs at least MIN_GROUP_SIZE patients, otherwise the
    # statistics are not trustworthy.
    big_enough_diseases = []
    for disease_name in disease_counts.index:
        if disease_counts[disease_name] >= MIN_GROUP_SIZE:
            big_enough_diseases.append(disease_name)

    can_analyse = []
    for one_disease in disease:
        if one_disease in big_enough_diseases:
            can_analyse.append(True)
        else:
            can_analyse.append(False)

    print("Keeping", len(big_enough_diseases), "groups with at least",
          MIN_GROUP_SIZE, "patients")
    print("That is", sum(can_analyse), "of", len(disease), "patients")

    # ---- Safety check ----
    # This dataset should have about 17 disease groups. If we only
    # found two or three, we almost certainly read the wrong line from
    # the file (tissue instead of diagnosis). Better to stop now than
    # to run for an hour and produce meaningless answers.
    if len(disease_counts) < 5:
        print("")
        print("WARNING: only", len(disease_counts), "different diagnoses found.")
        print("This dataset should have about 18. The values are:")
        for one_name in disease_counts.index:
            print("  ", one_name)
        print("If these look like tissue types rather than leukaemia")
        print("names, the wrong line was read from the data file.")
        print("")

    groups = pd.DataFrame()
    groups["sample"] = gene_data.columns
    groups["group"] = disease.values
    groups["group_source"] = "diagnosis"          # written down on purpose
    groups["tree_cluster"] = clusters["tree_cluster"].values
    groups["can_analyse"] = can_analyse
    groups.to_csv(RESULT_FOLDER + "/step7_groups.csv", index=False)

    # ---- Check our clustering against the real diagnoses ----
    # The clustering never saw these labels, so this is a fair test.
    disease_as_numbers = disease.astype("category").cat.codes.values
    score = adjusted_rand_score(disease_as_numbers,
                                clusters["tree_cluster"].values)
    print("Do our clusters match the real diagnoses? ARI =", round(score, 3))

    # NMI is a second agreement measure. Like ARI, 1.0 means perfect
    # agreement and 0 means none, but it is based on shared
    # information rather than counting pairs. Reporting both is
    # standard, because they can disagree slightly.
    nmi = normalized_mutual_info_score(disease_as_numbers,
                                       clusters["tree_cluster"].values)
    print("Same question measured a second way (NMI) =", round(nmi, 3))

    # A table showing which cluster contains which diseases
    comparison = pd.crosstab(disease, clusters["tree_cluster"])
    comparison.to_csv(RESULT_FOLDER + "/step7_cluster_vs_diagnosis.csv")

    # ---- How PURE is each cluster? ----
    # Purity asks: of everybody in this cluster, what fraction share
    # the single most common diagnosis? A pure cluster of 1.0 holds
    # one disease only. A low value means the cluster mixes diseases.

    cluster_numbers = []
    cluster_sizes = []
    main_diseases = []
    purities = []

    for one_cluster in comparison.columns:
        column = comparison[one_cluster]
        total_here = column.sum()
        biggest = column.max()
        main_disease = column.idxmax()

        cluster_numbers.append(one_cluster)
        cluster_sizes.append(total_here)
        main_diseases.append(main_disease)
        purities.append(biggest / total_here)

    purity_table = pd.DataFrame()
    purity_table["cluster"] = cluster_numbers
    purity_table["n_patients"] = cluster_sizes
    purity_table["main_diagnosis"] = main_diseases
    purity_table["purity"] = purities
    purity_table.to_csv(RESULT_FOLDER + "/step7_cluster_purity.csv",
                        index=False)

    print("How pure is each cluster?")
    for row in purity_table.itertuples(index=False):
        print("  cluster", row.cluster, "-", row.n_patients, "patients,",
              "mostly", str(row.main_diagnosis)[0:28],
              "- purity", round(row.purity, 2))

    # Save the agreement numbers for the report
    agreement_table = pd.DataFrame()
    agreement_table["measure"] = ["adjusted_rand_index",
                                  "normalized_mutual_info",
                                  "average_purity"]
    agreement_table["value"] = [score, nmi, purity_table["purity"].mean()]
    agreement_table.to_csv(RESULT_FOLDER + "/step7_agreement.csv", index=False)

    return groups


def check_for_batch_effects(pca_table, patients):
    # THE BATCH EFFECT CHECK. Do not skip this one.
    #
    # PC1 is the single biggest source of variation in the data. We
    # HOPE it is biology (the disease). But it could easily be
    # something technical instead: which machine ran the chip, which
    # hospital sent the sample, which tissue was taken.
    #
    # If a technical detail explains PC1 better than the diagnosis
    # does, our "biological findings" may partly be measuring
    # equipment. That must be reported, not hidden.
    #
    # THE TEST. For each fact we know about the patients, we ask:
    # do the different values of this fact have different PC scores?
    # The Kruskal-Wallis test answers that without assuming the
    # numbers follow a bell curve. Then we turn its H statistic into
    # eta-squared, a share from 0 to 1: "how much of this PC does
    # this fact explain?"

    print("Checking whether anything technical explains the data...")

    patients = patients.reindex(pca_table["sample"].values)

    result_rows = []

    for column_name in patients.columns:
        values = patients[column_name].astype(str)

        # Only test facts with a sensible number of different values.
        # One value explains nothing; hundreds of values (like a
        # patient ID) would fit anything by accident.
        how_many_different = len(set(values))
        if how_many_different < 2 or how_many_different > 60:
            continue

        for pc_number in range(1, 6):
            pc_column = "PC" + str(pc_number)
            scores = pca_table[pc_column].values

            # Split the PC scores into one list per value of this fact
            score_groups = []
            for one_value in sorted(set(values)):
                picked = scores[values.values == one_value]
                if len(picked) >= 3:
                    score_groups.append(picked)

            if len(score_groups) < 2:
                continue

            h_value, p_value = stats.kruskal(*score_groups)

            # eta-squared: what share of this PC does the fact explain?
            n = len(scores)
            n_groups = len(score_groups)
            share = (h_value - n_groups + 1) / (n - n_groups)
            if share < 0:
                share = 0.0

            result_rows.append({"fact": column_name,
                                "pc": pc_number,
                                "share_explained": share,
                                "p_value": p_value})

    if len(result_rows) == 0:
        print("Not enough patient facts to test. Skipping.")
        return

    table = pd.DataFrame(result_rows)
    table = table.sort_values(["pc", "share_explained"], ascending=[True, False])
    table.to_csv(RESULT_FOLDER + "/step7_batch_check.csv", index=False)

    print("What explains PC1, the biggest source of variation:")
    pc1_only = table[table["pc"] == 1]
    for row in pc1_only.head(5).itertuples(index=False):
        print("  ", row.fact, "explains", round(100 * row.share_explained, 1),
              "% of PC1")

    # ---- The verdict ----
    best_fact = pc1_only["fact"].iloc[0]
    if best_fact == "disease":
        print("GOOD: the diagnosis explains PC1 better than anything")
        print("technical. That is what we want to see.")
    else:
        print("WARNING: '" + best_fact + "' explains PC1 better than the")
        print("diagnosis does. This may be a batch effect. Report it")
        print("honestly in your write-up and discuss the confound.")


# =============================================================
# STEP 8 - FIND GENES THAT DIFFER BETWEEN GROUPS
# =============================================================

def correct_pvalues(pvalues):
    # We test ~15,000 genes at once. If we accept p < 0.05, then by
    # pure luck about 750 genes look "significant" when they are not.
    #
    # The Benjamini-Hochberg method makes the p-values stricter to
    # account for this. The result is called an FDR.
    #
    # How it works:
    #   1. Sort the p-values from small to big.
    #   2. Multiply each by (total tests / its rank).
    #   3. Walk backwards making sure the list never goes down.

    pvalues = np.asarray(pvalues, dtype="float64")
    total = len(pvalues)

    order_small_to_big = np.argsort(pvalues)
    adjusted = np.zeros(total)

    smallest_so_far = 1.0

    # Walk from the biggest p-value back to the smallest
    for i in range(total - 1, -1, -1):
        position = order_small_to_big[i]
        rank = i + 1

        new_value = pvalues[position] * total / rank

        if new_value < smallest_so_far:
            smallest_so_far = new_value

        adjusted[position] = smallest_so_far

    # p-values can never be above 1
    adjusted = np.clip(adjusted, 0, 1)
    return adjusted


def find_different_genes(gene_data, groups):
    # For every disease group we ask, for every gene:
    # "Is this gene more active in THIS disease than in all the
    #  other diseases put together?"
    #
    # This is called a one-versus-rest comparison.
    #
    # We use Welch's t-test. A normal t-test assumes both sides have
    # the same spread; Welch's does not. Our groups have very
    # different sizes (28 to 448 patients), so Welch's is safer.

    print("Comparing each disease group against all the others...")

    group_names = groups["group"].astype(str).values
    can_analyse = groups["can_analyse"].values

    # The list of groups we are allowed to test
    groups_to_test = []
    for i in range(len(group_names)):
        if can_analyse[i] == True:
            if group_names[i] not in groups_to_test:
                groups_to_test.append(group_names[i])
    groups_to_test.sort()

    all_genes = gene_data.index.values
    values = gene_data.values.astype("float64")

    # ---- Warning about the numbers in this dataset ----
    # This data was already squashed to the range 0 to 1 before it
    # was shared. So the difference between two averages is NOT a
    # "fold change". It is just a difference on a 0-to-1 scale, and
    # it is always small. Never call it log2 fold change.
    smallest = float(values.min())
    biggest = float(values.max())
    print("Data range:", round(smallest, 2), "to", round(biggest, 2))
    if biggest <= 1.01:
        print("NOTE: data is on a 0-1 scale, so differences are NOT fold changes")

    all_results = []

    for one_group in groups_to_test:
        # Which patients are in this group, and which are not?
        in_group = []
        for one_name in group_names:
            if one_name == one_group:
                in_group.append(True)
            else:
                in_group.append(False)
        in_group = np.array(in_group)

        inside = values[:, in_group]
        outside = values[:, in_group == False]

        # Welch's t-test for every gene at once
        t_value, p_value = stats.ttest_ind(inside, outside,
                                           axis=1, equal_var=False)

        # If a test failed we get NaN. Treat it as "not significant".
        p_value = np.nan_to_num(p_value, nan=1.0)

        average_inside = inside.mean(axis=1)
        average_outside = outside.mean(axis=1)
        difference = average_inside - average_outside

        one_result = pd.DataFrame()

        # CAREFUL: put a LIST column in first, before any single values.
        # A brand new pd.DataFrame() has zero rows. If you put a single
        # value like "CLL" into it first, pandas has no rows to put it
        # in, so it makes an empty column. Adding the gene list
        # afterwards creates the rows, but "CLL" is already gone and
        # you silently get a column full of blanks.
        one_result["gene"] = all_genes

        # Now the table has rows, so single values fill every row
        one_result["group"] = one_group
        one_result["patients_in_group"] = int(in_group.sum())

        one_result["average_inside"] = average_inside
        one_result["average_outside"] = average_outside
        one_result["difference"] = difference
        one_result["t_value"] = np.nan_to_num(t_value)
        one_result["p_value"] = p_value
        one_result["fdr"] = correct_pvalues(p_value)

        # Safety check: catch the blank-column mistake straight away
        # instead of finding out five steps later.
        if one_result["group"].isna().any():
            print("ERROR: the group column came out blank. Stopping.")
            raise SystemExit(1)

        all_results.append(one_result)

    # Stack all the group results into one big table
    results = pd.concat(all_results, ignore_index=True)

    # ---- Choose the size cut-off from the data itself ----
    # We cannot use the usual "fold change of 1" rule because our
    # numbers are on a 0-1 scale and never get that big. Instead we
    # keep the top 1% biggest differences we actually saw.
    size_cutoff = float(np.percentile(np.abs(results["difference"]), 99))
    print("Difference cut-off (top 1% of what we saw):", round(size_cutoff, 4))

    # Save the cut-off so later steps use the SAME number
    cutoff_table = pd.DataFrame()
    cutoff_table["cutoff"] = [size_cutoff]
    cutoff_table["fdr"] = [FDR_CUTOFF]
    cutoff_table["rule"] = ["99th percentile of observed differences"]
    cutoff_table.to_csv(RESULT_FOLDER + "/step8_cutoff.csv", index=False)

    results.to_csv(RESULT_FOLDER + "/step8_all_genes.csv", index=False)

    # ---- Count the winners in each group, UP and DOWN ----
    # Down-regulated genes matter too: a gene switched OFF in one
    # disease is just as useful a marker as one switched on.
    summary_group = []
    summary_n = []
    summary_up = []
    summary_down = []

    for one_group in groups_to_test:
        one_part = results[results["group"] == one_group]

        is_significant = one_part["fdr"] < FDR_CUTOFF
        is_big_and_up = one_part["difference"] >= size_cutoff
        is_big_and_down = one_part["difference"] <= -size_cutoff

        number_up = int((is_significant & is_big_and_up).sum())
        number_down = int((is_significant & is_big_and_down).sum())

        summary_group.append(one_group)
        summary_n.append(int(one_part["patients_in_group"].iloc[0]))
        summary_up.append(number_up)
        summary_down.append(number_down)

        print("  ", one_group[0:30], "-> up:", number_up,
              " down:", number_down)

    summary = pd.DataFrame()
    summary["group"] = summary_group
    summary["patients"] = summary_n
    summary["genes_up"] = summary_up
    summary["genes_down"] = summary_down
    summary.to_csv(RESULT_FOLDER + "/step8_summary.csv", index=False)

    # ---- Make a short list: the best 20 genes per disease ----
    # The full table has 250,000 rows. Nobody reads that. The short
    # list is what goes in the report.
    #
    # We rank by combining two things into one score:
    #   how sure we are (-log10 of the p-value)
    #   how big the difference is
    # A gene needs BOTH to score well.

    is_significant = results["fdr"] < FDR_CUTOFF
    is_big_and_up = results["difference"] >= size_cutoff
    winners = results[is_significant & is_big_and_up].copy()

    if len(winners) == 0:
        print("No genes passed. Try a lower cut-off.")
        return results, size_cutoff

    sureness = -np.log10(winners["p_value"] + 1e-300)
    winners["rank_score"] = sureness * winners["difference"]

    short_list = []
    for one_group in groups_to_test:
        just_this_group = winners[winners["group"] == one_group]
        just_this_group = just_this_group.sort_values("rank_score",
                                                     ascending=False)
        short_list.append(just_this_group.head(20))

    short_list = pd.concat(short_list, ignore_index=True)
    short_list.to_csv(RESULT_FOLDER + "/step8_shortlist.csv", index=False)
    print("Short list:", len(short_list), "genes saved for the report")

    # ---- Save a smaller file for the volcano plot ----
    # 250,000 dots is slow to draw and looks like a grey blob. We keep
    # every significant gene, plus a random 10% of the boring ones so
    # the background cloud still looks right.
    boring = results[is_significant == False]
    boring_sample = boring.sample(frac=0.1, random_state=SEED)

    volcano = pd.concat([results[is_significant], boring_sample])
    volcano = volcano[["group", "gene", "difference", "p_value", "fdr"]]
    volcano.to_csv(RESULT_FOLDER + "/step8_volcano.csv", index=False)

    # ---- Save a small table for the heatmap picture ----
    # The top few genes per disease, and a sample of patients from
    # each disease. Rows are z-scored so colours are comparable.
    heatmap_genes = []
    for one_group in groups_to_test:
        just_this_group = short_list[short_list["group"] == one_group]
        for one_gene in just_this_group["gene"].head(8):
            if one_gene not in heatmap_genes:
                heatmap_genes.append(one_gene)

    random_maker = np.random.default_rng(SEED)
    chosen_patients = []
    for one_group in groups_to_test:
        positions = np.where(group_names == one_group)[0]
        how_many = min(25, len(positions))
        picked = random_maker.choice(positions, how_many, replace=False)
        for one_position in picked:
            chosen_patients.append(one_position)
    chosen_patients = np.sort(np.array(chosen_patients))

    small = gene_data.loc[heatmap_genes].iloc[:, chosen_patients]

    # z-score each row so one loud gene does not wash out the colours
    row_average = small.mean(axis=1)
    row_spread = small.std(axis=1).replace(0, 1)
    small_z = small.sub(row_average, axis=0).div(row_spread, axis=0)
    small_z.round(4).to_csv(RESULT_FOLDER + "/step8_heatmap.csv")

    patient_labels = pd.DataFrame()
    patient_labels["sample"] = gene_data.columns[chosen_patients]
    patient_labels["group"] = group_names[chosen_patients]
    patient_labels.to_csv(RESULT_FOLDER + "/step8_heatmap_patients.csv",
                          index=False)

    print("Heatmap data saved:", len(heatmap_genes), "genes x",
          len(chosen_patients), "patients")

    return results, size_cutoff


# =============================================================
# STEP 9 - PICK CANDIDATE BIOMARKERS
# =============================================================

def best_group_excluding_own(gene, own_group, best_group, second_group):
    # Which disease is the RUNNER UP for this gene?
    #
    # If this gene's strongest disease is the one we are looking at,
    # the runner up is the second strongest. Otherwise the strongest
    # disease is itself the runner up (some other disease beat us).

    if best_group[gene] == own_group:
        return second_group[gene]

    return best_group[gene]


def find_biomarkers(results, gene_data, size_cutoff):
    # A biomarker is a gene we could use as a TEST for one disease.
    #
    # A good biomarker needs four things:
    #   1. significant  - the difference is real (low FDR)
    #   2. big          - the difference is large enough to measure
    #   3. specific     - it is high in ONE disease only. A gene that
    #                     is high in six diseases cannot tell them apart.
    #   4. expressed    - the gene is active enough that a real lab
    #                     test would see a signal
    #
    # Specificity = my difference / the next best disease's difference.
    # A big ratio means "much higher here than anywhere else".

    print("Looking for candidate biomarkers...")

    # Keep only significant and big
    is_significant = results["fdr"] < FDR_CUTOFF
    is_big = results["difference"] >= size_cutoff
    good = results[is_significant & is_big].copy()

    print(len(good), "gene-disease pairs are significant and big")

    if len(good) == 0:
        print("Nothing passed. Try a lower cut-off.")
        return None

    # ---- Work out the best and second best disease for each gene ----
    best_difference = {}
    best_group = {}
    second_difference = {}
    second_group = {}

    for row in results.itertuples(index=False):
        gene = row.gene
        value = row.difference

        if gene not in best_difference:
            best_difference[gene] = value
            best_group[gene] = row.group
            second_difference[gene] = -999.0
            second_group[gene] = "none"
        elif value > best_difference[gene]:
            # the old best becomes the new second best
            second_difference[gene] = best_difference[gene]
            second_group[gene] = best_group[gene]
            best_difference[gene] = value
            best_group[gene] = row.group
        elif value > second_difference[gene]:
            second_difference[gene] = value
            second_group[gene] = row.group

    # ---- Score each candidate ----
    gene_average = gene_data.mean(axis=1)
    expression_cutoff = float(np.median(gene_data.values))

    specificity_list = []
    is_own_best_list = []
    expression_list = []
    was_floored_list = []
    next_best_group_list = []

    for row in good.itertuples(index=False):
        gene = row.gene

        # Is this gene's strongest disease the one we are looking at?
        if best_group[gene] == row.group:
            is_own_best_list.append(True)
        else:
            is_own_best_list.append(False)

        # The next best difference. We put a floor of 0.1 on it
        # because dividing by a tiny number gives a silly huge ratio.
        next_best = abs(second_difference[gene])

        if next_best < SPECIFICITY_FLOOR:
            next_best = SPECIFICITY_FLOOR
            was_floored_list.append(True)
        else:
            was_floored_list.append(False)

        specificity_list.append(row.difference / next_best)
        expression_list.append(gene_average[gene])
        next_best_group_list.append(best_group_excluding_own(
            gene, row.group, best_group, second_group))

    good["specificity"] = specificity_list
    good["is_own_best"] = is_own_best_list
    good["mean_expression"] = expression_list
    good["specificity_floored"] = was_floored_list
    good["next_best_disease"] = next_best_group_list

    # ---- HONESTY CHECK on the specificity floor ----
    #
    # Specificity = my difference / the next best disease's difference.
    # If the next best is tiny we would divide by almost nothing and
    # get a silly huge number, so we floor the bottom at 0.1.
    #
    # BUT on this dataset the whole significant range is only about
    # 0.2 to 0.6. So 0.1 is NOT small here, and for many genes the
    # floor IS the denominator. When that happens, "specificity" is
    # really just the difference divided by a constant - it stops
    # measuring specificity at all.
    #
    # So we count how often that happened and warn if it is common.
    how_many_floored = sum(was_floored_list)
    percent_floored = 100.0 * how_many_floored / len(good)

    print("Specificity floor used for", how_many_floored, "of", len(good),
          "pairs (", round(percent_floored), "% )")

    if percent_floored >= 25.0:
        print("WARNING: the floor is doing the work for many genes, so")
        print("'specificity' is not really measuring specificity here.")
        print("Try SPECIFICITY_FLOOR = 0.02, or ignore the specificity")
        print("part of the score. Say which you did in your write-up.")

    # ---- Apply the last two filters ----
    keep_rows = []
    for row in good.itertuples(index=False):
        if row.is_own_best == False:
            keep_rows.append(False)
        elif row.specificity < SPECIFICITY_CUTOFF:
            keep_rows.append(False)
        elif row.mean_expression <= expression_cutoff:
            keep_rows.append(False)
        else:
            keep_rows.append(True)

    candidates = good[keep_rows].copy()
    print(len(candidates), "candidates are specific enough and active enough")

    if len(candidates) == 0:
        print("No candidates left. Try SPECIFICITY_CUTOFF = 1.2")
        return None

    # ---- Give each candidate a score out of 1 ----
    # We combine three things. Each is first squashed to 0-1 so they
    # can be added together fairly.
    significance_part = scale_to_0_1(-np.log10(candidates["p_value"] + 1e-300))
    size_part = scale_to_0_1(candidates["difference"])
    specificity_part = scale_to_0_1(np.clip(candidates["specificity"], 0, 10))

    candidates["score"] = (0.4 * significance_part +
                           0.3 * size_part +
                           0.3 * specificity_part)

    # ---- Keep only the best few PER DISEASE ----
    #
    # WHY THIS MATTERS A LOT. CLL has 448 patients, so its genes get
    # the strongest statistics and would fill the whole top of the
    # list. A biomarker panel made only of CLL genes can detect one
    # disease out of seventeen - useless.
    #
    # So we take the best CANDIDATES_PER_DISEASE from each disease
    # first, and only then sort them all together. Every disease gets
    # a fair chance to contribute.

    fair_list = []
    all_diseases = candidates["group"].unique()

    for one_disease in all_diseases:
        just_this_one = candidates[candidates["group"] == one_disease]
        just_this_one = just_this_one.sort_values("score", ascending=False)
        fair_list.append(just_this_one.head(CANDIDATES_PER_DISEASE))

    candidates = pd.concat(fair_list, ignore_index=True)
    candidates = candidates.sort_values("score", ascending=False)

    print("After keeping the best", CANDIDATES_PER_DISEASE, "per disease:",
          len(candidates), "candidates from", len(all_diseases), "diseases")

    candidates.to_csv(RESULT_FOLDER + "/step9_candidates.csv", index=False)

    # ---- One summary row per disease, for the report ----
    summary_disease = []
    summary_count = []
    summary_best = []

    for one_disease in candidates["group"].unique():
        just_this_one = candidates[candidates["group"] == one_disease]
        summary_disease.append(one_disease)
        summary_count.append(len(just_this_one))
        summary_best.append(just_this_one["gene"].iloc[0])

    summary = pd.DataFrame()
    summary["disease"] = summary_disease
    summary["how_many_candidates"] = summary_count
    summary["best_gene"] = summary_best
    summary.to_csv(RESULT_FOLDER + "/step9_summary.csv", index=False)

    print("Candidates found for each disease:")
    for row in summary.itertuples(index=False):
        print("  ", row.disease[0:35], "->", row.how_many_candidates,
              "( best:", row.best_gene, ")")

    print("Top 5 candidate biomarkers overall:")
    top_five = candidates.head(5)
    for row in top_five.itertuples(index=False):
        print("  ", row.gene, "for", row.group[0:25],
              " score =", round(row.score, 3))

    print("NOTE: these are ideas to test, not proven findings yet.")
    return candidates


def scale_to_0_1(numbers):
    # Squash any list of numbers into the range 0 to 1, so that
    # different measurements can be compared and added together.

    numbers = np.asarray(numbers, dtype="float64")

    smallest = numbers.min()
    biggest = numbers.max()
    spread = biggest - smallest

    if spread == 0:
        # Everything is the same value, so give everything 1
        return np.ones(len(numbers))

    scaled = (numbers - smallest) / spread
    return scaled


# =============================================================
# STEP 10 - MACHINE LEARNING AND SHAP
# =============================================================

class PickGenesForEveryDisease(SelectorMixin, BaseEstimator):
    # A gene picker that gives EVERY disease a fair say.
    #
    # THE PROBLEM WITH THE NORMAL PICKER
    # The usual SelectKBest asks each gene one question: "how well do
    # you separate all 17 diseases at once?" A gene that perfectly
    # marks T-ALL (174 patients) but says nothing about the other 16
    # gets a mediocre score, because most of what it is being judged
    # on is variation it knows nothing about.
    #
    # The original project found this threw away BCL11B, GATA3,
    # PRKCQ, TOX and CTSW - all textbook T-cell genes, exactly the
    # T-ALL markers we were hunting for.
    #
    # THE FIX
    # Ask the question once per disease instead: "how well do you
    # separate THIS disease from all the others?" Keep each disease's
    # top genes, then use everything that any disease voted for.
    #
    # This is written as a class because scikit-learn Pipelines need
    # objects with fit() and a support mask. It is the one class in
    # this program, and it exists only so the gene picking stays
    # inside the pipeline and cannot leak the test patients.

    def __init__(self, genes_per_disease=150):
        self.genes_per_disease = genes_per_disease

    def fit(self, features, labels):
        features = np.asarray(features)
        labels = np.asarray(labels)

        self.n_features_in_ = features.shape[1]

        # Start with "keep nothing", then let each disease vote genes in
        keep_mask = np.zeros(features.shape[1], dtype=bool)

        how_many = int(min(self.genes_per_disease, features.shape[1]))
        all_diseases = np.unique(labels)

        for one_disease in all_diseases:
            # Turn the 17-way question into a yes/no question:
            # "is this patient in THIS disease, yes or no?"
            is_this_disease = (labels == one_disease).astype(int)

            scores, ignore_p = f_classif(features, is_this_disease)
            scores = np.nan_to_num(scores, nan=0.0, posinf=0.0, neginf=0.0)

            # Find the best genes for this one disease
            best_positions = np.argsort(scores)[::-1][0:how_many]

            for one_position in best_positions:
                keep_mask[one_position] = True

        self.support_mask_ = keep_mask
        return self

    def _get_support_mask(self):
        return self.support_mask_


def build_one_model(model_name):
    # Build one model as a PIPELINE.
    #
    # A Pipeline is a to-do list done in order. Here the list is:
    #   1. pick_genes  - keep only the most useful genes
    #   2. scale       - (logistic regression only) put genes on one scale
    #   3. classifier  - actually learn the diagnosis
    #
    # WHY THE PIPELINE MATTERS. If we picked the useful genes by
    # looking at ALL the patients first, that choice would already
    # have been shaped by the hidden test patients, and our accuracy
    # score would be a lie. This is called DATA LEAKAGE. Inside a
    # pipeline, gene picking only ever sees the training patients.

    # Which gene picker? See PickGenesForEveryDisease above for why
    # "per_disease" usually works better on this data.
    if HOW_TO_PICK_GENES == "per_disease":
        gene_picker = PickGenesForEveryDisease(
            genes_per_disease=GENES_PER_DISEASE)
    else:
        gene_picker = SelectKBest(f_classif, k=GENES_FOR_MODEL)

    if model_name == "random_forest":
        # A random forest is 400 decision trees that vote.
        # class_weight="balanced_subsample" stops it ignoring the
        # small diseases just because they have fewer patients.
        model = Pipeline([
            ("pick_genes", gene_picker),
            ("classifier", RandomForestClassifier(n_estimators=400,
                                                  min_samples_leaf=2,
                                                  class_weight="balanced_subsample",
                                                  n_jobs=-1,
                                                  random_state=SEED))
        ])
        return model

    # Logistic regression is one big equation. It needs every gene on
    # the same scale first, or the loud genes dominate the equation.
    model = Pipeline([
        ("pick_genes", gene_picker),
        ("scale", StandardScaler()),
        ("classifier", LogisticRegression(max_iter=2000,
                                          class_weight="balanced",
                                          random_state=SEED))
    ])
    return model


def tidy_shap_shape(values, n_patients, n_genes, n_classes):
    # THIS FUNCTION EXISTS BECAUSE OF A REAL TRAP.
    #
    # SHAP gives back its answer in a different SHAPE depending on
    # which version of the shap library you have installed:
    #
    #   old versions  -> a LIST of grids, one grid per diagnosis
    #   new versions  -> one 3D block: patients x genes x diagnoses
    #
    # If we guess the shape wrongly, the numbers still look perfectly
    # reasonable - they are just attached to the WRONG GENES. That is
    # much worse than a crash, because you would never notice.
    #
    # So instead of guessing, we look at the actual sizes and match
    # them against what we know: how many patients, genes and
    # diagnoses we have. If nothing matches, we stop.

    # Old style: a list of grids. Stack them into one block.
    if isinstance(values, list):
        values = np.stack(values, axis=-1)

    values = np.asarray(values)

    # Only two diagnoses: SHAP gives a flat grid. Add a third side.
    if values.ndim == 2:
        return values[:, :, None]

    if values.ndim == 3:
        if values.shape == (n_patients, n_genes, n_classes):
            return values                                   # already right

        if values.shape == (n_classes, n_patients, n_genes):
            # np.transpose reorders the sides. (1,2,0) means
            # "new side 0 is old side 1, new side 1 is old side 2..."
            return np.transpose(values, (1, 2, 0))

        if values.shape == (n_patients, n_classes, n_genes):
            return np.transpose(values, (0, 2, 1))

    print("ERROR: I do not recognise the shape SHAP returned.")
    print("  I got:", values.shape)
    print("  I expected some order of:", (n_patients, n_genes, n_classes))
    print("Stopping, because guessing would attach numbers to the")
    print("wrong genes without telling you.")
    raise SystemExit(1)


def get_shap_importance(model_name, model, test_features, gene_names):
    # WHAT IS SHAP?
    #
    # The model can predict a diagnosis, but it cannot tell us why.
    # SHAP opens the box. For every patient and every gene it works
    # out: "how much did this gene push the answer towards this
    # diagnosis?" A big push means the gene mattered.
    #
    # This is what makes the project EXPLAINABLE. A model we cannot
    # explain is useless for biology, because we could not say which
    # genes to study next.
    #
    # We only explain SHAP_PATIENTS patients, not all of them,
    # because SHAP is slow: its work grows with
    # patients x genes x diagnoses.

    import shap

    # Which genes survived the pipeline gene picking?
    chosen_mask = model.named_steps["pick_genes"].get_support()
    chosen_genes = gene_names[chosen_mask]

    classifier = model.named_steps["classifier"]
    class_names = classifier.classes_

    how_many = min(SHAP_PATIENTS, len(test_features))
    small_test = test_features[0:how_many]

    # model[:-1] means "all the pipeline steps EXCEPT the last one".
    # So this runs the gene picking (and scaling) but stops before
    # the classifier, giving us exactly the numbers the classifier
    # itself sees.
    prepared = model[:-1].transform(small_test)

    print("  explaining", how_many, "patients,", len(chosen_genes),
          "genes,", len(class_names), "diagnoses")

    if model_name == "random_forest":
        # TreeExplainer is the fast exact method for tree models
        explainer = shap.TreeExplainer(classifier)
        raw_values = explainer.shap_values(prepared)
    else:
        # LinearExplainer is for equation-style models. It needs a
        # "background" of typical patients to compare against.
        background = prepared[0:min(100, len(prepared))]
        explainer = shap.LinearExplainer(classifier, background)
        raw_values = explainer.shap_values(prepared)

    values = tidy_shap_shape(raw_values, prepared.shape[0],
                             len(chosen_genes), len(class_names))

    # We only care HOW MUCH a gene pushed, not which way, so we take
    # the size of each push and ignore the plus or minus sign.
    push_size = np.abs(values)

    # Average over patients -> a grid of genes x diagnoses
    average_per_class = push_size.mean(axis=0)

    # Average over diagnoses too -> one score per gene
    overall = average_per_class.mean(axis=1)

    global_table = pd.DataFrame()
    global_table["gene"] = chosen_genes
    global_table["importance"] = overall
    global_table["method"] = "shap"

    # ---- Also keep the per-diagnosis scores ----
    #
    # WHY BOTHER? A gene that is brilliant for ONE rare leukaemia
    # only helps one of the seventeen predictions. Averaged over all
    # seventeen it looks unimportant. Keeping the per-diagnosis
    # scores stops us throwing away exactly the specific markers we
    # are hunting for.

    all_rows = []

    for class_number in range(len(class_names)):
        one_class = pd.DataFrame()
        one_class["gene"] = chosen_genes
        one_class["diagnosis"] = class_names[class_number]
        one_class["shap_score"] = average_per_class[:, class_number]

        # Rank within this diagnosis: 1 = most important
        one_class = one_class.sort_values("shap_score", ascending=False)
        one_class["rank_for_this_diagnosis"] = range(1, len(one_class) + 1)

        all_rows.append(one_class)

    per_class_table = pd.concat(all_rows, ignore_index=True)

    return global_table, per_class_table


def get_simple_importance(model_name, model, gene_names):
    # A fallback for when SHAP is switched off or fails.
    #
    # Every model can already report a rough "how much did I use this
    # gene" score. It is quicker but cruder than SHAP: it says a gene
    # was used, not how it changed each prediction.

    chosen_mask = model.named_steps["pick_genes"].get_support()
    chosen_genes = gene_names[chosen_mask]
    classifier = model.named_steps["classifier"]

    table = pd.DataFrame()
    table["gene"] = chosen_genes

    if model_name == "random_forest":
        table["importance"] = classifier.feature_importances_
        table["method"] = "built_in_tree_score"
    else:
        # For the equation model, a bigger number in the equation
        # means the gene mattered more. We take the size and ignore
        # the sign, then average over the diagnoses.
        weights = np.abs(classifier.coef_)
        table["importance"] = weights.mean(axis=0)
        table["method"] = "built_in_equation_weight"

    return table


def save_one_model(model_name, model):
    # Save the trained model to a file so it can be reused later
    # without retraining - for example to test it on patients from a
    # completely different hospital.

    try:
        import joblib
        file_path = MODEL_FOLDER + "/" + model_name + ".joblib"
        joblib.dump(model, file_path)
        print("  saved the model to", file_path)
    except Exception as problem:
        print("  could not save the model:", problem)


def check_candidates_against_shap(candidates, per_class_table):
    # AN HONEST COMPARISON WORTH REPORTING.
    #
    # Step 9 found biomarkers using statistics, one gene at a time.
    # SHAP found the genes the model actually leaned on. Do they
    # agree? Often much less than people expect.
    #
    # The reason is REDUNDANCY. If ten genes all switch on together
    # in T-ALL, the model only needs one or two of them. The other
    # eight get a low SHAP score even though they are perfectly good
    # markers biologically. Statistics tests each gene alone, so it
    # keeps all ten.
    #
    # This is not a contradiction. Both methods found the same
    # biology and picked different representatives from it. Say so
    # in your write-up rather than tuning until it disappears.

    if candidates is None or per_class_table is None:
        return

    print("")
    print("Do the Step 9 biomarkers match what the model relied on?")

    how_many_checked = 0
    how_many_in_top_50 = 0

    for row in candidates.drop_duplicates("gene").itertuples(index=False):
        # Find this gene's rank for its own disease
        matching = per_class_table[
            (per_class_table["gene"] == row.gene) &
            (per_class_table["diagnosis"] == row.group)]

        if len(matching) == 0:
            continue        # the gene never survived gene picking

        how_many_checked = how_many_checked + 1
        its_rank = matching["rank_for_this_diagnosis"].iloc[0]

        if its_rank <= 50:
            how_many_in_top_50 = how_many_in_top_50 + 1

    print("  ", how_many_checked, "candidates reached the model")
    print("  ", how_many_in_top_50, "of those are in the model's top 50")
    print("  The rest are probably redundant, not wrong.")


def report_which_candidates_survived(model, gene_names):
    # Say plainly which Step 9 biomarkers the gene picking threw away.
    #
    # This is about being fair when we compare the two methods later.
    # If a gene never made it into the model, "the model did not use
    # it" tells us nothing about the gene - it only tells us the gene
    # picker dropped it.

    file_path = RESULT_FOLDER + "/step9_candidates.csv"

    if os.path.exists(file_path) == False:
        return

    candidates = pd.read_csv(file_path)

    kept_mask = model.named_steps["pick_genes"].get_support()
    kept_genes = []
    for one_gene in gene_names[kept_mask]:
        kept_genes.append(str(one_gene))

    survived = []
    dropped = []

    for one_gene in candidates["gene"].unique():
        if str(one_gene) in kept_genes:
            survived.append(one_gene)
        else:
            dropped.append(one_gene)

    print("  Step 9 candidates:", len(survived), "survived gene picking,",
          len(dropped), "dropped")

    if len(dropped) > 0:
        shown = dropped[0:10]
        text = ""
        for one_gene in shown:
            text = text + str(one_gene) + " "
        print("    dropped:", text)
        print("    these cannot appear in any importance list, so do not")
        print("    count them as 'the model disagreed'")


def choose_the_model(scores):
    # DECIDE which model's gene ranking feeds Step 11, and say why.
    #
    # This used to be a name typed at the top of the file. That is an
    # assumption dressed up as a setting: the number that would
    # justify it (how well each model scored) is not known until
    # after both models have run.
    #
    # So now we look at the scores first and print the reasoning.

    print("")

    # ---- Case 1: the user named a model, so just use it ----
    if MAIN_MODEL != "best":
        if MAIN_MODEL not in list(scores["model"]):
            print("ERROR: MAIN_MODEL is set to '" + str(MAIN_MODEL) + "'")
            print("but that is not one of the models we trained.")
            print("Use one of:", MODEL_NAMES, 'or "best"')
            raise SystemExit(1)

        print("Step 11 will use:", MAIN_MODEL, "(you chose it by name)")

        # Still say whether the numbers agree with that choice, so a
        # forced setting cannot quietly hide a worse model.
        sorted_scores = scores.sort_values(BEST_MODEL_MEASURE,
                                           ascending=False)
        would_have_picked = sorted_scores["model"].iloc[0]

        if would_have_picked != MAIN_MODEL:
            print("NOTE: on", BEST_MODEL_MEASURE, "the better model is")
            print("  actually", would_have_picked + ".",
                  "You are overriding that on purpose.")
            print("  Say why in your write-up, or set MAIN_MODEL = \"best\".")

        chosen = MAIN_MODEL

    else:
        # ---- Case 2: let the scores decide ----
        if BEST_MODEL_MEASURE not in scores.columns:
            print("ERROR: BEST_MODEL_MEASURE is '" + BEST_MODEL_MEASURE + "'")
            print("but that column does not exist. Choose one of:")
            for one_column in scores.columns:
                print("  ", one_column)
            raise SystemExit(1)

        sorted_scores = scores.sort_values(BEST_MODEL_MEASURE,
                                           ascending=False)
        chosen = sorted_scores["model"].iloc[0]
        winning_score = float(sorted_scores[BEST_MODEL_MEASURE].iloc[0])

        print("Step 11 will use:", chosen)
        print("Reason: best", BEST_MODEL_MEASURE, "=",
              round(winning_score, 4))

        # ---- Is the win actually meaningful? ----
        #
        # Cross-validation gives us a wobble as well as an average. If
        # the two models are closer together than the wobble, the
        # "winner" is noise and we should not pretend otherwise.
        if len(sorted_scores) > 1:
            runner_up = sorted_scores["model"].iloc[1]
            runner_up_score = float(sorted_scores[BEST_MODEL_MEASURE].iloc[1])
            gap = winning_score - runner_up_score

            print("Runner up:", runner_up, "=", round(runner_up_score, 4),
                  "( gap", round(gap, 4), ")")

            if "cv_wobble" in scores.columns:
                biggest_wobble = float(scores["cv_wobble"].max())

                if gap == 0.0:
                    # An exact tie. Whichever model we take is an
                    # arbitrary pick, and pretending otherwise would
                    # be inventing a result.
                    print("")
                    print("EXACT TIE on this measure. The pick is arbitrary.")
                    print("Try BEST_MODEL_MEASURE = \"test_auc\", or just")
                    print("report both models and say they tied.")

                elif gap <= biggest_wobble:
                    # Note: <= not <. If the gap only equals the
                    # wobble it is still not a real difference.
                    print("")
                    print("CAREFUL: the gap (" + str(round(gap, 4)) + ") is")
                    print("no bigger than the cross-validation wobble ("
                          + str(round(biggest_wobble, 4)) + ").")
                    print("So the two models are not really distinguishable")
                    print("on this measure. Report both, and do not claim one")
                    print("model is better than the other.")
                else:
                    print("The gap is bigger than the wobble ("
                          + str(round(biggest_wobble, 4))
                          + "), so this is a real difference.")

    # ---- Write the decision down ----
    # A choice that only exists in the terminal scrollback is a choice
    # nobody can check later.
    decision = pd.DataFrame()
    decision["item"] = ["model_used_for_step_11", "how_it_was_chosen",
                        "measure_used"]
    decision["value"] = [chosen,
                         "by score" if MAIN_MODEL == "best" else "set by hand",
                         BEST_MODEL_MEASURE]
    decision.to_csv(RESULT_FOLDER + "/step10_chosen_model.csv", index=False)

    return chosen


def run_machine_learning(gene_data, groups):
    # Now we ask a computer to LEARN the diagnosis from the genes.
    # If it can, then the gene patterns really do carry the
    # information a doctor uses.
    #
    # Vocabulary:
    #   feature  = one piece of information the model uses (one gene)
    #   label    = the answer we want it to predict (the diagnosis)
    #   training = the patients the model learns from
    #   test     = patients kept hidden, used to check it honestly

    print("Training the machine learning models...")

    can_analyse = groups["can_analyse"].values
    all_labels = groups["group"].astype(str).values

    # Models want patients as rows, so we flip the table
    all_features = gene_data.values.T

    features = all_features[can_analyse]
    labels = all_labels[can_analyse]
    gene_names = gene_data.index.values

    print(features.shape[0], "patients,", features.shape[1], "genes,",
          len(set(labels)), "diagnoses")

    # ---- Split into training and test ----
    # stratify keeps the same mix of diseases in both halves.
    train_features, test_features, train_labels, test_labels = train_test_split(
        features, labels, test_size=0.3, stratify=labels, random_state=SEED)

    print("Training patients:", len(train_labels))
    print("Test patients (kept hidden):", len(test_labels))

    score_rows = []
    all_importance = []
    all_per_class = []

    # Keep each model's answers, so we can choose AFTER we have seen
    # how they all scored. Choosing first and scoring later would be
    # picking the winner before the race.
    importance_by_model = {}
    per_class_by_model = {}

    # ---- Train each model in turn ----
    for model_name in MODEL_NAMES:
        print("")
        print("--- model:", model_name, "---")

        model = build_one_model(model_name)

        # ---- CROSS-VALIDATION on the training patients ----
        #
        # One train/test split gives one number, and that number
        # partly depends on which patients happened to land in the
        # test half. Cross-validation gives a range instead.
        #
        # It splits the TRAINING patients into 5 parts, then trains 5
        # times, each time holding out a different part. Five scores
        # tell us the average AND how much it wobbles.
        #
        # The hidden test patients are not involved at all here. They
        # stay untouched until the very end.

        print("  cross-validating on the training patients...")

        splitter = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True,
                                   random_state=SEED)

        cv_scores = cross_val_score(model, train_features, train_labels,
                                    cv=splitter,
                                    scoring="balanced_accuracy")

        cv_average = float(cv_scores.mean())
        cv_wobble = float(cv_scores.std())

        print("  cross-validation score:", round(cv_average, 3),
              "plus or minus", round(cv_wobble, 3))

        print("  learning on all the training patients...")
        model.fit(train_features, train_labels)

        predictions = model.predict(test_features)
        probabilities = model.predict_proba(test_features)

        # Balanced accuracy treats small diseases as importantly as
        # big ones. Plain accuracy would look great just by always
        # guessing the most common disease.
        accuracy = balanced_accuracy_score(test_labels, predictions)

        # AUC asks how well the model separates the groups.
        # 1.0 is perfect, 0.5 is guessing.
        #
        # sklearn wants a different shape depending on how many
        # groups there are, so we check first.
        number_of_diagnoses = len(model.classes_)

        if number_of_diagnoses == 2:
            second_column = probabilities[:, 1]
            auc = roc_auc_score(test_labels, second_column)
        else:
            auc = roc_auc_score(test_labels, probabilities,
                                multi_class="ovr", average="macro",
                                labels=model.classes_)

        print("  balanced accuracy:", round(accuracy, 3))
        print("  AUC:", round(auc, 3))

        score_rows.append({"model": model_name,
                           "cv_balanced_accuracy": round(cv_average, 4),
                           "cv_wobble": round(cv_wobble, 4),
                           "test_balanced_accuracy": round(accuracy, 4),
                           "test_auc": round(auc, 4),
                           "n_train": len(train_labels),
                           "n_test": len(test_labels),
                           "n_diagnoses": number_of_diagnoses,
                           "gene_picker": HOW_TO_PICK_GENES,
                           "label_source": "diagnosis"})

        # ---- Save the confusion matrix ----
        # This grid shows which diseases got mixed up with which.
        # Rows = the true answer, columns = what the model guessed.
        grid = confusion_matrix(test_labels, predictions,
                                labels=model.classes_)
        confusion = pd.DataFrame(grid,
                                 index=model.classes_,
                                 columns=model.classes_)
        confusion.to_csv(RESULT_FOLDER + "/step10_confusion_" +
                         model_name + ".csv")

        save_one_model(model_name, model)

        # ---- Which Step 9 candidates even reached the model? ----
        # A gene thrown out by gene picking cannot possibly appear in
        # any importance ranking. Counting it as "not confirmed" would
        # be unfair, so we say plainly which ones never got a chance.
        report_which_candidates_survived(model, gene_names)

        # ---- Work out which genes mattered ----
        per_class_table = None
        importance = None

        if USE_SHAP:
            print("  running SHAP...")
            try:
                importance, per_class_table = get_shap_importance(
                    model_name, model, test_features, gene_names)
                print("  SHAP finished")
            except Exception as problem:
                # Print the REAL problem, not a vague guess, so you
                # can actually fix it.
                print("  SHAP failed:", type(problem).__name__, "-", problem)
                print("  falling back to the built-in score")
                importance = None

        if importance is None:
            importance = get_simple_importance(model_name, model, gene_names)

        importance["model"] = model_name
        importance = importance.sort_values("importance", ascending=False)
        all_importance.append(importance)

        if per_class_table is not None:
            per_class_table["model"] = model_name
            all_per_class.append(per_class_table)

        importance_by_model[model_name] = importance
        per_class_by_model[model_name] = per_class_table

        print("  top 5 genes for this model:")
        for row in importance.head(5).itertuples(index=False):
            print("    ", row.gene, round(row.importance, 6))

    # ---- Save everything ----
    scores = pd.DataFrame(score_rows)
    scores.to_csv(RESULT_FOLDER + "/step10_model_scores.csv", index=False)

    everything = pd.concat(all_importance, ignore_index=True)
    everything.to_csv(RESULT_FOLDER + "/step10_gene_importance.csv",
                      index=False)

    if len(all_per_class) > 0:
        per_class_all = pd.concat(all_per_class, ignore_index=True)
        per_class_all.to_csv(RESULT_FOLDER + "/step10_shap_per_diagnosis.csv",
                             index=False)
        print("")
        print("Per-diagnosis SHAP saved:", len(per_class_all), "rows")

    print("")
    print("Comparing the two models:")
    print(scores.to_string(index=False))

    # ---- Now choose which model's gene ranking goes to Step 11 ----
    chosen_model = choose_the_model(scores)

    main_importance = importance_by_model[chosen_model]
    main_per_class = per_class_by_model[chosen_model]

    return main_importance, main_per_class



# =============================================================
# STEP 11 - COMBINE BOTH KINDS OF EVIDENCE
# =============================================================

def make_consensus(candidates, importance):
    # We now have two separate opinions about which genes matter:
    #
    #   Statistics (Step 9): this gene differs between diseases
    #   Machine learning (Step 10): the model needed this gene
    #
    # A gene supported by BOTH is much stronger evidence than a gene
    # supported by only one. We combine them with a harmonic mean,
    # which punishes genes that are strong on one side and weak on
    # the other. A normal average would hide that weakness.

    print("Combining statistics and machine learning...")

    if candidates is None:
        print("No candidates from step 9, skipping")
        return None

    # Ranks: 1 = most important
    importance = importance.sort_values("importance", ascending=False)
    importance_rank = {}
    rank_number = 1
    for row in importance.itertuples(index=False):
        importance_rank[row.gene] = rank_number
        rank_number = rank_number + 1

    # Keep one row per gene (its best one)
    candidates = candidates.drop_duplicates("gene")

    shared_genes = []
    for row in candidates.itertuples(index=False):
        if row.gene in importance_rank:
            shared_genes.append(row)

    print(len(shared_genes), "genes have BOTH kinds of support")

    if len(shared_genes) == 0:
        print("No overlap found. This happens and is worth reporting.")
        return None

    gene_list = []
    disease_list = []
    statistics_score_list = []
    model_rank_list = []

    for row in shared_genes:
        gene_list.append(row.gene)
        disease_list.append(row.group)
        statistics_score_list.append(row.score)
        model_rank_list.append(importance_rank[row.gene])

    consensus = pd.DataFrame()
    consensus["gene"] = gene_list
    consensus["disease"] = disease_list
    consensus["statistics_score"] = statistics_score_list
    consensus["model_rank"] = model_rank_list

    # Turn both into 0-1 scores where bigger is better.
    # For rank, small is good, so we flip it.
    statistics_part = scale_to_0_1(consensus["statistics_score"])
    model_part = scale_to_0_1(-consensus["model_rank"])

    # Harmonic mean = 2ab / (a + b)
    top = 2.0 * statistics_part * model_part
    bottom = statistics_part + model_part
    bottom = np.where(bottom == 0, 1e-9, bottom)   # never divide by zero

    consensus["consensus_score"] = top / bottom
    consensus = consensus.sort_values("consensus_score", ascending=False)
    consensus.to_csv(RESULT_FOLDER + "/step11_consensus.csv", index=False)

    print("Best genes supported by both methods:")
    top_ten = consensus.head(10)
    for row in top_ten.itertuples(index=False):
        print("  ", row.gene, "for", row.disease[0:25],
              " score =", round(row.consensus_score, 3))

    # Show the spread of diseases, not just the top of the list
    print("Diseases covered by the consensus list:")
    disease_counts = consensus["disease"].value_counts()
    for one_disease in disease_counts.index:
        print("  ", one_disease[0:35], "->", disease_counts[one_disease],
              "genes")

    return consensus


def test_each_biomarker(consensus, gene_data, groups):
    # DOES EACH BIOMARKER ACTUALLY WORK AS A TEST?
    #
    # Everything so far said "this gene looks different between
    # diseases". That is not the same as "I could use this gene as a
    # diagnostic test". So now we try it as a test and score it.
    #
    # AUC (area under the ROC curve) is the score:
    #   1.0  = perfect, never wrong
    #   0.5  = useless, same as flipping a coin
    #   0.8+ = genuinely useful
    #
    # THE CRITICAL PART: we score it ONLY on the hidden test
    # patients. Scoring a gene on the same patients we used to pick it
    # would be like marking your own exam - the number would look
    # great and mean nothing.

    print("Testing each biomarker on the HIDDEN patients only...")

    if consensus is None:
        print("No consensus genes, skipping")
        return None

    can_analyse = groups["can_analyse"].values
    all_labels = groups["group"].astype(str).values

    # Rebuild the SAME split Step 10 used. Same seed, same test size,
    # so the same patients are hidden. If these ever disagree, the
    # "hidden" patients would not really be hidden.
    all_positions = np.arange(gene_data.shape[1])
    usable_positions = all_positions[can_analyse]

    train_positions, test_positions = train_test_split(
        usable_positions, test_size=0.3,
        stratify=all_labels[can_analyse], random_state=SEED)

    print("Scoring on", len(test_positions), "hidden patients")

    auc_list = []
    curve_pieces = []

    # A fixed set of x-axis points, so all the curves can be drawn
    # on the same axes in R
    x_axis_points = np.linspace(0, 1, 200)

    for row in consensus.itertuples(index=False):
        one_gene = row.gene

        if one_gene not in gene_data.index:
            auc_list.append(np.nan)
            continue

        # This gene's expression in the hidden patients
        gene_values = gene_data.loc[one_gene].values[test_positions]

        # The right answer: is this patient in the gene's disease?
        true_answer = (all_labels[test_positions] == str(row.disease))
        true_answer = true_answer.astype(int)

        # We need both kinds of patient to score anything
        how_many_yes = int(true_answer.sum())
        if how_many_yes < 3 or how_many_yes == len(true_answer):
            auc_list.append(np.nan)
            continue

        score = roc_auc_score(true_answer, gene_values)
        auc_list.append(score)

        # Save the curve for the first few genes, for the picture
        if len(curve_pieces) < 12:
            false_rate, true_rate, ignore = roc_curve(true_answer, gene_values)

            one_curve = pd.DataFrame()
            one_curve["gene"] = one_gene
            one_curve["disease"] = str(row.disease)
            one_curve["auc"] = score
            one_curve["false_positive_rate"] = x_axis_points

            # Fill in the curve at our fixed x points
            one_curve["true_positive_rate"] = np.interp(x_axis_points,
                                                        false_rate, true_rate)
            curve_pieces.append(one_curve)

    consensus = consensus.copy()
    consensus["auc_on_hidden_patients"] = auc_list
    consensus.to_csv(RESULT_FOLDER + "/step11_consensus.csv", index=False)

    if len(curve_pieces) > 0:
        all_curves = pd.concat(curve_pieces, ignore_index=True)
        all_curves.to_csv(RESULT_FOLDER + "/step11_roc_curves.csv", index=False)

    middle_auc = float(np.nanmedian(np.array(auc_list, dtype="float64")))
    print("Middle (median) single-gene AUC:", round(middle_auc, 3))

    print("Best single-gene tests:")
    best = consensus.sort_values("auc_on_hidden_patients", ascending=False)
    for row in best.head(8).itertuples(index=False):
        print("  ", row.gene, "for", row.disease[0:28],
              " AUC =", round(row.auc_on_hidden_patients, 3))

    return consensus, train_positions, test_positions


def test_the_whole_panel(consensus, gene_data, groups,
                         train_positions, test_positions):
    # DOES A SMALL PANEL OF GENES WORK ON ITS OWN?
    #
    # One gene can only really answer one yes/no question. A real
    # diagnostic test uses a PANEL of several genes together.
    #
    # So we take the best PANEL_SIZE consensus genes, train a fresh
    # model using ONLY those genes, and score it on the hidden
    # patients. If a 25-gene panel gets close to the 14,781-gene
    # model, that is a genuinely useful result: a lab could actually
    # build a 25-gene test, but not a 14,781-gene one.

    print("Testing a small panel of genes on its own...")

    if consensus is None:
        return None

    all_labels = groups["group"].astype(str).values

    # Pick the panel genes, skipping any not in our table
    panel_genes = []
    for one_gene in consensus["gene"].head(PANEL_SIZE):
        if one_gene in gene_data.index:
            panel_genes.append(one_gene)

    if len(panel_genes) < 2:
        print("Not enough panel genes, skipping")
        return None

    print("Panel size:", len(panel_genes), "genes")

    panel_values = gene_data.loc[panel_genes].values.T

    panel_model = RandomForestClassifier(n_estimators=400,
                                         min_samples_leaf=2,
                                         class_weight="balanced_subsample",
                                         n_jobs=-1,
                                         random_state=SEED)

    panel_model.fit(panel_values[train_positions],
                    all_labels[train_positions])

    guesses = panel_model.predict(panel_values[test_positions])
    chances = panel_model.predict_proba(panel_values[test_positions])

    panel_accuracy = balanced_accuracy_score(all_labels[test_positions],
                                             guesses)

    if len(panel_model.classes_) == 2:
        panel_auc = roc_auc_score(all_labels[test_positions], chances[:, 1])
    else:
        panel_auc = roc_auc_score(all_labels[test_positions], chances,
                                  multi_class="ovr", average="macro",
                                  labels=panel_model.classes_)

    print("Panel balanced accuracy:", round(panel_accuracy, 3))
    print("Panel AUC:", round(panel_auc, 3))

    panel_table = pd.DataFrame()
    panel_table["item"] = ["panel_size", "balanced_accuracy", "auc", "genes"]
    panel_table["value"] = [len(panel_genes),
                            round(panel_accuracy, 4),
                            round(panel_auc, 4),
                            ";".join(panel_genes)]
    panel_table.to_csv(RESULT_FOLDER + "/step11_panel.csv", index=False)

    # ---- HOW STRONG IS OUR EVIDENCE, REALLY? ----
    #
    # There are three levels of proof, from weakest to strongest:
    #
    #   Level 1  single genes, hidden patients, SAME hospital cohort
    #   Level 2  a gene panel, hidden patients, SAME cohort
    #   Level 3  a completely DIFFERENT cohort of patients
    #
    # Only level 3 is validation in the proper sense. Levels 1 and 2
    # still share all the same equipment, labs and processing quirks
    # as the data we trained on. Say this plainly in the write-up
    # rather than letting a good AUC speak for itself.

    middle_auc = float(np.nanmedian(consensus["auc_on_hidden_patients"]))

    levels = pd.DataFrame()
    levels["level"] = [1, 2, 3]
    levels["what_was_tested"] = ["single genes, hidden patients",
                                 "gene panel, hidden patients",
                                 "a different cohort entirely"]
    levels["score"] = [round(middle_auc, 4), round(panel_auc, 4), np.nan]
    levels["is_real_validation"] = [False, False, True]
    levels.to_csv(RESULT_FOLDER + "/step11_evidence_levels.csv", index=False)

    print("")
    print("Level 3 (a different cohort) is EMPTY. Until you run")
    print("simple_external_validation.py, nothing here is validated")
    print("outside this one dataset. That is a limitation to report.")

    return panel_table



# =============================================================
# STEP 12 - GENE NETWORK (WHICH GENES WORK TOGETHER?)
# =============================================================

def build_gene_network(gene_data, groups):
    # Genes do not work alone. They switch on and off in teams.
    #
    # If two genes always go up and down together across patients,
    # they are probably part of the same biological process. We draw
    # that as a NETWORK: genes are dots, and we draw a line between
    # two dots when they move together strongly.
    #
    # "Move together" is measured by CORRELATION:
    #   +1 = always up together
    #    0 = unrelated
    #   -1 = one goes up when the other goes down
    # We draw a line when the strength is at least NETWORK_STRENGTH.
    #
    # ==========================================================
    # READ THIS BEFORE BELIEVING THE PICTURE
    # ==========================================================
    # Our patients have 17 different diseases. Any two genes that are
    # both high in CLL will look strongly correlated across the whole
    # cohort - even if they have nothing to do with each other
    # biologically. They are not partners; they just both happen to
    # mark the same disease.
    #
    # So we build TWO networks:
    #
    #   "everyone"   correlation across all patients. Full of huge
    #                blobs that just re-describe the diseases.
    #   "within"     correlation worked out separately INSIDE each
    #                cluster, then averaged. This cancels out the
    #                between-disease effect and leaves gene teams that
    #                hold together regardless of disease.
    #
    # The "within" network is the one to interpret. Report both and
    # say which you used.

    print("Building the gene network...")

    # ---- Use only the most variable genes ----
    # A network of 14,781 genes has 100 million possible lines. It
    # would not fit in memory and could not be drawn.
    values = gene_data.values.astype("float64")

    middle = np.median(values, axis=1, keepdims=True)
    mad = np.median(np.abs(values - middle), axis=1)

    order_big_to_small = np.argsort(mad)[::-1]
    how_many = min(NETWORK_GENES, gene_data.shape[0])
    chosen = np.sort(order_big_to_small[0:how_many])

    small = gene_data.iloc[chosen]
    gene_names = small.index.values
    small_values = small.values.astype("float64")

    print("Network over the", len(gene_names), "most variable genes")

    # ---- Network 1: correlation across everyone ----
    everyone_grid = np.corrcoef(small_values)
    everyone_grid = np.nan_to_num(everyone_grid, nan=0.0)

    everyone_lines = find_strong_pairs(everyone_grid, gene_names,
                                       NETWORK_STRENGTH)
    everyone_lines.to_csv(RESULT_FOLDER + "/step12_network_everyone.csv",
                          index=False)
    print("Across everyone:", len(everyone_lines), "strong gene pairs")

    # ---- Network 2: correlation inside each cluster, then averaged ----
    cluster_labels = groups["tree_cluster"].values

    total_grid = np.zeros((len(gene_names), len(gene_names)))
    total_weight = 0.0

    for one_cluster in np.unique(cluster_labels):
        in_this_cluster = (cluster_labels == one_cluster)
        how_many_patients = int(in_this_cluster.sum())

        # Too few patients gives a meaningless correlation
        if how_many_patients < 10:
            continue

        just_this_cluster = small_values[:, in_this_cluster]
        one_grid = np.corrcoef(just_this_cluster)
        one_grid = np.nan_to_num(one_grid, nan=0.0)

        # Bigger clusters count for more
        total_grid = total_grid + one_grid * how_many_patients
        total_weight = total_weight + how_many_patients

    if total_weight > 0:
        within_grid = total_grid / total_weight
    else:
        within_grid = everyone_grid

    within_lines = find_strong_pairs(within_grid, gene_names,
                                     NETWORK_STRENGTH)
    within_lines.to_csv(RESULT_FOLDER + "/step12_network_within.csv",
                        index=False)

    print("Within clusters:", len(within_lines), "strong gene pairs")

    if len(everyone_lines) > 0:
        shrunk = 100.0 * len(within_lines) / len(everyone_lines)
        print("That is", round(shrunk), "% of the first number.")
        print("The drop is the fake disease-driven correlation going away.")

    # ---- Find the gene teams (modules) ----
    # We cluster the GENES this time, not the patients. Genes that
    # correlate strongly end up in the same team.
    distance = 1.0 - np.abs(everyone_grid)
    np.fill_diagonal(distance, 0.0)

    gene_tree = linkage(squareform(distance, checks=False), method="average")
    team_numbers = fcluster(gene_tree, t=NETWORK_TEAMS, criterion="maxclust")

    # How many lines does each gene have? A gene with many is a "hub".
    everyone_count = count_lines_per_gene(everyone_lines, gene_names)
    within_count = count_lines_per_gene(within_lines, gene_names)

    dots = pd.DataFrame()
    dots["gene"] = gene_names
    dots["team"] = team_numbers
    dots["lines_everyone"] = everyone_count
    dots["lines_within"] = within_count
    dots["variability"] = mad[chosen]
    dots["average_expression"] = small_values.mean(axis=1)
    dots = dots.sort_values("lines_everyone", ascending=False)
    dots.to_csv(RESULT_FOLDER + "/step12_network_genes.csv", index=False)

    print("Gene teams found:")
    for one_team in sorted(set(team_numbers)):
        just_this_team = dots[dots["team"] == one_team]
        hub_gene = just_this_team["gene"].iloc[0]
        print("   team", one_team, "-", len(just_this_team), "genes, hub =",
              hub_gene)

    biggest_team = 0
    for one_team in set(team_numbers):
        size = int((team_numbers == one_team).sum())
        if size > biggest_team:
            biggest_team = size

    if biggest_team > 0.5 * len(gene_names):
        print("")
        print("WARNING: one team holds over half the genes. That is the")
        print("disease-domination problem in the comment above. Use the")
        print("'within' network for your biology, not this one.")

    # ---- A note on sex genes ----
    # Some gene teams are really just "was the donor male or female".
    # TXLNGY sits on the Y chromosome, ZNF711 on the X. If they turn
    # up as hubs, that team is about donor sex, not about leukaemia.
    for sex_gene in ["TXLNGY", "ZNF711", "XIST", "RPS4Y1", "DDX3Y"]:
        if sex_gene in list(gene_names):
            print("NOTE:", sex_gene, "is in the network. It is a sex")
            print("  chromosome gene, so ignore it biologically.")

    return dots


def find_strong_pairs(grid, gene_names, strength):
    # Turn a correlation grid into a list of "gene A - gene B" lines.
    #
    # We only take the top-right half of the grid. The grid is a
    # mirror image of itself, so taking everything would list every
    # pair twice.

    top_half = np.triu_indices(len(gene_names), k=1)
    all_strengths = grid[top_half]

    # Keep only the strong ones. We use abs() because a strong
    # NEGATIVE correlation is just as interesting as a positive one.
    is_strong = np.abs(all_strengths) >= strength

    lines = pd.DataFrame()
    lines["gene_a"] = gene_names[top_half[0][is_strong]]
    lines["gene_b"] = gene_names[top_half[1][is_strong]]
    lines["strength"] = np.round(all_strengths[is_strong], 4)

    direction = []
    for one_strength in all_strengths[is_strong]:
        if one_strength > 0:
            direction.append("same_way")
        else:
            direction.append("opposite_ways")
    lines["direction"] = direction

    return lines


def count_lines_per_gene(lines, gene_names):
    # How many lines does each gene have?

    counter = {}
    for one_gene in gene_names:
        counter[one_gene] = 0

    for row in lines.itertuples(index=False):
        counter[row.gene_a] = counter[row.gene_a] + 1
        counter[row.gene_b] = counter[row.gene_b] + 1

    answer = []
    for one_gene in gene_names:
        answer.append(counter[one_gene])

    return answer



# =============================================================
# STEP 13 - GO AND KEGG (WHAT DO THESE GENES ACTUALLY DO?)
# =============================================================

def read_gene_sets(file_path):
    # Read one "gene set library" file.
    #
    # WHAT IS A GENE SET? Biologists have spent decades writing down
    # lists like:
    #
    #   "T cell activation" -> CD3D, CD3E, LCK, ZAP70, ...
    #   "DNA repair"        -> BRCA1, BRCA2, ATM, TP53, ...
    #
    # GO (Gene Ontology) lists what genes DO. KEGG lists PATHWAYS -
    # chains of genes that work in sequence. If our disease genes fill
    # up one of these lists, that tells us what is going wrong in the
    # disease.
    #
    # The file format is one set per line, separated by tabs:
    #   set name <tab> [description] <tab> GENE <tab> GENE <tab> ...
    #
    # Annoyingly the description is sometimes there and sometimes not.
    # A description contains spaces; a gene name never does. So we
    # skip any piece with a space in it.

    all_sets = {}

    open_file = open(file_path, "r", encoding="utf-8", errors="replace")

    for line in open_file:
        pieces = line.strip().split("\t")

        if len(pieces) < 2:
            continue

        set_name = pieces[0].strip()
        genes_in_set = set()

        for one_piece in pieces[1:]:
            one_piece = one_piece.strip()

            if one_piece == "":
                continue

            # A description or a web link, not a gene
            if " " in one_piece:
                continue
            if one_piece.lower().startswith("http"):
                continue

            # Sometimes a weight is stuck on like "CD3D,1.0"
            clean_gene = one_piece.split(",")[0].strip().upper()

            if clean_gene != "":
                genes_in_set.add(clean_gene)

        if set_name != "" and len(genes_in_set) >= 2:
            all_sets[set_name] = genes_in_set

    open_file.close()
    return all_sets


def test_one_gene_set(my_genes, background_genes, one_set):
    # TEST 1: THE COUNTING TEST (over-representation)
    #
    # The question: "I have 200 disease genes. 15 of them are in the
    # 'T cell activation' list. Is 15 more than I would expect by
    # pure chance?"
    #
    # We work out the expected number, then use the hypergeometric
    # test for the p-value. That test is exactly the maths of drawing
    # coloured balls from a bag without putting them back.
    #
    # ==========================================================
    # THE BACKGROUND IS THE PART PEOPLE GET WRONG
    # ==========================================================
    # "By chance" needs a pool to compare against. We must use OUR
    # 14,781 measured genes, NOT all 20,000 human genes.
    #
    # Why: our chip only measures some genes, and we filtered out the
    # switched-off ones. If we compared against the whole genome, then
    # simply "being measurable" would look like enrichment. This is
    # the single most common way people accidentally invent results.

    # Only count genes we actually measured
    my_genes_in_pool = set(my_genes) & set(background_genes)
    set_genes_in_pool = one_set & set(background_genes)

    pool_size = len(background_genes)
    how_many_mine = len(my_genes_in_pool)
    how_many_in_set = len(set_genes_in_pool)

    overlap = my_genes_in_pool & set_genes_in_pool
    how_many_overlap = len(overlap)

    # Too small to say anything
    if how_many_overlap < 2:
        return None

    expected = how_many_mine * how_many_in_set / pool_size

    if expected > 0:
        fold = how_many_overlap / expected
    else:
        fold = np.nan

    # sf(x - 1, ...) means "the chance of getting x OR MORE"
    p_value = stats.hypergeom.sf(how_many_overlap - 1, pool_size,
                                 how_many_in_set, how_many_mine)

    gene_names_found = ";".join(sorted(overlap)[0:30])

    answer = {}
    answer["set_name"] = "placeholder"
    answer["genes_overlapping"] = how_many_overlap
    answer["genes_in_set"] = how_many_in_set
    answer["genes_i_gave"] = how_many_mine
    answer["expected_by_chance"] = round(expected, 3)
    answer["fold_more_than_expected"] = fold
    answer["p_value"] = p_value
    answer["which_genes"] = gene_names_found
    return answer


def test_one_gene_set_by_rank(all_genes, all_scores, one_set):
    # TEST 2: THE RANKING TEST (no cut-off needed)
    #
    # The counting test has a weakness: it needs us to draw a line and
    # say "these 200 genes are my disease genes". A gene just below
    # the line is thrown away completely.
    #
    # This test avoids the line. We line up ALL 14,781 genes from most
    # changed to least changed, then ask:
    #
    #   "are the 'T cell activation' genes bunched up near the top of
    #    my list, or scattered evenly through it?"
    #
    # Bunched near the top = the whole process is shifted, even if no
    # single gene passed our cut-off. This is a Mann-Whitney U test on
    # the positions in the list.
    #
    # Running BOTH tests is the point. A set that both tests agree on
    # is much stronger evidence than one only found by either.

    gene_positions = {}
    for i in range(len(all_genes)):
        gene_positions[str(all_genes[i]).upper()] = i

    # Turn the scores into positions in the queue (1st, 2nd, 3rd...)
    positions_in_queue = stats.rankdata(np.asarray(all_scores, dtype=float))

    # Where in the queue are this set's genes?
    my_positions = []
    for one_gene in one_set:
        if one_gene in gene_positions:
            my_positions.append(positions_in_queue[gene_positions[one_gene]])

    set_size = len(my_positions)
    total_genes = len(all_genes)

    if set_size < ENRICHMENT_MIN_SET or set_size >= total_genes:
        return None

    # The Mann-Whitney U statistic
    sum_of_positions = float(np.sum(my_positions))
    u_value = sum_of_positions - set_size * (set_size + 1) / 2.0

    # What U would we expect if the genes were scattered randomly?
    expected_u = set_size * (total_genes - set_size) / 2.0
    wobble = np.sqrt(set_size * (total_genes - set_size) *
                     (total_genes + 1) / 12.0)

    if wobble == 0:
        return None

    # How many wobbles away from expected are we?
    z_value = (u_value - expected_u) / wobble

    # Turn that into a p-value (times 2 because we look both ways)
    p_value = 2.0 * stats.norm.sf(abs(z_value))

    if z_value > 0:
        direction = "shifted_up"
    else:
        direction = "shifted_down"

    answer = {}
    answer["set_name"] = "placeholder"
    answer["genes_in_set"] = set_size
    answer["auc"] = u_value / (set_size * (total_genes - set_size))
    answer["z_value"] = z_value
    answer["direction"] = direction
    answer["p_value"] = p_value
    return answer


def tidy_set_name(raw_name):
    # Gene set names arrive with an ID stuck on the end, like
    #   "T cell activation (GO:0042110)"
    # We split that into a readable label and the ID.

    label = raw_name
    set_id = ""

    if "(GO:" in raw_name:
        pieces = raw_name.split("(GO:")
        label = pieces[0].strip()
        set_id = "GO:" + pieces[1].replace(")", "").strip()

    return label, set_id


def do_enrichment(results, gene_data, size_cutoff):
    # Run both tests, for every disease, against every library.

    print("Looking up what the disease genes actually do...")

    if os.path.exists(GENESET_FOLDER) == False:
        print("No gene set folder at:", GENESET_FOLDER)
        print("Get the files first:  python3 simple_get_genesets.py")
        return None

    # Find the library files
    library_files = []
    for one_name in sorted(os.listdir(GENESET_FOLDER)):
        if one_name.endswith(".txt") or one_name.endswith(".gmt"):
            library_files.append(one_name)

    if len(library_files) == 0:
        print("The gene set folder is empty.")
        print("Get the files first:  python3 simple_get_genesets.py")
        return None

    # ---- The background pool: our measured, switched-on genes ----
    background_genes = []
    for one_gene in gene_data.index:
        background_genes.append(str(one_gene).upper())

    print("Background pool:", len(background_genes), "measured genes")
    print("(NOT all 20,000 human genes - see the comment above)")

    all_diseases = sorted(results["group"].unique())

    counting_results = []
    ranking_results = []
    summary_rows = []

    for one_file in library_files:
        full_path = GENESET_FOLDER + "/" + one_file
        library_name = one_file.replace(".txt", "").replace(".gmt", "")

        gene_sets = read_gene_sets(full_path)

        if len(gene_sets) == 0:
            print(library_name, "- could not read any gene sets, skipping")
            continue

        # How much of this library do we even measure?
        library_genes = set()
        for one_set_name in gene_sets:
            library_genes = library_genes | gene_sets[one_set_name]

        overlap_with_us = len(library_genes & set(background_genes))

        print("")
        print(library_name, "-", len(gene_sets), "gene sets,",
              overlap_with_us, "of its genes measured by us")

        if overlap_with_us < 500:
            print("  WARNING: poor overlap. Are the gene names the right kind?")

        for one_disease in all_diseases:
            just_this_disease = results[results["group"] == one_disease]

            # My disease genes = significant AND clearly up
            is_significant = just_this_disease["fdr"] < FDR_CUTOFF
            is_big = just_this_disease["difference"] >= size_cutoff
            my_genes = []
            for one_gene in just_this_disease[is_significant & is_big]["gene"]:
                my_genes.append(str(one_gene).upper())

            # ---- Test 1: counting ----
            counting_rows = []
            for one_set_name in gene_sets:
                one_set = gene_sets[one_set_name]

                # Skip sets that are too small or absurdly broad
                size_in_pool = len(one_set & set(background_genes))
                if size_in_pool < ENRICHMENT_MIN_SET:
                    continue
                if size_in_pool > ENRICHMENT_MAX_SET:
                    continue

                one_answer = test_one_gene_set(my_genes, background_genes,
                                               one_set)
                if one_answer is None:
                    continue

                one_answer["set_name"] = one_set_name
                counting_rows.append(one_answer)

            how_many_counting = 0
            if len(counting_rows) > 0:
                counting_table = pd.DataFrame(counting_rows)
                counting_table["fdr"] = correct_pvalues(
                    counting_table["p_value"].values)
                counting_table["library"] = library_name
                counting_table["disease"] = one_disease
                counting_table["test"] = "counting"

                labels = []
                ids = []
                for one_name in counting_table["set_name"]:
                    one_label, one_id = tidy_set_name(one_name)
                    labels.append(one_label)
                    ids.append(one_id)
                counting_table["label"] = labels
                counting_table["set_id"] = ids

                how_many_counting = int(
                    (counting_table["fdr"] < FDR_CUTOFF).sum())
                counting_results.append(counting_table)

            # ---- Test 2: ranking ----
            how_many_ranking = 0
            ranking_rows = []

            for one_set_name in gene_sets:
                one_answer = test_one_gene_set_by_rank(
                    just_this_disease["gene"].values,
                    just_this_disease["t_value"].values,
                    gene_sets[one_set_name])

                if one_answer is None:
                    continue

                one_answer["set_name"] = one_set_name
                ranking_rows.append(one_answer)

            if len(ranking_rows) > 0:
                ranking_table = pd.DataFrame(ranking_rows)
                ranking_table["fdr"] = correct_pvalues(
                    ranking_table["p_value"].values)
                ranking_table["library"] = library_name
                ranking_table["disease"] = one_disease
                ranking_table["test"] = "ranking"

                labels = []
                ids = []
                for one_name in ranking_table["set_name"]:
                    one_label, one_id = tidy_set_name(one_name)
                    labels.append(one_label)
                    ids.append(one_id)
                ranking_table["label"] = labels
                ranking_table["set_id"] = ids

                how_many_ranking = int(
                    (ranking_table["fdr"] < FDR_CUTOFF).sum())
                ranking_results.append(ranking_table)

            summary_rows.append({"library": library_name,
                                 "disease": one_disease,
                                 "genes_i_gave": len(my_genes),
                                 "counting_test_hits": how_many_counting,
                                 "ranking_test_hits": how_many_ranking})

            print("   ", one_disease[0:28], "-> gave", len(my_genes),
                  "genes, counting:", how_many_counting,
                  " ranking:", how_many_ranking)

    if len(counting_results) == 0 and len(ranking_results) == 0:
        print("No enrichment results at all.")
        return None

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(RESULT_FOLDER + "/step13_enrichment_summary.csv",
                   index=False)

    # ---- Save the full results ----
    keep_columns = ["library", "disease", "test", "label", "set_id",
                    "genes_in_set", "p_value", "fdr"]

    best_pieces = []

    if len(counting_results) > 0:
        everything = pd.concat(counting_results, ignore_index=True)
        everything.to_csv(RESULT_FOLDER + "/step13_counting_test.csv",
                          index=False)

        good = everything[everything["fdr"] < FDR_CUTOFF].copy()
        good["effect"] = good["fold_more_than_expected"]

        for one_disease in good["disease"].unique():
            just_one = good[good["disease"] == one_disease]
            just_one = just_one.sort_values("p_value")
            best_pieces.append(just_one.head(ENRICHMENT_TOP)[
                keep_columns + ["effect"]])

    if len(ranking_results) > 0:
        everything = pd.concat(ranking_results, ignore_index=True)
        everything.to_csv(RESULT_FOLDER + "/step13_ranking_test.csv",
                          index=False)

        good = everything[everything["fdr"] < FDR_CUTOFF]
        good = good[good["direction"] == "shifted_up"].copy()
        good["effect"] = good["auc"]

        for one_disease in good["disease"].unique():
            just_one = good[good["disease"] == one_disease]
            just_one = just_one.sort_values("p_value")
            best_pieces.append(just_one.head(ENRICHMENT_TOP)[
                keep_columns + ["effect"]])

    if len(best_pieces) == 0:
        print("Nothing was significant in either test.")
        return None

    best = pd.concat(best_pieces, ignore_index=True)
    best.to_csv(RESULT_FOLDER + "/step13_enrichment_best.csv", index=False)

    # ---- Which results did BOTH tests agree on? ----
    # These are the ones to lead with in the write-up. Agreement
    # between two different tests is much stronger than either alone.
    print("")
    print("Results BOTH tests agreed on (report these first):")

    how_many_agreed = 0

    for one_disease in best["disease"].unique():
        just_one = best[best["disease"] == one_disease]

        for one_label in just_one["label"].unique():
            same_label = just_one[just_one["label"] == one_label]

            if len(same_label["test"].unique()) > 1:
                how_many_agreed = how_many_agreed + 1
                if how_many_agreed <= 10:
                    print("  ", one_disease[0:25], "->", one_label[0:40])

    print("  (", how_many_agreed, "in total )")

    if how_many_agreed == 0:
        print("  none - the two tests found different things. Say so.")

    return best



# =============================================================
# STEP 14 - WRITE THE FINAL REPORT
# =============================================================

def read_result_file(file_name):
    # Read one of our own result files, or None if it is not there.
    # Returning None instead of crashing means a missing file only
    # costs us one section of the report, not the whole thing.

    full_path = RESULT_FOLDER + "/" + file_name

    if os.path.exists(full_path) == False:
        return None

    try:
        table = pd.read_csv(full_path)
    except Exception:
        return None

    if len(table) == 0:
        return None

    return table


def make_markdown_table(table, columns=None, how_many_rows=15):
    # Turn a pandas table into markdown text you can paste into a
    # document. Markdown tables look like:
    #
    #   | gene | score |
    #   |---|---|
    #   | CD3D | 0.97 |

    if table is None:
        return "_not available - that step has not been run_\n"

    if columns is not None:
        # Only keep columns that really exist
        safe_columns = []
        for one_column in columns:
            if one_column in table.columns:
                safe_columns.append(one_column)
        if len(safe_columns) == 0:
            return "_not available_\n"
        small = table[safe_columns]
    else:
        small = table

    small = small.head(how_many_rows)

    lines = []

    # The header row
    header = "| "
    for one_column in small.columns:
        header = header + str(one_column) + " | "
    lines.append(header)

    # The dashes under the header
    dashes = "|"
    for one_column in small.columns:
        dashes = dashes + "---|"
    lines.append(dashes)

    # The data rows
    for row in small.itertuples(index=False):
        one_line = "| "
        for one_value in row:
            if isinstance(one_value, float):
                one_line = one_line + str(round(one_value, 4)) + " | "
            else:
                one_line = one_line + str(one_value) + " | "
        lines.append(one_line)

    return "\n".join(lines) + "\n"


def write_final_report():
    # Put every number we calculated into one document.
    #
    # WHY BOTHER? Copying numbers by hand from 20 CSV files into a
    # thesis is slow and it is where typing mistakes creep in. This
    # fills in every number automatically.
    #
    # Where the report says TODO, that is a place where a human has to
    # write the biology. A computer can count genes; it cannot tell
    # you what they mean.

    print("Writing the final report...")

    lines = []

    lines.append("# Leukaemia Biomarker Report")
    lines.append("")
    lines.append("Dataset: GSE13159 (MILE study)")
    lines.append("")
    lines.append("Numbers below are filled in automatically.")
    lines.append("Every **TODO** is a place where you write the biology.")
    lines.append("")

    # ---- 1. The data ----
    lines.append("## 1. The data and cleaning")
    lines.append("")

    missing = read_result_file("step2_missing_values.csv")
    if missing is not None:
        lines.append(make_markdown_table(missing))

    quality = read_result_file("step2_quality.csv")
    if quality is not None:
        average_similarity = round(float(quality["mean_correlation"].mean()), 3)
        how_many_odd = int(quality["outlier"].sum())

        lines.append("")
        lines.append("- Average similarity between chips: "
                     + str(average_similarity))
        lines.append("- Unusual chips: " + str(how_many_odd)
                     + " (flagged, NOT removed)")

    lines.append("")
    lines.append("- The values in this dataset run from 0 to 1. They were")
    lines.append("  rescaled before being shared, so differences between")
    lines.append("  group averages are NOT log2 fold changes.")
    lines.append("")
    lines.append("**TODO** one paragraph on why this cohort suits the question.")
    lines.append("")

    # ---- 2. Clustering ----
    lines.append("## 2. Grouping patients by their genes")
    lines.append("")

    cluster_summary = read_result_file("step4_summary.csv")
    if cluster_summary is not None:
        lines.append(make_markdown_table(cluster_summary))

    lines.append("")
    lines.append("Two ways of choosing the number of groups:")
    lines.append("")

    k_choice = read_result_file("step4_k_choice.csv")
    if k_choice is not None:
        lines.append(make_markdown_table(k_choice))

    lines.append("")
    lines.append("**TODO** the two methods may disagree. Say which you used")
    lines.append("and why. That disagreement is a finding, not a problem.")
    lines.append("")

    comparison = read_result_file("step7_cluster_vs_diagnosis.csv")
    if comparison is not None:
        lines.append("Which clusters hold which diseases:")
        lines.append("")
        lines.append(make_markdown_table(comparison, how_many_rows=25))
        lines.append("")
        lines.append("**TODO** which clusters are clean, which mix two")
        lines.append("diseases, and whether those mixes make biological sense.")
        lines.append("")

    # ---- 3. Gene network ----
    network = read_result_file("step12_network_genes.csv")
    if network is not None:
        lines.append("## 3. Gene network")
        lines.append("")
        lines.append("Genes with the most connections (the hubs):")
        lines.append("")
        lines.append(make_markdown_table(
            network, ["gene", "team", "lines_everyone", "lines_within"], 15))
        lines.append("")
        lines.append("**TODO** say which network you interpret. The")
        lines.append("'within clusters' one is the meaningful one, because")
        lines.append("the other is dominated by disease identity.")
        lines.append("")

    # ---- 4. Differential expression ----
    lines.append("## 4. Genes that differ between diseases")
    lines.append("")

    de_summary = read_result_file("step8_summary.csv")
    if de_summary is not None:
        lines.append(make_markdown_table(de_summary, None, 25))

    cutoff = read_result_file("step8_cutoff.csv")
    if cutoff is not None:
        lines.append("")
        lines.append("Difference cut-off used: "
                     + str(round(float(cutoff["cutoff"].iloc[0]), 4)))
        lines.append("(the top 1% of the differences we actually saw)")

    lines.append("")
    lines.append("**TODO** comment on group size versus gene counts. Small")
    lines.append("groups have less statistical power, but a very uniform")
    lines.append("disease can still give tiny p-values.")
    lines.append("")

    # ---- 5. Enrichment ----
    enrichment = read_result_file("step13_enrichment_best.csv")
    if enrichment is not None:
        lines.append("## 5. What those genes do (GO and KEGG)")
        lines.append("")
        lines.append(make_markdown_table(
            enrichment,
            ["library", "disease", "test", "label", "effect", "fdr"], 20))
        lines.append("")
        lines.append("**TODO** lead with anything BOTH tests agreed on.")
        lines.append("State that the background was our own measured genes,")
        lines.append("not the whole genome.")
        lines.append("")

    # ---- 6. Biomarkers ----
    lines.append("## 6. Candidate biomarkers")
    lines.append("")

    candidate_summary = read_result_file("step9_summary.csv")
    if candidate_summary is not None:
        lines.append("How many candidates each disease contributed:")
        lines.append("")
        lines.append(make_markdown_table(candidate_summary, None, 20))
        lines.append("")

    candidates = read_result_file("step9_candidates.csv")
    if candidates is not None:
        lines.append(make_markdown_table(
            candidates,
            ["gene", "group", "difference", "fdr", "specificity", "score"], 20))
        lines.append("")

    # ---- 7. Machine learning ----
    lines.append("## 7. Machine learning")
    lines.append("")

    model_scores = read_result_file("step10_model_scores.csv")
    if model_scores is not None:
        lines.append(make_markdown_table(model_scores))
        lines.append("")
        lines.append("The labels were the doctors' diagnoses, NOT our own")
        lines.append("clusters, so these numbers are not circular. Gene")
        lines.append("picking sat inside the pipeline and never saw the")
        lines.append("hidden test patients.")
        lines.append("")

    # Which model did Step 10 actually end up using?
    chosen_table = read_result_file("step10_chosen_model.csv")
    if chosen_table is not None:
        which_model = str(chosen_table["value"].iloc[0])
        lines.append("Model used for the gene ranking: " + which_model)
        lines.append("Chosen: " + str(chosen_table["value"].iloc[1]))
        lines.append("")
    else:
        which_model = MODEL_NAMES[0]

    importance = read_result_file("step10_gene_importance.csv")
    if importance is not None:
        main_only = importance[importance["model"] == which_model]
        lines.append("Genes the " + which_model + " leaned on most:")
        lines.append("")
        lines.append(make_markdown_table(main_only, ["gene", "importance"], 15))
        lines.append("")

    # ---- 8. Consensus and testing ----
    lines.append("## 8. Consensus biomarkers and how well they work")
    lines.append("")

    consensus = read_result_file("step11_consensus.csv")
    if consensus is not None:
        lines.append(str(len(consensus)) + " genes have BOTH statistical and")
        lines.append("machine learning support.")
        lines.append("")
        lines.append(make_markdown_table(
            consensus,
            ["gene", "disease", "model_rank", "consensus_score",
             "auc_on_hidden_patients"], 20))
        lines.append("")

    panel = read_result_file("step11_panel.csv")
    if panel is not None:
        lines.append("The small gene panel on its own:")
        lines.append("")
        lines.append(make_markdown_table(panel))
        lines.append("")

    levels = read_result_file("step11_evidence_levels.csv")
    if levels is not None:
        lines.append("How strong is the evidence?")
        lines.append("")
        lines.append(make_markdown_table(levels))
        lines.append("")

    external = read_result_file("external_results.csv")
    if external is not None:
        lines.append("### A different cohort")
        lines.append("")
        lines.append(make_markdown_table(external))
        lines.append("")
    else:
        lines.append("**No different cohort was tested.** Everything above")
        lines.append("comes from this one dataset. That is a limitation, not")
        lines.append("a failure, but it must be said plainly: nothing here")
        lines.append("is validated outside GSE13159.")
        lines.append("")

    lines.append("**TODO** for each top gene, one or two sentences of")
    lines.append("literature background. Genes already known to mark that")
    lines.append("disease are your proof the pipeline works. Unexpected")
    lines.append("genes are hypotheses, NOT discoveries.")
    lines.append("")

    # ---- 9. Limitations (mostly pre-written) ----
    lines.append("## 9. Limitations")
    lines.append("")
    lines.append("- We compared each disease against ALL the others lumped")
    lines.append("  together. That reference group is a mixture, which can")
    lines.append("  inflate significance for a gene missing from just one")
    lines.append("  big competing disease.")
    lines.append("- Where several probes measured one gene, we kept the")
    lines.append("  strongest. A different rule would give a slightly")
    lines.append("  different gene list.")
    lines.append("- These are bulk tissue samples: a difference may reflect")
    lines.append("  a change in which CELLS are present, rather than genes")
    lines.append("  changing inside one cell type.")
    lines.append("- We measured whether technical factors explain PC1, but")
    lines.append("  we did not correct for them.")
    lines.append("- The values are on a 0-to-1 scale, so effect sizes cannot")
    lines.append("  be reported as fold changes.")
    lines.append("")
    lines.append("**TODO** add anything else you ran into.")
    lines.append("")

    lines.append("## Reproducing this")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 simple_pipeline.py")
    lines.append("Rscript simple_figures.R")
    lines.append("```")
    lines.append("")
    lines.append("Random seed " + str(SEED) + " throughout.")
    lines.append("")

    report_text = "\n".join(lines)

    file_path = RESULT_FOLDER + "/FINAL_REPORT.md"
    open(file_path, "w", encoding="utf-8").write(report_text)

    how_many_todo = report_text.count("TODO")

    print("Report written to:", file_path)
    print(how_many_todo, "TODO markers left for you to write")


# =============================================================
# MAIN PROGRAM - read this to see the whole story
# =============================================================

def mymain():
    print("=" * 60)
    print("SIMPLE LEUKAEMIA GENE EXPRESSION PIPELINE")
    print("=" * 60)

    make_folders()
    check_the_data_file()

    # ---- Load ----
    print("\n--- STEP 1: LOAD THE DATA ---")
    lines_to_skip = find_table_start()
    patients = load_patient_information()
    data = load_gene_data(lines_to_skip)

    # ---- Clean ----
    print("\n--- STEP 2: CLEAN THE DATA ---")
    data = fill_missing_values(data)
    data = remove_control_probes(data)
    data = normalize_data(data)
    check_quality(data)

    # Keep a copy of the full clean table. Step 6 needs ALL the
    # genes, not just the most variable ones.
    clean_data = data.copy()

    # ---- Filter and scale ----
    print("\n--- STEP 3: PICK THE USEFUL PROBES ---")
    data = filter_probes(data)
    scaled_data = zscore_data(data)

    # ---- Group patients without using the diagnosis ----
    print("\n--- STEP 4: CLUSTER THE PATIENTS ---")
    clusters = do_clustering(scaled_data)

    print("\n--- STEP 5: PCA ---")
    pca_table = do_pca(scaled_data, clusters)

    # ---- Move to gene level ----
    print("\n--- STEP 6: PROBES TO GENES ---")
    gene_data = map_probes_to_genes(clean_data)

    # ---- Bring in the doctors' diagnoses ----
    print("\n--- STEP 7: MAKE PATIENT GROUPS ---")
    groups = make_groups(gene_data, patients, clusters)

    # Was PC1 driven by biology, or by something technical?
    check_for_batch_effects(pca_table, patients)

    # ---- Statistics ----
    print("\n--- STEP 8: FIND GENES THAT DIFFER ---")
    results, size_cutoff = find_different_genes(gene_data, groups)

    print("\n--- STEP 9: PICK BIOMARKERS ---")
    candidates = find_biomarkers(results, gene_data, size_cutoff)

    # ---- Machine learning ----
    print("\n--- STEP 10: MACHINE LEARNING AND SHAP ---")
    importance, per_class = run_machine_learning(gene_data, groups)

    # Is the model leaning on the same genes the statistics found?
    check_candidates_against_shap(candidates, per_class)

    # ---- Combine the two kinds of evidence, then TEST it ----
    print("\n--- STEP 11: COMBINE AND TEST THE EVIDENCE ---")
    consensus = make_consensus(candidates, importance)

    tested = test_each_biomarker(consensus, gene_data, groups)

    if tested is not None:
        consensus, train_positions, test_positions = tested
        test_the_whole_panel(consensus, gene_data, groups,
                             train_positions, test_positions)

    # ---- Which genes work together? ----
    print("\n--- STEP 12: GENE NETWORK ---")
    build_gene_network(gene_data, groups)

    # ---- What do the disease genes actually do? ----
    print("\n--- STEP 13: GO AND KEGG ---")
    do_enrichment(results, gene_data, size_cutoff)

    # ---- Put every number in one document ----
    print("\n--- STEP 14: FINAL REPORT ---")
    write_final_report()

    print("\n" + "=" * 60)
    print("FINISHED. All answers are in the '" + RESULT_FOLDER + "' folder.")
    print("Report: " + RESULT_FOLDER + "/FINAL_REPORT.md")
    print("Now run: Rscript simple_figures.R  to draw the pictures")
    print("")
    print("For the strongest kind of proof, test on a different cohort:")
    print("  python3 simple_external_validation.py")
    print("=" * 60)


if __name__ == "__main__":
    mymain()
