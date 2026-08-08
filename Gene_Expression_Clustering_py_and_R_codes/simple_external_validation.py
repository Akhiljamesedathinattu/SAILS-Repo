#!/usr/bin/env python3

# =============================================================
# TEST THE MODEL ON A COMPLETELY DIFFERENT COHORT
#
# THIS IS THE ONLY STEP THAT COUNTS AS REAL VALIDATION.
#
# Everything in simple_pipeline.py used ONE dataset. Even the
# "hidden" test patients came from the same hospitals, the same
# machines, the same batch of chips, processed the same way. A model
# can learn those local quirks and still score well.
#
# The real test is: take the frozen, already-trained model, show it
# patients from a different study it has never seen, and see if it
# still works.
#
# WHAT YOU NEED FIRST
#   1. Run simple_pipeline.py, so models/ contains a saved model
#   2. Download another GPL570 leukaemia dataset from GEO as a
#      series matrix .gz file, and put it in the raw folder
#
# THEN EDIT the OTHER_DATA_FILE line below and run:
#   python3 simple_external_validation.py
#
# EXPECT THE SCORE TO DROP. A different cohort means different
# equipment and different handling. A drop is normal and is itself
# a result worth reporting. Do not hide it.
# =============================================================

import gzip
import os

import numpy as np
import pandas as pd

from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.metrics import confusion_matrix


# ---- Settings ----

BASE_FOLDER = "/home/sails/SAILS-Repo/Gene_Expression_Clustering"

RAW_FOLDER = BASE_FOLDER + "/raw"
RESULT_FOLDER = BASE_FOLDER + "/results"
MODEL_FOLDER = BASE_FOLDER + "/models"

# CHANGE THIS to the other dataset you downloaded
OTHER_DATA_FILE = RAW_FOLDER + "/GSE13164_series_matrix.txt.gz"

ANNOTATION_FILE = RAW_FOLDER + "/GPL570_annot.csv"

WHICH_MODEL = "random_forest"


def check_everything_is_ready():
    # Look for all three things we need before we start, so a missing
    # file gives a clear message instead of a crash halfway through.

    model_path = MODEL_FOLDER + "/" + WHICH_MODEL + ".joblib"

    if os.path.exists(model_path) == False:
        print("ERROR: no saved model at", model_path)
        print("Run simple_pipeline.py first.")
        raise SystemExit(1)

    if os.path.exists(OTHER_DATA_FILE) == False:
        print("ERROR: cannot find the other dataset:")
        print("  " + OTHER_DATA_FILE)
        print("")
        print("Download another GPL570 leukaemia series matrix from GEO,")
        print("put it in the raw folder, and edit OTHER_DATA_FILE at the")
        print("top of this script.")
        raise SystemExit(1)

    if os.path.exists(ANNOTATION_FILE) == False:
        print("ERROR: no annotation file at", ANNOTATION_FILE)
        print("Run: Rscript simple_make_annotation.R")
        raise SystemExit(1)

    print("Everything needed is present")
    return model_path


def load_other_dataset():
    # Read the other cohort's series matrix. This is the same job as
    # Step 1 of the pipeline, so the code looks familiar.

    print("Reading the other cohort...")

    # ---- Find where the numbers start ----
    line_number = 0
    found = False

    open_file = gzip.open(OTHER_DATA_FILE, "rt", errors="replace")
    for line in open_file:
        if line.startswith("!series_matrix_table_begin"):
            found = True
            break
        line_number = line_number + 1
    open_file.close()

    if found == False:
        print("ERROR: that file is not a GEO series matrix.")
        raise SystemExit(1)

    lines_to_skip = line_number + 1

    # ---- Get the patient IDs and diagnoses ----
    patient_ids = []
    all_facts = {}

    open_file = gzip.open(OTHER_DATA_FILE, "rt", errors="replace")
    for line in open_file:
        if line.startswith("!series_matrix_table_begin"):
            break

        pieces = line.strip().split("\t")

        if pieces[0] == "!Sample_geo_accession":
            for one_piece in pieces[1:]:
                patient_ids.append(one_piece.replace('"', ""))

        if pieces[0] == "!Sample_characteristics_ch1":
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
            if this_label == "":
                this_label = "unnamed_" + str(len(all_facts) + 1)
            all_facts[this_label] = these_values
    open_file.close()

    # Pick the diagnosis line, the same way the main pipeline does
    chosen_label = ""
    for one_label in all_facts:
        if "leukemia" in one_label.lower():
            chosen_label = one_label
            break

    diseases = None
    if chosen_label != "":
        print("Found diagnoses under '" + chosen_label + "'")
        diseases = all_facts[chosen_label]
    else:
        print("No diagnosis field found. We can still make predictions,")
        print("but we cannot score them without the right answers.")

    # ---- Read the numbers ----
    data = pd.read_csv(OTHER_DATA_FILE, sep="\t",
                       skiprows=lines_to_skip, index_col=0,
                       low_memory=False)

    good_rows = []
    for probe_name in data.index:
        if str(probe_name).startswith("!"):
            good_rows.append(False)
        else:
            good_rows.append(True)
    data = data[good_rows]
    data = data.astype("float32")

    clean_names = []
    for one_name in data.columns:
        clean_names.append(str(one_name).replace('"', ""))
    data.columns = clean_names

    print("Other cohort:", data.shape[0], "probes x", data.shape[1],
          "patients")

    patients = None
    if diseases is not None:
        if len(diseases) == len(patient_ids):
            patients = pd.DataFrame()
            patients["sample"] = patient_ids
            patients["disease"] = diseases
            patients = patients.set_index("sample")

    return data, patients


def prepare_the_other_data(data, our_gene_names):
    # Put the other cohort through the SAME preparation the model was
    # trained on. A model trained on prepared data cannot understand
    # raw data - the numbers would be on a different scale entirely.

    print("Preparing the other cohort the same way...")

    # ---- Remove control probes ----
    keep_rows = []
    for probe_name in data.index:
        if str(probe_name).upper().startswith("AFFX"):
            keep_rows.append(False)
        else:
            keep_rows.append(True)
    data = data[keep_rows]

    # ---- Turn probes into genes ----
    annotation = pd.read_csv(ANNOTATION_FILE)
    probe_to_gene = {}
    for row in annotation.itertuples(index=False):
        probe_to_gene[str(row.probe)] = str(row.symbol)

    average_value = data.mean(axis=1)

    best_probe = {}
    best_average = {}

    for probe_name in data.index:
        probe_name = str(probe_name)
        if probe_name not in probe_to_gene:
            continue

        gene_name = probe_to_gene[probe_name]
        this_average = average_value[probe_name]

        if gene_name not in best_average:
            best_probe[gene_name] = probe_name
            best_average[gene_name] = this_average
        elif this_average > best_average[gene_name]:
            best_probe[gene_name] = probe_name
            best_average[gene_name] = this_average

    chosen_probes = []
    chosen_genes = []
    for gene_name in best_probe:
        chosen_genes.append(gene_name)
        chosen_probes.append(best_probe[gene_name])

    gene_data = data.loc[chosen_probes]
    gene_data.index = chosen_genes

    print("Other cohort has", gene_data.shape[0], "named genes")

    # ---- Line the genes up with what the model expects ----
    #
    # The model expects exactly our gene list, in exactly our order.
    # Some of our genes will be missing from the other cohort.
    #
    # HONEST WARNING: for a missing gene we have to put in SOMETHING.
    # We use the average, which is the least biased filler available,
    # but it pushes predictions towards whatever is most common. The
    # more genes are missing, the less the result means.

    how_many_present = 0
    for one_gene in our_gene_names:
        if one_gene in gene_data.index:
            how_many_present = how_many_present + 1

    percent_present = 100.0 * how_many_present / len(our_gene_names)

    print(how_many_present, "of", len(our_gene_names),
          "of our genes are present (", round(percent_present, 1), "% )")

    if percent_present < 50.0:
        print("")
        print("WARNING: fewer than half our genes are here. Most of what")
        print("the model sees will be filled-in averages, so treat the")
        print("result with real caution.")

    lined_up = gene_data.reindex(our_gene_names)

    # Fill each missing gene with its own average across the patients
    # we DO have. If a gene is missing entirely, use 0.5 (the middle
    # of the 0-to-1 scale this data uses).
    for one_gene in lined_up.index:
        one_row = lined_up.loc[one_gene]
        if one_row.isna().all():
            lined_up.loc[one_gene] = 0.5
        elif one_row.isna().any():
            lined_up.loc[one_gene] = one_row.fillna(one_row.mean())

    return lined_up, percent_present


def score_the_predictions(model, guesses, chances, patients, sample_names):
    # Compare the model's guesses with the real diagnoses.

    if patients is None:
        print("No real diagnoses available, so no score can be worked out.")
        return None

    real_answers = patients.reindex(sample_names)["disease"].astype(str).values

    known_diseases = list(model.classes_)

    # ---- Only score patients whose disease the model was taught ----
    # The other study may use different names, or include diseases
    # ours never saw. Guessing at those would not be fair either way.
    can_score = []
    for one_answer in real_answers:
        if one_answer in known_diseases:
            can_score.append(True)
        else:
            can_score.append(False)
    can_score = np.array(can_score)

    print(int(can_score.sum()), "of", len(real_answers),
          "patients have a disease our model was taught")

    if int(can_score.sum()) < 10:
        print("")
        print("Too few comparable patients to score anything. The two")
        print("studies probably name their diseases differently. Compare")
        print("the two name lists by hand and rename them to match.")
        print("")
        print("Our model knows these names:")
        for one_name in known_diseases:
            print("  ", one_name)
        return None

    scorable_answers = real_answers[can_score]
    scorable_guesses = guesses[can_score]
    scorable_chances = chances[can_score]

    accuracy = balanced_accuracy_score(scorable_answers, scorable_guesses)

    try:
        if len(known_diseases) == 2:
            auc = roc_auc_score(scorable_answers, scorable_chances[:, 1])
        else:
            auc = roc_auc_score(scorable_answers, scorable_chances,
                                multi_class="ovr", average="macro",
                                labels=known_diseases)
    except ValueError:
        # This happens when some disease has no patients here at all
        auc = np.nan

    print("")
    print("ON THE DIFFERENT COHORT:")
    print("  balanced accuracy:", round(accuracy, 3))
    print("  AUC:", round(auc, 3))

    # ---- Save the confusion matrix ----
    all_names = sorted(set(scorable_answers) | set(scorable_guesses))
    grid = confusion_matrix(scorable_answers, scorable_guesses,
                            labels=all_names)
    confusion = pd.DataFrame(grid, index=all_names, columns=all_names)
    confusion.to_csv(RESULT_FOLDER + "/external_confusion.csv")

    return accuracy, auc, int(can_score.sum())


def compare_with_our_own_score(accuracy, auc):
    # Show the drop between our own hidden patients and the other
    # cohort. The drop IS the honest measure of how far these findings
    # travel beyond one dataset.

    our_scores_file = RESULT_FOLDER + "/step10_model_scores.csv"

    if os.path.exists(our_scores_file) == False:
        return

    our_scores = pd.read_csv(our_scores_file)
    just_our_model = our_scores[our_scores["model"] == WHICH_MODEL]

    if len(just_our_model) == 0:
        return

    if "test_auc" in just_our_model.columns:
        our_auc = float(just_our_model["test_auc"].iloc[0])
    else:
        return

    print("")
    print("THE COMPARISON THAT MATTERS:")
    print("  our own hidden patients : AUC", round(our_auc, 3))
    print("  a different cohort      : AUC", round(auc, 3))
    print("")
    print("Report BOTH numbers. The gap between them is how much of")
    print("our result was about leukaemia, and how much was about")
    print("this one dataset's equipment and handling.")


def mymain():
    print("=" * 60)
    print("TESTING ON A DIFFERENT COHORT")
    print("=" * 60)

    model_path = check_everything_is_ready()

    import joblib
    model = joblib.load(model_path)

    print("Loaded the", WHICH_MODEL, "model")
    print("It was taught", len(model.classes_), "diseases")

    # Which genes does the model expect, in which order?
    gene_list_file = RESULT_FOLDER + "/step6_gene_list.csv"
    if os.path.exists(gene_list_file) == False:
        print("ERROR: cannot find", gene_list_file)
        print("Run simple_pipeline.py first.")
        raise SystemExit(1)

    our_gene_names = pd.read_csv(gene_list_file)["gene"].astype(str).values

    data, patients = load_other_dataset()

    lined_up, percent_present = prepare_the_other_data(data, our_gene_names)

    # Models want patients as rows
    features = lined_up.values.T

    print("")
    print("Making predictions...")
    guesses = model.predict(features)
    chances = model.predict_proba(features)

    # ---- Save every prediction ----
    predictions = pd.DataFrame()
    predictions["sample"] = lined_up.columns
    predictions["predicted_disease"] = guesses
    predictions["how_confident"] = chances.max(axis=1)

    if patients is not None:
        real_answers = patients.reindex(lined_up.columns)["disease"]
        predictions["real_disease"] = real_answers.values

    predictions.to_csv(RESULT_FOLDER + "/external_predictions.csv",
                       index=False)
    print("Predictions saved to results/external_predictions.csv")

    # ---- Score them ----
    scored = score_the_predictions(model, guesses, chances, patients,
                                   lined_up.columns)

    if scored is None:
        print("")
        print("Finished, but with no score.")
        return

    accuracy, auc, how_many_scored = scored

    results = pd.DataFrame()
    results["item"] = ["cohort_file", "model", "patients_scored",
                       "percent_our_genes_present",
                       "balanced_accuracy", "auc"]
    results["value"] = [OTHER_DATA_FILE, WHICH_MODEL, how_many_scored,
                        round(percent_present, 2),
                        round(accuracy, 4), round(auc, 4)]
    results.to_csv(RESULT_FOLDER + "/external_results.csv", index=False)

    compare_with_our_own_score(accuracy, auc)

    print("")
    print("=" * 60)
    print("Saved to results/external_results.csv")
    print("Re-run simple_pipeline.py's report step to include this,")
    print("or just paste the numbers into your write-up.")
    print("=" * 60)


if __name__ == "__main__":
    mymain()
