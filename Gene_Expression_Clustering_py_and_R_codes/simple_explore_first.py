#!/usr/bin/env python3

# =============================================================
# STEP ZERO - LOOK AT THE DATA BEFORE ANALYSING IT
#
# WHAT THIS SCRIPT IS FOR
#
# Every analysis needs settings: how many genes to keep, how big a
# difference counts as big, how strong a correlation counts as
# strong. Those numbers are usually copied from somebody else's
# project. That is a bad habit, because the right number depends on
# the dataset in front of you.
#
# This script measures the dataset and works the numbers out from
# what it finds. At the end it prints a settings block you can paste
# into simple_pipeline.py, with the reason for each number.
#
# RUN THIS FIRST, BEFORE simple_pipeline.py:
#     python3 simple_explore_first.py
#
# It only READS the data. It changes nothing and trains nothing, so
# it is safe to run as many times as you like.
#
# HONEST LIMIT: not every number can be measured. Some are
# conventions the whole field agreed on (like FDR < 0.05) and some
# can only be decided by running the analysis (like the number of
# patient clusters). This script labels which is which, so you know
# what you can defend with data and what you are just following.
# =============================================================

import gzip
import os

import numpy as np
import pandas as pd

from scipy import stats


# ---- The only two things you must set ----

BASE_FOLDER = "/home/sails/SAILS-Repo/Gene_Expression_Clustering"
DATA_FILE = BASE_FOLDER + "/raw/GSE13159_series_matrix.txt.gz"

# How many probes to sample when we look at the numbers. We do not
# need all 54,000 to work out what the data looks like, and reading
# all of them takes minutes.
HOW_MANY_TO_PEEK = 3000


# A place to collect every recommendation as we go, so we can print
# them all together at the end.
FINDINGS = []


def note(setting_name, value, reason, is_measured=True):
    # Record one recommendation, with WHY.
    #
    # A number without a reason is just as unjustified as a number
    # copied from someone else, so we never store one without the
    # other.

    FINDINGS.append({"setting": setting_name,
                     "value": value,
                     "reason": reason,
                     "measured_from_your_data": is_measured})


# =============================================================
# PART 1 - WHAT IS IN THE FILE AT ALL?
# =============================================================

def read_the_header():
    # Read only the "!" lines at the top. This tells us how many
    # patients there are and what we know about them, without touching
    # the huge number table.

    print("=" * 62)
    print("PART 1: WHAT IS IN THE FILE?")
    print("=" * 62)

    if os.path.exists(DATA_FILE) == False:
        print("ERROR: cannot find", DATA_FILE)
        raise SystemExit(1)

    size_in_gb = os.path.getsize(DATA_FILE) / 1024 / 1024 / 1024
    print("File size:", round(size_in_gb, 2), "GB")

    header_lines = 0
    found_the_table = False

    patient_ids = []
    all_facts = {}
    series_title = ""

    open_file = gzip.open(DATA_FILE, "rt", errors="replace")

    for line in open_file:
        if line.startswith("!series_matrix_table_begin"):
            found_the_table = True
            break

        header_lines = header_lines + 1

        pieces = line.strip().split("\t")
        label = pieces[0]

        if label == "!Series_title":
            if len(pieces) > 1:
                series_title = pieces[1].replace('"', "")

        if label == "!Sample_geo_accession":
            for one_piece in pieces[1:]:
                patient_ids.append(one_piece.replace('"', ""))

        # Every other "!Sample_" line is a fact about the patients.
        # We keep them ALL, because we do not yet know which one holds
        # the diagnosis.
        if label.startswith("!Sample_"):
            field_name = label.replace("!Sample_", "")

            values = []
            inner_label = ""

            for one_piece in pieces[1:]:
                clean_text = one_piece.replace('"', "").strip()

                # Values often look like "leukemia class: CLL"
                if ":" in clean_text:
                    two_parts = clean_text.split(":", 1)
                    inner_label = two_parts[0].strip()
                    values.append(two_parts[1].strip())
                else:
                    values.append(clean_text)

            # Prefer the inner label, it is more descriptive
            if inner_label != "":
                field_name = inner_label

            # The same line name can appear several times, so do not
            # overwrite an earlier one
            while field_name in all_facts:
                field_name = field_name + "_2"

            all_facts[field_name] = values

    open_file.close()

    if found_the_table == False:
        print("ERROR: no '!series_matrix_table_begin' line.")
        print("This is not a GEO series matrix, or the download stopped")
        print("part way through. Check the file size against GEO.")
        raise SystemExit(1)

    print("Title:", series_title[0:70])
    print("Patients:", len(patient_ids))
    print("Header lines to skip:", header_lines + 1)

    return header_lines + 1, patient_ids, all_facts


def count_the_probes(lines_to_skip):
    # Count the rows of the number table without loading them.
    # We read the file line by line and just count, which uses almost
    # no memory however big the file is.

    print("Counting probes (reading the file once)...")

    how_many = 0
    line_number = 0

    open_file = gzip.open(DATA_FILE, "rt", errors="replace")

    for line in open_file:
        line_number = line_number + 1

        # Skip the header and the column-names row
        if line_number <= lines_to_skip + 1:
            continue

        # The end marker
        if line.startswith("!"):
            break

        how_many = how_many + 1

    open_file.close()

    print("Probes:", how_many)
    return how_many


# =============================================================
# PART 2 - WHICH FACT ABOUT THE PATIENTS IS THE DIAGNOSIS?
# =============================================================

def find_the_group_label(all_facts, how_many_patients):
    # WHY THIS COMES BEFORE EVERYTHING ELSE
    #
    # Every statistic in the pipeline compares one group of patients
    # against another. If we pick the wrong field to define the
    # groups, every result afterwards is meaningless - and it will
    # still look perfectly normal, which is worse.
    #
    # HOW WE JUDGE A FIELD. A useful grouping field has:
    #   - more than one value (a field that is the same for everybody
    #     tells us nothing)
    #   - not a different value for every patient (that is an ID, not
    #     a group)
    #   - enough patients per value to do statistics on
    #
    # We score every field and show them ranked, but we do NOT decide
    # silently. You have to look at the names and confirm.

    print("")
    print("=" * 62)
    print("PART 2: WHICH FIELD DEFINES THE PATIENT GROUPS?")
    print("=" * 62)

    rows = []

    for field_name in all_facts:
        values = all_facts[field_name]

        if len(values) != how_many_patients:
            continue

        counts = pd.Series(values).value_counts()
        how_many_different = len(counts)

        # How many groups would be big enough to analyse at 20+?
        big_enough = 0
        for one_value in counts.index:
            if counts[one_value] >= 20:
                big_enough = big_enough + 1

        # Judge what kind of field this is
        if how_many_different == 1:
            verdict = "useless - same for everybody"
        elif how_many_different == how_many_patients:
            verdict = "an ID, not a group"
        elif how_many_different > how_many_patients * 0.5:
            verdict = "too many values to be a group"
        elif big_enough < 2:
            verdict = "groups too small to test"
        elif how_many_different == 2:
            verdict = "possible, but only two groups"
        else:
            verdict = "GOOD CANDIDATE"

        rows.append({"field": field_name,
                     "different_values": how_many_different,
                     "groups_with_20_or_more": big_enough,
                     "smallest_group": int(counts.min()),
                     "largest_group": int(counts.max()),
                     "example": str(counts.index[0])[0:40],
                     "verdict": verdict})

    table = pd.DataFrame(rows)
    table = table.sort_values("groups_with_20_or_more", ascending=False)

    print("")
    print(table.to_string(index=False))

    good_ones = table[table["verdict"] == "GOOD CANDIDATE"]

    print("")

    if len(good_ones) == 0:
        print("No field looks like a usable grouping. Look at the table")
        print("above yourself - you may need to combine two fields, or")
        print("this dataset may not have diagnoses attached.")
        note("LABEL_FIELD", "UNKNOWN - look at the table above",
             "no field had 2+ groups of 20 or more patients")
        return None, table

    best_field = good_ones["field"].iloc[0]

    print("Best candidate:", best_field)
    print("  ", int(good_ones["different_values"].iloc[0]), "different values,",
          int(good_ones["groups_with_20_or_more"].iloc[0]),
          "of them big enough to test")

    # ---- THE TRAP WORTH NAMING OUT LOUD ----
    if len(good_ones) > 1:
        print("")
        print("CHECK THIS YOURSELF. More than one field could work:")
        for one_field in good_ones["field"]:
            print("   ", one_field)
        print("Fields like 'tissue' or 'sample type' describe WHERE the")
        print("sample came from, not what the patient has. Picking one of")
        print("those by mistake makes the whole analysis compare bone")
        print("marrow against blood. Read the names and be sure.")

    note("LABEL_FIELD", best_field,
         "the field with the most groups of 20+ patients - CONFIRM BY EYE")

    return best_field, table


def decide_minimum_group_size(all_facts, best_field):
    # How many patients does a group need before we trust its
    # statistics?
    #
    # There is no magic number, but there is a real trade-off, and we
    # can show it with the actual data: raise the threshold and you
    # get more reliable tests but fewer diseases; lower it and you
    # test more diseases less reliably.
    #
    # A Welch t-test gets shaky below roughly 10 per group. Below 20
    # a single unusual patient can drive the result.

    print("")
    print("How many patients should a group need?")

    if best_field is None:
        note("MIN_GROUP_SIZE", 20, "convention - could not check your groups",
             is_measured=False)
        return 20

    counts = pd.Series(all_facts[best_field]).value_counts()

    print("")
    print("  threshold | groups kept | patients kept")

    options = [10, 15, 20, 25, 30, 50]
    results = {}

    for one_threshold in options:
        groups_kept = 0
        patients_kept = 0

        for one_value in counts.index:
            if counts[one_value] >= one_threshold:
                groups_kept = groups_kept + 1
                patients_kept = patients_kept + int(counts[one_value])

        results[one_threshold] = (groups_kept, patients_kept)
        print("     ", str(one_threshold).rjust(5), "|",
              str(groups_kept).rjust(11), "|", str(patients_kept).rjust(13))

    total_patients = int(counts.sum())

    # ---- Pick the highest threshold that still keeps almost everyone ----
    # "Almost everyone" = 95% of patients. Dropping a handful of
    # patients from tiny groups costs little; dropping a fifth of the
    # cohort costs a lot.
    chosen = 10

    for one_threshold in options:
        groups_kept, patients_kept = results[one_threshold]

        if groups_kept < 2:
            continue

        if patients_kept >= 0.95 * total_patients:
            chosen = one_threshold

    groups_kept, patients_kept = results[chosen]
    lost = total_patients - patients_kept

    print("")
    print("Suggested:", chosen)
    print("  keeps", groups_kept, "groups and", patients_kept, "patients")
    print("  loses", lost, "patients from groups too small to test")

    note("MIN_GROUP_SIZE", chosen,
         "highest threshold that still keeps 95% of patients (" +
         str(groups_kept) + " groups, " + str(lost) + " patients dropped)")

    return chosen


# =============================================================
# PART 3 - WHAT DO THE NUMBERS ACTUALLY LOOK LIKE?
# =============================================================

def peek_at_the_numbers(lines_to_skip, how_many_probes):
    # Load a SAMPLE of the number table.
    #
    # IMPORTANT DETAIL: we take probes spread evenly through the file,
    # not the first few thousand. On an Affymetrix chip the first rows
    # are AFFX control probes, whose values are not typical of real
    # genes. Judging the whole dataset from them would be wrong.

    print("")
    print("=" * 62)
    print("PART 3: WHAT DO THE NUMBERS LOOK LIKE?")
    print("=" * 62)

    print("Sampling", HOW_MANY_TO_PEEK, "probes spread through the file...")

    # Which rows do we want? Every Nth one.
    rows_to_skip = list(range(lines_to_skip))

    if how_many_probes > HOW_MANY_TO_PEEK:
        step = how_many_probes // HOW_MANY_TO_PEEK

        first_data_row = lines_to_skip + 1
        for row_number in range(how_many_probes):
            if row_number % step != 0:
                rows_to_skip.append(first_data_row + row_number)

    data = pd.read_csv(DATA_FILE, sep="\t", index_col=0,
                       skiprows=rows_to_skip, comment="!",
                       low_memory=False)

    data = data.apply(pd.to_numeric, errors="coerce")

    # Drop the control probes before looking at the values
    keep_rows = []
    how_many_control = 0

    for probe_name in data.index:
        if str(probe_name).upper().startswith("AFFX"):
            keep_rows.append(False)
            how_many_control = how_many_control + 1
        else:
            keep_rows.append(True)

    if how_many_control > 0:
        print("Ignoring", how_many_control, "AFFX control probes")
        print("  (their values are not typical of real genes)")
        data = data[keep_rows]

    print("Working with", data.shape[0], "probes x", data.shape[1],
          "patients")

    return data


def work_out_the_scale(data):
    # WHAT SCALE ARE THESE NUMBERS ON? This one decision changes how
    # you report every effect size in your write-up.
    #
    # Four possibilities:
    #
    #   0 to 1        already rescaled by whoever shared it.
    #                 Differences are NOT fold changes and cannot be
    #                 turned back into fold changes.
    #   log2          about 0 to 20. Differences ARE log2 fold changes,
    #                 so "log2FC >= 1" is meaningful.
    #   linear        thousands. Needs a log2 transform first, or the
    #                 brightest probes dominate every distance.
    #   something else  stop and look at it by hand.

    print("")
    print("Working out the scale...")

    values = data.values

    smallest = float(np.nanmin(values))
    biggest = float(np.nanmax(values))
    middle = float(np.nanmedian(values))

    print("  smallest:", round(smallest, 4))
    print("  largest:", round(biggest, 4))
    print("  middle:", round(middle, 4))

    if biggest <= 1.01 and smallest >= -0.01:
        scale = "already_0_to_1"
        needs_log2 = False
        can_use_fold_change = False
        print("")
        print("SCALE: already rescaled to 0-1.")
        print("  Do NOT apply log2 - it is already been squashed.")
        print("  Do NOT call differences 'fold changes'. You cannot")
        print("  recover fold changes from this file at all; you would")
        print("  need the original raw CEL files.")

    elif biggest <= 30.0:
        scale = "log2"
        needs_log2 = False
        can_use_fold_change = True
        print("")
        print("SCALE: log2 expression.")
        print("  Do NOT apply log2 again.")
        print("  Differences between group averages ARE log2 fold")
        print("  changes, so the usual 'log2FC >= 1' rule works here.")

    elif biggest > 50.0:
        scale = "linear"
        needs_log2 = True
        can_use_fold_change = True
        print("")
        print("SCALE: linear intensities.")
        print("  You MUST apply log2(x + 1) before any distance,")
        print("  correlation or t-test. Without it the brightest probes")
        print("  swamp everything else.")
        print("  After the log2, differences are log2 fold changes.")

    else:
        scale = "unclear"
        needs_log2 = False
        can_use_fold_change = False
        print("")
        print("SCALE: unclear - the largest value sits between the log2")
        print("  and linear ranges. Do not guess. Plot a histogram of")
        print("  the values and decide by eye before going further.")

    note("value_scale", scale, "largest value seen = " + str(round(biggest, 3)))
    note("apply_log2_transform", needs_log2,
         "follows from the scale above")
    note("effect_sizes_are_fold_changes", can_use_fold_change,
         "follows from the scale above")

    return scale, can_use_fold_change


def check_how_it_was_normalised(data):
    # HAS SOMEONE ALREADY NORMALISED THIS, AND HOW?
    #
    # Two different things both squash data, and they have very
    # different consequences:
    #
    #   quantile normalised   every patient forced to the same
    #                         distribution shape. Good: patients are
    #                         directly comparable.
    #   min-max per patient   every patient stretched to exactly 0-1
    #                         individually. Awkward: a value now
    #                         depends on that patient's own brightest
    #                         and dimmest probe, so the same number
    #                         means different things in different
    #                         patients.
    #
    # We tell them apart by looking at each patient's own minimum and
    # maximum. If EVERY patient runs exactly 0 to 1, it is min-max.

    print("")
    print("Checking whether it has already been normalised...")

    values = data.values

    each_patient_min = np.nanmin(values, axis=0)
    each_patient_max = np.nanmax(values, axis=0)

    spread_of_maxes = float(np.nanmax(each_patient_max) -
                           np.nanmin(each_patient_max))

    print("  patient maximums vary by:", round(spread_of_maxes, 6))

    all_start_at_zero = bool(np.allclose(each_patient_min,
                                        np.nanmin(each_patient_min),
                                        atol=1e-6))
    all_end_at_one = bool(np.allclose(each_patient_max,
                                      np.nanmax(each_patient_max),
                                      atol=1e-6))

    if all_start_at_zero and all_end_at_one and spread_of_maxes < 0.01:
        print("")
        print("Every patient spans exactly the same range.")
        print("This is MIN-MAX SCALING PER PATIENT.")
        print("  Say so in your Methods section. It means each value")
        print("  depends on that patient's own extremes, so effect")
        print("  sizes are bounded and small.")
        note("already_normalised", "min_max_per_sample",
             "every patient has an identical minimum and maximum")
        return "min_max_per_sample"

    if spread_of_maxes < 0.01 * max(abs(float(np.nanmax(values))), 1e-9):
        print("")
        print("All patients share almost the same maximum, so the arrays")
        print("look QUANTILE NORMALISED already.")
        print("  Running quantile normalisation again is harmless but")
        print("  pointless. Keep it in for safety and say you did.")
        note("already_normalised", "quantile",
             "all patients share nearly the same maximum")
        return "quantile"

    print("")
    print("Patients differ in range, so this looks RAW.")
    print("  You definitely need quantile normalisation, especially if")
    print("  the samples came from more than one laboratory.")
    note("already_normalised", "no",
         "patient ranges differ by " + str(round(spread_of_maxes, 3)))
    return "raw"


def check_missing_values(data):
    # How much is missing, and does that change what we should do?

    print("")
    print("Checking for missing values...")

    how_many = int(data.isna().sum().sum())
    total = data.shape[0] * data.shape[1]
    percent = 100.0 * how_many / total

    print("  missing:", how_many, "of", total,
          "(", round(percent, 4), "% )")

    if percent == 0:
        print("  Nothing missing. No filling needed.")
        note("how_to_fill_gaps", "nothing_needed", "0% missing")

    elif percent < 1.0:
        print("  Very few. Filling each gap with that probe's average")
        print("  is safe and will not change any conclusion.")
        note("how_to_fill_gaps", "probe_average",
             str(round(percent, 4)) + "% missing - too few to matter")

    elif percent < 5.0:
        print("  Noticeable. The probe average is still reasonable, but")
        print("  check whether the gaps are concentrated in a few")
        print("  patients - if so, consider dropping those patients.")
        note("how_to_fill_gaps", "probe_average_but_check",
             str(round(percent, 2)) + "% missing - check if clustered")

    else:
        print("  A LOT is missing. Filling this much invents data.")
        print("  Find out why before continuing: a failed batch? a")
        print("  merge of two incompatible platforms?")
        note("how_to_fill_gaps", "INVESTIGATE_FIRST",
             str(round(percent, 2)) + "% missing - too much to fill blindly")

    return percent


def check_chip_quality(data):
    # Are any patients wildly unlike the rest?
    #
    # We also use the spread of the correlations to choose the outlier
    # threshold, rather than copying "5 standard deviations" from
    # somewhere. If the correlations are tightly bunched, 5 is very
    # strict; if they are spread out, 5 may catch nothing.

    print("")
    print("Checking how similar the patients are to each other...")

    values = data.values.astype("float64")

    # Fill gaps temporarily, or correlation returns nothing
    for i in range(values.shape[0]):
        one_row = values[i]
        if np.isnan(one_row).any():
            values[i] = np.nan_to_num(one_row, nan=np.nanmean(one_row))

    grid = np.corrcoef(values.T)
    np.fill_diagonal(grid, np.nan)

    average_similarity = np.nanmean(grid, axis=1)

    typical = float(np.nanmean(average_similarity))
    spread = float(np.nanstd(average_similarity))
    worst = float(np.nanmin(average_similarity))

    print("  typical similarity:", round(typical, 3))
    print("  spread:", round(spread, 4))
    print("  worst patient:", round(worst, 3))

    if spread > 0:
        worst_in_spreads = (worst - typical) / spread
    else:
        worst_in_spreads = 0.0

    print("  the worst patient is", round(abs(worst_in_spreads), 1),
          "spreads below typical")

    # ---- Choose a threshold that flags a sensible number ----
    # We want a line that catches genuine oddities without throwing
    # away rare diseases that look different for real reasons.
    print("")
    print("  threshold | patients flagged")

    chosen = 5.0

    for one_threshold in [3.0, 4.0, 5.0, 6.0]:
        flagged = 0
        for one_value in average_similarity:
            if spread > 0 and (one_value - typical) / spread < -one_threshold:
                flagged = flagged + 1

        percent_flagged = 100.0 * flagged / len(average_similarity)
        print("       ", str(one_threshold).rjust(5), "|", str(flagged).rjust(7),
              "(", round(percent_flagged, 2), "% )")

        # Prefer the loosest threshold that still flags under 2%.
        # More than that and we are probably flagging biology.
        if percent_flagged <= 2.0:
            chosen = one_threshold
            break

    print("")
    print("Suggested:", chosen, "spreads below typical")
    print("  FLAG them, do not delete them. A rare leukaemia genuinely")
    print("  looks different from everything else, and deleting it")
    print("  would remove exactly what you are looking for.")

    if typical < 0.3:
        print("")
        print("  NOTE: typical similarity is low. Either these patients")
        print("  are very diverse, or something is wrong with the")
        print("  normalisation. Worth a look before continuing.")

    note("OUTLIER_THRESHOLD", chosen,
         "loosest cutoff flagging under 2% of patients")

    return average_similarity


# =============================================================
# PART 4 - HOW MANY GENES SHOULD WE KEEP?
# =============================================================

def choose_the_probe_filters(data, how_many_probes):
    # TWO FILTERS, BOTH DERIVED FROM THIS DATASET.
    #
    # Filter 1: throw away probes that are switched off in everybody.
    #   The usual rule is "keep probes whose 95th percentile is above
    #   the middle of the whole matrix". That is data-derived already,
    #   so we just report what the number comes out as here.
    #
    # Filter 2: keep the N most variable probes. N is normally just
    #   guessed at 5000. We can do better: look at where the
    #   variability curve flattens out. Past that point the probes
    #   being added are barely more variable than the ones already in,
    #   so they add noise rather than signal.

    print("")
    print("=" * 62)
    print("PART 4: HOW MANY GENES TO KEEP?")
    print("=" * 62)

    values = data.values
    values = np.nan_to_num(values, nan=float(np.nanmedian(values)))

    # ---- Filter 1 ----
    high_value = np.percentile(values, 95, axis=1)
    threshold = float(np.median(values))

    how_many_pass = 0
    for one_value in high_value:
        if one_value > threshold:
            how_many_pass = how_many_pass + 1

    fraction_passing = how_many_pass / len(high_value)

    print("")
    print("Filter 1: is the probe switched on anywhere?")
    print("  threshold (middle of the matrix):", round(threshold, 4))
    print("  probes passing:", round(100 * fraction_passing, 1), "%")
    print("  on the full file that is about",
          int(fraction_passing * how_many_probes), "probes")

    if fraction_passing > 0.98:
        print("")
        print("  NOTE: almost every probe passes, so this filter is not")
        print("  actually removing anything. That happens when all the")
        print("  probes sit in a similar range. Harmless, but do not")
        print("  claim in your write-up that you filtered out silent")
        print("  probes when you effectively did not.")
    elif fraction_passing < 0.20:
        print("")
        print("  NOTE: this filter removes over 80% of probes. Check that")
        print("  is what you want before continuing - it is a lot to")
        print("  throw away on one rule.")

    note("expression_filter_threshold", round(threshold, 4),
         "the middle value of your matrix - " +
         str(round(100 * fraction_passing, 1)) + "% of probes pass")

    # ---- Filter 2: find the elbow ----
    middle_of_each = np.median(values, axis=1, keepdims=True)
    mad = np.median(np.abs(values - middle_of_each), axis=1)

    # Only rank the probes that passed filter 1
    passed = mad[high_value > threshold]
    sorted_mad = np.sort(passed)[::-1]

    print("")
    print("Filter 2: how many of the most variable probes to keep?")
    print("")
    print("  Variability at each position in the ranking:")

    # Scale positions in our sample up to the full file
    scale_up = how_many_probes / len(high_value)

    for fraction in [0.02, 0.05, 0.10, 0.20, 0.40, 0.80]:
        position = int(fraction * len(sorted_mad))
        if position >= len(sorted_mad):
            continue

        equivalent = int(position * scale_up)
        print("    top", str(round(100 * fraction, 0)) + "% (about",
              equivalent, "probes ) - variability",
              round(float(sorted_mad[position]), 4))

    # THE ELBOW. We find where the curve stops falling steeply.
    # Concretely: the point where variability has dropped to a fifth
    # of the most variable probe. Past there the probes are nearly
    # flat and contribute mostly noise.
    top_variability = float(sorted_mad[0])
    elbow_level = top_variability * 0.2

    elbow_position = len(sorted_mad)
    found_an_elbow = False

    for i in range(len(sorted_mad)):
        if sorted_mad[i] < elbow_level:
            elbow_position = i
            found_an_elbow = True
            break

    suggested = int(elbow_position * scale_up)

    # Keep the answer sane. Too few and clustering has nothing to work
    # with; too many and noise takes over.
    if suggested < 1000:
        suggested = 1000
        reason_extra = " (raised to the 1000 minimum)"
    elif suggested > 10000:
        suggested = 10000
        reason_extra = " (capped at 10000)"
    else:
        reason_extra = ""

    # Round to something readable
    suggested = int(round(suggested / 500.0) * 500)

    print("")

    if found_an_elbow:
        print("The variability curve drops off at about position", suggested)
        print("Suggested TOP_PROBES:", suggested)
        print("  Past this point probes are less than a fifth as variable")
        print("  as the most variable one, so they add noise not signal.")

        note("TOP_PROBES", suggested,
             "where variability drops below 20% of the maximum" + reason_extra)
    else:
        # BE HONEST. If variability declines gently all the way down,
        # there is no natural place to cut, and pretending we found one
        # would be inventing a justification.
        print("NO NATURAL CUT-OFF FOUND.")
        print("  Variability declines gently all the way down - it never")
        print("  drops to a fifth of the maximum. So there is no elbow in")
        print("  the curve and no number here is 'the right one'.")
        print("")
        print("  What that means in practice: your choice of how many")
        print("  probes to keep is a trade-off you make, not a fact the")
        print("  data gives you. Look at the table above and pick a")
        print("  number, then say in your Methods that you chose it and")
        print("  why - do not claim it came from the data.")
        print("")
        print("  A reasonable default is the top 25%, about",
              int(0.25 * how_many_probes), "probes.")

        suggested = int(round(0.25 * how_many_probes / 500.0) * 500)
        if suggested < 1000:
            suggested = 1000
        if suggested > 10000:
            suggested = 10000

        note("TOP_PROBES", suggested,
             "NO elbow in the variability curve - this is the top 25%, "
             "a choice not a measurement", is_measured=False)

    return suggested


# =============================================================
# PART 5 - HOW BIG IS A "BIG" DIFFERENCE HERE?
# =============================================================

def measure_real_effect_sizes(data, all_facts, best_field,
                              can_use_fold_change):
    # THE MOST IMPORTANT MEASUREMENT IN THIS SCRIPT.
    #
    # Every guide says "keep genes with a fold change above 2". On a
    # 0-to-1 rescaled dataset that rule selects exactly nothing, and
    # you get an empty results table with no explanation.
    #
    # So instead of trusting a rule, we actually run a few
    # comparisons on the sampled probes and look at the distribution
    # of differences we get. Then we pick a cutoff from that.

    print("")
    print("=" * 62)
    print("PART 5: HOW BIG IS A 'BIG' DIFFERENCE IN THIS DATASET?")
    print("=" * 62)

    if best_field is None:
        print("No grouping field, so this cannot be measured.")
        note("EFFECT_CUTOFF", "unknown", "no grouping field found",
             is_measured=False)
        return None

    labels = np.array(all_facts[best_field])

    if len(labels) != data.shape[1]:
        print("The label list and the number table disagree in length.")
        print("  labels:", len(labels), " patients in table:", data.shape[1])
        print("Skipping this part.")
        return None

    values = data.values.astype("float64")
    values = np.nan_to_num(values, nan=float(np.nanmedian(values)))

    counts = pd.Series(labels).value_counts()
    testable = []
    for one_value in counts.index:
        if counts[one_value] >= 20:
            testable.append(one_value)

    print("Running one-against-the-rest on", len(testable), "groups,")
    print("using the", data.shape[0], "sampled probes...")

    all_differences = []
    best_per_gene = {}

    for one_group in testable:
        in_group = (labels == one_group)

        inside = values[:, in_group]
        outside = values[:, in_group == False]

        difference = inside.mean(axis=1) - outside.mean(axis=1)

        for i in range(len(difference)):
            all_differences.append(difference[i])

            # Track each gene's best and second best, for the
            # specificity floor further down
            gene_name = str(data.index[i])

            if gene_name not in best_per_gene:
                best_per_gene[gene_name] = []
            best_per_gene[gene_name].append(abs(difference[i]))

    all_differences = np.abs(np.array(all_differences))

    print("")
    print("The differences we actually saw:")
    print("  middle (median):", round(float(np.median(all_differences)), 4))
    print("  90th percentile:", round(float(np.percentile(all_differences, 90)), 4))
    print("  99th percentile:", round(float(np.percentile(all_differences, 99)), 4))
    print("  largest:", round(float(all_differences.max()), 4))

    # ---- Would the textbook rule work? ----
    print("")
    if can_use_fold_change:
        how_many_above_one = int((all_differences >= 1.0).sum())
        percent = 100.0 * how_many_above_one / len(all_differences)

        print("The usual 'log2FC >= 1' rule would select",
              round(percent, 3), "% of comparisons.")

        if percent < 0.01:
            print("  That is almost nothing. Use a percentile instead.")
            use_percentile = True
        else:
            print("  That works on this data. You may use log2FC >= 1")
            print("  and report it as a fold change, which reviewers")
            print("  will find familiar.")
            use_percentile = False
    else:
        print("The 'log2FC >= 1' rule CANNOT be used here: the largest")
        print("difference in the whole dataset is only",
              round(float(all_differences.max()), 4) , ".")
        print("  A cutoff of 1.0 would return an empty table and you")
        print("  would have no idea why. Use a percentile.")
        use_percentile = True

    if use_percentile:
        cutoff = float(np.percentile(all_differences, 99))
        print("")
        print("Suggested cutoff:", round(cutoff, 4),
              "(the top 1% of what we saw)")
        print("  This says 'unusually large FOR THIS DATASET', which is")
        print("  an honest claim. Say in your Methods that the cutoff is")
        print("  a percentile of the observed distribution, and why.")

        note("EFFECT_CUTOFF_PERCENTILE", 99,
             "top 1% of observed differences = about " + str(round(cutoff, 4)))
        note("EFFECT_CUTOFF_VALUE", round(cutoff, 4),
             "measured from " + str(len(all_differences)) + " comparisons")
    else:
        cutoff = 1.0
        note("EFFECT_CUTOFF_VALUE", 1.0,
             "the conventional log2FC >= 1 works on this scale")

    # ---- The specificity floor ----
    #
    # Specificity divides a gene's own effect by the runner-up's. When
    # the runner-up is near zero, that division explodes, so we floor
    # the bottom. But the floor must be SMALL compared to real effect
    # sizes, or the floor becomes the answer.
    print("")
    print("Choosing the specificity floor...")

    runner_ups = []
    for gene_name in best_per_gene:
        effects = sorted(best_per_gene[gene_name], reverse=True)
        if len(effects) >= 2:
            runner_ups.append(effects[1])

    if len(runner_ups) > 0:
        runner_ups = np.array(runner_ups)
        typical_runner_up = float(np.median(runner_ups))

        print("  typical runner-up effect:", round(typical_runner_up, 4))

        # The floor should sit well below a typical runner-up, or it
        # replaces the real denominator most of the time. A tenth of
        # the typical runner-up is a defensible choice.
        floor = typical_runner_up / 10.0
        floor = float(round(floor, 4))

        if floor < 0.001:
            floor = 0.001

        # How often would the usual 0.1 have kicked in?
        how_often_default = 100.0 * float((runner_ups < 0.1).mean())
        how_often_ours = 100.0 * float((runner_ups < floor).mean())

        print("  a floor of 0.1 would replace the real value",
              round(how_often_default), "% of the time")
        print("  a floor of", floor, "would do so",
              round(how_often_ours), "% of the time")

        if how_often_default > 25:
            print("")
            print("  So the common default of 0.1 is TOO BIG for this")
            print("  dataset. It would be doing the work instead of the")
            print("  data, and 'specificity' would stop meaning")
            print("  specificity.")

        print("")
        print("Suggested SPECIFICITY_FLOOR:", floor)

        note("SPECIFICITY_FLOOR", floor,
             "a tenth of the typical runner-up effect (" +
             str(round(typical_runner_up, 4)) + ")")

    return cutoff


# =============================================================
# PART 6 - HOW STRONG IS A "STRONG" CORRELATION HERE?
# =============================================================

def choose_the_network_cutoff(data):
    # The gene network draws a line between two genes when they move
    # together strongly. "Strongly" is normally set to 0.7 because
    # that is what everybody uses.
    #
    # But 0.7 might catch a million pairs in one dataset and none in
    # another. What we actually want is a cutoff that leaves a network
    # you can look at: sparse enough to read, dense enough to have
    # structure. So we measure the distribution and pick from it.

    print("")
    print("=" * 62)
    print("PART 6: HOW STRONG IS A 'STRONG' CORRELATION HERE?")
    print("=" * 62)

    values = data.values.astype("float64")
    values = np.nan_to_num(values, nan=float(np.nanmedian(values)))

    # Use the most variable few hundred, which is what the network
    # step will use
    middle_of_each = np.median(values, axis=1, keepdims=True)
    mad = np.median(np.abs(values - middle_of_each), axis=1)

    how_many = min(300, values.shape[0])
    top = np.argsort(mad)[::-1][0:how_many]

    grid = np.corrcoef(values[top])
    grid = np.nan_to_num(grid, nan=0.0)

    top_half = np.triu_indices(how_many, k=1)
    strengths = np.abs(grid[top_half])

    total_pairs = len(strengths)

    print("Looking at", how_many, "genes =", total_pairs, "possible pairs")
    print("")
    print("  cutoff | pairs kept | how dense the network would be")

    chosen = None

    for one_cutoff in [0.5, 0.6, 0.7, 0.8, 0.9]:
        kept = int((strengths >= one_cutoff).sum())
        density = 100.0 * kept / total_pairs

        print("    ", one_cutoff, "|", str(kept).rjust(10), "|",
              round(density, 2), "%")

        # A readable network has roughly 1-5% of possible lines drawn.
        # Below that it is dust; above it is a hairball.
        if chosen is None and density <= 5.0:
            chosen = one_cutoff

    if chosen is None:
        chosen = 0.9

    kept = int((strengths >= chosen).sum())

    print("")
    print("Suggested NETWORK_STRENGTH:", chosen)
    print("  gives about", kept, "lines, which is readable")

    # A sudden collapse between two neighbouring cutoffs means the
    # answer is very sensitive to where you put the line, which is
    # worth knowing before you quote a network as a finding.
    at_point_six = int((strengths >= 0.6).sum())
    at_point_seven = int((strengths >= 0.7).sum())

    if at_point_six > 10 * max(at_point_seven, 1):
        print("")
        print("  CAREFUL: the number of lines collapses sharply between")
        print("  0.6 and 0.7 (" + str(at_point_six), "down to",
              str(at_point_seven) + ").")
        print("  So your network depends heavily on exactly where you put")
        print("  the line. Report the cutoff you used, and ideally show")
        print("  that your conclusion holds at a neighbouring value too.")

    if kept < 20:
        print("")
        print("  WARNING: even the loosest cutoff gives very few lines.")
        print("  Either these genes really are independent, or your")
        print("  sample of probes is too small to see the structure.")
        print("  The full pipeline uses all the genes, so it will find")
        print("  more than this preview does.")

    note("NETWORK_STRENGTH", chosen,
         "loosest cutoff keeping the network under 5% dense (about " +
         str(kept) + " lines)")

    return chosen


# =============================================================
# PART 7 - WHAT WILL THIS COST TO RUN?
# =============================================================

def estimate_the_cost(how_many_probes, how_many_patients):
    # Before you start a run that takes hours, know what it will cost.
    # These are rough estimates from the sizes, not measurements, and
    # they are labelled as such.

    print("")
    print("=" * 62)
    print("PART 7: WHAT WILL THIS COST TO RUN?")
    print("=" * 62)

    # The full matrix as float32
    matrix_gb = how_many_probes * how_many_patients * 4 / 1024 / 1024 / 1024

    print("The full matrix in memory:", round(matrix_gb, 2), "GB")
    print("  You need roughly three times that free, because several")
    print("  steps hold a copy while they work:",
          round(matrix_gb * 3, 1), "GB")

    # The patient-by-patient grids used by clustering
    pair_grid_mb = how_many_patients * how_many_patients * 4 / 1024 / 1024
    print("A patient-by-patient grid:", round(pair_grid_mb, 1), "MB each")

    print("")

    if how_many_patients > 3000:
        repeats = 10
        print("With", how_many_patients, "patients the stability test is")
        print("expensive. Suggest CONSENSUS_REPEATS = 10 to start.")
    else:
        repeats = 20
        print("Suggest CONSENSUS_REPEATS = 20 (the usual choice).")

    note("CONSENSUS_REPEATS", repeats,
         "based on " + str(how_many_patients) + " patients")

    if how_many_patients > 1500:
        shap_patients = 300
    else:
        shap_patients = min(300, int(0.3 * how_many_patients))

    print("Suggest SHAP_PATIENTS =", shap_patients)
    print("  SHAP is the slowest single thing in the pipeline. Its work")
    print("  grows with patients x genes x diseases, so explaining a")
    print("  few hundred patients is enough to rank genes reliably.")

    note("SHAP_PATIENTS", shap_patients,
         "enough to rank genes without an unreasonable wait")

    if matrix_gb > 4:
        print("")
        print("WARNING: this matrix is large. If the run dies with a")
        print("MemoryError, process the data in chunks or use a machine")
        print("with more RAM.")


# =============================================================
# PART 8 - WRITE THE SETTINGS OUT
# =============================================================

def write_the_settings():
    # Print everything we worked out, and save it as a file you can
    # paste into simple_pipeline.py.
    #
    # Every line carries its reason. That is the whole point: when
    # somebody asks "why 4300 probes?", the answer is in the file, not
    # in somebody's memory.

    print("")
    print("=" * 62)
    print("WHAT I FOUND, AND WHY")
    print("=" * 62)

    table = pd.DataFrame(FINDINGS)

    print("")
    for row in table.itertuples(index=False):
        if row.measured_from_your_data:
            tag = "[measured]"
        else:
            tag = "[convention]"

        print(tag, row.setting, "=", row.value)
        print("      because:", row.reason)

    table.to_csv(BASE_FOLDER + "/results/step0_settings_found.csv",
                 index=False)

    # ---- Build the paste-in block ----
    lines = []
    lines.append("# ============================================")
    lines.append("# SETTINGS WORKED OUT FROM YOUR OWN DATA")
    lines.append("# by simple_explore_first.py")
    lines.append("#")
    lines.append("# Paste these over the matching lines in")
    lines.append("# simple_pipeline.py. Each one has its reason next")
    lines.append("# to it, so you can defend it or change it.")
    lines.append("# ============================================")
    lines.append("")

    # Only the ones the pipeline actually has a setting for
    pipeline_settings = ["MIN_GROUP_SIZE", "TOP_PROBES",
                         "EFFECT_CUTOFF_PERCENTILE", "SPECIFICITY_FLOOR",
                         "NETWORK_STRENGTH", "CONSENSUS_REPEATS",
                         "SHAP_PATIENTS", "OUTLIER_THRESHOLD"]

    for row in table.itertuples(index=False):
        if row.setting in pipeline_settings:
            lines.append("# " + row.reason)
            lines.append(str(row.setting) + " = " + str(row.value))
            lines.append("")

    lines.append("# ---- Things you must decide yourself ----")
    lines.append("#")

    for row in table.itertuples(index=False):
        if row.setting not in pipeline_settings:
            lines.append("# " + str(row.setting) + " = " + str(row.value))
            lines.append("#     " + row.reason)

    lines.append("#")
    lines.append("# ---- Things this script CANNOT decide ----")
    lines.append("#")
    lines.append("# FDR_CUTOFF = 0.05")
    lines.append("#     A convention, not a measurement. The whole field")
    lines.append("#     uses 0.05. Keep it unless you have a reason.")
    lines.append("#")
    lines.append("# FORCE_K / number of clusters")
    lines.append("#     Cannot be known before the analysis runs. The")
    lines.append("#     pipeline measures stability and tightness and")
    lines.append("#     tells you what it found. Leave FORCE_K = 0.")
    lines.append("#")
    lines.append("# SPECIFICITY_CUTOFF = 1.5")
    lines.append("#     How one-disease-only a marker must be. This is a")
    lines.append("#     judgement about what you want, not a fact about")
    lines.append("#     the data. Raise it for stricter markers.")
    lines.append("#")
    lines.append("# MAIN_MODEL")
    lines.append('#     Leave it on "best". The pipeline compares the')
    lines.append("#     models and picks using their scores.")

    settings_text = "\n".join(lines)

    file_path = BASE_FOLDER + "/results/step0_my_settings.txt"
    open(file_path, "w", encoding="utf-8").write(settings_text)

    print("")
    print("=" * 62)
    print("Saved to:")
    print("  " + file_path)
    print("  " + BASE_FOLDER + "/results/step0_settings_found.csv")
    print("")
    print("Open the .txt file, read the reasons, and paste the settings")
    print("into simple_pipeline.py. Then run the pipeline.")
    print("=" * 62)


def mymain():
    print("=" * 62)
    print("LOOKING AT THE DATA BEFORE ANALYSING IT")
    print("=" * 62)
    print("")
    print("This reads your dataset and works out the settings from")
    print("what it finds, instead of copying somebody else's numbers.")
    print("Nothing is changed and nothing is trained.")
    print("")

    os.makedirs(BASE_FOLDER + "/results", exist_ok=True)

    # Part 1: the file
    lines_to_skip, patient_ids, all_facts = read_the_header()
    how_many_probes = count_the_probes(lines_to_skip)

    # Part 2: the groups
    best_field, field_table = find_the_group_label(all_facts,
                                                  len(patient_ids))
    field_table.to_csv(BASE_FOLDER + "/results/step0_metadata_fields.csv",
                       index=False)
    decide_minimum_group_size(all_facts, best_field)

    # Part 3: the numbers
    data = peek_at_the_numbers(lines_to_skip, how_many_probes)
    scale, can_use_fold_change = work_out_the_scale(data)
    check_how_it_was_normalised(data)
    check_missing_values(data)
    check_chip_quality(data)

    # Part 4: how many genes
    choose_the_probe_filters(data, how_many_probes)

    # Part 5: effect sizes
    measure_real_effect_sizes(data, all_facts, best_field,
                              can_use_fold_change)

    # Part 6: correlation
    choose_the_network_cutoff(data)

    # Part 7: cost
    estimate_the_cost(how_many_probes, len(patient_ids))

    # Part 8: write it down
    write_the_settings()


if __name__ == "__main__":
    mymain()
