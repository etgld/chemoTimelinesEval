import argparse
import json
import logging
import os
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Union, cast

import dateutil.parser

from docker_output_to_timeline import TimelineDict

VERSION = "v20240305"

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
"""
    Set logger level as INFO for per patient-level evaluation results. 
    Set DEBUG for more detailed results.
"""


parser = argparse.ArgumentParser(
    description="Evaluate predicted output against gold annotations"
)

parser.add_argument(
    "--patient_level_metrics",
    required=False,
    help="Output metrics by patient ID",
    action="store_true",
)
parser.add_argument(
    "--debug",
    required=False,
    help="Output false positives and negatives organized by patient ID",
    action="store_true",
)
parser.add_argument("--gold_path", required=True, help="A gold annotation json file")
parser.add_argument("--pred_path", required=True, help="A predicted output json file")
parser.add_argument(
    "--all_id_path",
    required=True,
    help="Path to file with list of ids, delimited by new line characters",
)
parser.add_argument(
    "--gold_id_path",
    required=False,
    help="(Only for test evaluation) Path to file with list of gold annotated ids, delimited by new line characters",
)

parser.add_argument(
    "--strict", action="store_true", help="do strict eval", default=False
)
parser.add_argument(
    "--relaxed_to",
    help="Type 'year' to only evaluate year, 'month' to evaluate year and month, "
    "or 'day' to evaluate year-month-day",
    choices=["day", "month", "year"],
)
# ideally will standardize this to the tuple case
TimelineTuple = Union[Tuple[str, str, str], List[str]]
TimelineTuples = List[TimelineTuple]
# from "sensitivity and specificity"
# https://en.wikipedia.org/wiki/Sensitivity_and_specificity
SSMatrix = Tuple[Dict[str, int], Dict[str, int], Dict[str, int]]
MetricMatrix = Tuple[Dict[str, float], Dict[str, float], Dict[str, float]]
DebugDict = Dict[str, Dict[str, TimelineTuples]]


class Chemo:
    def __init__(
        self,
        text: str,
        first_start: Optional[datetime] = None,
        last_end: Optional[datetime] = None,
        cui: Optional[str] = None,
    ) -> None:
        self.text: str = text
        self.first_start: Optional[datetime] = first_start
        self.last_end: Optional[datetime] = last_end
        self.cui: Optional[str] = cui

    def __str__(self) -> str:
        return "\t".join(
            [self.text if self.text else "Null", self.cui if self.cui else "Null"]
        )


def datetime_normalize(time: str) -> datetime:
    time = time.replace("w", "W")
    if "W" in time:
        return datetime.strptime(time + "-1", "%Y-W%W-%w")
    return dateutil.parser.parse(time)


def read_all_patients(data_path: str) -> TimelineDict:
    # Note that all key/value pairs of JSON are always of the type str.
    # https://docs.python.org/3/library/json.html
    with open(data_path, "r") as fr:
        all_patient_timelines = json.load(fr)
    return all_patient_timelines


def relaxed_rel_eval(
    incorrect: List[TimelineTuple],
    missing: List[TimelineTuple],
    preds: List[TimelineTuple],
    golds: List[TimelineTuple],
) -> Tuple[List[TimelineTuple], List[TimelineTuple]]:
    not_truly_incorrect: List[TimelineTuple] = []
    not_truly_missing: List[TimelineTuple] = []
    for ptup in incorrect:
        is_not_truly_incorrect = False
        chemo, rel, timex = ptup
        # Basically we think contains-1 can be replaced by begins-on/ends-on,
        # and begins-on/ends-on can be replaced by contains-1.
        if rel in ["begins-on", "ends-on"]:
            if [chemo, "contains-1", timex] in golds:
                is_not_truly_incorrect = True
        elif rel == "contains-1":
            if [chemo, "begins-on", timex] in golds or [
                chemo,
                "ends-on",
                timex,
            ] in golds:
                is_not_truly_incorrect = True
        if is_not_truly_incorrect:
            not_truly_incorrect.append(ptup)

    for gtup in missing:
        is_not_truly_missing = False
        chemo, rel, timex = gtup

        if rel in ["begins-on", "ends-on"]:
            if [chemo, "contains-1", timex] in preds:
                is_not_truly_missing = True
        elif rel == "contains-1":
            if [chemo, "begins-on", timex] in preds or [
                chemo,
                "ends-on",
                timex,
            ] in preds:
                is_not_truly_missing = True
        if is_not_truly_missing:
            not_truly_missing.append(gtup)
    return not_truly_incorrect, not_truly_missing


def relaxed_within_range_eval(
    incorrect: List[TimelineTuple],
    missing: List[TimelineTuple],
    gold_chemos,
    pred_chemos,
) -> Tuple[List[TimelineTuple], List[TimelineTuple]]:
    """
    incorrect: false positive,
    missing: false negative,
    gold_chemos and pred_chemos, basically use Chemo object to get the start and end dates for each chemo
    """
    not_truly_incorrect = []
    not_truly_missing = []
    for ptup in incorrect:
        is_not_truly_incorrect = False
        source, rel, raw_target = ptup

        target = datetime_normalize(raw_target)

        if source in gold_chemos:
            gold_start, gold_end = (
                gold_chemos[source].first_start,
                gold_chemos[source].last_end,
            )
            if not gold_start or not gold_end:
                continue
            if rel in ["ENDS-ON", "ends-on"]:
                # The end date predicted by the system (target), is before the gold start date, then it's wrong.
                if target <= gold_start:
                    continue
            if rel in ["BEGINS-ON", "begins-on"]:
                # The start date predicted by the system (target) is after the gold end date, then it's wrong.
                if target >= gold_end:
                    continue
            # If the predicted date is in between the gold start and end date, i.e. in the correct range,
            # we consider this is correct, not truly false positive.
            if gold_start <= target <= gold_end:
                is_not_truly_incorrect = True
        if is_not_truly_incorrect:
            not_truly_incorrect.append(ptup)

    for gtup in missing:
        is_not_truly_missing = False
        source, rel, raw_target = gtup

        target = datetime_normalize(raw_target)

        if source in pred_chemos:
            pred_start, pred_end = (
                pred_chemos[source].first_start,
                pred_chemos[source].last_end,
            )
            if not pred_start or not pred_end:
                continue
            if rel in ["ENDS-ON", "ends-on"]:
                if target <= pred_start:
                    continue
            if rel in ["BEGINS-ON", "begins-on"]:
                if target >= pred_end:
                    continue
            # This is saying, for example, <taxol, contains-1, 2011-03-01> is missing in predictions, i.e. is false negative,
            # however, we can find <taxol, begins-on, 2011-01-01> and <taxol, ends-on, 2011-05-31> in predictions,
            # we consider <taxol, contains-1, 2011-03-01> is not false negative, because the system predicted the
            # correct range that covers the gold timeline.
            if pred_start <= target <= pred_end:
                is_not_truly_missing = True
        if is_not_truly_missing:
            not_truly_missing.append(gtup)
    return not_truly_incorrect, not_truly_missing


def group_chemo_dates(golds: List[TimelineTuple]) -> Dict[str, Chemo]:
    gold_group_by_start_end: Dict[str, Dict[str, List[datetime]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for tup in golds:
        source, label, raw_target = tup
        if label.upper() not in ["BEGINS-ON", "ENDS-ON"]:
            continue
        target = datetime_normalize(raw_target)

        gold_group_by_start_end[source][label].append(target)
    all_gold_chemos = {}
    # all_gold_chemos: maps from the text of this chemo to the chemo event obj
    for chemo, labels in gold_group_by_start_end.items():
        first_date = (
            min(labels["BEGINS-ON".lower()]) if "BEGINS-ON".lower() in labels else None
        )
        last_date = (
            max(labels["ENDS-ON".lower()]) if "ENDS-ON".lower() in labels else None
        )
        # For each chemo, find the earliest start date, and the latest end date,
        # use them to get the span of this chemo, the span will be used when doing relaxed evaluation.
        chemo_event = Chemo(text=chemo, first_start=first_date, last_end=last_date)
        all_gold_chemos[chemo] = chemo_event
    return all_gold_chemos


def normalize_to_month_and_year(
    golds: List[TimelineTuple],
) -> Tuple[List[TimelineTuple], List[TimelineTuple]]:
    month_only_pairs = []
    year_only_pairs = []
    for tup in golds:
        source, label, raw_target = tup
        target = datetime_normalize(raw_target)
        year = target.year
        month = target.month
        if month < 10:
            normalized_month = str(year) + "-0" + str(month)
        else:
            normalized_month = str(year) + "-" + str(month)
        normalized_year = str(year)

        month_pair = [source, label, normalized_month]
        year_pair = [source, label, normalized_year]

        if month_pair not in month_only_pairs:
            month_only_pairs.append(month_pair)
        if year_pair not in year_only_pairs:
            year_only_pairs.append(year_pair)

    has_more_specific_chemos_month = summarize_timeline(
        cast(List[TimelineTuple], month_only_pairs)
    )
    has_more_specific_chemos_year = summarize_timeline(
        cast(List[TimelineTuple], year_only_pairs)
    )
    month_only_pairs = [
        tup for tup in month_only_pairs if tup not in has_more_specific_chemos_month
    ]
    year_only_pairs = [
        tup for tup in year_only_pairs if tup not in has_more_specific_chemos_year
    ]
    return cast(List[TimelineTuple], month_only_pairs), cast(
        List[TimelineTuple], year_only_pairs
    )


def summarize_timeline(timelines: List[TimelineTuple]) -> List[TimelineTuple]:
    """
    This is to postprocess timelines one more time after we normalized original timeline
    to year only or month only timelines. What it does is: if we have a generic chemo mention,
    e.g. <chemotherapy, contains-1, 2011-01>, or <chemoradiation, contains-1, 2011-01>, we want to see
    if we can have more specific chemo mention happened on the same date with the same label,
    e.g. <Taxol, contains-1, 2011-01>. If we find a more specific chemo mention,
    we would ignore the generic chemo mention, only add <Taxol, contains-1, 2011-01> to the timeline.
    """
    date_rel_to_chemo: Dict[str, Dict[str, List[str]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for tup in timelines:
        chemo, rel, date = tup
        date_rel_to_chemo[date][rel].append(chemo)

    has_more_specific_chemos = []
    for date, rel_chemos in date_rel_to_chemo.items():
        for rel, chemos in rel_chemos.items():
            for chemo in chemos:
                # chemo.startswith("chemo") is how we check if this is a generic chemo mention
                if chemo.startswith("chemo"):
                    if len(date_rel_to_chemo[date][rel]) > 1:
                        has_more_specific_chemos.append(
                            cast(TimelineTuple, [chemo, rel, date])
                        )
    return has_more_specific_chemos


def compute_f(true_pos, false_pos, false_neg):
    precision = len(true_pos) / (len(true_pos) + len(false_pos))
    recall = len(true_pos) / (len(true_pos) + len(false_neg))
    f1 = 2 * (precision * recall) / (precision + recall)
    return f1


def strict_eval(
    gold: List[TimelineTuple], pred: List[TimelineTuple]
) -> Tuple[List[TimelineTuple], List[TimelineTuple], List[TimelineTuple]]:
    true_pos = [prediction for prediction in pred if prediction in gold]
    false_pos = [prediction for prediction in pred if prediction not in gold]
    false_neg = [correct for correct in gold if correct not in pred]
    not_truly_fp = fp_fn_single_count(false_pos, false_neg)
    false_pos = [pred for pred in false_pos if pred not in not_truly_fp]
    return true_pos, false_pos, false_neg


def fp_fn_single_count(false_pos, false_neg):
    """
    What it does here is: let's say in pred we have <Taxol, BEGINS-ON, 2011-01-01>,
    in gold we have <Taxol, CONTAINS-1, 2011-01-01>, then <Taxol, BEGINS-ON, 2011-01-01> would be false positive,
    <Taxol, CONTAINS-1, 2011-01-01> would be false negative, that means, the same mistake is counted twice,
    once in fp, once in fn. So, here, we want to make sure, we count <Taxol, CONTAINS-1, 2011-01-01> as false negative,
    and don't count <Taxol, BEGINS-ON, 2011-01-01> as false positive.
    """
    not_truly_fp = []
    # false_neg_tracker: (chemo, timex) to label
    false_neg_tracker = {(item[0], item[-1]): item[1] for item in false_neg}
    for ptup in false_pos:
        # E.g. we check in <Taxol, 2011-01-01> is already in false negative.
        if (ptup[0], ptup[-1]) in false_neg_tracker:
            not_truly_fp.append(ptup)
    return not_truly_fp


def relaxed_eval(
    gold: List[TimelineTuple],
    gold_chemo: Dict[str, Chemo],
    pred: List[TimelineTuple],
    pred_chemo: Dict[str, Chemo],
) -> Tuple[
    List[TimelineTuple],
    List[TimelineTuple],
    List[TimelineTuple],
    List[TimelineTuple],
    List[TimelineTuple],
    List[TimelineTuple],
    List[TimelineTuple],
]:
    true_pos = [prediction for prediction in pred if prediction in gold]
    false_pos = [prediction for prediction in pred if prediction not in gold]
    false_neg = [correct for correct in gold if correct not in pred]
    not_truly_fp_with_range, not_truly_fn_with_range = relaxed_within_range_eval(
        incorrect=false_pos,
        missing=false_neg,
        gold_chemos=gold_chemo,
        pred_chemos=pred_chemo,
    )
    not_truly_fp_with_label, not_truly_fn_with_label = relaxed_rel_eval(
        incorrect=false_pos, missing=false_neg, preds=pred, golds=gold
    )

    not_truly_fp_as_label_single_count = fp_fn_single_count(false_pos, false_neg)

    truly_tp = true_pos
    truly_fp, truly_fn = [], []
    for tup in false_pos:
        if tup in not_truly_fp_with_range or tup in not_truly_fp_with_label:
            # Add this one to true positive if it's not considered as true fp
            truly_tp.append(tup)
        elif tup in not_truly_fp_as_label_single_count:
            continue
        else:
            truly_fp.append(tup)
    for tup in false_neg:
        if tup in not_truly_fn_with_range or tup in not_truly_fn_with_label:
            continue
        truly_fn.append(tup)
    return (
        truly_tp,
        truly_fp,
        truly_fn,
        not_truly_fp_with_range,
        not_truly_fn_with_range,
        not_truly_fp_with_label,
        not_truly_fn_with_label,
    )


def evaluation(
    gold: List[TimelineTuple], pred: List[TimelineTuple], strict: bool, relaxed_to: str
):
    # Get the earliest start and latest end dates for each chemo
    all_gold_chemos = group_chemo_dates(gold)
    all_pred_chemos = group_chemo_dates(pred)

    gold_month_timeline, gold_year_timeline = normalize_to_month_and_year(gold)
    pred_month_timeline, pred_year_timeline = normalize_to_month_and_year(pred)

    rmv_from_fp_range: List[TimelineTuple] = []
    rmv_from_fn_range: List[TimelineTuple] = []
    rmv_from_fp_label: List[TimelineTuple] = []
    rmv_from_fn_label: List[TimelineTuple] = []
    if strict:
        true_pos, false_pos, false_neg = strict_eval(gold, pred)
    else:
        assert relaxed_to, "For relaxed evaluation, please specify --relaxed_to flag"
        if relaxed_to == "day":
            # rmv_from_fp_range: remove from false positive because it's in the right range;
            # rmv_from_fp_label: remove from false positive because the label is consider correct.
            (
                true_pos,
                false_pos,
                false_neg,
                rmv_from_fp_range,
                rmv_from_fn_range,
                rmv_from_fp_label,
                rmv_from_fn_label,
            ) = relaxed_eval(gold, all_gold_chemos, pred, all_pred_chemos)
        elif relaxed_to == "month":
            (
                true_pos,
                false_pos,
                false_neg,
                rmv_from_fp_range,
                rmv_from_fn_range,
                rmv_from_fp_label,
                rmv_from_fn_label,
            ) = relaxed_eval(
                gold_month_timeline,
                all_gold_chemos,
                pred_month_timeline,
                all_pred_chemos,
            )
        elif relaxed_to == "year":
            (
                true_pos,
                false_pos,
                false_neg,
                rmv_from_fp_range,
                rmv_from_fn_range,
                rmv_from_fp_label,
                rmv_from_fn_label,
            ) = relaxed_eval(
                gold_year_timeline, all_gold_chemos, pred_year_timeline, all_pred_chemos
            )
        else:
            raise ValueError("--relaxed_to must be one of 'day', 'month', or 'year'")

    logger.debug(f"true_pos... {len(true_pos)}")
    for item in true_pos:
        logger.debug(f"{item}")
        pass
    logger.debug(f"false_pos... {len(false_pos)}")
    for item in false_pos:
        logger.debug(f"{item}")
        pass
    logger.debug(f"false_neg... {len(false_neg)}")
    for item in false_neg:
        logger.debug(f"{item}")
        pass

    if not strict:
        logger.debug(f"removed from false_pos_range... {len(rmv_from_fp_range)}")
        logger.debug(f"{rmv_from_fp_range}")
        logger.debug(f"removed from false_pos_label... {len(rmv_from_fp_label)}")
        logger.debug(f"{rmv_from_fp_label}")

        logger.debug(f"removed from false_neg_range... {len(rmv_from_fn_range)}")
        logger.debug(f"{rmv_from_fn_range}")
        logger.debug(f"removed from false_neg_label... {len(rmv_from_fn_label)}")
        logger.debug(f"{rmv_from_fn_label}")

    if len(true_pos) + len(false_neg) == 0:
        precision, recall, f1 = 0.0, 0.0, 0.0
    elif len(true_pos) + len(false_pos) == 0:
        precision, recall, f1 = 0.0, 0.0, 0.0
    else:
        precision = len(true_pos) / (len(true_pos) + len(false_pos))
        recall = len(true_pos) / (len(true_pos) + len(false_neg))
        if precision + recall:
            f1 = 2 * (precision * recall) / (precision + recall)
        else:
            f1 = 0.0

    logger.info(f"precision: {precision}")
    logger.info(f"recall: {recall}")
    logger.info(f"f1: {f1}")

    return true_pos, false_pos, false_neg, precision, recall, f1


def read_files(
    gold_id_path: str,
    all_id_path: str,
    gold_timeline: str | TimelineDict,
    pred_timeline: str | TimelineDict,
) -> Tuple[TimelineDict, TimelineDict, List[str], List[str]]:
    gold_all_patient = (
        read_all_patients(gold_timeline)
        if isinstance(gold_timeline, str)
        else gold_timeline
    )
    pred_all_patient = (
        read_all_patients(pred_timeline)
        if isinstance(pred_timeline, str)
        else pred_timeline
    )

    with open(all_id_path, mode="rt") as fp:
        all_ids = [line.splitlines()[0] for line in fp.readlines()]

    # Sanity check
    print(gold_all_patient)
    if len(all_ids) == 0 or len(all_ids) != len(pred_all_patient):
        # raise ValueError(f"Malformated or some patients are missing in prediction file. all ids {sorted(all_ids)} vs {sorted(pred_all_patient)}")
        # print(f"before {pred_all_patient}")
        for empty_patient in set(all_ids) - pred_all_patient.keys():
            pred_all_patient[empty_patient] = []
        # print(f"after {pred_all_patient}")
    if not gold_id_path:
        if len(all_ids) != len(gold_all_patient):
            raise ValueError(
                "Error in gold annotated ids file. Check the content of file in gold_id_path. Length of all_ids: %s, gold_all_patient: %s"
                % (len(all_ids), len(gold_all_patient))
            )

    # Only for orgnizers: for test dataset - screening silver datasets
    if gold_id_path:
        if not os.path.exists(gold_id_path):
            raise ValueError("Error in gold annotated ids file path")
        with open(gold_id_path, mode="rt") as fp:
            gold_ids = [line.splitlines()[0] for line in fp.readlines()]

        if (
            len(gold_ids) == 0
            or len(gold_ids) != len(gold_all_patient)
            or len(gold_ids) > len(pred_all_patient)
        ):
            raise ValueError(
                "Error in gold annotated ids file. Check the content of file in gold_id_path"
            )
    else:
        gold_ids = list(gold_all_patient.keys())

    return pred_all_patient, gold_all_patient, all_ids, gold_ids


def micro_average_metrics(
    all_true_pos: Dict[str, int],
    all_false_pos: Dict[str, int],
    all_false_neg: Dict[str, int],
) -> float:
    # Micro average metrics
    logger.info(
        f"tp, fp, fn over all patients: {sum(all_true_pos.values())}, {sum(all_false_pos.values())}, {sum(all_false_neg.values())}"
    )
    if len(all_true_pos) + len(all_false_pos) != 0:
        micro_precision = sum(all_true_pos.values()) / (
            sum(all_true_pos.values()) + sum(all_false_pos.values())
        )
    else:
        micro_precision = 0
    if len(all_true_pos) + len(all_false_neg) != 0:
        micro_recall = sum(all_true_pos.values()) / (
            sum(all_true_pos.values()) + sum(all_false_neg.values())
        )
    else:
        micro_recall = 0
    if micro_precision + micro_recall:
        micro_f1 = (
            2 * (micro_precision * micro_recall) / (micro_precision + micro_recall)
        )
    else:
        micro_f1 = 0

    print("Micro average metrics")
    print("Micro precision:", micro_precision)
    print("Micro recall:", micro_recall)
    print("Micro f1:", micro_f1)
    print()

    return micro_f1


def macro_average_metrics(
    local_precision: Dict[str, float],
    local_recall: Dict[str, float],
    local_f1: Dict[str, float],
    local_relations: Dict[str, int],
) -> Tuple[float, float]:
    # Macro average metrics
    print("Macro average metrics")

    print("Type A evaluation: including the notes with no true relations")
    total_patients = len(local_f1)
    assert (
        len(local_precision) == len(local_recall) == len(local_f1) == total_patients
    ), "total_patients should be equal to the number of local metrics"

    type_a_macro_prec = sum(local_precision.values()) / total_patients
    type_a_macro_rec = sum(local_recall.values()) / total_patients
    type_a_macro_f1 = sum(local_f1.values()) / total_patients
    print("[Type A] Macro precision:", type_a_macro_prec)
    print("[Type A] Macro recall:", type_a_macro_rec)
    print("[Type A] Macro f1: ", type_a_macro_f1)
    print()

    print("Type B evaluation: excluding the notes with no true relations")
    type_b_prec = [
        score
        for pat_id, score in local_precision.items()
        if local_relations[pat_id] != 0
    ]
    type_b_rec = [
        score for pat_id, score in local_recall.items() if local_relations[pat_id] != 0
    ]
    type_b_f1 = [
        score for pat_id, score in local_f1.items() if local_relations[pat_id] != 0
    ]
    assert (
        len(type_b_prec) == len(type_b_rec) == len(type_b_f1)
    ), "The number of local metrics should be the same."

    type_b_macro_prec = sum(type_b_prec) / len(type_b_prec)
    type_b_macro_rec = sum(type_b_rec) / len(type_b_rec)
    type_b_macro_f1 = sum(type_b_f1) / len(type_b_f1)
    print("[Type B] Macro precision:", type_b_macro_prec)
    print("[Type B] Macro recall:", type_b_macro_rec)
    print("[Type B] Macro f1: ", type_b_macro_f1)
    print()

    return type_a_macro_f1, type_b_macro_f1


def results_and_instances(
    pred_all_patient: TimelineDict,
    gold_all_patient: TimelineDict,
    gold_ids: List[str],
    strict: bool,
    relaxed_to: str,
    debug: bool,
) -> Tuple[SSMatrix, MetricMatrix, DebugDict, Dict[str, int]]:
    all_true_pos, all_false_pos, all_false_neg = {}, {}, {}
    local_relations = {}  # Key = patient ID, Value = number of timeline in the patient
    local_f1, local_precision, local_recall = (
        {},
        {},
        {},
    )  # Key = patient ID, Value = score for the patient
    fn_fp_debug: DebugDict = defaultdict(
        lambda: defaultdict(list)
    )  # patient_id -> <FP | FN> -> List[<chemo, timex, rel>]
    for pred_patient, pred_timeline in pred_all_patient.items():
        # pred_patient: patient ID; pred_timeline: a list of <chemo, label, timex>
        if pred_patient not in gold_ids:
            continue

        if pred_patient not in gold_all_patient:
            raise ValueError(
                "The given patient ID: '%s' does not exist in the gold annotated dataset."
                % pred_patient
            )
        gold_timeline = gold_all_patient.get(pred_patient, [])

        # Handling patients without timeline. Models are expected to output "no timeline".
        """
            For type A evaluation, metric includes the patients with no gold timelines:
            If models successfully make no prediction, they will achieve scores of 
            1.0, 1.0, and 1.0 for the given patient in local precision, recall, and F1, respectively.
        """
        if len(gold_timeline) == 0:
            true_pos: List[TimelineTuple] = []
            false_pos: List[TimelineTuple] = []
            false_neg: List[TimelineTuple] = []
            if len(pred_timeline) == 0:
                true_pos = []
                false_pos = []
                false_neg = []
                p, r, f_score = 1, 1, 1
            else:
                true_pos = []
                false_pos = cast(List[TimelineTuple], pred_timeline)
                false_neg = []
                p, r, f_score = 0, 0, 0
            local_relations[pred_patient] = 0
        else:
            logger.info(f"pred_patient ID: {pred_patient}")
            true_pos, false_pos, false_neg, p, r, f_score = evaluation(
                gold=cast(List[TimelineTuple], gold_timeline),
                pred=cast(List[TimelineTuple], pred_timeline),
                strict=strict,
                relaxed_to=relaxed_to,
            )
            local_relations[pred_patient] = len(gold_timeline)

        all_true_pos[pred_patient] = len(true_pos)
        all_false_pos[pred_patient] = len(false_pos)
        all_false_neg[pred_patient] = len(false_neg)
        if debug:
            fn_fp_debug[pred_patient]["false_positive"].extend(false_pos)
            fn_fp_debug[pred_patient]["false_negative"].extend(false_neg)
        local_precision[pred_patient] = float(p)
        local_recall[pred_patient] = float(r)
        local_f1[pred_patient] = float(f_score)
    ss_matrix = all_true_pos, all_false_pos, all_false_neg
    metric_matrix = local_f1, local_precision, local_recall
    return ss_matrix, metric_matrix, fn_fp_debug, local_relations


def evaluate_and_log(
    gold_id_path: str,
    all_id_path: str,
    gold_timeline: str | TimelineDict,
    pred_timeline: str | TimelineDict,
    strict: bool,
    relaxed_to: str,
    debug: bool,
    patient_level_metrics: bool,
) -> Optional[DebugDict]:
    final_out: Optional[DebugDict] = None
    print(f"Evaluation code for ChemoTimelines Shared Task. Version: {VERSION}")
    print("Reading from files...")
    pred_all_patient, gold_all_patient, all_ids, gold_ids = read_files(
        gold_id_path=gold_id_path,
        all_id_path=all_id_path,
        gold_timeline=gold_timeline,
        pred_timeline=pred_timeline,
    )
    print(
        f"predicted output: {len(pred_all_patient)}, gold annotation: {len(gold_all_patient)}, all ids: {len(all_ids)}\n"
    )
    ss_matrix, metrics_matrix, fn_fp_debug, local_relations = results_and_instances(
        pred_all_patient=pred_all_patient,
        gold_all_patient=gold_all_patient,
        gold_ids=gold_ids,
        strict=strict,
        relaxed_to=relaxed_to,
        debug=debug,
    )
    all_true_pos, all_false_pos, all_false_neg = ss_matrix
    local_f1, local_precision, local_recall = metrics_matrix
    _ = micro_average_metrics(
        all_true_pos=all_true_pos,
        all_false_pos=all_false_pos,
        all_false_neg=all_false_neg,
    )

    if patient_level_metrics:
        with open("patient_level_metrics.txt", mode="wt") as patient_metrics_file:
            for patient_id, precision in local_precision.items():
                recall = local_recall[patient_id]
                f1 = local_f1[patient_id]
                patient_metrics_file.write(f"{patient_id} micro average metrics\n\n")
                patient_metrics_file.write(f"Micro precision: {precision}\n")
                patient_metrics_file.write(f"Micro recall: {recall}\n")
                patient_metrics_file.write(f"Micro f1: {f1}\n\n")
    if debug:
        with open("patient_level_debug.json", mode="wt") as patient_debug_file:
            json.dump(fn_fp_debug, patient_debug_file, indent=4)
        final_out = fn_fp_debug
    type_a_macro_f1, type_b_macro_f1 = macro_average_metrics(
        local_precision=local_precision,
        local_recall=local_recall,
        local_f1=local_f1,
        local_relations=local_relations,
    )

    print(
        "Official Score: Arithmetic mean of two types of Macro F1, type A and B, "
        + "in 'relaxed to month' setting will be used for the leaderboard. "
    )
    if not (strict) and relaxed_to == "month":
        print(f"Official Score: {(type_a_macro_f1 + type_b_macro_f1) / 2}")
    else:
        print(
            "To see the official score, please run without --strict flag, and set 'relaxed to month' setting by --relaxed_to=month "
        )

    print("Evaluation completed!")
    return final_out


def driver(args: argparse.Namespace) -> None:
    logger.info(args)
    evaluate_and_log(
        args.gold_id_path,
        args.all_id_path,
        args.gold_path,
        args.pred_path,
        args.strict,
        args.relaxed_to,
        args.debug,
        args.patient_level_metrics,
    )


def main() -> None:
    args = parser.parse_args()
    driver(args)


if __name__ == "__main__":
    main()
