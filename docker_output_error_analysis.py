import argparse
import json
from datetime import datetime
from typing import Dict, List, Tuple, Union
from collections import Counter

import dateutil.parser
import pandas as pd
from tabulate import tabulate
from enum import Enum

parser = argparse.ArgumentParser(description="")

parser.add_argument("--docker_tsv", type=str)
parser.add_argument("--error_json", type=str)
parser.add_argument("--output_dir", type=str)
parser.add_argument("--strict", action="store_true", help="do strict eval")
parser.add_argument(
    "--relaxed_to",
    help="Type 'year' to only evaluate year, 'month' to evaluate year and month, "
    "or 'day' to evaluate year-month-day",
    choices=["day", "month", "year"],
)
# chemo, rel, timex, filename
instance = List[str]
label_to_hierarchy = {
    "begins-on": 1,
    "ends-on": 1,
    "contains": 2,
    "contains-1": 2,
    "before": 3,
}

source_header = ["chemo", "tlink", "normed timex", "mode"]
pred_header = ["chemo", "tlink", "normed timex", "note name"]
eval_modes = {"strict", "day", "month", "year"}


class ErrorType(Enum):
    FALSE_POSITIVE = 1
    FALSE_NEGATIVE = 2


class ErrorCause(Enum):
    TLINK = 1
    CHEMO = 2
    TIMEX = 3


class ErrorDebug:
    def __init__(
        self, instances: List[instance], error_type: ErrorType, error_cause: ErrorCause
    ) -> None:
        self.source_instance = instances[0]
        self.pred_instances = instances[1:]
        self.error_type = error_type
        self.error_cause = error_cause

    def generate_report(self) -> str:
        source_table = tabulate(
            [source_header, self.source_instance], headers="firstrow"
        )
        pred_table = (
            tabulate([pred_header, *self.pred_instances], headers="firstrow")
            if len(self.pred_instances) > 0
            else ""
        )
        return f"\n\nSource Instance:\n\n{source_table}\n\nPredicted Instances:\n\n{pred_table}"

    def __str__(self) -> str:
        report = self.generate_report()
        return report

    def get_error_cause(self ) -> ErrorCause:
        return self.error_cause

def raw_normalize(time: str) -> str:
    return time.lower().split("t")[0]


def datetime_normalize(time: str) -> datetime:
    time = raw_normalize(time)
    if "W" in time:
        return datetime.strptime(time + "-1", "%Y-W%W-%w")
    return dateutil.parser.parse(time)


def strict_eval(time1: str, time2: str) -> bool:
    return raw_normalize(time1) == raw_normalize(time2)


def relaxed_day_eval(time1: str, time2: str) -> bool:
    norm_time1 = datetime_normalize(time1)
    norm_time2 = datetime_normalize(time2)
    # as far as I can tell this is how it works
    # in eval_timeline.py
    return norm_time1 == norm_time2


def relaxed_month_eval(time1: str, time2: str) -> bool:
    norm_time1 = datetime_normalize(time1)
    norm_time2 = datetime_normalize(time2)
    return (norm_time1.year, norm_time1.month) == (norm_time2.year, norm_time2.month)


def relaxed_year_eval(time1: str, time2: str) -> bool:
    norm_time1 = datetime_normalize(time1)
    norm_time2 = datetime_normalize(time2)
    return norm_time1.year == norm_time2.year


mode_to_eval = {
    "strict": strict_eval,
    "relaxed_day": relaxed_day_eval,
    "relaxed_month": relaxed_month_eval,
    "relaxed_year": relaxed_year_eval,
}


def compatible_time(time1: str, time2: str, eval_mode: str) -> bool:
    if eval_mode in mode_to_eval:
        return mode_to_eval[eval_mode](time1, time2)
    else:
        ValueError(f"You picked a non existant mode: {eval_mode}")
        return False


# Pandas has really dumb syntax
def compatible_chemos(chemo1: str, chemo2: str):
    return chemo1.lower() == chemo2.lower()


def collect_error_events(
    patient_id: str,
    eval_mode: str,
    error_dict: Dict[str, List[List[str]]],
    docker_df: pd.DataFrame,
) -> Tuple[List[ErrorDebug], List[ErrorDebug]]:
    patient_df = docker_df[docker_df["patient_id"].isin([patient_id])]
    fp_events = collect_fp_events(error_dict["false_positive"], eval_mode, patient_df)
    fn_events = collect_fn_events(error_dict["false_negative"], eval_mode, patient_df)
    return fp_events, fn_events


def preimage_and_cause(
    error_type: ErrorType,
    chemo_text: str,
    tlink: str,
    normed_timex: str,
    eval_mode: str,
    patient_df: pd.DataFrame,
) -> Tuple[List[instance], ErrorCause]:
    def row_chemo_compatible(row_chemo: str) -> bool:
        return compatible_chemos(row_chemo, chemo_text)

    chemo_matches = patient_df.loc[patient_df["chemo_text"].apply(row_chemo_compatible)]
    if len(chemo_matches) == 0:
        return [], ErrorCause.CHEMO

    def row_timex_compatible(row_timex: str) -> bool:
        return compatible_time(row_timex, normed_timex, eval_mode)

    timex_chemo_matches = chemo_matches.loc[
        chemo_matches["normed_timex"].apply(row_timex_compatible)
    ]

    if len(timex_chemo_matches) == 0:
        return [], ErrorCause.TIMEX

    if error_type == ErrorType.FALSE_NEGATIVE:
        return (
            timex_chemo_matches[
                ["chemo_text", "tlink", "normed_timex", "note_name"]
            ].values.tolist(),
            ErrorCause.TLINK,
        )

    def row_tlink_compatible(row_tlink: str) -> bool:
        return row_tlink.lower() == tlink.lower()

    full_matches = timex_chemo_matches.loc[
        timex_chemo_matches["tlink"].apply(row_tlink_compatible)
    ]
    return (
        full_matches[
            ["chemo_text", "tlink", "normed_timex", "note_name"]
        ].values.tolist(),
        ErrorCause.TLINK,
    )


def collect_fp_events(
    false_positives: List[List[str]], eval_mode: str, patient_df: pd.DataFrame
) -> List[ErrorDebug]:
    def fp_instance(timelines_event: List[str]) -> ErrorDebug:
        # since the resulting tlink is predetermined
        # given the conflict resolution rules
        chemo_text, tlink, normed_timex = timelines_event
        summarization_preimage, error_cause = preimage_and_cause(
            ErrorType.FALSE_POSITIVE,
            chemo_text,
            tlink,
            normed_timex,
            eval_mode,
            patient_df,
        )
        return ErrorDebug(
            [[*timelines_event, eval_mode], *summarization_preimage],
            ErrorType.FALSE_POSITIVE,
            error_cause,
        )

    return [fp_instance(event) for event in false_positives]


def collect_fn_events(
    false_negatives: List[List[str]], eval_mode: str, patient_df: pd.DataFrame
) -> List[ErrorDebug]:
    def fn_instance(timelines_event: List[str]) -> ErrorDebug:
        # since the resulting tlink is predetermined
        # given the conflict resolution rules
        chemo_text, tlink, normed_timex = timelines_event
        summarization_preimage, error_cause = preimage_and_cause(
            ErrorType.FALSE_POSITIVE,
            chemo_text,
            tlink,
            normed_timex,
            eval_mode,
            patient_df,
        )
        return ErrorDebug(
            [[*timelines_event, eval_mode], *summarization_preimage],
            ErrorType.FALSE_NEGATIVE,
            error_cause,
        )

    return [fn_instance(event) for event in false_negatives]


def write_patient_error_reports(
    patient_id: str,
    fp_events: List[ErrorDebug],
    fn_events: List[ErrorDebug],
    output_dir: str,
):
    fn = output_dir + "/" + patient_id + "_error_analysis.txt"

    fp_error_causes = Counter(map(ErrorDebug.get_error_cause, fp_events))
    fn_error_causes = Counter(map(ErrorDebug.get_error_cause, fn_events))

    fp_str = (
        "\n\nFalse Positives\n\n" + "\n".join(map(str, fp_events))
        if len(fp_events) > 0
        else ""
    )
    fn_str = (
        "\n\nFalse Negatives\n\n" + "\n".join(map(str, fn_events))
        if len(fn_events) > 0
        else ""
    )
    with open(fn, "wt") as fn_out:
        fn_out.write(fp_str)
        fn_out.write(fn_str)


def write_instances(docker_tsv: str, error_json: str, output_dir: str, eval_mode: str):
    docker_df = pd.read_csv(docker_tsv, sep="\t")
    with open(error_json, mode="rt") as json_f:
        patient_to_errors = json.load(json_f)
    for patient_id, error_dict in patient_to_errors.items():
        fp_events, fn_events = collect_error_events(
            patient_id, eval_mode, error_dict, docker_df
        )
        write_patient_error_reports(patient_id, fp_events, fn_events, output_dir)


def main():
    args = parser.parse_args()
    if args.strict:
        eval_mode = "strict"
    else:
        eval_mode = args.relaxed_to
    write_instances(args.docker_tsv, args.error_json, args.output_dir, eval_mode)


if __name__ == "__main__":
    main()
