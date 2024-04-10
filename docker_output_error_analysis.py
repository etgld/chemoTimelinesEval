import argparse
import json
from collections import Counter
from datetime import datetime
from enum import Enum
from functools import reduce
from typing import Dict, List, Tuple, cast

import dateutil.parser
import pandas as pd
from tabulate import tabulate

from eval_timeline import DebugDict, TimelineTuple, TimelineTuples

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
pred_header = ["chemo", "tlink", "normed timex", "note name", "classifier input"]
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
        if len(self.pred_instances) > 0:
            print(self.source_instance)
            print(self.pred_instances)
            print(pred_table)
        return f"\n\n{self.error_cause} Instance in Summary:\n\n{source_table}\n\nInstances from Docker Output:\n\n{pred_table}"

    def __str__(self) -> str:
        report = self.generate_report()
        return report

    def get_error_cause(self) -> ErrorCause:
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


def compatible_chemos(chemo1: str, chemo2: str):
    return chemo1.lower() == chemo2.lower()


def collect_error_events(
    patient_id: str,
    eval_mode: str,
    error_dict: Dict[str, TimelineTuples],
    docker_df: pd.DataFrame,
) -> Tuple[List[ErrorDebug], List[ErrorDebug]]:
    patient_df = docker_df.loc[docker_df["patient_id"] == patient_id]
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
    def row_chemo_compatible(pandas_row: pd.Series) -> bool:
        row_chemo = cast(str, pandas_row.chemo_text)
        return compatible_chemos(row_chemo, chemo_text)

    chemo_matches = patient_df.loc[patient_df.apply(row_chemo_compatible, axis=1)]
    if len(chemo_matches) == 0:
        return [], ErrorCause.CHEMO

    def row_timex_compatible(pandas_row: pd.Series) -> bool:
        row_timex = cast(str, pandas_row.normed_timex)
        return compatible_time(row_timex, normed_timex, eval_mode)

    timex_chemo_matches = chemo_matches.loc[
        chemo_matches.apply(row_timex_compatible, axis=1)
    ]

    if len(timex_chemo_matches) == 0:
        return [], ErrorCause.TIMEX

    if error_type == ErrorType.FALSE_NEGATIVE:
        return (
            timex_chemo_matches[
                ["chemo_text", "tlink", "normed_timex", "note_name", "tlink_inst"]
            ].values.tolist(),
            ErrorCause.TLINK,
        )

    def row_tlink_compatible(pandas_row: pd.Series) -> bool:
        row_tlink = cast(str, pandas_row.tlink)
        return row_tlink.lower() == tlink.lower()

    full_matches = timex_chemo_matches.loc[
        timex_chemo_matches.apply(row_tlink_compatible, axis=1)
    ]
    return (
        full_matches[
            ["chemo_text", "tlink", "normed_timex", "note_name", "tlink_inst"]
        ].values.tolist(),
        ErrorCause.TLINK,
    )


def collect_fp_events(
    false_positives: TimelineTuples, eval_mode: str, patient_df: pd.DataFrame
) -> List[ErrorDebug]:
    def fp_instance(timelines_event: TimelineTuple) -> ErrorDebug:
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
    false_negatives: TimelineTuples, eval_mode: str, patient_df: pd.DataFrame
) -> List[ErrorDebug]:
    def fn_instance(timelines_event: TimelineTuple) -> ErrorDebug:
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


def get_error_cause_count(events: List[ErrorDebug]) -> Counter[ErrorCause]:
    return Counter(map(ErrorDebug.get_error_cause, events))


def get_type_raw_total(
    error_type: ErrorType,
    patient_error_dict: Dict[str, Dict[ErrorType, Counter[ErrorCause]]],
) -> int:
    return sum(v[error_type].total() for v in patient_error_dict.values())


def get_type_cause_totals(
    error_type: ErrorType,
    patient_error_dict: Dict[str, Dict[ErrorType, Counter[ErrorCause]]],
) -> Counter[ErrorCause]:
    return cast(
        Counter,
        reduce(
            Counter.__add__,
            (v[error_type] for v in patient_error_dict.values()),
        ),
    )


def patients_by_type(
    error_type: ErrorType,
    patient_error_dict: Dict[str, Dict[ErrorType, Counter[ErrorCause]]],
) -> Counter[str]:
    return Counter({k: v[error_type].total() for k, v in patient_error_dict.items()})


def patients_by_cause(
    error_cause: ErrorCause,
    patient_error_dict: Dict[str, Dict[ErrorType, Counter[ErrorCause]]],
) -> Counter[str]:
    def cause_total(type_dict: Dict[ErrorType, Counter[ErrorCause]]) -> int:
        total_counter = cast(
            Counter,
            reduce(Counter.__add__, type_dict.values()),
        )
        return total_counter[error_cause]

    return Counter({k: cause_total(v) for k, v in patient_error_dict.items()})


def write_summaries(
    patient_error_dict: Dict[str, Dict[ErrorType, Counter[ErrorCause]]], output_dir: str
) -> None:
    fp_total: int = get_type_raw_total(ErrorType.FALSE_POSITIVE, patient_error_dict)
    fn_total: int = get_type_raw_total(ErrorType.FALSE_NEGATIVE, patient_error_dict)
    fp_causes: Counter[ErrorCause] = get_type_cause_totals(
        ErrorType.FALSE_POSITIVE, patient_error_dict
    )
    fn_causes: Counter[ErrorCause] = get_type_cause_totals(
        ErrorType.FALSE_NEGATIVE, patient_error_dict
    )
    total_causes: Counter[ErrorCause] = fp_causes + fn_causes
    cause_to_patient_count: Dict[ErrorCause, Counter[str]] = {
        error_cause: patients_by_cause(error_cause, patient_error_dict)
        for error_cause in ErrorCause
    }
    type_to_patient_count: Dict[ErrorType, Counter[str]] = {
        error_type: patients_by_type(error_type, patient_error_dict)
        for error_type in ErrorType
    }

    def write_patient_count_dict(patient_count_dict, writer, top_k):
        for error_type, patient_counts in patient_count_dict.items():
            writer.write(f"\n\n\nTop {top_k} Patients for {error_type}\n")
            for patient_id, count in patient_counts.most_common(top_k):
                if count > 0:
                    writer.write(f"{patient_id}\t{count}\n")

    with open(output_dir + "/metrics.txt", mode="wt") as outfile:
        outfile.write(
            f"Total False Positives: {fp_total}\tTotal False Negatives: {fn_total}\n\n\n"
        )
        outfile.write("ErrorCause\tTotal\n")
        for error_cause, count in total_causes.most_common():
            outfile.write(f"{error_cause}\t{count}\n")
        write_patient_count_dict(cause_to_patient_count, outfile, 10)
        write_patient_count_dict(type_to_patient_count, outfile, 10)


def write_patient_error_reports(
    patient_id: str,
    fp_events: List[ErrorDebug],
    fn_events: List[ErrorDebug],
    output_dir: str,
) -> None:
    fn = output_dir + "/" + patient_id + "_error_analysis.txt"

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
    with open(fn, mode="wt") as fn_out:
        fn_out.write(fp_str)
        fn_out.write(fn_str)


def write_instances_and_summaries(
    docker_tsv: str, error_dict: str | DebugDict, output_dir: str, eval_mode: str
) -> None:
    docker_df = pd.read_csv(docker_tsv, delimiter="\t")
    patient_error_dict: Dict[str, Dict[ErrorType, Counter[ErrorCause]]] = {}
    if isinstance(error_dict, str):
        with open(error_dict, mode="rt") as json_f:
            patient_to_errors = cast(
                Dict[str, Dict[str, TimelineTuples]], json.load(json_f)
            )
    else:
        patient_to_errors = error_dict
    for patient_id, error_type_dict in patient_to_errors.items():
        fp_events, fn_events = collect_error_events(
            patient_id, eval_mode, error_type_dict, cast(pd.DataFrame, docker_df)
        )
        write_patient_error_reports(patient_id, fp_events, fn_events, output_dir)
        patient_error_dict[patient_id] = {
            ErrorType.FALSE_NEGATIVE: get_error_cause_count(fn_events),
            ErrorType.FALSE_POSITIVE: get_error_cause_count(fp_events),
        }
    write_summaries(patient_error_dict, output_dir)


def main():
    args = parser.parse_args()
    if args.strict:
        eval_mode = "strict"
    else:
        eval_mode = args.relaxed_to
    write_instances_and_summaries(
        args.docker_tsv, args.error_json, args.output_dir, eval_mode
    )


if __name__ == "__main__":
    main()
