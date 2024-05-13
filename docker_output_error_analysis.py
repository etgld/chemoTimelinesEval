import argparse
import json
import os
import random
import textwrap
from collections import Counter
from datetime import datetime
from enum import Enum
from functools import reduce
from itertools import groupby
from math import ceil
from operator import itemgetter
from typing import Dict, Iterable, List, Sequence, Tuple, Union, cast

import dateutil.parser
import pandas as pd
from more_itertools import flatten, partition
from scipy.stats import norm
from tabulate import tabulate

from docker_output_to_timeline import clean_timex, keep_normalized_timex
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
CLASSIFIER_INST_WIDTH = 60
# chemo, rel, timex, filename
instance = List[str]
label_to_hierarchy = {
    "begins-on": 1,
    "ends-on": 1,
    "contains": 2,
    "contains-1": 2,
    "before": 3,
    "after": 3,
    "none": 4,
}

source_header = ["chemo", "tlink", "normed timex", "evaluation mode"]
pred_header = ["chemo", "tlink", "normed timex", "note", "DCT", "instance"]
pandas_columns = [
    "ERROR TYPE",
    "ERROR CAUSE",
    "CHEMO",
    "NORMED TIMEX",
    "TLINK",
    "TEXT",
    "NOTE",
    "DCT",
    "CHEMO ERROR PRESENT",
    "CHEMO ERROR COMMENTS",
    "TLINK ERROR PRESENT",
    "TLINK ERROR COMMENTS",
    "TIMEX DETECTION ERROR PRESENT",
    "TIMEX DETECTION ERROR COMMENTS",
    "TIMEX NORMALIZATION ERROR PRESENT",
    "TIMEX NORMALIZATION ERROR COMMENTS",
    "SUMMARIZATION ERROR PRESENT",
    "SUMMARIZATION ERROR COMMENTS",
    "ANALYSIS ERROR PRESENT",
    "ANALYSIS ERROR COMMENTS",
    "ANNOTATION ERROR PRESENT",
    "ANNOTATION ERROR COMMENTS",
]
eval_modes = {"strict", "day", "month", "year"}
docker_output_columns = [
    "chemo_text",
    "tlink",
    "normed_timex",
    "note_name",
    "DCT",
    "tlink_inst",
]


def sample_size(
    population_size: int,
    confidence_level: float = 0.95,
    margin_of_error: float = 0.05,
    estimated_proportion: float = 0.5,
) -> int:
    z_score = norm.ppf(1 - (1 - confidence_level) / 2)
    repeated_numerator = pow(z_score, 2.0) * (
        estimated_proportion * (1 - estimated_proportion)
    )
    margin_squared = pow(margin_of_error, 2.0)
    unlimited_sample_size = repeated_numerator / margin_squared
    finite_case_denominator = 1 + (
        repeated_numerator / (population_size * margin_squared)
    )
    float_finite_sample_size = unlimited_sample_size / finite_case_denominator
    return ceil(float_finite_sample_size)


def normalize_str(string: str) -> str:
    return string[0].capitalize() + string[1:].lower()


class ErrorType(Enum):
    def __str__(self) -> str:
        return " ".join(normalize_str(s) for s in self.name.split("_"))

    FALSE_POSITIVE = 1
    FALSE_NEGATIVE = 2


class ErrorCause(Enum):
    def __str__(self) -> str:
        return normalize_str(self.name)

    TLINK = 1
    CHEMO = 2
    TIMEX = 3


class ErrorDebug:
    def __init__(
        self, instances: List[instance], error_type: ErrorType, error_cause: ErrorCause
    ) -> None:
        self.source_instance = instances[0]
        source_tlink = self.source_instance[1]
        # sort by label precedence then note name
        sorted_insts = ErrorDebug._format_and_sort_instances(instances[1:])
        tlink_non_matches, tlink_matches = partition(
            lambda inst: inst[1] == source_tlink, sorted_insts
        )
        self.tlink_matches = list(tlink_matches)
        self.chemo_and_timex_matches = list(tlink_non_matches)
        self.error_type = error_type
        self.error_cause = error_cause

    @staticmethod
    def _format_and_sort_instances(instances: Iterable[instance]) -> List[instance]:
        def format_date(raw_date: str) -> str:
            year = raw_date[0:4]
            month = raw_date[4:6]
            day = raw_date[6:]
            return year + "_" + month + "_" + day

        def wrap_tlink_inst(inst: instance) -> instance:
            return [
                *inst[:-2],
                format_date(str(inst[-2])),
                "\n".join(textwrap.wrap(inst[-1], width=CLASSIFIER_INST_WIDTH)),
            ]

        wrapped_instances = (wrap_tlink_inst(inst) for inst in instances)
        return sorted(
            wrapped_instances, key=lambda inst: (label_to_hierarchy[inst[1]], inst[3])
        )

    def generate_report(self) -> str:
        source_table = tabulate(
            [source_header, self.source_instance], headers="firstrow"
        )
        tlink_table = (
            tabulate([pred_header, *self.tlink_matches], headers="firstrow")
            if len(self.tlink_matches) > 0
            else ""
        )
        chemo_and_timex_table = (
            tabulate([pred_header, *self.chemo_and_timex_matches], headers="firstrow")
            if len(self.chemo_and_timex_matches) > 0
            else ""
        )
        core_instance = f"\n\n{str(self.error_type)} {str(self.error_cause)} Instance:\n\n{source_table}\n\n"
        tlink_table_str = (
            f"TLINK, Chemo, and Timex Matches from Docker Output:\n\n{tlink_table}\n\n"
            if len(tlink_table) > 0
            else ""
        )
        chemo_and_timex_table_str = (
            f"Chemo and Timex Only Matches from Docker Output:\n\n{chemo_and_timex_table}"
            if len(chemo_and_timex_table) > 0
            else ""
        )
        if len(tlink_table_str) > 0 or len(chemo_and_timex_table_str) > 0:
            suffix = tlink_table_str + chemo_and_timex_table_str
        else:
            suffix = "No matches found in Docker output"
        return core_instance + suffix

    def __str__(self) -> str:
        report = self.generate_report()
        return report

    def get_error_cause(self) -> ErrorCause:
        return self.error_cause

    @staticmethod
    def get_row(full_instance: instance) -> List[str]:
        (chemo, tlink, timex, note, DCT, text_instance) = full_instance
        pandas_row = ["", "", chemo, timex, tlink, text_instance, note, DCT]
        pandas_row += [""] * (len(pandas_columns) - len(pandas_row))
        return pandas_row

    def to_pandas(self) -> pd.DataFrame:
        (
            summarized_chemo,
            summarized_tlink,
            summarized_timex,
            _,  # since we're always using month
        ) = self.source_instance
        first_row = [
            self.error_type,
            self.error_cause,
            summarized_chemo,
            summarized_timex,
            summarized_tlink,
            "SUMMARIZED INSTANCE",
        ]
        first_row += [""] * (len(pandas_columns) - len(first_row))
        remaining_rows = (
            ErrorDebug.get_row(instance)
            for instance in flatten((self.tlink_matches, self.chemo_and_timex_matches))
        )
        return pd.DataFrame.from_records(
            flatten(([first_row], remaining_rows)), columns=pandas_columns
        )


def raw_normalize(time: str) -> str:
    date_only = time.lower().split("t")[0]
    week_fixed = date_only.replace("w", "W")
    return week_fixed


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
    "day": relaxed_day_eval,
    "month": relaxed_month_eval,
    "year": relaxed_year_eval,
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
    fp_events = false_events_by_type(
        ErrorType.FALSE_POSITIVE, error_dict["false_positive"], eval_mode, patient_df
    )
    fn_events = false_events_by_type(
        ErrorType.FALSE_NEGATIVE, error_dict["false_negative"], eval_mode, patient_df
    )
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

    chemo_matches["normed_timex"] = chemo_matches.apply(clean_timex, axis=1)

    chemo_matches_with_potential_timexes = chemo_matches.loc[
        chemo_matches.apply(keep_normalized_timex, axis=1)
    ]

    def row_timex_compatible(pandas_row: pd.Series) -> bool:
        row_timex = cast(str, pandas_row.normed_timex)
        return compatible_time(row_timex, normed_timex, eval_mode)

    if len(chemo_matches_with_potential_timexes) == 0:
        return [], ErrorCause.TIMEX

    timex_chemo_matches = chemo_matches_with_potential_timexes.loc[
        chemo_matches_with_potential_timexes.apply(row_timex_compatible, axis=1)
    ]
    if len(timex_chemo_matches) == 0:
        return [], ErrorCause.TIMEX

    if error_type == ErrorType.FALSE_NEGATIVE:
        return (
            timex_chemo_matches[docker_output_columns].values.tolist(),
            ErrorCause.TLINK,
        )

    def row_tlink_compatible(pandas_row: pd.Series) -> bool:
        row_tlink = cast(str, pandas_row.tlink)
        return row_tlink.lower() == tlink.lower()

    full_matches = timex_chemo_matches.loc[
        timex_chemo_matches.apply(row_tlink_compatible, axis=1)
    ]
    return (
        full_matches[docker_output_columns].values.tolist(),
        ErrorCause.TLINK,
    )


def false_events_by_type(
    error_type: ErrorType,
    false_of_type: TimelineTuples,
    eval_mode: str,
    patient_df: pd.DataFrame,
) -> List[ErrorDebug]:
    def false_of_type_instance(timelines_event: TimelineTuple) -> ErrorDebug:
        chemo_text, tlink, normed_timex = timelines_event
        summarization_preimage, error_cause = preimage_and_cause(
            error_type,
            chemo_text,
            tlink,
            normed_timex,
            eval_mode,
            patient_df,
        )
        return ErrorDebug(
            [[*timelines_event, eval_mode], *summarization_preimage],
            error_type,
            error_cause,
        )

    return [false_of_type_instance(event) for event in false_of_type]


def get_error_cause_count(events: Iterable[ErrorDebug]) -> Counter[ErrorCause]:
    return Counter(ErrorDebug.get_error_cause(event) for event in events)


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
    patient_error_dict: Dict[str, Dict[ErrorType, Counter[ErrorCause]]],
    output_dir: str,
    original_fps: int,
    original_fns: int,
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
            writer.write(f"\n\n\nTop {top_k} Patients for {error_type} errors\n")
            for patient_id, count in patient_counts.most_common(top_k):
                if count > 0:
                    writer.write(f"{patient_id}\t{count}\n")

    with open(output_dir + "/error_metrics.txt", mode="wt") as outfile:
        outfile.write(
            f"{fp_total} Sampled False Positives Out of {original_fps} \t{fn_total} Sampled False Negatives Out of {original_fns}\n\n\n"
        )
        outfile.write("ErrorCause Ranked\n")
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
        "\n\nFalse Positives\n\n" + "\n".join(str(fp) for fp in fp_events)
        if len(fp_events) > 0
        else ""
    )
    fn_str = (
        "\n\nFalse Negatives\n\n" + "\n".join(str(fn) for fn in fn_events)
        if len(fn_events) > 0
        else ""
    )
    with open(fn, mode="wt") as fn_out:
        fn_out.write(fp_str)
        fn_out.write(fn_str)


def write_to_excel(
    patient_to_errors_of_type: Dict[str, List[ErrorDebug]],
    error_type: str,
    output_dir: str,
) -> None:
    for patient_id, errors in patient_to_errors_of_type.items():
        pt_error_type_frames = (error.to_pandas() for error in errors)
        full_dataframe = pd.concat(pt_error_type_frames, ignore_index=False)
        fn = f"{patient_id}_{error_type}.xlsx"
        final_path = os.path.join(output_dir, fn)
        full_dataframe.to_excel(
            final_path
        )

def write_instances_and_summaries(
    docker_tsv: str,
    error_dict: Union[str, DebugDict],
    output_dir: str,
    eval_mode: str,
) -> None:
    docker_df = pd.read_csv(docker_tsv, delimiter="\t")
    patient_error_dict: Dict[str, Dict[ErrorType, Counter[ErrorCause]]] = {}
    patient_to_false_positives: Dict[str, List[ErrorDebug]] = {}
    patient_to_false_negatives: Dict[str, List[ErrorDebug]] = {}
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
        patient_to_false_positives[patient_id] = fp_events
        patient_to_false_negatives[patient_id] = fn_events

    def dict_to_pairs(
        d: Dict[str, List[ErrorDebug]]
    ) -> Sequence[Tuple[str, ErrorDebug]]:
        return [(key, value) for key, values in d.items() for value in values]

    def pairs_to_dict(
        i: Iterable[Tuple[str, ErrorDebug]]
    ) -> Dict[str, List[ErrorDebug]]:
        def group_to_list(g):
            return list(map(itemgetter(1), g))

        return {
            patient_id: group_to_list(group)
            for patient_id, group in groupby(
                sorted(i, key=itemgetter(0)), key=itemgetter(0)
            )
        }

    total_fp = sum(map(len, patient_to_false_positives.values()))
    total_fn = sum(map(len, patient_to_false_negatives.values()))
    fp_sample_size = sample_size(total_fp)
    fn_sample_size = sample_size(total_fn)
    _sampled_fps = random.sample(
        dict_to_pairs(patient_to_false_positives), k=fp_sample_size
    )
    _sampled_fns = random.sample(
        dict_to_pairs(patient_to_false_negatives), k=fn_sample_size
    )
    patient_to_sampled_fp = pairs_to_dict(_sampled_fps)
    patient_to_sampled_fn = pairs_to_dict(_sampled_fns)
    for patient_id in {*patient_to_sampled_fp.keys(), *patient_to_sampled_fn.keys()}:
        sampled_fps = patient_to_sampled_fp.get(patient_id, [])
        sampled_fns = patient_to_sampled_fn.get(patient_id, [])
        write_patient_error_reports(patient_id, sampled_fps, sampled_fns, output_dir)
        patient_error_dict[patient_id] = {
            ErrorType.FALSE_NEGATIVE: get_error_cause_count(sampled_fns),
            ErrorType.FALSE_POSITIVE: get_error_cause_count(sampled_fps),
        }
    write_summaries(
        patient_error_dict,
        output_dir,
        total_fp,
        total_fn,
    )


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
