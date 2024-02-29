import argparse
import json
from datetime import datetime
from typing import Dict, List, Tuple, Union

import dateutil.parser
import pandas as pd
from tabulate import tabulate

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


class FPDebug:
    def __init__(self, instances: List[instance]) -> None:
        self.source_instance = instances[0]
        self.pred_instances = instances[1:]

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


class FNDebug:
    def __init__(self, instances: List[instance]) -> None:
        self.source_instance = instances[0]
        self.pred_instances = instances[1:]

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


DebugEvent = Union[FPDebug, FNDebug]


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
) -> Tuple[List[FPDebug], List[FNDebug]]:
    patient_df = docker_df[docker_df["patient_id"].isin([patient_id])]
    fp_events = collect_fp_events(error_dict["false_positive"], eval_mode, patient_df)
    fn_events = collect_fn_events(error_dict["false_negative"], eval_mode, patient_df)
    return fp_events, fn_events


def collect_fp_events(
    false_positives: List[List[str]], eval_mode: str, patient_df: pd.DataFrame
) -> List[FPDebug]:
    def fp_instance(timelines_event: List[str]) -> FPDebug:
        # since the resulting tlink is predetermined
        # given the conflict resolution rules
        chemo_text, tlink, normed_timex = timelines_event
        summarization_preimage = patient_df.loc[
            (patient_df["chemo_text"].apply(lambda t: compatible_chemos(t, chemo_text)))
            & (patient_df["tlink"] == tlink.lower())
            & (
                patient_df["normed_timex"].apply(
                    lambda t: compatible_time(t, normed_timex, eval_mode)
                )
            )
        ][["chemo_text", "tlink", "normed_timex", "note_name"]].values.tolist()
        return FPDebug([[*timelines_event, eval_mode], *summarization_preimage])

    return [fp_instance(event) for event in false_positives]


def collect_fn_events(
    false_negatives: List[List[str]], eval_mode: str, patient_df: pd.DataFrame
) -> List[FNDebug]:
    def fn_instance(timelines_event: List[str]) -> FNDebug:
        # since the resulting tlink is predetermined
        # given the conflict resolution rules
        chemo_text, tlink, normed_timex = timelines_event
        summarization_preimage = patient_df.loc[
            (patient_df["chemo_text"].apply(lambda t: compatible_chemos(t, chemo_text)))
            # need to turn off none-filtering here
            # & (patient_df["tlink"] != "none")
            & (
                patient_df["normed_timex"].apply(
                    lambda t: compatible_time(t, normed_timex, eval_mode)
                )
            )
        ][["chemo_text", "tlink", "normed_timex", "note_name"]].values.tolist()
        return FNDebug([[*timelines_event, eval_mode], *summarization_preimage])

    return [fn_instance(event) for event in false_negatives]


def write_patient_error_reports(
    patient_id: str, fp_events: List[FPDebug], fn_events: List[FNDebug], output_dir: str
):
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
