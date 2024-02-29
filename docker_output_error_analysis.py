import pandas as pd
import argparse
import json
from typing import Tuple, List, Dict, Union
from datetime import datetime
import dateutil.parser

parser = argparse.ArgumentParser(description="")

parser.add_argument("--docker_tsv", type=str)
parser.add_argument("--error_json", type=str)
parser.add_argument("--output_dir", type=str)
parser.add_argument(
    "--eval_mode",
    choices=["strict", "day", "month", "year"],
)
# chemo, rel, timex, filename
instance = Tuple[str, str, str, str]
label_to_hierarchy = {
    "begins-on": 1,
    "ends-on": 1,
    "contains": 2,
    "contains-1": 2,
    "before": 3,
}


class FPDebug:
    def __init__(self, pred_instances: List[instance]) -> None:
        self.pred_instances = pred_instances

    def generate_report(self) -> str:
        def inst2str(inst: instance) -> str:
            return "\t".join(inst)

        return "\n".join(map(inst2str, self.pred_instances))

    def __str__(self) -> str:
        report = self.generate_report()
        return report


class FNDebug:
    def __init__(self, pred_instances: List[instance]) -> None:
        self.pred_instances = pred_instances

    def generate_report(self) -> str:
        def inst2str(inst: instance) -> str:
            return "\t".join(inst)

        return "\n".join(map(inst2str, self.pred_instances))

    def __str__(self) -> str:
        report = self.generate_report()
        return report


DebugEvent = Union[FPDebug, FNDebug]


def raw_normalize(time: str) -> datetime:
    time = time.replace("w", "W")
    if "W" in time:
        return datetime.strptime(time + "-1", "%Y-W%W-%w")
    return dateutil.parser.parse(time)


def strict_eval(time1: str, time2: str) -> bool:
    return time1 == time2


def relaxed_day_eval(time1: str, time2: str) -> bool:
    norm_time1 = raw_normalize(time1)
    norm_time2 = raw_normalize(time2)
    # as far as I can tell this is how it works
    # in eval_timeline.py
    return norm_time1 == norm_time2


def relaxed_month_eval(time1: str, time2: str) -> bool:
    norm_time1 = raw_normalize(time1)
    norm_time2 = raw_normalize(time2)
    return (norm_time1.year, norm_time1.month) == (norm_time2.year, norm_time2.month)


def relaxed_year_eval(time1: str, time2: str) -> bool:
    norm_time1 = raw_normalize(time1)
    norm_time2 = raw_normalize(time2)
    return norm_time1.year == norm_time2.year


def compatible_time(time1: str, time2: str, eval_mode: str) -> bool:
    if eval_mode in mode_to_eval:
        return mode_to_eval[eval_mode](time1, time2)
    else:
        ValueError(f"You picked a non existant mode: {eval_mode}")
        return False


mode_to_eval = {
    "strict": strict_eval,
    "relaxed_day": relaxed_day_eval,
    "relaxed_month": relaxed_month_eval,
    "relaxed_year": relaxed_year_eval,
}


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
            (patient_df["chemo_text"].lower() == chemo_text.lower())
            & (patient_df["tlink"].lower() == tlink.lower())
            & (
                patient_df["normed_timex"].apply(
                    lambda t: compatible_time(t, normed_timex, eval_mode)
                )
            )
        ][["chemo_text", "normed_timex", "tlink", "note_name"]].values.tolist()
        return FPDebug(summarization_preimage)

    return [fp_instance(event) for event in false_positives]


def collect_fn_events(
    false_negatives: List[List[str]], eval_mode: str, patient_df: pd.DataFrame
) -> List[FNDebug]:
    def fn_instance(timelines_event: List[str]) -> FNDebug:
        # since the resulting tlink is predetermined
        # given the conflict resolution rules
        chemo_text, _, normed_timex = timelines_event
        summarization_preimage = patient_df.loc[
            (patient_df["chemo_text"] == chemo_text)
            # need to turn off none-filtering here
            # & (patient_df["tlink"] != "none")
            & (
                patient_df["normed_timex"].apply(
                    lambda t: compatible_time(t, normed_timex, eval_mode)
                )
            )
        ][["chemo_text", "normed_timex", "tlink", "note_name"]].values.tolist()
        
        return FNDebug(summarization_preimage)

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
    write_instances(args.docker_tsv, args.error_json, args.output_dir, args.eval_mode)


if __name__ == "__main__":
    main()
