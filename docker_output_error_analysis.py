import pandas as pd
import argparse
import json
from typing import Tuple, List, Dict, Union

parser = argparse.ArgumentParser(description="")

parser.add_argument("--docker_tsv", type=str)
parser.add_argument("--error_json", type=str)
parser.add_argument("--output_dir", type=str)

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
    def __init__(self, event_instances: List[instance]) -> None:
        self.instances = event_instances

    def generate_report(self) -> str:
        return ""

    def __str__(self) -> str:
        report = self.generate_report()
        return report


class FNDebug:
    def __init__(self, event_instances: List[instance]) -> None:
        self.instances = event_instances

    def generate_report(self) -> str:
        return ""

    def __str__(self) -> str:
        report = self.generate_report()
        return report


DebugEvent = Union[FPDebug, FNDebug]


def collect_error_events(
    patient_id: str, error_dict: Dict[str, List[List[str]]], docker_df: pd.DataFrame
) -> List[DebugEvent]:
    return []


def write_patient_error_reports(patient_id: str, error_events: List[DebugEvent]):
    pass


def write_instances(docker_tsv: str, error_json: str, output_dtr: str):
    docker_df = pd.read_csv(docker_tsv, sep="\t")
    with open(error_json, mode="rt") as json_f:
        patient_to_errors = json.load(json_f)
    for patient_id, error_dict in patient_to_errors.items():
        error_events = collect_error_events(patient_id, error_dict, docker_df)
        write_patient_error_reports(patient_id, error_events)


def main():
    args = parser.parse_args()
    write_instances(args.docker_tsv, args.error_json, args.output_dir)


if __name__ == "__main__":
    main()
