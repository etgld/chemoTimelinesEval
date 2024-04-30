"""
Convert predictions of run_glue.py to event-timex pairs and summarize to timelines.
"""

import argparse
import json
import os
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Set, Tuple, cast

import pandas as pd

parser = argparse.ArgumentParser(description="")

parser.add_argument("--input_tsv", type=str)
parser.add_argument("--impute_relative", action="store_true")
parser.add_argument("--cancer_type", choices=["ovarian", "breast", "melanoma"])
parser.add_argument("--output_dir", type=str)

CHEMO_MENTIONS = {
    "chemotherapy",
    "chemo",
    "chem",
    "chemo therapy",
    "chemo-radiation",
    "chemo-rt",
    "chemoembolization",
    "chemorad",
    "chemoirradiation",
    "chemort",
    "chemotherapeutic",
    "chemotherap",
    "chemotherapies",
    "chemotherapeutic",
    "chemotherapy's",
    "chemotheray",
    "chemoradiation",
}

label_to_hierarchy = {
    "begins-on": 1,
    "ends-on": 1,
    "contains": 2,
    "contains-1": 2,
    "before": 3,
    "after": 3
}

NORMALIZED_TIMEXES_TO_SKIP = {"Luz 5", "P2000D"}
TimelineDict = Dict[str, List[List[str]]]
TimelineTuples = List[Tuple[str, str, str, str, str]]


def rank_labels(labels: Iterable[str]) -> str:
    label_rankings = [(label, label_to_hierarchy[label]) for label in labels]
    label_rankings = sorted(label_rankings, key=lambda x: x[1])
    return label_rankings[0][0]


def deduplicate(timeline_tuples: TimelineTuples) -> TimelineDict:
    merged_rows: Dict[str, Dict[Tuple[str, str], Set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    chemo_date_map: Dict[str, Dict[Tuple[str, str], List[str]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for row in timeline_tuples:
        source_id, source_text, rel, target_id, target_text = row
        source_text = source_text.lower()
        target_text = target_text.lower()
        # Taking care of timex like this: 2013-10-30T06:47
        if "t" in target_text:
            target_text = target_text.split("t")[0]
        if "p" == target_text[0] and "d" in target_text[-1]:
            continue
        note_id = source_id.strip().split("@")[-2]
        patient_id = note_id.split("_")[0]

        merged_rows[patient_id][(source_text, rel)].add(target_text)

        chemo_date_map[patient_id][(target_text, rel)].append(source_text)

    deduplicated = defaultdict(list)
    for patient, treatments in merged_rows.items():
        one_patient_timelines = []
        chemos_same_day_rel = chemo_date_map[patient]
        for k, v in treatments.items():
            for target in v:
                if k[0] in CHEMO_MENTIONS:
                    has_specific_chemo = False
                    if (target, k[1]) in chemos_same_day_rel:
                        for medication in chemos_same_day_rel[(target, k[1])]:
                            if medication not in CHEMO_MENTIONS:
                                has_specific_chemo = True
                    if not has_specific_chemo:
                        if [k[0], k[1], target] not in one_patient_timelines:
                            one_patient_timelines.append([k[0], k[1], target])
                else:
                    if [k[0], k[1], target] not in one_patient_timelines:
                        one_patient_timelines.append([k[0], k[1], target])
        deduplicated[patient] = one_patient_timelines

    return deduplicated


def conflict_resolution(timelines: TimelineDict) -> TimelineDict:
    resolved_timelines = defaultdict(list)
    for patient, treatments in timelines.items():
        source_target_to_rel = defaultdict(list)
        for tup in treatments:
            s, r, t = tup
            source_target_to_rel[(s, t)].append(r)
        for pair, labels in source_target_to_rel.items():
            if len(labels) > 1:
                more_specific_lbl = rank_labels(labels)
                resolved_timelines[patient].append(
                    [pair[0], more_specific_lbl, pair[1]]
                )
            else:
                resolved_timelines[patient].append([pair[0], labels[0], pair[1]])
    return resolved_timelines


def write_to_output(data: TimelineDict, outfile_name: str) -> None:
    with open(outfile_name, "wt", encoding="utf-8") as fw:
        json.dump(data, fw)


def keep_normalized_timex(pandas_row: pd.Series) -> bool:
    normalized_timex = cast(str, pandas_row.normed_timex)
    components = normalized_timex.split("-")
    valid_prefix = components[0].isnumeric()
    valid_suffix = True
    if len(components) > 1:
        # avoiding seasonal expressions
        # valid_suffix = components[1].lower() not in {"sp", "su", "fa", "wi"}
        valid_suffix = components[1].isnumeric() or (
            "w" in components[1].lower() and components[1][1:].isnumeric()
        )
    return valid_prefix and valid_suffix


def clean_timex(pandas_row: pd.Series) -> str:
    normalized_timex = cast(str, pandas_row.normed_timex)
    components = normalized_timex.split("-")
    valid_prefix = components[0].isnumeric()
    valid_suffix = True
    if len(components) > 1:
        # avoiding seasonal expressions
        # valid_suffix = components[1].lower() not in {"sp", "su", "fa", "wi"}
        valid_suffix = components[1].isnumeric() or (
            "w" in components[1].lower() and components[1][1:].isnumeric()
        )
    if valid_prefix and not valid_suffix:
        return components[0]
    return normalized_timex


def impute_relative_timexes(dataframe: pd.DataFrame) -> pd.DataFrame:
    dataframe.loc[
        dataframe["normed_timex"] == "PRESENT_REF", "normed_timex"
    ] = dataframe.loc[dataframe["normed_timex"] == "PRESENT_REF", "DCT"]
    return dataframe


# not implementing prune by modality and
# prune by polarity since that's currently happening
# upstream to save processing time.
# you can turn that off in
# timeline_delegator.py in the Docker
def convert_docker_output(
    docker_tsv_output_path: str, impute_relative: bool
) -> Tuple[TimelineTuples, Set[str]]:
    docker_output_dataframe = cast(
        pd.DataFrame, pd.read_csv(docker_tsv_output_path, delimiter="\t")
    )

    no_none_tlinks = docker_output_dataframe.loc[
        docker_output_dataframe["tlink"] != "none"
    ]

    normed_timexes_with_tlinks = no_none_tlinks.loc[
        no_none_tlinks["normed_timex"] != "none"
    ]

    if impute_relative:
        normed_timexes_with_tlinks = impute_relative_timexes(normed_timexes_with_tlinks)

    normed_timexes_with_tlinks["normed_timex"] = normed_timexes_with_tlinks.apply(
        clean_timex, axis=1
    )

    acceptable_normed_timexes_with_tlinks = normed_timexes_with_tlinks.loc[
        normed_timexes_with_tlinks.apply(keep_normalized_timex, axis=1)
    ]

    no_discovery_pt_ids = set(
        cast(Iterable[str], docker_output_dataframe["patient_id"])
    ) - set(cast(Iterable[str], acceptable_normed_timexes_with_tlinks["patient_id"]))

    timeline_tups = acceptable_normed_timexes_with_tlinks[
        [
            "chemo_annotation_id",
            "chemo_text",
            "tlink",
            "timex_annotation_id",
            "normed_timex",
        ]
    ].values.tolist()

    return timeline_tups, no_discovery_pt_ids


def convert_resolve_write(
    input_tsv: str,
    cancer_type: Optional[str] = None,
    output_dir: Optional[str] = None,
    impute_relative: bool = False,
) -> Optional[TimelineDict]:
    timelines_tups, no_discovery_pt_ids = convert_docker_output(
        input_tsv, impute_relative
    )

    timelines_deduplicated = deduplicate(timelines_tups)
    resolved_timelines = conflict_resolution(timelines_deduplicated)

    # dumbest hack I've written so far this year but
    for patient_id in no_discovery_pt_ids:
        resolved_timelines[patient_id] = []

    if output_dir is None and cancer_type is None:
        return resolved_timelines

    assert cancer_type is not None and output_dir is not None

    outfile_name = cancer_type + "_timelines.json"
    if impute_relative:
        outfile_name += "_impute_relative"
    write_to_output(
        resolved_timelines,
        os.path.join(output_dir, outfile_name),
    )
    return None  # so mypy doesn't compain about the Optional


def main() -> None:
    args = parser.parse_args()
    convert_resolve_write(
        args.input_tsv, args.cancer_type, args.output_dir, args.impute_relative
    )


if __name__ == "__main__":
    main()
