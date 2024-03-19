import argparse

from .docker_output_error_analysis import write_instances_and_summaries
from .docker_output_to_timeline import TimelineDict, convert_resolve_write
from .eval_timeline import (DebugDict, evaluate_and_log,
                            read_all_patients)

parser = argparse.ArgumentParser(description="Knitting spaghetti")

parser.add_argument("--input_tsv", type=str)
parser.add_argument("--gold_path", required=True, help="A gold annotation json file")
parser.add_argument("--cancer_type", choices=["ovarian", "breast", "melanoma"])
parser.add_argument("--output_dir", type=str)
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
parser.add_argument(
    "--strict", action="store_true", help="do strict eval", type=bool, default=False
)
parser.add_argument(
    "--relaxed_to",
    help="Type 'year' to only evaluate year, 'month' to evaluate year and month, "
    "or 'day' to evaluate year-month-day",
    type=str,
    choices=["day", "month", "year"],
)


def driver(args: argparse.Namespace) -> None:
    pred_timeline: TimelineDict | None = convert_resolve_write(
        input_tsv=args.input_tsv, cancer_type=args.cancer_type
    )
    if pred_timeline is None:
        return
    gold_timeline: TimelineDict = read_all_patients(args.gold_path)
    error_dict: DebugDict | None = evaluate_and_log(
        gold_id_path=args.gold_id_path,
        all_id_path=args.all_id_path,
        gold_timeline=gold_timeline,
        pred_timeline=pred_timeline,
        strict=args.strict,
        relaxed_to=args.relaxed_to,
        debug=args.debug,
        patient_level_metrics=args.patient_level_metrics,
    )
    if args.strict:
        eval_mode = "strict"
    else:
        eval_mode = args.relaxed_to
    if error_dict is not None:
        write_instances_and_summaries(
            docker_tsv=args.input_tsv,
            error_dict=error_dict,
            output_dir=args.output_dir,
            eval_mode=eval_mode,
        )


def main() -> None:
    args = parser.parse_args()
    driver(args)


if __name__ == "__main__":
    main()
