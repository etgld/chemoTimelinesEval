import argparse
import os
from collections import deque
import pandas as pd

parser = argparse.ArgumentParser(description="")

parser.add_argument("--input_dir", type=str)


# Run with:
#  (312) etg@laptop:~/Repos/chemoTimelinesEval$ python error_analysis_totals.py --input_dir ~/r/DeepPhe_Phase2/treatment_timeline/2023_pilot/test_error_analysis/spreadsheets/Melanoma_Test_Month_Spreadsheets
# Output is:
# ANALYSIS ERROR PRESENT                0
# ANNOTATION ERROR PRESENT              4
# CHEMO ERROR PRESENT                   0
# SUMMARIZATION ERROR PRESENT           0
# TIMEX DETECTION ERROR PRESENT         0
# TIMEX NORMALIZATION ERROR PRESENT     0
# TLINK ERROR PRESENT                  20
# dtype: int64
# total unsummarized preimages: 24
# total summarized predictions: 20
def concatenated_frame(input_dir: str) -> pd.DataFrame:
    frames = deque()
    for root, dirs, files in os.walk(input_dir):
        if "False_Positives.xlsx" in files:
            frames.append(pd.read_excel(os.path.join(root, "False_Positives.xlsx")))
    concatenated_df = pd.concat(frames, ignore_index=False)
    return concatenated_df


def get_counts(input_dir: str) -> None:
    concatenated_df = concatenated_frame(input_dir)
    preimage_error_causes = sorted(
        column_name
        for column_name in concatenated_df.columns
        if column_name.endswith("PRESENT")
    )
    preimages_df = concatenated_df.loc[
        concatenated_df["TEXT"] != "SUMMARIZED INSTANCE"
    ][preimage_error_causes]
    summaries_df = concatenated_df.loc[concatenated_df["TEXT"] == "SUMMARIZED INSTANCE"]
    print(preimages_df.count())
    # the - 1 is since len gets the columns too
    print(f"total unsummarized preimages: {len(preimages_df) - 1}")
    print(f"total summarized predictions: {len(summaries_df) - 1}")


def main():
    args = parser.parse_args()
    get_counts(args.input_dir)


if __name__ == "__main__":
    main()
