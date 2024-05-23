import argparse
import os
from collections import deque
import pandas as pd

parser = argparse.ArgumentParser(description="")

parser.add_argument("--input_dir", type=str)


def concatenated_frame(input_dir: str) -> pd.DataFrame:
    frames = deque()
    for root, dirs, files in os.walk(input_dir):
        if "False_Positives.xlsx" in files:
            frames.append(pd.read_excel(os.path.join(root, "False_Positives.xlsx")))
    concatenated_df = pd.concat(frames, ignore_index=False)
    return concatenated_df


def get_counts(input_dir: str) -> None:
    concatenated_df = concatenated_frame(input_dir)


def main():
    args = parser.parse_args()
    get_counts(args.input_dir)


if __name__ == "__main__":
    main()
