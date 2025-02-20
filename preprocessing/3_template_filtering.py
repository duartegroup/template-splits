"""
Filters USPTO and Pistachio by the number of template occurrences.

Note: Both datasets must have their templates extracted between
      2_cleaning.py and this script, with templates saved in a
      column named 'template' in data/raw/2_{dataset}.csv. This
      is done via the LocalTemplate extractor algorithm in
      LocalRetro: https://github.com/kaist-amsg/LocalRetro.
"""

import pandas as pd
import os
import argparse

def load_df(dataset_name):
    df = pd.read_csv(f"data/raw/2_{dataset_name}.csv", index_col="dataset_id")
    if 'template' not in df.columns:
        raise FileNotFoundError("No template column found in csv. Must be extracted via LocalTemplate.")
    return df


def filter_templates(df, threshold=5):
    print(len(df))
    g = df.groupby('template')
    filtered_df = g.filter(lambda x: len(x) > threshold)
    print(filtered_df["template"].value_counts())
    print(len(filtered_df))
    return filtered_df


def save_df(df, dataset_name):
    data_dir = f"data/processed/{dataset_name}_retro"
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    df.to_csv(f"{data_dir}/{dataset_name}_retro.csv")
    print(f"Saved to {data_dir}/{dataset_name}_retro.csv")


def filter_df(dataset_name, threshold=5):
    df = load_df(dataset_name)
    df = filter_templates(df)
    df = save_df(df, dataset_name)


def main():
    filter_df("uspto", 5)
    if not args["uspto_only"]:
        filter_df("pistachio", 20)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--uspto_only", action="store_true", help="Only run USPTO code.")
    args = vars(parser.parse_args())
    main(args)
    print("done!")