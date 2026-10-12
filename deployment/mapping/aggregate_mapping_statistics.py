#!/usr/bin/env python3
"""Aggregate per-dataset mapping_statistics.tsv into a study-level table.

Scans ``<study-output>/*/mapping/mapping_statistics.tsv`` (written by the
MappingStatistics process) and concatenates them into a single table with one
row per dataset, used to study mapping recall and recall recovery across a
study.
"""

import argparse
import glob
import os
import sys

import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description="Aggregate per-dataset mapping statistics for a study.")
    parser.add_argument("--study-output", required=True, help="Study output directory containing dataset folders.")
    parser.add_argument("--output", default=None, help="Output TSV (default: <study-output>/mapping_statistics_all.tsv).")
    return parser.parse_args()


def main():
    args = get_args()
    pattern = os.path.join(args.study_output, "*", "mapping", "mapping_statistics.tsv")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"No mapping_statistics.tsv found under {args.study_output}/*/mapping/", file=sys.stderr)
        return 1

    frames = [pd.read_csv(path, sep="\t") for path in files]
    combined = pd.concat(frames, ignore_index=True)
    combined = combined.sort_values("data_set").reset_index(drop=True)

    output = args.output or os.path.join(args.study_output, "mapping_statistics_all.tsv")
    combined.to_csv(output, sep="\t", index=False)
    print(f"Aggregated {len(files)} dataset(s) -> {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
