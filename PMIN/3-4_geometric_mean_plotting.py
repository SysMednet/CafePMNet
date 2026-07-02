import os
import numpy as np
import pandas as pd
import argparse
import time
from utils import (plot_threshold_selection, calculate_geometric_mean, save_table, find_best_threshold)

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=("Calculate geometric mean scores from R² and joint average RSS results.")
    )

    parser.add_argument(
        "-r",
        "--r2",
        required=True,
        help="R² result file"
    )

    parser.add_argument(
        "-s",
        "--rss",
        required=True,
        help="joint average RSS result file"
    )

    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output result file. A plot with the same basename will also be generated."
    )

    return parser.parse_args()

def main():

    args = parse_args()

    base = os.path.splitext(args.output_file)[0]

    plot_file = f"{base}.png"

    print("Loading R² table...")

    df_r2 = pd.read_csv(args.r2,sep="\t")

    print("Loading joint average RSS table...")

    df_rss = pd.read_csv(args.rss,sep="\t")

    required_r2 = {"threshold","R2"}

    required_rss = {"threshold","joint average RSS"}

    if not required_r2.issubset(df_r2.columns):

        raise ValueError("R² file must contain columns: threshold, R2")

    if not required_rss.issubset(df_rss.columns):

        raise ValueError("joint avgRSS file must contain columns: threshold, joint average RSS")

    print("Calculating geometric mean...")

    result_df = calculate_geometric_mean(df_r2,df_rss, "CafePMNets")

    save_table(result_df,args.output_file)

    best_threshold, best_gm = (find_best_threshold(result_df))

    print("\nThreshold selection summary")

    print(f"Best threshold = {best_threshold}")

    print(f"Geometric mean = {best_gm:.4f}")

    print("Generating threshold plot...")

    plot_threshold_selection(result_df,best_threshold,plot_file, "meta-z threshold")

    print("Done.")

if __name__ == "__main__":

    start_time = time.time()

    main()

    elapsed_time = (time.time() - start_time)

    print(f"Elapsed time: "f"{elapsed_time:.2f} seconds")