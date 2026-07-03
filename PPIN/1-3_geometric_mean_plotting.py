import os
import time
import argparse
import pandas as pd
from utils import (calculate_geometric_mean, plot_threshold_selection, find_best_threshold,save_table)

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=("Calculate geometric mean scores from R² and AvgRSS results.")
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
        help="Average RSS result file"
    )

    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output result file. A plot with the same basename will also be generated."
    )

    return parser.parse_args()

def main():
    """
    Calculate geometric mean scores for threshold selection.

    Geometric mean is calculated from:

        GM = sqrt(R² × AvgRSS)

    The optimal threshold is determined
    as the threshold with the highest
    geometric mean score
    """

    args = parse_args()

    base = os.path.splitext(args.output_file)[0]

    plot_file = f"{base}.png"

    print("Loading R² table...")

    df_r2 = pd.read_csv(args.r2,sep="\t")

    print("Loading AvgRSS table...")

    df_rss = pd.read_csv(args.rss,sep="\t")

    # --------------------------------------------------
    # Check threshold columns
    # --------------------------------------------------

    required_r2 = {"threshold", "R2"}
    required_rss = {"threshold", "AvgRSS"}

    if not required_r2.issubset(df_r2.columns):
        raise ValueError("R² file must contain columns: PCC threshold, R2")

    if not required_rss.issubset(df_rss.columns):
        raise ValueError("AvgRSS file must contain columns: PCC threshold, AvgRSS")

    # -----------------------------------------------------------------------------------------------------
    # Calculate geometric mean
    # Please note that the optimal threshold may not be the best choice due to the few edges in PPIN
    # Please refer to Supplementary Table
    # -----------------------------------------------------------------------------------------------------

    print("Calculating geometric mean...")

    result_df  = calculate_geometric_mean(df_r2,df_rss, "PPIN")

    save_table(result_df,args.output_file)

    best_threshold, best_gm = (find_best_threshold(result_df))

    print("\nThreshold selection summary")
    print(f"Best threshold = {best_threshold}")
    print(f"Geometric mean = {best_gm:.4f}")

    edge_row = df_r2[df_r2["threshold"] == best_threshold]

    if not edge_row.empty:

        best_edge_count = int(edge_row["PPIN edge count"].iloc[0])

        print(f"Edge count = {best_edge_count:,}")

    else:

        print("Edge count = N/A")

    print("Please note that the amount of edge may not sufficient if best threshold is selected")

    print("Generating threshold plot...")

    plot_threshold_selection(result_df,best_threshold,plot_file,"PPIN", "PCC threshold", show_best_line = False)

    print("Done.")

if __name__ == "__main__":

    start_time = time.time()

    main()

    elapsed_time = (time.time() - start_time)

    print(f"Elapsed time: "f"{elapsed_time:.2f} seconds")
