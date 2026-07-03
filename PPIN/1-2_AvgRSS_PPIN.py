import os
import time
import argparse
import numpy as np
import pandas as pd
from utils import (get_thresholds, load_exp_ppi, load_rss_cache, load_pcc, calculate_avg_rss, save_table)

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=("Calculate average RSS values under different PCC thresholds.")
    )

    parser.add_argument(
        "-e",
        "--expPPI",
        required=True,
        help="Experimentally validated PPI file"
    )

    parser.add_argument(
        "-b",
        "--bp",
        required=True,
        help="BP score file"
    )

    parser.add_argument(
        "-c",
        "--cc",
        required=True,
        help="CC score file"
    )

    parser.add_argument(
        "-p",
        "--pcc",
        required=True,
        help="Co-expressed gene pairs file"
    )

    parser.add_argument(
        "--threshold_start",
        type=float,
        default=0.1
    )

    parser.add_argument(
        "--threshold_end",
        type=float,
        default=0.9
    )

    parser.add_argument(
        "--threshold_step",
        type=float,
        default=0.1
    )

    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output file (.txt or .csv)"
    )

    return parser.parse_args()

def main():
    """
    Calculate average RSS values for threshold-dependent PPINs.

    Experimentally validated PPIs are filtered according to a series of
    user-defined PCC thresholds. For each threshold, the average
    Relative Specificity Similarity (RSS) score of retained PPIs is
    calculated.

    Output
    ------
    A table containing the AvgRSS value calculated for each evaluated PCC threshold.
    """
    args = parse_args()

    thresholds = get_thresholds(args.threshold_start,args.threshold_end,args.threshold_step)

    print("Loading experimentally validated PPIs...")

    ppi_dict = load_exp_ppi(args.expPPI)

    print("Loading BP and CC scores...")

    rss_cache = load_rss_cache(args.bp,args.cc)

    sample_name = os.path.splitext(os.path.basename(args.pcc))[0]

    print("Loading PCC file...")

    print(f"Processing {sample_name} ...")

    pcc_dict = load_pcc(args.pcc)

    result = []

    best_threshold = None
    best_rss = -np.inf

    for threshold in thresholds:

        threshold = round(threshold,1)

        avg_rss, edge_count = calculate_avg_rss(ppi_dict,pcc_dict,rss_cache,threshold)

        result.append({"threshold": threshold, "PPIN edge count": edge_count,"AvgRSS": avg_rss})

        if (not np.isnan(avg_rss)and avg_rss > best_rss):
            best_rss = avg_rss
            best_threshold = threshold

    if best_threshold is not None:
        print(f"best threshold = {best_threshold}, "f"AvgRSS = {best_rss:.4f}")
    else:
        print("No valid threshold found.")

    final_df = pd.DataFrame(result)

    print("Saving results...")

    output_dir = os.path.dirname(args.output_file)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    save_table(final_df,args.output_file)

if __name__ == "__main__":

    start_time = time.time()

    main()

    elapsed_time = (time.time()- start_time)

    print("Finished.")
    print(f"Elapsed time: "f"{elapsed_time:.2f} seconds")
