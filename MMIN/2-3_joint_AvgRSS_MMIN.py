import math
import os
import numpy as np
import pandas as pd
import argparse
import time
from utils import (load_modules, load_hyperP_file, load_rss_cache, build_mmi_rss_cache, get_meta_z_thresholds, save_table)

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=("Calculate joint average RSS values under different meta-z score thresholds.")
    )

    parser.add_argument(
        "-k",
        "--experiment",
        required=True,
        help="Experimental PPI file"
    )

    parser.add_argument(
        "-e",
        "--expPPI",
        required=True,
        help="Hypergeometric enrichment results from experimental PPI network"
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
        "-m",
        "--module",
        required=True,
        help="Module file"
    )

    parser.add_argument(
        "-p",
        "--coexp",
        required=True,
        help="Hypergeometric enrichment results from co-expressed gene pairs"
    )

    parser.add_argument(
        "--threshold_start",
        type=float,
        default=1.0
    )

    parser.add_argument(
        "--threshold_end",
        type=float,
        default=9.0
    )

    parser.add_argument(
        "--threshold_step",
        type=float,
        default=0.5
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output file (.txt or .csv)"
    )

    return parser.parse_args()

def main():

    args = parse_args()

    # --------------------------------------------------------
    # Load data
    # --------------------------------------------------------

    print("Loading module information ...")

    module_proteins = load_modules(args.module)

    print("Loading experimental PPI MMI results ...")

    exp_mmi = load_hyperP_file(args.expPPI)

    print("Loading co-expression MMI results ...")

    cor_mmi = load_hyperP_file(args.coexp)

    # --------------------------------------------------------
    # Build RSS cache
    # --------------------------------------------------------

    print("Loading RSS scores ...")

    ppi_rss_cache = load_rss_cache(args.bp,args.cc)

    print("Building MMI RSS cache ...")

    mmi_rss_cache = build_mmi_rss_cache(exp_mmi.keys(),module_proteins,ppi_rss_cache)

    # --------------------------------------------------------
    # Precompute meta-z scores
    # --------------------------------------------------------

    meta_z_cache = {}

    for mmi in exp_mmi:

        if mmi not in cor_mmi:
            continue

        meta_z_cache[mmi] = (exp_mmi[mmi] +cor_mmi[mmi]) / math.sqrt(2)

    # --------------------------------------------------------
    # Evaluate AvgRSS
    # --------------------------------------------------------

    thresholds = get_meta_z_thresholds(args.threshold_start,args.threshold_end,args.threshold_step)

    results = []

    best_threshold = None
    best_avg_rss = -np.inf

    for threshold in thresholds:

        score_sum = 0

        edge_count = 0

        for mmi, meta_z in meta_z_cache.items():

            if meta_z < threshold:
                continue

            if mmi not in mmi_rss_cache:
                continue

            score_sum += mmi_rss_cache[mmi]

            edge_count += 1

        if edge_count == 0:

            avg_rss = np.nan

        else:

            avg_rss = score_sum / edge_count

        results.append(
            {
                "threshold": threshold,
                "MMIN edge count": edge_count,
                "joint average RSS": avg_rss
            }
        )

        if not np.isnan(avg_rss):

            if avg_rss > best_avg_rss:

                best_avg_rss = avg_rss

                best_threshold = threshold

    # --------------------------------------------------------
    # Report best threshold
    # --------------------------------------------------------

    if best_threshold is None:

        print("No valid threshold found.")

    else:

        print(f"Best threshold = {best_threshold}, joint average RSS = {best_avg_rss:.4f}")

    # --------------------------------------------------------
    # Save results
    # --------------------------------------------------------

    final_df = pd.DataFrame(results)

    print("Saving results ...")

    output_dir = os.path.dirname(os.path.abspath(args.output))

    if output_dir:

        os.makedirs(output_dir,exist_ok=True)

    save_table(final_df,args.output)

if __name__ == "__main__":

    start_time = time.time()

    main()

    elapsed_time = (time.time()- start_time)

    print("Finished.")
    print(f"Elapsed time: "f"{elapsed_time:.2f} seconds")
