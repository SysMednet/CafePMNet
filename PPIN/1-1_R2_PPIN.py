import os
import time
import argparse
import numpy as np
import pandas as pd
from utils import (get_thresholds, load_exp_ppi, load_pcc, scalefree_cal, save_table)

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Calculate scale-free R² for a PPIN under different PCC thresholds."
    )

    parser.add_argument(
        "-e",
        "--expPPI",
        required=True,
        help="Experimentally validated PPI file"
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
    Calculate scale-free network fit (R²) under different PCC thresholds.

    Experimentally validated PPIs are filtered using a series of PCC
    thresholds to construct threshold-dependent PPINs. For each threshold,
    the coefficient of determination (R²) of the log-log degree
    distribution is calculated to evaluate network scale-free topology.

    Output
    ------
    A table containing the R² value calculated for each evaluated PCC threshold.
    """
    args = parse_args()

    thresholds = get_thresholds(args.threshold_start,args.threshold_end,args.threshold_step)
    
    print("Loading experimentally validated PPIs...")

    ppi_dict = load_exp_ppi(args.expPPI)

    print("Loading PCC file...")

    sample_name = os.path.splitext(os.path.basename(args.pcc))[0]

    print(f"Processing {sample_name} ...")

    pcc_dict = load_pcc(args.pcc)

    result = []

    best_threshold = None
    best_r2 = -np.inf

    for threshold in thresholds:

        threshold = round(threshold,1)

        selected_ppi = []

        for ppi in ppi_dict:

            if ppi not in pcc_dict:
                continue

            if abs(pcc_dict[ppi]) < threshold:
                continue

            selected_ppi.append(ppi)

        edge_count = len(selected_ppi)

        try:
            scale_free_exp, r2 = scalefree_cal(selected_ppi)
        except Exception:
            r2 = np.nan

        result.append({"threshold": threshold,"PPIN edge count":edge_count, "R2":r2})

        if not np.isnan(r2):

            if r2 > best_r2:

                best_r2 = r2
                best_threshold = threshold

    if best_threshold is None:
        print("No valid threshold found.")
    else:
        print(f"Best threshold = {best_threshold}, "f"R² = {best_r2:.4f}")

    final_df = pd.DataFrame(result)

    print("Saving results...")
    
    save_table(final_df,args.output_file)

if __name__ == "__main__":

    start_time = time.time()

    main()

    elapsed_time = (time.time()- start_time)

    print("Finished.")
    print(f"Elapsed time: "f"{elapsed_time:.2f} seconds")