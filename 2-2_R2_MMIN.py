import math
import os
import numpy as np
import pandas as pd
import time
import argparse
from utils import (load_hyperP_file, scalefree_cal, save_table, get_meta_z_thresholds) 

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Calculate scale-free R² for a MMIN under different PCC thresholds."
    )

    parser.add_argument(
        "-e",
        "--expPPI",
        required=True,
        help="Hypergeometric enrichment results from experimental PPI network"
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
        default=1.0,
        help="Starting meta-z threshold (default: 1.0)"
    )

    parser.add_argument(
        "--threshold_end",
        type=float,
        default=9.0,
        help="Ending meta-z threshold (default: 9.0)"
    )

    parser.add_argument(
        "--threshold_step",
        type=float,
        default=0.5,
        help="Meta-z threshold increment (default: 0.5)"
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
    # Load MMI hypergeometric enrichment results
    # --------------------------------------------------------

    print("Loading experimental PPI MMI results ...")

    exp_mmi = load_hyperP_file(args.expPPI)

    print("Loading co-expressed MMI results ...")

    cor_mmi = load_hyperP_file(args.coexp)

    sample_name1 = os.path.splitext(os.path.basename(args.expPPI))[0]

    sample_name2 = os.path.splitext(os.path.basename(args.coexp))[0]

    print(f"Processing {sample_name1} and {sample_name2} ...")

    # --------------------------------------------------------
    # Evaluate scale-free topology under different thresholds
    # --------------------------------------------------------

    thresholds = get_meta_z_thresholds(args.threshold_start, args.threshold_end, args.threshold_step)

    results = []

    best_threshold = None
    best_r2 = -np.inf

    for threshold in thresholds:

        selected_mmi = []
        for mmi in exp_mmi:
            if mmi not in cor_mmi:
                continue
            meta_z = (exp_mmi[mmi]+ cor_mmi[mmi]) / math.sqrt(2)
            if meta_z >= threshold:
                selected_mmi.append(mmi)
        final_network = selected_mmi 

        try:

            scale_free_exp, r2 = scalefree_cal(final_network)

        except Exception:

            scale_free_exp = np.nan

            r2 = np.nan

        results.append({"threshold": threshold,"MMIN edge count":len(final_network), "R2":r2})

        if not np.isnan(r2):

            if r2 > best_r2:

                best_r2 = r2
                best_threshold = threshold

    if best_threshold is None:
        print("No valid threshold found.")
    else:
        print(f"Best threshold = {best_threshold}, "f"R² = {best_r2:.4f}")

    # --------------------------------------------------------
    # Save results
    # --------------------------------------------------------

    final_df = pd.DataFrame(results)

    print("Saving results...")

    output_dir = os.path.dirname(os.path.abspath(args.output))

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    save_table(final_df,args.output)

if __name__ == "__main__":

    start_time = time.time()

    main()

    elapsed_time = (time.time()- start_time)

    print("Finished.")
    print(f"Elapsed time: "f"{elapsed_time:.2f} seconds")