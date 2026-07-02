import math
import os
import numpy as np
import pandas as pd
import time
import argparse
from utils import (scalefree_cal, load_hyperP_file, get_meta_z_thresholds, save_table)


def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Calculate scale-free R² for CafePMNets under different meta-z thresholds."
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
        "-pp",
        "--ppin",
        required=True,
        help="PPIN file"
    )

    parser.add_argument(
        "-mm",
        "--mmin",
        required=True,
        help="MMIN file"
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
        help="Ending meta-z threshold (default: 10.0)"
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


def load_base_network(ppin_file,mmin_file):
    """
    Load PPIN and MMIN edges.

    Parameters
    ----------
    ppin_file : str

    mmin_file : str

    Returns
    -------
    list
        Existing network edges.
    """

    network_edges = set()

    print("Loading PPIN...")

    with open(ppin_file) as f:

        for line in f:

            network_edges.add(line.rstrip())

    print("Loading MMIN...")

    with open(mmin_file) as f:

        for line in f:

            network_edges.add(line.rstrip())

    print(f"Loaded {len(network_edges):,} base network edges")

    return list(network_edges)

def main():

    args = parse_args()

    # --------------------------------------------------------
    # Load backbone network
    # --------------------------------------------------------

    print("Loading PPIN backbone network ...")

    network_edges = load_base_network(args.ppin,args.mmin)

    # --------------------------------------------------------
    # Load PMI hypergeometric enrichment results
    # --------------------------------------------------------

    print("Loading experimental PPI PMI results ...")

    exp_pmi = load_hyperP_file(args.expPPI)

    print("Loading co-expressed PMI results ...")

    cor_pmi = load_hyperP_file(args.coexp)

    sample_name1 = os.path.splitext(os.path.basename(args.expPPI))[0]

    sample_name2 = os.path.splitext(os.path.basename(args.coexp))[0]

    print(f"Processing {sample_name1} and {sample_name2} ...")

    print("Building meta-z cache...")

    sqrt2 = math.sqrt(2)

    meta_z_cache = {}

    for pmi in exp_pmi:

        if pmi not in cor_pmi:
            continue

        meta_z_cache[pmi] = (exp_pmi[pmi] +cor_pmi[pmi]) / sqrt2

    # --------------------------------------------------------
    # Calculate scale-free R² under different thresholds
    # --------------------------------------------------------

    thresholds = get_meta_z_thresholds(args.threshold_start,args.threshold_end,args.threshold_step)

    results = []

    best_threshold = None
    best_r2 = -np.inf

    for threshold in thresholds:

        selected_pmi = []

        for pmi, meta_z in meta_z_cache.items():

            if meta_z >= threshold:

                selected_pmi.append(pmi)

        # ----------------------------------------------------
        # CafePMNets = PPIN + MMIN + selected PMI
        # ----------------------------------------------------

        final_network = network_edges + selected_pmi

        try:

            scale_free_exp, r2 = scalefree_cal(final_network)

        except Exception:

            scale_free_exp = np.nan
            r2 = np.nan

        results.append(
            {
                "threshold": threshold,
                "CafePMNets edge count": len(final_network),
                "R2": r2
            }
        )

        if not np.isnan(r2):

            if r2 > best_r2:

                best_r2 = r2
                best_threshold = threshold

    # --------------------------------------------------------
    # Report best threshold
    # --------------------------------------------------------

    if best_threshold is None:

        print("No valid threshold found.")

    else:

        print(f"Best threshold = {best_threshold}, R² = {best_r2:.4f}")

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