import argparse
import math
import os
import numpy as np
import pandas as pd
import time

from utils import (zscore_from_pvalue, pair_key,save_table,load_rss_cache, load_modules, build_mmi_rss_cache, get_meta_z_thresholds, protein_module_key)

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=("Calculate joint average RSS values under different meta-z score thresholds.")
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

def load_pmi_hyperP_file(filename):

    pmi_dict = {}

    with open(filename) as f:

        next(f)

        for line in f:

            fields = line.rstrip().split("\t")

            key = protein_module_key(
                fields[0],
                fields[1]
            )

            pmi_dict[key] = zscore_from_pvalue(
                float(fields[-1])
            )

    return pmi_dict

def calculate_ppin_rss(ppin_file,ppi_rss_cache):
    """
    Calculate RSS contribution from PPIN.

    Returns
    -------
    score_sum : float

    edge_count : int
    """

    score_sum = 0

    edge_count = 0

    with open(ppin_file) as f:

        for line in f:

            edge = line.rstrip()

            if edge not in ppi_rss_cache:
                continue

            score_sum += ppi_rss_cache[edge]

            edge_count += 1

    return score_sum, edge_count

def calculate_pmi_rss(protein,module,module_proteins,ppi_rss_cache):
    """
    Calculate RSS of one PMI edge.

    Parameters
    ----------
    protein : str

    module : str

    module_proteins : dict

    ppi_rss_cache : dict

    Returns
    -------
    average_rss : float or None

    valid : bool
    """
    participating = set()

    score_sum = 0

    pair_count = 0

    for module_protein in module_proteins[module]:
        edge = pair_key(protein,module_protein)
        if edge not in ppi_rss_cache:
            continue
        
        participating.add(module_protein)

        score_sum += ppi_rss_cache[edge]

        pair_count += 1

    if len(participating) < 2:
        return None, False

    if pair_count == 0:
        return None, False

    return score_sum / pair_count, True

def build_pmi_rss_cache(exp_pmi,module_proteins,ppi_rss_cache):
    """
    Precompute RSS scores for all valid PMIs.

    Parameters
    ----------
    exp_pmi : dict

    module_proteins : dict

    ppi_rss_cache : dict

    Returns
    -------
    dict
        PMI -> AvgRSS
    """

    pmi_rss_cache = {}

    print("Building PMI RSS cache...")

    for pmi in exp_pmi:

        protein, module = pmi.split("\t")

        avg_rss, valid = calculate_pmi_rss(protein,module,module_proteins,ppi_rss_cache)

        if valid:

            pmi_rss_cache[pmi] = avg_rss

    return pmi_rss_cache

def main():

    args = parse_args()

    # --------------------------------------------------------
    # Load modules
    # --------------------------------------------------------

    print("Loading modules...")

    module_proteins = load_modules(args.module)

    print(f"Loaded {len(module_proteins):,} modules")

    # --------------------------------------------------------
    # Load RSS scores
    # --------------------------------------------------------

    print("Loading RSS cache...")

    ppi_rss_cache = load_rss_cache(bp_file=args.bp,cc_file=args.cc)

    print(f"Cached RSS scores for {len(ppi_rss_cache):,} PPIs")

    # --------------------------------------------------------
    # Calculate fixed PPIN RSS
    # --------------------------------------------------------

    print("Calculating PPIN RSS...")

    score_ppin, count_ppin = calculate_ppin_rss(args.ppin,ppi_rss_cache)

    print(f"PPIN:{count_ppin:,} edges, RSS sum = {score_ppin:.4f}")

    # --------------------------------------------------------
    # Load PMI hypergeometric enrichment results
    # --------------------------------------------------------

    print("Loading experimental PPI PMI results...")

    exp_pmi = load_pmi_hyperP_file(args.expPPI)

    print("Loading co-expressed PMI results...")

    cor_pmi = load_pmi_hyperP_file(args.coexp)

    # --------------------------------------------------------
    # Load MMIN network
    # --------------------------------------------------------

    print("Loading MMIN network...")

    mmin_edges = []

    with open(args.mmin) as f:

        for line in f:

            mmin_edges.append(line.rstrip())

    print(f"Loaded {len(mmin_edges):,} MMIN edges")

    # --------------------------------------------------------
    # Build MMI RSS cache
    # --------------------------------------------------------

    print("Building MMI RSS cache...")

    exp_mmi = {}

    for mmi in mmin_edges:

        exp_mmi[mmi] = 1

    mmi_rss_cache = build_mmi_rss_cache(exp_mmi,module_proteins,ppi_rss_cache)

    score_mmin = sum(mmi_rss_cache.values())

    count_mmin = len(mmi_rss_cache)

    # --------------------------------------------------------
    # Build PMI RSS cache
    # --------------------------------------------------------

    print("Building PMI RSS cache...")

    pmi_rss_cache = build_pmi_rss_cache(exp_pmi,module_proteins,ppi_rss_cache)

    # --------------------------------------------------------
    # Evaluate AvgRSS under different thresholds
    # --------------------------------------------------------

    thresholds = get_meta_z_thresholds(args.threshold_start,args.threshold_end,args.threshold_step)

    results = []

    best_threshold = None

    best_rss = -np.inf

    for threshold in thresholds:

        score_pmi = 0

        count_pmi = 0

        for pmi in exp_pmi:

            if pmi not in cor_pmi:
                continue

            meta_z = (exp_pmi[pmi] + cor_pmi[pmi]) / math.sqrt(2)

            if (meta_z < threshold) or (pmi not in pmi_rss_cache):
                continue

            score_pmi += pmi_rss_cache[pmi]

            count_pmi += 1

        total_score = (score_ppin + score_mmin + score_pmi)

        total_count = (count_ppin + count_mmin + count_pmi)

        if total_count == 0:

            avg_rss = np.nan

        else:

            avg_rss = total_score / total_count

        results.append(
            {
                "threshold": threshold,
                "CafePMNets edge count": total_count,
                "joint average RSS": avg_rss
            }
        )

        if not np.isnan(avg_rss):

            if avg_rss > best_rss:

                best_rss = avg_rss

                best_threshold = threshold

    # --------------------------------------------------------
    # Report best threshold
    # --------------------------------------------------------

    if best_threshold is None:

        print("No valid threshold found.")

    else:

        print(f"Best threshold = {best_threshold}, joint average RSS = {best_rss:.4f}")

    # --------------------------------------------------------
    # Save results
    # --------------------------------------------------------

    final_df = pd.DataFrame(results)

    print("Saving results...")

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
