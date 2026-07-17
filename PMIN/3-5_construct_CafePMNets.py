import math
import os
import time
import argparse
from utils import (zscore_from_pvalue)

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Construct CafePMNets using meta-z score threshold."
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
        help="PPIN network"
    )

    parser.add_argument(
        "-mm",
        "--mmin",
        required=True,
        help="MMIN network"
    )

    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        required=True,
        help="Meta-z score threshold"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output CafePMNets file"
    )

    return parser.parse_args()

def load_pmi_zscore(filename):

    """
    Load PMI hypergeometric enrichment results
    and convert p-values to Z-scores.

    Parameters
    ----------
    filename : str

    Returns
    -------
    dict
        PMI -> Z-score
    """

    pmi_zscore = {}

    with open(filename) as f:

        next(f)

        for line in f:

            fields = line.rstrip().split("\t")

            interaction =  f"{fields[0]}\t{fields[1]}"

            pmi_zscore[interaction] = zscore_from_pvalue(float(fields[-1]))

    return pmi_zscore

def build_pmin(ppin_file, mmin_file,exp_pmi,cor_pmi,threshold,output_file):

    """
    Construct CafePMNet.

    CafePMNet consists of:

        PPIN
      + MMIN
      + PMI with meta-z >= threshold

    Parameters
    ----------
    ppin_file : str

    mmin_file : str

    exp_pmi : dict
        Experimental PMI Z-scores.

    cor_pmi : dict
        Co-expression PMI Z-scores.

    threshold : float

    output_file : str
    """

    print("Building meta-z cache...")

    sqrt2 = math.sqrt(2)

    meta_z_cache = {}

    for pmi in exp_pmi:

        if pmi not in cor_pmi:
            continue

        meta_z_cache[pmi] = (exp_pmi[pmi] + cor_pmi[pmi]) / sqrt2

    edge_count = 0

    with open(output_file, "w") as out:

        # ----------------------------------------------------
        # PPIN
        # ----------------------------------------------------

        with open(ppin_file) as f:

            for line in f:

                out.write(line)

                edge_count += 1

        print(f"PPIN edges : {edge_count:,}")

        # ----------------------------------------------------
        # MMIN
        # ----------------------------------------------------

        with open(mmin_file) as f:

            for line in f:

                out.write(line)

                edge_count += 1

        print(f"PPIN + MMIN edges : {edge_count:,}")

        # ----------------------------------------------------
        # PMI
        # ----------------------------------------------------

        selected_count = 0

        for pmi, meta_z in meta_z_cache.items():

            if meta_z < threshold:
                continue

            out.write(pmi + "\n")

            edge_count += 1

            selected_count += 1

    print(f"Selected PMI edges : {selected_count:,}")

    print(f"Total network edges : {edge_count:,}")

def main():

    args = parse_args()

    # --------------------------------------------------------
    # Load experimental PMI
    # --------------------------------------------------------

    print("Loading experimental PMI results...")

    exp_pmi = load_pmi_zscore(args.expPPI)

    # --------------------------------------------------------
    # Load co-expression PMI
    # --------------------------------------------------------

    print("Loading co-expressed PMI results...")

    cor_pmi = load_pmi_zscore(args.coexp)

    sample_name1 = os.path.splitext(os.path.basename(args.expPPI))[0]

    sample_name2 = os.path.splitext(os.path.basename(args.coexp))[0]

    print(f"Processing {sample_name1} and {sample_name2} ...")

    # --------------------------------------------------------
    # Build PMIN
    # --------------------------------------------------------

    build_pmin(
        ppin_file=args.ppin,
        mmin_file=args.mmin,
        exp_pmi=exp_pmi,
        cor_pmi=cor_pmi,
        threshold=args.threshold,
        output_file=args.output
    )

    print("Saving network...")

    output_dir = os.path.dirname(os.path.abspath(args.output))

    if output_dir:

        os.makedirs(output_dir, exist_ok=True)

    print("Done.")

if __name__ == "__main__":

    start_time = time.time()

    main()

    elapsed_time = (time.time() - start_time)

    print(f"Elapsed time: "f"{elapsed_time:.2f} seconds")
