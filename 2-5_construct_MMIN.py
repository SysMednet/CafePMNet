import math
import os
import time
import argparse
import re
from utils import zscore_from_pvalue

def parse_args():

    parser = argparse.ArgumentParser(
        description=("Construct MMINs using experimentally validated PPIs and co-expressed gene pairs.")
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
        "-o",
        "--output",
        required=True,
        help="Output file"
    )

    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        required=True,
        help="meta-z threshold"
    )

    return parser.parse_args()


def load_hyperP_file2(filename):

    hyperP_file = {}

    with open(filename) as f:

        next(f)

        for line in f:

            fields = line.rstrip().split("\t")

            pair_keys = f"{fields[0]}\t{fields[1]}"

            hyperP_file[pair_keys] = zscore_from_pvalue(float(fields[-1]))

    return hyperP_file

def build_mmin(output_file, true_file, cor_file, threshold):
    """
    Construct a Module-Module Interaction Network (MMIN)
    using experimentally validated PPIs and co-expression
    derived hypergeometric enrichment results.

    Parameters
    ----------
    output_file : str
        Output MMIN file.

    true_file : str
        Experimental PPI hypergeometric enrichment result.

    cor_file : str
        Co-expression hypergeometric enrichment result.

    threshold : float
        Meta-z threshold.

    Returns
    -------
    node_count : int
        Number of modules in the MMIN.

    edge_count : int
        Number of MMIs in the MMIN.
    """

    print("Loading experimental PPI MMI results...")

    true_z = load_hyperP_file2(true_file)

    print("Loading co-expression MMI results...")

    cor_z = load_hyperP_file2(cor_file)

    sqrt2 = math.sqrt(2)

    nodes = set()
    
    output_dir = os.path.dirname(os.path.abspath(output_file))

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    selected_mmis = []

    # --------------------------------------------------------
    # Select MMIs passing meta-z threshold
    # --------------------------------------------------------

    for mmi, true_score in true_z.items():

        cor_score = cor_z.get(mmi)

        if cor_score is None:
            continue

        meta_z = (true_score + cor_score) / sqrt2

        if meta_z < threshold:
            continue
        
        selected_mmis.append(mmi)

    with open(output_file, "w") as out:

        for mmi in selected_mmis:

            module1, module2 =  mmi.split("\t")

            nodes.add(module1)

            nodes.add(module2)

            out.write(f"{module1}\t{module2}\n")

    node_count = len(nodes)

    edge_count = len(selected_mmis)

    print(f"{output_file}: {node_count:,} modules, {edge_count:,} module-module interactions,")

def main():

    args = parse_args()

    print("Building MMIN ...")

    build_mmin(output_file=args.output,true_file=args.expPPI,cor_file=args.coexp,threshold=args.threshold)

    print("Finished.")

if __name__ == "__main__":

    start_time = time.time()

    main()

    elapsed_time = (time.time()- start_time)

    print("Finished.")
    print(f"Elapsed time: "f"{elapsed_time:.2f} seconds")