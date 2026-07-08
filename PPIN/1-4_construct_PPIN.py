# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 03:22:14 2022

@author: junmao
"""
import time
import os
from utils import (pair_key, load_pcc)

import argparse

def parse_args():

    parser = argparse.ArgumentParser(
        description=("Construct PPINs using experimentally validated PPIs and co-expressed gene pairs.")
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
        "-o",
        "--output_file",
        required=True,
        help="Output file (.txt or .csv)"
    )

    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        required=True,
        help="PCC threshold within [0,1]"
    )

    return parser.parse_args()

def load_ppi_network(ppi_file):
    """
    Load experimentally validated PPIs.
    """

    ppi_set = set()

    with open(ppi_file, "r") as f:

        next(f)

        for line in f:

            fields = line.rstrip().split("\t")

            ppi_set.add(pair_key(fields[0], fields[1]))

    return ppi_set

def build_ppin(sample, ppi_set, pcc_file, output_file, threshold):
    """
    Construct PPIN for a single cancer.
    """
    pcc_dict = load_pcc(pcc_file)

    edge_count = 0

    node_set = set()

    output_dir = os.path.dirname(os.path.abspath(output_file))

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(output_file, "w") as out:

        for ppi in ppi_set:

            if ppi not in pcc_dict:
                continue

            if abs(pcc_dict[ppi]) < threshold:
                continue

            out.write(ppi + "\n")

            edge_count += 1 

            p1, p2 = ppi.split("\t")

            node_set.add(p1)
            node_set.add(p2)

    node_count = len(node_set)

    print(f"{sample}: "f"{edge_count:,} edges retained, "f"{node_count:,} nodes retained")


def main():
    """
    Construct a PPIN from experimentally validated
    PPIs and PCC values using a user-defined PCC threshold.

    PPIs are retained if:

        1. Present in experimentally validated PPI dataset
        2. Present in cancer-specific PCC file
        3. |PCC| >= user-defined threshold

    Outputs one PPIN edge list per cancer type.
    """
    args = parse_args()

    threshold = args.threshold

    if not (0 <= threshold <= 1):
        raise ValueError("Threshold must be within [0,1]")

    print("Loading experimentally validated PPIs...")

    ppi_set = load_ppi_network(args.expPPI)

    print(f"Total PPIs: {len(ppi_set):,}")

    sample = os.path.splitext(os.path.basename(args.pcc))[0]

    print(f"Processing {sample}")

    output_dir = os.path.dirname(args.output_file)

    if output_dir:

        os.makedirs(output_dir, exist_ok=True)

    build_ppin(sample=sample,ppi_set=ppi_set,pcc_file=args.pcc,output_file=args.output_file, threshold=threshold)

if __name__ == "__main__":
    start_time = time.time()

    main()

    elapsed_time = time.time() - start_time

    print("Finished.")
    print(f"Elapsed time: {elapsed_time:.2f} seconds")

