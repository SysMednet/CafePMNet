import sys
import os
import time
import argparse
from utils import (hypergeometric_pvalue, load_modules, calculate_global_statistics, protein_module_key, pair_key, load_exp_ppi, load_pcc)

# ============================================================
# Arguments
# ============================================================
def parse_args():

    parser = argparse.ArgumentParser(
        description=("Calculate module-module interaction hypergeometric enrichment statistics.")
    )

    parser.add_argument(
        "--mode",
        required=True,
        choices=["expPPI", "coexpressed"],
        help=("Reference interaction source. 'expPPI' uses experimentally validated PPIs; 'coexpressed' uses co-expressed network。")
    )

    parser.add_argument(
        "-m",
        "--module",
        required=True,
        help=("Module file.")
    )

    parser.add_argument(
        "-n",
        "--network",
        required=True,
        help=("Interaction network file. expPPI mode: experimental PPI file； coexpressed mode: PCC file.")
    )

    parser.add_argument(
        "-k",
        "--ref",
        required=True,
        help=("Reference genome file.")
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help=("Absolute PCC threshold.Required when mode=coexpressed.")
    )

    parser.add_argument(
        "-o",
        "--output_file",
        required = True,
        help=("Output file")
    )

    return parser.parse_args()

# ============================================================
# Load reference proteome
# ============================================================
def load_reference_genome(filename):

    whole_genome = []
    
    with open(filename,"r") as f:

        for line in f:

            protein = line.rstrip().split("\t")[0]

            whole_genome.append(protein)

    return whole_genome

# ============================================================
# Calculate PMI statistics
# ============================================================
def calculate_pmi_statistics(module_proteins,interaction_dict,reference_proteins,network_type,threshold=None):
    """
    Calculate protein-module interaction statistics.

    Parameters
    ----------
    module_proteins : dict
        {module_id : [protein1, protein2, ...]}

    interaction_dict : dict
        Experimental PPI network or co-expressed network.

    reference_proteins : iterable
        Background reference genome.

    network_type : {"expPPI", "coexpressed"}

    threshold : float, optional
        PCC threshold used when network_type="coexpressed".

    Returns
    -------
    pair_true_count : dict

    pair_total_count : dict

    pair_participating_count : dict
    """

    if network_type not in {"expPPI", "coexpressed"}:
        raise ValueError(f"Unsupported network_type: {network_type}")

    if network_type == "coexpressed" and threshold is None:
        raise ValueError("threshold is required when network_type='coexpressed'")

    pair_true_count = {}
    pair_total_count = {}
    pair_participating_count = {}

    print("Calculating protein-module statistics...")

    for protein in reference_proteins:

        for module, members in module_proteins.items():

            # ------------------------------------------------
            # Skip if protein already belongs to module
            # ------------------------------------------------

            if protein in members:
                continue

            key = protein_module_key(protein, module)

            true_count = 0

            total_count = 0

            participating_members = set()

            for member in members:

                edge = pair_key(protein, member)

                # --------------------------------------------
                # Coexpressed PMI
                # --------------------------------------------

                if network_type == "coexpressed":

                    if (edge in interaction_dict and interaction_dict[edge] >= threshold):

                        true_count += 1

                        participating_members.add(member)

                    total_count += 1

                # --------------------------------------------
                # Experimental PPI PMI
                # --------------------------------------------

                elif network_type == "expPPI":

                    if edge in interaction_dict:

                        true_count += 1

                        participating_members.add(member)

                    total_count += 1

            pair_true_count[key] = true_count

            pair_total_count[key] = total_count

            pair_participating_count[key] = len(participating_members)

    return (pair_true_count,pair_total_count,pair_participating_count)

# ============================================================
# Output
# ============================================================

def write_pmi_results(output_file,pair_true_count,pair_total_count,pair_participating_count,M,N,mode):

    """
    Write PMI hypergeometric results.

    Parameters
    ----------
    output_file : str

    pair_true_count : dict

    pair_total_count : dict

    pair_participating_count : dict

    M : int
        Global interaction count.

    N : int
        Global possible interaction count.

    mode : str
        expPPI or coexpressed.
    """

    if mode == "expPPI":
        interaction_label = "(m)experimental PPI gene pairs count"
    elif mode == "coexpressed":
        interaction_label = "(m)co-expressed gene pairs count"
    else:
        raise ValueError(f"Unsupported mode: {mode}")
    
    header = [
        "protein",
        "module",
        "module participate count",
        interaction_label,
        "(n)all gene pairs count",
        "M",
        "N",
        "hyper p-value"
    ]

    output_dir = os.path.dirname(os.path.abspath(output_file))

    if output_dir:

        os.makedirs(output_dir, exist_ok=True)

    
    with open(output_file, "w") as out:

        out.write("\t".join(header) + "\n")

        for pmi in sorted(pair_true_count):

            # ----------------------------------------
            # Original rule:
            # At least 2 module proteins participate
            # ----------------------------------------

            if pair_participating_count[pmi] < 2:
                continue

            p1, m1 = pmi.split("\t")

            hyper_p = hypergeometric_pvalue(pair_true_count[pmi],pair_total_count[pmi],M,N)

            output_line = [
                p1,
                m1,
                str(pair_participating_count[pmi]),
                str(pair_true_count[pmi]),
                str(pair_total_count[pmi]),
                str(M),
                str(N),
                str(hyper_p)
            ]

            out.write("\t".join(output_line) + "\n")

    print(f"Results written to: {output_file}")

def main():
    """
    Calculate protein-module interaction (PMI)
    hypergeometric enrichment statistics.

    For each protein-module pair:

        1. Count observed interactions (m)
        2. Count possible interactions (n)
        3. Calculate hypergeometric enrichment p-value

    Supports:

        - expPPI
        - coexpressed
    """

    args = parse_args()

    print("Loading modules...")

    module_proteins = load_modules(args.module)

    print(f"Loaded {len(module_proteins):,} modules")

    print("Loading reference proteome...")

    reference_proteins = load_reference_genome(args.ref)

    print(f"Loaded {len(reference_proteins):,} proteins")

    # ----------------------------------------------------
    # Load interaction network
    # ----------------------------------------------------

    if args.mode == "expPPI":

        print("Loading experimental PPIs...")

        interaction_dict = load_exp_ppi(args.network)

        threshold = None

    elif args.mode == "coexpressed":

        print("Loading PCC network...")

        interaction_dict = load_pcc(args.network)

        threshold = args.threshold

    else:

        raise ValueError(f"Unsupported mode: {args.mode}")

    print(f"Loaded {len(interaction_dict):,} interactions")

    # ----------------------------------------------------
    # Calculate PMI statistics
    # ----------------------------------------------------

    (
        pair_true_count,
        pair_total_count,
        pair_participating_count
    ) = calculate_pmi_statistics(
        module_proteins=module_proteins,
        interaction_dict=interaction_dict,
        reference_proteins=reference_proteins,
        network_type=args.mode,
        threshold=threshold
    )

    # ----------------------------------------------------
    # Calculate global M and N
    # ----------------------------------------------------

    M, N = calculate_global_statistics(pair_true_count,pair_total_count)

    print("Global statistics:")

    print(f"M = {M:,}")

    print(f"N = {N:,}")

    # ----------------------------------------------------
    # Write results
    # ----------------------------------------------------

    print("Writing results...")

    write_pmi_results(
        output_file=args.output_file,
        pair_true_count=pair_true_count,
        pair_total_count=pair_total_count,
        pair_participating_count=pair_participating_count,
        M=M,
        N=N,
        mode=args.mode
    )

    print("Finished.")

if __name__ == "__main__":

    start_time = time.time()

    main()

    elapsed_time = (time.time()- start_time)

    print(f"Elapsed time: "f"{elapsed_time:.2f} seconds")