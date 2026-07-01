import os
import time
import argparse
import warnings
from utils import (load_modules,load_exp_ppi,load_pcc,hypergeometric_pvalue,calculate_global_statistics)


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
        help=("Interaction network file. expPPI mode: experimental PPI file； coexpressed mode: Co-expressed gene pairs file.")
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help=("Absolute PCC threshold.Required when mode=correlation.")
    )

    parser.add_argument(
        "-o",
        "--output_file",
        required = True,
        help=("Output file")
    )

    return parser.parse_args()

def build_module_pairs(module_proteins):

    module_pairs = {}

    module_ids = list(module_proteins.keys())

    for i in range(len(module_ids)):
        m1 = module_ids[i]
        for j in range(i + 1, len(module_ids)):
            m2 = module_ids[j]
            unique1 = sorted(set(module_proteins[m1]) - set(module_proteins[m2]))
            unique2 = sorted(set(module_proteins[m2]) - set(module_proteins[m1]))

            if len(unique1) >= 2 and len(unique2) >= 2:
                module_pairs[f"{m1}\t{m2}"] = (unique1, unique2)

    return module_pairs

def is_true_interaction(protein1,protein2,mode,interaction_dict,threshold=None):
    """
    Determine whether a protein pair is a true interaction.
    """

    edge_key = "\t".join(sorted([protein1, protein2]))

    if mode == "expPPI":
        return edge_key in interaction_dict
    elif mode == "coexpressed":
        return (edge_key in interaction_dict and abs(interaction_dict[edge_key]) >= threshold)

    return False

def calculate_module_statistics(module_pairs,interaction_dict,mode,threshold=None):
    """
    Calculate module-module interaction statistics.
    Parameters
    ----------
    module_pairs : dict
        Output from build_module_pairs().

    interaction_dict : dict
        Experimental PPI dictionary or PCC dictionary.

    mode : str
        "expPPI" or "correlation".

    threshold : float, optional
        PCC threshold when mode="correlation".

    Returns
    -------
    tuple
        (
            pair_true_count,
            pair_total_count,
            module1_participating_count,
            module2_participating_count
        )
    """
    pair_true_count = {}
    pair_total_count = {}

    module1_participating_count = {}
    module2_participating_count = {}

    print("Calculating module-module statistics...")

    for module_pair, (group1, group2) in module_pairs.items():

        true_count = 0
        total_count = 0

        participating1 = set()
        participating2 = set()

        for protein1 in group1:
            for protein2 in group2:
                if is_true_interaction(protein1,protein2,mode,interaction_dict,threshold):
                    true_count += 1
                    participating1.add(protein1)
                    participating2.add(protein2)

                total_count += 1

        pair_true_count[module_pair] = true_count
        pair_total_count[module_pair] = total_count
        module1_participating_count[module_pair] = len(participating1)
        module2_participating_count[module_pair] = len(participating2)

    return (pair_true_count,pair_total_count,module1_participating_count,module2_participating_count)

def write_mmi_results(output_file,pair_true_count,pair_total_count,module1_participating_count,module2_participating_count,M,N,mode):

    """
    Write MMI hypergeometric results.

    Parameters
    ----------
    output_file : str

    pair_true_count : dict

    pair_total_count : dict

    module1_participating_count : dict

    module2_participating_count : dict

    M : int
        Global interaction count.

    N : int
        Global possible interaction count.

    mode : str
        expPPI or coexpressed.
    """

    if mode == "expPPI":
        interaction_label = "(m)experimental PPI gene pairs count"
    else:
        interaction_label = "(m)co-expressed gene pairs count"

    header = [
        "module 1",
        "module 2",
        "module 1 participate count",
        "module 2 participate count",
        interaction_label,
        "(n)all gene pairs count",
        "M",
        "N",
        "hyper p-value"
    ]
    
    with open(output_file, "w") as out:

        out.write("\t".join(header) + "\n")

        for module_pair in pair_true_count:

            if (module1_participating_count[module_pair] < 2 or module2_participating_count[module_pair] < 2):
                continue

            m1, m2 = module_pair.split("\t")

            hyper_p = hypergeometric_pvalue(pair_true_count[module_pair],pair_total_count[module_pair],M,N)

            output_line = [
                m1,
                m2,
                str(module1_participating_count[module_pair]),
                str(module2_participating_count[module_pair]),
                str(pair_true_count[module_pair]),
                str(pair_total_count[module_pair]),
                str(M),
                str(N),
                str(hyper_p)
            ]

            out.write("\t".join(output_line) + "\n")

    print(f"Results written to: {output_file}")

def main():
    """
    Calculate module-module interaction (MMI) hypergeometric enrichment statistics.

    Modes
    -----
    expPPI:
        Use experimentally validated PPIs.

    correlation:
        Use co-expressed protein pairs filtered by PCC threshold.
    """

    args = parse_args()

    output_dir = os.path.dirname(args.output_file)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    print("Loading modules...")

    module_proteins = load_modules(args.module)

    module_pairs = build_module_pairs(module_proteins)

    print(f"Module pairs: {len(module_pairs):,}")

    # --------------------------------------------------
    # Load interaction network
    # --------------------------------------------------

    if args.mode == "expPPI":

        print("Loading experimentally validated PPIs...")

        if args.threshold is not None:
            warnings.warn(
                "--threshold is ignored when mode='expPPI'.",
                UserWarning
            )
        interaction_dict = load_exp_ppi(args.network)

        threshold = None

    else:

        print("Loading co-expressed network...")

        if args.threshold is None:
            raise ValueError("--threshold is required when mode=correlation")

        if not (0 <= args.threshold <= 1):
            raise ValueError("Threshold must be within [0,1]")

        interaction_dict = load_pcc(args.network)

        threshold = args.threshold

    # --------------------------------------------------
    # Calculate module statistics
    # --------------------------------------------------

    print("Calculating module-module statistics...")

    (pair_true_count,pair_total_count,module1_participating_count,module2_participating_count
    ) = calculate_module_statistics(
        module_pairs=module_pairs,
        interaction_dict=interaction_dict,
        mode=args.mode,
        threshold=threshold
    )

    # --------------------------------------------------
    # Calculate global statistics
    # --------------------------------------------------

    M, N = calculate_global_statistics(pair_true_count,pair_total_count)

    print(f"Global statistics: M={M:,}, N={N:,}")

    print(f"Writing results to {args.output_file} ...")

    write_mmi_results(
    output_file=args.output_file,
    pair_true_count=pair_true_count,
    pair_total_count=pair_total_count,
    module1_participating_count=module1_participating_count,
    module2_participating_count=module2_participating_count,
    M=M,
    N=N,
    mode=args.mode)

if __name__ == "__main__":

    start_time = time.time()

    main()

    elapsed_time = (time.time()- start_time)

    print("Finished.")
    print(f"Elapsed time: "f"{elapsed_time:.2f} seconds")
