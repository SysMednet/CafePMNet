# utils.py
import os
import math
import numpy as np
import pandas as pd
import networkx as nx
from scipy.optimize import curve_fit
from scipy.stats import hypergeom

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

def load_exp_ppi(filename):
    """
    Loading experimental PPI files
    Format required: Protein1\tProtein2 
    (delimited by tab)
    """
    ppi_dict = {}

    with open(filename) as f:

        next(f)

        for line in f:

            fields = line.rstrip().split("\t")

            p1 = fields[0]
            p2 = fields[1]

            edge = "\t".join(sorted([p1, p2]))

            ppi_dict[edge] = 1

    return ppi_dict

def load_pcc(filename):
    """
    Loading Pearson correlation coefficient (PCC) files of each cancer
    Format required: Protein1\tProtein2\tPCC
    (delimited by tab)
    """
    pcc_dict = {}

    with open(filename) as f:

        for line in f:

            fields = line.rstrip().split("\t")

            p1 = fields[0]
            p2 = fields[1]

            edge = "\t".join(sorted([p1, p2]))

            pcc_dict[edge] = float(fields[-1])

    return pcc_dict

def load_rss_cache(bp_file, cc_file):
    """
    Precompute RSS scores for all protein pairs
    with both BP and CC scores.

    Returns
    -------
    dict
        protein pair -> RSS
    """

    bp_dict = {}
    cc_dict = {}

    with open(bp_file) as f:

        for line in f:

            fields = line.rstrip().split("\t")

            edge = pair_key(fields[0],fields[1])

            bp_dict[edge] = float(fields[-1])

    with open(cc_file) as f:

        for line in f:

            fields = line.rstrip().split("\t")

            edge = pair_key(fields[0],fields[1])

            cc_dict[edge] = float(fields[-1])

    rss_cache = {}

    common_edges = set(bp_dict) & set(cc_dict)

    for edge in common_edges:

        rss_cache[edge] = rss_score(edge,bp_dict,cc_dict)

    return rss_cache

def rss_score(pair_key_str, bp_dict, cc_dict):
    """
    RSS = sqrt(BP × CC)
    """
    return math.sqrt(bp_dict[pair_key_str] * cc_dict[pair_key_str])

def get_thresholds(start,end,step):
    """
    Validate PCC thresholds.

    Parameters
    ----------
    thresholds : list of float

    Returns
    -------
    list of float
        Sorted unique threshold values.
    """

    if not (0 <= start <= 1):
        raise ValueError(f"threshold_start must be within [0,1], got {start}")

    if not (0 <= end <= 1):
        raise ValueError(f"threshold_end must be within [0,1], got {end}")

    if start > end:
        raise ValueError("threshold_start must be <= threshold_end")

    if step <= 0:
        raise ValueError("threshold_step must be > 0")

    thresholds = np.arange(start,end + step / 2,step)

    return [round(float(t), 6) for t in thresholds]

def pair_key(a, b):
    """
    Generate order-independent pair key.
    """

    return "\t".join(sorted([a, b]))

def scalefree_cal(edge_l):
    """
    Estimate the scale-free property of a network from its degree distribution.

    The input edge list is converted into a NetworkX graph. The degree
    distribution P(k) is calculated and transformed into log-log space:

        log(P(k)) = a + b * log(k)

    A linear regression is then fitted using scipy.optimize.curve_fit().
    The returned scaling exponent is:

        gamma = -b

    where b is the fitted slope in log-log space.

    The coefficient of determination (R²) is also calculated to quantify
    how well the degree distribution follows a power-law relationship.
    """
    def func(X,a,b):
        return a+b*X
    net = nx.parse_edgelist(edge_l,delimiter = '\t',nodetype = str)
    degree_d = {}
    for _ , d in nx.degree(net) :
       degree_d[d] = degree_d.get(d, 0) + 1
    x = np.array(sorted(degree_d.keys()))
    y = np.array([degree_d[d] for d in x])
    x , y = np.log(x),  np.log(y)
    popt,_ = curve_fit(func,x,y)
    _ , b = popt[0] , popt[1]
    residuals = y - func(x, *popt)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    return[-b,r_squared]

def calculate_avg_rss(ppi_dict,pcc_dict,rss_cache,threshold):
    """
    Calculate average RSS of PPIs passing
    a PCC threshold.
    """
    score_sum = 0
    count = 0

    for ppi in ppi_dict:

        if (ppi not in pcc_dict or abs(pcc_dict[ppi]) < threshold or ppi not in rss_cache):
            continue

        score_sum += rss_cache[ppi]

        count += 1

    if count == 0:
        return np.nan, count

    return score_sum / count, count

def calculate_geometric_mean(df_r2,df_rss, network_type):
    """
    Calculate geometric mean.

    GM = sqrt(R² × AvgRSS)
    """

    if network_type == "PPIN":
        rss_column = "AvgRSS"
    else:
        rss_column = "joint average RSS"
    
    result_df = pd.merge(df_r2[["threshold", "R2"]],df_rss[["threshold", rss_column]],on="threshold",how="inner")

    if len(result_df) != len(df_r2) or len(result_df) != len(df_rss):
        raise ValueError("R² and AvgRSS files contain different threshold values.")

    result_df["GM"] = np.sqrt(result_df["R2"] *result_df[rss_column])

    return result_df

def prepare_plot_dataframe(df_geo,df_rss,df_r2):
    """
    Convert dataframes into long format for boxplot visualization.
    """
    df_geo_plot = df_geo.copy()
    df_geo_plot["method"] = "Geometric Mean"

    geo_melt = pd.melt(df_geo_plot,id_vars=["cancer", "method"],var_name="threshold",value_name="score")

    df_rss_plot = df_rss.copy()
    df_rss_plot["method"] = "RSS"

    rss_melt = pd.melt(df_rss_plot,id_vars=["cancer", "method"],var_name="threshold",value_name="score")

    df_r2_plot = df_r2.copy()
    df_r2_plot["method"] = "R²"

    r2_melt = pd.melt(df_r2_plot,id_vars=["cancer", "method"],var_name="threshold",value_name="score")

    return pd.concat([geo_melt, rss_melt, r2_melt],ignore_index=True)

def plot_threshold_selection(result_df, best_threshold, output_file, network_type, threshold_name = "Threshold", show_best_line = True):
    """
    Plot R², AvgRSS and Geometric Mean across thresholds.

    Parameters
    ----------
    result_df : pandas.DataFrame

    best_threshold : float

    output_file : str

    threshold_name : str
        Label for x-axis.

    show_best_line : bool, default=True
        Whether to draw the vertical line indicating
        the optimal threshold.
    """
    if network_type == "PPIN":
        rss_column = "AvgRSS"
    else:
        rss_column = "joint average RSS"

    plt.figure(figsize=(8, 6), dpi=300)

    plt.plot(
        result_df["threshold"],
        result_df["R2"],
        marker="o",
        linewidth=2,
        label="R²"
    )

    plt.plot(
        result_df["threshold"],
        result_df[rss_column],
        marker="o",
        linewidth=2,
        label="AvgRSS"
    )

    plt.plot(
        result_df["threshold"],
        result_df["GM"],
        marker="o",
        linewidth=2,
        label="Geometric Mean"
    )

    if show_best_line:

        plt.axvline(
            x=best_threshold,
            linestyle="--",
            linewidth=1.5,
            label=f"Best threshold = {best_threshold}"
        )

    plt.xlabel(threshold_name)

    plt.ylabel("Score")

    plt.legend()

    plt.tight_layout()

    plt.savefig(output_file,dpi=300,bbox_inches="tight")

    plt.close()

def find_best_threshold(result_df):
    """
    Identify threshold with highest mean geometric mean score.
    """
    idx = result_df["GM"].idxmax()

    best_threshold = result_df.loc[idx,"threshold"]

    best_gm = result_df.loc[idx,"GM"]

    return best_threshold, best_gm

def save_table(df, output_file):
    """
    Save dataframe as TSV (.txt) or CSV (.csv)
    according to file extension.
    """

    ext = os.path.splitext(output_file)[1].lower()

    if ext == ".txt":
        sep = "\t"

    elif ext == ".csv":
        sep = ","

    else:
        raise ValueError(f"Unsupported output format: {ext}. ""Use .txt or .csv")

    df.to_csv(output_file,sep=sep,index=False)

def hypergeometric_pvalue(x, n, M, N):
    """
    x : observed successes
    n : sample size
    M : total successes in population
    N : total population size
    """
    return hypergeom.sf(x - 1, N, M, n)

def calculate_global_statistics(pair_true_count,pair_total_count):
    """
    Calculate global M and N.

    Returns
    -------
    tuple
        (M, N)
    """

    M = sum(pair_true_count.values())

    N = sum(pair_total_count.values())

    return M, N

def load_modules(module_file):
    """
    Load module definitions.

    Parameters
    ----------
    module_file : str
        Module definition file.

    Returns
    -------
    dict{module_id : [protein1, protein2, ...]}
    """

    module_proteins = {}

    with open(module_file, "r") as f:

        for line in f:

            fields = line.rstrip().split("\t")

            module_id = fields[0]

            module_proteins[module_id] = fields[2].split()

    return module_proteins
