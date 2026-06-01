# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header

from tqdm import tqdm
from pathlib import Path
from collections import Counter
from openbabel import openbabel, pybel
from matplotlib.ticker import MultipleLocator

import matplotlib.pyplot as plt
import matplotlib.lines as ln
import scipy.stats as stats
import pandas as pd
import numpy as np

INPUT_FILES = {
        "tmQM-RDF": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev"),
        "tmQM-RDF-1Ksel": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "data", "v2025dev", "%s", "1k_selection.csv"), # %s is tm_mode (earlyTM/lateTM)
        "details": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", "v2025dev", "details.csv"),
        "ligands_misc_info": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", "v60k", "ligands_misc_info.csv"),
        "ligands_fingerprints": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", "v60k", "ligands_fingerprints.csv"),
        "scaffolds": os.path.join(ROOT_DIR, "computational", "reconstruction", "intermediate", "scaffold", "%s", "test"), # %s is tm_mode
        "recon_results": os.path.join(ROOT_DIR, "computational", "reconstruction", "results", "%s-latmod-s_10_12", "ranks", "%s", "ranks.txt") # %s, %s is tm_mode and filter (n/hd/c)
        #"recon_results": os.path.join(ROOT_DIR, "experiments", "frequency_based_reconstruction", "outputs", "%s-latmod-s_10_12", "ranks", "%s", "ranks.txt") # %s, %s is tm_mode and filter (n/hd/c)
    }

OUTPUT_FILES = {
        "results": os.path.join(ROOT_DIR, "computational", "reconstruction", "results", "figures")
    }

DATASET_SHORT_ID = "latmod" # Ligand ATtachment MODes

# Range of sizes of interes (extremes included)
SIZE_RANGE = [10, 12]

SETTINGS = {
    "metrics": ["cossim", "delgsim"],
    "variants": {
        "cossim": ["features_proxy", "features_semantic"],
        "delgsim": ["fixed_weights", "estimated_weights"]
    },
    "aggregations": ["max", "set_median"],
    "filters": ["n", "c", "hd"],
    "linestyles": ["solid", "dashed", "dotted"], 
    "filter_labels": ["No filter", "Same hap./den. orders", "Same charge"],
    "k": [1, 5, 10]
}

# %% Utility functions ----
def ta_hitATk(k, Ni_1, ni):
    """
    Tie-aware hit@k (see Saha et al., 2022)
    """
    if k <= Ni_1:
        return 0
    if k >= Ni_1 + ni:
        return 1
    return (k - Ni_1)/ni

def compute_ligand_frquencies(tm_mode):
    selection = pd.read_csv(INPUT_FILES["tmQM-RDF-1Ksel"] % tm_mode, index_col = 0)
    selection = list(selection.loc[selection["partition"] == "train", "TMC"])

    tmc_details = pd.read_csv(INPUT_FILES["details"], index_col = 0)
    mcs = list(set(tmc_details.loc[selection, "core"]))

    frequencies = {
        mc: [] for mc in ["*"] + mcs
    }

    for tmc in tqdm(selection):
        mc, ligands = tmc_details.loc[tmc, :]

        for l in ligands.strip().split(" "):
            frequencies["*"] += [l]
            frequencies[mc] += [l]

    return {mc: Counter(ligs) for mc, ligs in frequencies.items()}

def retrieve_ligand_filtering_info(ligands):
    misc_info = pd.read_csv(INPUT_FILES["ligands_misc_info"], sep = ";")
    descriptors = pd.read_csv(INPUT_FILES["ligands_fingerprints"], sep = ";")
    
    # Ex: [[1,2,3], [1]] -> [1, 3]
    list_to_hd_order = lambda y: sorted([len(z) for z in eval(y)])
    
    return {
            l : {
                    "n": None,
                    "c": descriptors.loc[descriptors["name"] == l, "charge"].iloc[0],
                    "hd": list_to_hd_order(misc_info.loc[misc_info["name"] == l, "smiles_metal_bond_node_idx_groups"].iloc[0])
                }
            for l in ligands
        }


def compute_binding_modalities(ligands):
    misc_info = pd.read_csv(INPUT_FILES["ligands_misc_info"], sep = ";")

    def binding_modality(lig_id):
        smiles = misc_info.loc[misc_info["name"] == lig_id, "smiles"].item()

        bindings = eval(misc_info.loc[misc_info["name"] == lig_id, "smiles_metal_bond_node_idx_groups"].item())
        mol = pybel.readstring("smi", smiles)
        atoms = [openbabel.GetSymbol(a.atomicnum) for a in mol.atoms]

        processed_bindings = []
        for bs in bindings:
            processed_bindings += [[atoms[i] for i in bs]]

        return "-".join(sorted(
                        [
                            ''.join(sorted([a.capitalize() for a in bs])) for bs in processed_bindings
                        ]
                    ))

    return {lig_id: binding_modality(lig_id) for lig_id in ligands}

def retrieve_bn_results(tm_mode, filter, blacklist = []):
    topk = dict()

    with open(INPUT_FILES["recon_results"] % (tm_mode, filter), "r") as f:
        rankings_raw = eval(f.read())
    
    #print(rankings_raw.keys())

    for metric in SETTINGS["metrics"]:
        topk.setdefault(metric, dict())
        for variant in SETTINGS["variants"][metric]:
            topk[metric].setdefault(variant, dict())
            for aggregation in SETTINGS["aggregations"]:
                ranks = [
                        rank_data for tmc, rank_data in rankings_raw[f"{metric}/{variant}/{aggregation}"].items()
                        if tmc not in blacklist
                    ]

                maxk = max([a + b for a, b in ranks])

                n = len(ranks)
                topk[metric][variant][aggregation] = [0]
                for k in range(1, maxk + 1):
                    topk[metric][variant][aggregation] += [sum([ta_hitATk(k, *rank) for rank in ranks])/n]

    return topk

def reconstruct_with_frequent_ligands(frequencies, tm_mode, filtering_info, ligands_modalities, count_upper_bound = np.inf):
    # Assemble dataset tag
    if SIZE_RANGE[1] > SIZE_RANGE[0]:
        size_tag = f"s_{SIZE_RANGE[0]}_{SIZE_RANGE[1]}"
    else:
        size_tag = f"s_{SIZE_RANGE[0]}"
    
    dataset_tag = f"{tm_mode}-{DATASET_SHORT_ID}-{size_tag}"

    tmc_details = pd.read_csv(INPUT_FILES["details"], index_col = 0)

    # Retrieve target scaffold
    target_scaffolds = [x for x in os.listdir(INPUT_FILES["scaffolds"] % tm_mode) if x.endswith(".ttl")]

    # Prepare rank mappings
    def rank(ligand, mc, filter):
        admissible = [
            l for l in frequencies[mc].keys()
            if filtering_info[l][filter] == filtering_info[ligand][filter]
        ]

        admissible_counts = {
            l: frequencies[mc][l] for l in admissible
        }

        scounts = sorted(set(admissible_counts.values()), reverse = True)

        c = {
            cnt: len(set([
                ligands_modalities[l] for l, x in admissible_counts.items() if x == cnt
            ]))
            for cnt in scounts
        }
        #c = Counter(admissible_counts.values())

        counts_to_level = {
            cnt: i + 1 for i, cnt in enumerate(scounts)
        }

        N = {
            i + 1: sum([c[x] for x in scounts if x >= scounts[i]])
            for i in range(len(c))
        }
        N[0] = 0

        correct_level = counts_to_level[frequencies[mc][ligand]]

        return (N[correct_level - 1], c[frequencies[mc][ligand]])

    # Compute rankings
    topk_table = dict()

    baseline_topk = dict()
    bn_topk = dict()

    for filter in SETTINGS["filters"]:
        blacklist = []
        tmc_rankings = dict()
        mcs = []
        for target in tqdm(target_scaffolds):
            tmc, gt = target.replace(".ttl", "").split("__")
            mc = tmc_details.loc[tmc, "core"]
            mcs += [mc]
            key = mc if gt in frequencies[mc] else "*"
            if gt in frequencies[key] and frequencies[key][gt] <= count_upper_bound:
                tmc_rankings[tmc] = rank(gt, key, filter)
            else:
                blacklist += [tmc]
                continue
            
        n = len(tmc_rankings)
        maxk = max([a + b for a, b in tmc_rankings.values()])

        bn_data = retrieve_bn_results(tm_mode, filter, blacklist)
        baseline_data = [0] + [
            sum([ta_hitATk(k, *rank) for rank in tmc_rankings.values()])/n
            for k in range(1, maxk + 1)
        ]

        baseline_topk[filter] = baseline_data
        bn_topk[filter] = bn_data

        topk_table.setdefault(filter, dict())
        topk_table[filter]["baseline"] = {
            k: baseline_topk[filter][k]
            for k in SETTINGS["k"]    
        }
        for metric in SETTINGS["metrics"]:
            topk_table[filter].setdefault(metric, dict())
            for variant in SETTINGS["variants"][metric]:
                topk_table[filter][metric].setdefault(variant, dict())   
                for aggregation in SETTINGS["aggregations"]:
                    topk_table[filter][metric][variant][aggregation] = {
                                        k: bn_data[metric][variant][aggregation][k]
                                        for k in SETTINGS["k"]
                                    }

    # Prepare plot data
    plot_data = []
    mins = []
    for filter, linestyle, label in zip(SETTINGS["filters"], SETTINGS["linestyles"], SETTINGS["filter_labels"]):
        y = bn_topk[filter]["cossim"]["features_proxy"]["max"]
        x = baseline_topk[filter]
        
        if len(x) != len(y):
            to_pad = x if len(x) < len(y) else y
            to_pad += [1.0] * abs(len(x) - len(y))

        x = np.array(x)
        y = np.array(y)
        #mins += [np.min(y - x)]

        topk_v_topk = {
            "args": [x, y],
            "kwargs": {
                "linestyle": linestyle,
                "color": "red" if tm_mode == "lateTM" else "blue",
                "label": label
            }
        }

        delta_v_k = {
            "args": [np.arange(len(x)), y - x],
            "kwargs": {
                "linestyle": linestyle,
                "color": "red" if tm_mode == "lateTM" else "blue",
                "label": label
            }
        }

        plot_data += [(topk_v_topk, delta_v_k)]

    return topk_table, plot_data

def plot_all(plot_data, count_upper_bound):
    # Plotting (BN Top-k acc. vs Baseline Top-k acc.)
    fig, ax = plt.subplots() # topk vs topk
    fig2, ax2 = plt.subplots() # (topk - topk) vs k
    mins = []

    axs = [ax, ax2]

    for selIdx in [0, 1]:
        for fltrIdx in range(len(SETTINGS["filters"])):
            for plotIdx in [0, 1]:
                axs[plotIdx].plot(*plot_data[selIdx][fltrIdx][plotIdx]["args"], **plot_data[selIdx][fltrIdx][plotIdx]["kwargs"])

    ax.plot([0,1], [0, 1], color = "grey", linewidth = 0.5, transform = ax.transAxes)
    ax2.plot([0, 1], [0, 0], color = "grey", linewidth = 0.5, transform = ax2.get_yaxis_transform())

    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax2.xaxis.set_minor_locator(MultipleLocator(10))
    ax2.xaxis.set_major_locator(MultipleLocator(50))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.1))

    def dummy_line(color, linestyle, label):
        return ln.Line2D([], [], color = color, linestyle = linestyle, label = label)
    dummy_handles = [
        dummy_line("blue", "solid", "earlyTM"), 
        dummy_line("black", "solid", "No filter"), 
        dummy_line("black", "dashed", "Same hap./den. orders"),
        dummy_line("black", "dotted", "Same charge"), 
        dummy_line("red", "solid", "lateTM"), 
        dummy_line("#0000", "solid", ""), 
        dummy_line("#0000", "solid", ""),  
        dummy_line("#0000", "solid", "")
    ]
    ax.legend(handles = dummy_handles, ncols = 2, loc = "lower right", fontsize = "x-small")
    ax2.legend(handles = dummy_handles, ncols = 2, loc = "upper right", fontsize = "x-small")

    #ax.legend(loc = "lower right", fontsize = "x-small")
    #ax2.legend(loc = "upper right", fontsize = "x-small")

    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    ax2.set_ylim(top = 1.02)

    ax.set_xlabel("Top-$k$ frequent ligand (top-$k$ accuracy)")
    ax.set_ylabel("BN-based scoring (top-$k$ accuracy)")
    ax2.set_xlabel("$k$")
    ax2.set_ylabel("$\\Delta_{\\mathrm{BN}-\\mathrm{baseline}}$Top$-k$ accuracy")

    ax.set_box_aspect(1)

    if not os.path.exists(OUTPUT_FILES["results"]):
        os.makedirs(OUTPUT_FILES["results"])

    fig.savefig(os.path.join(OUTPUT_FILES["results"], f"topkUB{count_upper_bound}.png"), format = "png", dpi = 200)
    fig.savefig(os.path.join(OUTPUT_FILES["results"], f"topkUB{count_upper_bound}.pdf"), format = "pdf")
    fig2.savefig(os.path.join(OUTPUT_FILES["results"], f"deltaUB{count_upper_bound}.png"), format = "png", dpi = 200)
    fig2.savefig(os.path.join(OUTPUT_FILES["results"], f"deltaUB{count_upper_bound}.pdf"), format = "pdf")

def nice_table(topk_table, add_to_clipboard = False):
    """
    Prints topk_table (plot_top_k_accuracy's output) in a latex-friendly format.
    
    - Arguments:
        - topk_table: a dictionary with keys [filter] -> ["baseline"] / [metric][variant][aggregation]
        - add_to_clipboard: should the resulting latex code be copied to the system's clipboard? Default: False
    """
    
    n_ks = len(SETTINGS["k"])
    
    py2tex = {
            "n": "No filter",
            "hd": "Hapticity/denticity order",
            "c": "Charge",
            "baseline": "Top$-k$ frequent ligand",
            "cossim": "$s_{cos}$",
            "delgsim": "$s_{DELG}$",
            "features_semantic": "$s_{cos;\,s}$",
            "features_proxy": "$s_{cos;\,p}$",
            "fixed_weights": "$s_{DELG;\,n}$",
            "estimated_weights": "$s_{DELG;\,l}$",
            "max": "$a_{max}$",
            "set_median": "$a_{median}$"
        }
    
    filters = ["n", "hd", "c"]
    variants = {
            "cossim": ["features_proxy", "features_semantic"],
            "delgsim": ["fixed_weights", "estimated_weights"]
        }
    
    entries = []
    
    for metric in ["baseline", "cossim", "delgsim"]:
        metric_entry = py2tex[metric]
        
        if metric == "baseline":
            row = [metric_entry, "\t", "\t"]
            for fltr in filters:
                try:
                    k_values = topk_table[fltr][metric].values()
                except Exception:
                    k_values = [np.nan] * n_ks
                
                row += [""] + ["%.3f" % acc for acc in k_values]
            
            entries += [row]
            continue

        for variant in variants[metric]:
            variant_entry = py2tex[variant]
            
            for aggregation in ["max", "set_median"]:
                aggregation_entry = py2tex[aggregation]
                
                row = [metric_entry, variant_entry, aggregation_entry]
                for fltr in filters:
                    try:
                        k_values = topk_table[fltr][metric][variant][aggregation].values()
                    except Exception:
                        k_values = [np.nan] * n_ks
                    
                    row += [""] + ["%.3f" % acc for acc in k_values]
                
                entries += [row]
                
                metric_entry = "\t"
                variant_entry = "\t"
    
    entries = np.array(entries, dtype = "<U50")
    
    print("RAW TABLE")
    print("\n".join(";".join(r) for r in entries))

    means = []
    
    for i in range(entries.shape[1]):
        try:
            x = entries[:,i].astype(float)
            
            means += [np.nanmean(x)]
            
            for j in np.where(x == np.nanmax(x))[0]:
                entries[j, i] = "$\\mathbf{%s}$" % entries[j, i].item()
                x[j] = -1
            
            for j in np.where(x == np.nanmax(x))[0]:
                entries[j, i] = "$\\emph{%s}$" % entries[j, i].item()
                
            for j in range(len(x)):
                if np.isnan(x[j]):
                    entries[j, i] = "--"
        except:
            continue
    
    n_headers = 3 + 4*n_ks
    
    table = ["\\toprule"]
    table += ["Metric & Variant & Aggregation & & \\multicolumn{%d}{c}{Filter}\\\\\n\\cmidrule(lr){5-%d}" % (n_headers - 4, n_headers)]
    table += [" & " * 3 + " & ".join(["& \\multicolumn{%d}{c}{%s}" % (n_ks, py2tex[f]) for f in filters]) + "\\\\"]
    table += ["".join(["\\cmidrule(lr){%d-%d}" % (5 + (n_ks + 1)*i, 5 + (n_ks + 1)*(i + 1) - 2) for i in range(n_ks)])]
    table += [" & " * 3 + " & ".join(["& \\multicolumn{%d}{c}{Top$-k$ accuracy}" % n_ks] * 3) + "\\\\"]
    table += ["".join(["\\cmidrule(lr){%d-%d}" % (5 + (n_ks + 1)*i, 5 + (n_ks + 1)*(i + 1) - 2) for i in range(n_ks)])]
    table += [" & " * 3 + " & ".join(([""] + ["$k = %d$" % k for k in SETTINGS["k"]])*3) + "\\\\"]
    table += ["\\midrule\n\\midrule"]
    
    table += [" & ".join(row) + "\\\\" for row in entries]
    
    table += ["\\bottomrule"]
    
    print("\n".join(table))
    
    if add_to_clipboard:
        try:
            import pyperclip
            
            pyperclip.copy("\n".join(table))
        except ModuleNotFoundError:
            print("--Warning---")
            print("pyperclip module not found! Unable to copy latex output to clipboard!")
    
    print("\n\nMeans:")
    print(" || ".join("%d) %.3f" % (i, j) for i, j in enumerate(means)))

# %% Main statement ----
if __name__ == "__main__":

    plot_data_inf = []
    plot_data_1 = []

    for tm_mode in ["lateTM", "earlyTM"]:
        frequencies = compute_ligand_frquencies(tm_mode)
        filtering_info = retrieve_ligand_filtering_info(frequencies["*"].keys())
        ligands_modalities = compute_binding_modalities(frequencies["*"].keys())

        topk_table, plot_data = reconstruct_with_frequent_ligands(frequencies, tm_mode, filtering_info, ligands_modalities)
        plot_data_inf += [plot_data]

        _, plot_data = reconstruct_with_frequent_ligands(frequencies, tm_mode, filtering_info, ligands_modalities, 5)
        plot_data_1 += [plot_data]

        print("")
        print("")

        nice_table(topk_table)
        print(list((k, topk_table[k]["baseline"]) for k in topk_table))

        print("")
        print("")

    plot_all(plot_data_inf, np.inf)
    plot_all(plot_data_1, 5)