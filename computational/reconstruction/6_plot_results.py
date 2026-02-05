# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header

import numpy as np
import matplotlib.pyplot as plt

INPUT_FILES = {
    "recon_results": os.path.join(ROOT_DIR, "computational", "reconstruction", "results", "%s", "ranks")
}

OUTPUT_FILES = {
    "top_k": os.path.join(ROOT_DIR, "computational", "reconstruction", "results", "%s", "figures", "top_k_accuracy")
}

TM_MODE = "earlyTM"
DATASET_SHORT_ID = "latmod" # Ligand ATtachment MODes

# Range of sizes of interes (extremes included)
SIZE_RANGE = [10, 12]

SETTINGS = {
    "rank_file_name": "ranks.txt",
    "plt_labels": {
            "cossim": {
                    "features_proxy": {
                            "max": "proxy - max",
                            "set_median": "proxy - set median",
                        },
                    "features_semantic": {
                            "max": "semantic - max",
                            "set_median": "semantic - set median",
                        }
                },
            "delgsim": {
                    "fixed_weights": {
                            "max": "fixed - max",
                            "set_median": "fixed - set median",
                        },
                    "estimated_weights": {
                            "max": "learned - max",
                            "set_median": "learned - set median",
                        }
                }     
        },
    "plt_colors": ["blue", "red"], # Each metric is a plot, each feature/weight type is a color
    "plt_styles": {"max": "-", "set_median": "--"}, # line style differentiate between aggregation functions
    "k": [1, 5, 10],
    "max_k": 90
}

# %% Utility functions

def plot_top_k_accuracy():
    # Assemble dataset tag
    if SIZE_RANGE[1] > SIZE_RANGE[0]:
        size_tag = f"s_{SIZE_RANGE[0]}_{SIZE_RANGE[1]}"
    else:
        size_tag = f"s_{SIZE_RANGE[0]}"
    
    dataset_tag = f"{TM_MODE}-{DATASET_SHORT_ID}-{size_tag}"
    
    topk_table = {}
    
    for fltr in os.listdir(INPUT_FILES["recon_results"] % dataset_tag):
        if not os.path.isdir(os.path.join(INPUT_FILES["recon_results"] % dataset_tag, fltr)):
            continue
        
        print(f"\n\nFILTER: {fltr}")
        
        # Retrieve rankings and parse into a dictionary with keys [metric][variant][aggregation]
        # Original rankings are stored in a dictionary where keys are [metric/variant/aggregation]
        with open(os.path.join(INPUT_FILES["recon_results"] % dataset_tag, fltr, SETTINGS["rank_file_name"]), "r") as f:
            ranks_raw = eval(f.read())
            
            ranks = {}
            
            for wp in ranks_raw.keys():
                entries = wp.split(os.sep)
                
                ranks.setdefault(entries[0], {}) # Add metric key
                ranks[entries[0]].setdefault(entries[1], {}) # Add variant key
                ranks[entries[0]][entries[1]].setdefault(entries[2], {}) # Add aggregation key
                ranks[entries[0]][entries[1]][entries[2]] = np.array(list(ranks_raw[wp].values())) # Add values
        
        topk_table[fltr] = {metric: {} for metric in ranks.keys()}
        
        # Plot
        if not os.path.exists(OUTPUT_FILES["top_k"] % dataset_tag):
            os.makedirs(OUTPUT_FILES["top_k"] % dataset_tag)
            
        for m, metric in enumerate(ranks.keys()):
            fig, ax = plt.subplots()
            fig.suptitle(metric)
            ax.set_ylim([-0.05,1.05])
    
            cols = {v: c for v, c in zip(ranks[metric].keys(), SETTINGS["plt_colors"])}
    
            topk_table[fltr][metric] = {variant: {} for variant in ranks[metric].keys()}
    
            print(f" - metric: {metric}")
            for variant in ranks[metric].keys():
                topk_table[fltr][metric][variant] = {aggregation: {} for aggregation in ranks[metric][variant].keys()}
                
                print(f"\t - variant: {variant}")
                for aggregation in ranks[metric][variant].keys():
                    print(f"\t\t - aggregation: {aggregation}")
    
                    # Top-k accuracy
                    y = np.sort(ranks[metric][variant][aggregation])
                    topk = lambda z: np.float64(0) if z < min(y) else (1 + np.max(np.where(y <= z)[0]))/len(y) # the 1 accounts for indexes starting at 0
                    
                    # Plotting
                    yy = np.concatenate( ( y, np.arange(y[-1] + 1, SETTINGS["max_k"] + 1) ) )
                    t = [0] # replicates the effect of ax.ecdf(y)
                    for z in yy[:-1]:
                        t += [topk(z), topk(z)]
                    t += [1]
                    
                    ax.plot(
                        np.repeat(yy, 2), 
                        t, 
                        color = cols[variant], 
                        ls = SETTINGS["plt_styles"][aggregation], 
                        label = SETTINGS["plt_labels"][metric][variant][aggregation]
                    )
                    
                    ax.set_xlabel("$k$")
                    
                    topk_table[fltr][metric][variant][aggregation] = {
                            k: topk(k)
                            for k in SETTINGS["k"]
                        }
            
            plt.show()
            
            fig.suptitle("Top-$k$ accuracy")
            fig.savefig(os.path.join(OUTPUT_FILES["top_k"] % dataset_tag, f"{metric}__filter_{fltr}__tka.png"), format = "png")
            fig.savefig(os.path.join(OUTPUT_FILES["top_k"] % dataset_tag, f"{metric}__filter_{fltr}__tka.pdf"), format = "pdf")
    
    return topk_table

def nice_table(topk_table, add_to_clipboard = False):
    """
    Prints topk_table (plot_top_k_accuracy's output) in a latex-friendly format.
    
    - Arguments:
        - topk_table: a dictionary with keys [filter][metric][variant][aggregation]
        - add_to_clipboard: should the resulting latex code be copied to the system's clipboard? Default: False
    """
    
    n_ks = len(SETTINGS["k"])
    
    py2tex = {
            "n": "No filter",
            "hd": "Hapticity/denticity order",
            "c": "Charge",
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
    
    for metric in ["cossim", "delgsim"]:
        metric_entry = py2tex[metric]
        
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
    
# %% Main statement
if __name__ == "__main__":
    topk_table = plot_top_k_accuracy()
    nice_table(topk_table, True)
