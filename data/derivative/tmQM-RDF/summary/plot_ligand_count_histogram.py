"""
This script plots the histogram of ligand counts. For simplicity, it relies on the preliminary information computed for tmQM-RDF-1Ksel.
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as pltpatch
import matplotlib.transforms as pltran

from collections import defaultdict
from tqdm import tqdm

DATASET_VERSION = "v2025dev"

INPUT_FILES = {
        "pubchem": os.path.join(ROOT_DIR, "data", "raw", "pubChem", "data", "PubChemElements_all.csv"),
        "details": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", DATASET_VERSION, "details.csv")
    }

OUTPUT_FILES = {
        "figures": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "summary", DATASET_VERSION, "figures")
    }

# %% Utility functions

def retrieve_ligand_counts():
    """
    Counts the ligands that appear in different scenarios (i.e., different datasets and different metal centres)
    
    Returns:
        - two structures made of two nested dictionaries indexed by: 
            metal centre (including the option 'all'), ligand
          the first structure counts all the occurrences of a ligand, whereas the second only counts the number of distinct 
          TMCs in which each ligand appears
    """
    
    # Retrieve TMC compositions
    details = pd.read_csv(INPUT_FILES["details"], index_col = 0)
    
    centres = set(details.loc[:, "core"])
    
    # Initialise working variables
    ligand_count = {c: defaultdict(int) for c in ["all"] + list(centres)}
    ligand_count_no_rep = {c: defaultdict(int) for c in ["all"] + list(centres)}

    # Count ligands within tmQM-RDF
    for tmc in tqdm(details.index):
        if details.loc[tmc, "core"] in centres:
            g = details.loc[tmc, :]
            
            ligands = g.loc["ligands"].split(" ")
            
            for lig in ligands:
                ligand_count["all"][lig] += 1
                ligand_count[g.loc["core"]][lig] += 1
                
            for lig in set(ligands):
                ligand_count_no_rep["all"][lig] += 1
                ligand_count_no_rep[g.loc["core"]][lig] += 1

    return ligand_count, ligand_count_no_rep

def print_quantiles_all_centres(ligand_count, ligand_count_no_rep):
    """
    Prints the quantiles of observed ligand counts (with and without multiplicity)
    
    Arguments:
        - ligand_count, ligand_count_no_rep: see output of retrieve_ligand_counts
    """
    
    alphas = np.concatenate([np.arange(0,0.9,0.1), np.arange(0.9, 1, 0.01)])
    
    q_rep = np.quantile(list(ligand_count["all"].values()), alphas)
    q_no_rep = np.quantile(list(ligand_count_no_rep["all"].values()), alphas)
    
    print("\t a   \t\t REP \t\t NO_REP")
    for q in zip(alphas, q_rep, q_no_rep):
        print(f"\t {q[0]:.2f} \t\t {q[1]:.5f} \t\t {q[2]:.5f}")

def plot_quantiles(ligand_count, ligand_count_no_rep, alpha = 0.95):
    """
    Plots the 95% quantiles of the empirical distributions of ligand counts (with and without multiplicity)
    
    Arguments:
        - ligand_count, ligand_count_no_rep: see output of retrieve_ligand_counts
        - alpha: the level of the desired quantiles (default: 0.95)
    """
    
    # Compute quantiles
    q_rep = {c: np.quantile(list(ligand_count[c].values()), alpha) for c in ligand_count}
    q_no_rep = {c: np.quantile(list(ligand_count_no_rep[c].values()), alpha) for c in ligand_count_no_rep}
    
    # Prepare element based colors
    pubchem = pd.read_csv(INPUT_FILES["pubchem"], sep = ",")
    centre_colors = {centre: "#" + pubchem.loc[pubchem["Symbol"] == centre, "CPKHexColor"].iloc[0] if centre != "all" else "#000000"
                     for centre in ligand_count.keys()}
    
    # Check that the colors extracted from PubChem are valid HEX colors, using a regex
    #   If a color is not valid, use a default color (pink)
    for centre, color in centre_colors.items():
        if not re.match(r"^#(?:[0-9a-fA-F]{2}){3}$", color):
            centre_colors[centre] = "#F78BB2"
            
    # Plot
    fig, ax = plt.subplots(figsize = (6.2, 4.8))
    
    x = range(len(ligand_count))
    ax.set_xlim([-1, max(x) + 1])
    max_y = max([max(q_rep.values()), max(q_no_rep.values())])
    ax.set_ylim([0, max_y + 1])
    
    sorted_centres = ["all"] + sorted(
        [x for x in ligand_count.keys() if x != "all"], 
        key = lambda el: pubchem.loc[pubchem["Symbol"] == el, "AtomicNumber"].iloc[0]
    )
    
    for i, c in enumerate(sorted_centres):
        legend_handles = [{}, {}]
        if c == "all":
            #legend_handles = [{"label": "Counts of ligand instance"}, {"label": "Counts of TMCs with copy of ligand"}]
            legend_handles = [{"label": "$q_{0.95}^{inst}$"}, {"label": "$q_{0.95}^{copy}$"}]
        
        ax.scatter(i, q_rep[c], edgecolors = centre_colors[c], facecolors = "none", marker = "^", **legend_handles[0])
        ax.scatter(i, q_no_rep[c], c = centre_colors[c], marker = "v", **legend_handles[1])
        ax.vlines(i, 0, max(q_rep[c], q_no_rep[c]), colors = centre_colors[c], alpha = 0.5, linestyle = "dashed")
    
    # Set x ticks
    ax.set_xticks(x, [l.capitalize() for l in sorted_centres], fontsize = "x-small")
    plt.setp(ax.xaxis.get_majorticklabels(), rotation = 60, rotation_mode = "anchor", ha = "right")

    # Shift x axis tick labels by -2 points in x direction and by 3 in y direction
    # Notice: matplotlip uses 75 points per inches
    dx = -2/72.; dy = 3/72. # shifts (in inches)
    offset = pltran.ScaledTranslation(dx, dy, fig.dpi_scale_trans)

    for label in ax.xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)
        
    # Set x/y labels
    ax.set_xlabel("Metal centres (ordered by atomic number)")
    ax.set_ylabel("95% quantile")

    # Legend
    ax.legend(fontsize = "x-small") 
    
    # Save plot
    if not os.path.exists(OUTPUT_FILES["figures"]):
        os.makedirs(OUTPUT_FILES["figures"])
        
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "ligand_count_quantiles.png"), format = "png", dpi = 200)
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "ligand_count_quantiles.pdf"), format = "pdf")
    
# %% Main statement
if __name__ == "__main__":
    ligand_count, ligand_count_no_rep = retrieve_ligand_counts()
    print_quantiles_all_centres(ligand_count, ligand_count_no_rep)
    plot_quantiles(ligand_count, ligand_count_no_rep)
