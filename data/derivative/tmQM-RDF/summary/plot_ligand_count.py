# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header
import sys
sys.path.append(ROOT_DIR)

from rdkit import Chem
from rdkit.Chem import Draw
from collections import defaultdict
from matplotlib.gridspec import GridSpec

import matplotlib.transforms as pltran
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cairosvg

DATASET_VERSION = "v2025dev"

TMQMG = {
        "v2025dev": "v74.637k"    
    }
TMQMGL = {
        "v2025dev": "v60k"    
    }

INPUT_FILES = {
        "viable": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", DATASET_VERSION, "viable_tmcs.txt"),
        "centres": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", DATASET_VERSION, "centres.csv"),
        "details": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", DATASET_VERSION, "details.csv"),
        "sizes": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", DATASET_VERSION, "tmc_sizes.txt"),
        "ligands": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", DATASET_VERSION, "%s", "ligands.csv"),
        "outliers": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", DATASET_VERSION, "outliers.txt"),
        "tmQMg_graphs": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg", TMQMG[DATASET_VERSION], "uNatQ_graphs"),
        "selection": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "data", DATASET_VERSION, "%s", "1k_selection.csv"),
        "ligands_fingerprints": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", TMQMGL[DATASET_VERSION], "ligands_fingerprints.csv"),
        "ligands_misc_info": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", TMQMGL[DATASET_VERSION], "ligands_misc_info.csv")
    }

OUTPUT_FILES = {
        "figures": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "summary", DATASET_VERSION, "figures")
    }

# Read PubChem information
#pubchem = pd.read_csv("../data/preprocessing/misc/PubChemElements_all.csv", sep = ",")

# Read selection info
# sel_info = pd.read_csv("../data/extract_selection_from_tmQM-RDF/tmc_info_for_selection/detail.csv", index_col = 0)
# sel_df = pd.read_csv("../data/extract_selection_from_tmQM-RDF/1k_selection.csv")

# %% Utility functions

def retrieve_ligand_counts():
    """
    Counts the ligands that appear in tmQM-RDF
    
    Returns:
        - two dictionaries, indexed by ligand, that count all the occurrences of a ligand (the first), and the number of distinct 
          TMCs in which each ligand appears (the second)
    """
    
    # Retrieve TMC compositions
    details = pd.read_csv(INPUT_FILES["details"], index_col = 0)
    
    # Initialise working variables

    ligand_count = defaultdict(int)
    ligand_count_no_rep = defaultdict(int)

    # Count linands    
    for tmc in details.index:
        g = details.loc[tmc, :]
        
        ligands = g.loc["ligands"].split(" ")
        
        for lig in ligands:
            ligand_count[lig] += 1
            
        for lig in set(ligands):
            ligand_count_no_rep[lig] += 1

    return ligand_count, ligand_count_no_rep

# %% Plotting functions

def plot_ligands(ligand_count, ligand_count_no_rep, n_lig_to_highlight = 10):
    """
    Plots the most frequent ligands in the datasets (tmQM-RDF, training, validation, test) by metal centre
    
    - Arguments:
        - ligand_count, ligand_count_no_rep: the outputs of retrieve_ligand_counts
        - n_lig_to_highlight: the n most frequent ligands to highlight in the plots
    """
    
    if not os.path.exists(OUTPUT_FILES["figures"]):
        os.makedirs(OUTPUT_FILES["figures"])
    
    # Acquire ligand info (size) from tmQMg-L
    ligand_fingerprints = pd.read_csv(INPUT_FILES["ligands_fingerprints"], sep = ";")
    
    ligand_names = np.array(ligand_fingerprints.loc[:, "name"])
    ligand_sizes = np.array(ligand_fingerprints.loc[:, "n_atoms"])
    
    # order ligands by size
    perm = np.argsort(ligand_sizes, stable = True)
    
    sorted_lig_sizes = ligand_sizes[perm]
    sorted_lig_names = ligand_names[perm]

    sorted_lig_count = [ligand_count[lig] for lig in sorted_lig_names]
    sorted_lig_count_no_rep = [ligand_count_no_rep[lig] for lig in sorted_lig_names]
    
    # Pick indexes of ligands to highlight
    frequent_ligands_idxs = np.argsort(sorted_lig_count)[-n_lig_to_highlight:]
    
    # Order indexes of ligands to highlight so that lexycographic order is locally preserved
    temp = []
    
    for i in sorted(list(set(sorted_lig_sizes[frequent_ligands_idxs]))):
        loc_freq_ligands = [j for j in frequent_ligands_idxs if sorted_lig_sizes[j] == i]
        temp += sorted(loc_freq_ligands, key = lambda j: sorted_lig_names[j])
        
    frequent_ligands_idxs = temp
        
    print([sorted_lig_names[i] for i in frequent_ligands_idxs])
    
    # ----------
    # Plot frequent ligands counts and structures (by centre)
    # ----------
    fig, ax = plt.subplots(figsize = (6.2, 4.8))
    # gs = GridSpec(1, 2, width_ratios = [1, 2])
    
    # axs = []
    # axs += [fig.add_subplot(gs[0])]
    # axs += [fig.add_subplot(gs[1])]
    
    # ----- Counts
    # ax = axs[0]
    
    # x position of bars
    x = np.arange(len(frequent_ligands_idxs))
    
    # graphical parameters
    barwidth = 0.17
    
    # Add counts
    
    y = [sorted_lig_count[i] for i in frequent_ligands_idxs]
    
    ax.bar(
            x, 
            y, 
            barwidth, 
            facecolor = "white",
            edgecolor = "orange",
            hatch = "///",
            label = "Instances"
        )
    
    # Add counts (no rep)
    y = [sorted_lig_count_no_rep[i] for i in frequent_ligands_idxs]
    
    ax.bar(
            x, 
            y, 
            barwidth, 
            color = "orange",
            linewidth = 1,
            label = "N. of TMCs with $\\geq\\ 1$ copy"
        )
        
    # Add axis labels
    ax.set_xlabel("Ligands (ordered by size)")
    ax.set_ylabel("Counts")
    
    # Add and rotate x axis tick labels
    ax.set_xticks(x, [sorted_lig_names[i] for i in frequent_ligands_idxs], fontsize = "x-small")
    plt.setp(ax.xaxis.get_majorticklabels(), rotation = 80, rotation_mode = "anchor", ha = "right")
    
    # Shift x axis tick labels by -2 points in x direction and by 3 in y direction
    # Notice: matplotlip uses 75 points per inches
    dx = -2/72.; dy = 3/72. # shifts (in inches)
    offset = pltran.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
    
    for label in ax.xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)
    
    # Add legend
    ax.legend(fontsize = "x-small")
    
    # # Save figure
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "lig_counts.pdf"), format = "pdf", bbox_inches = "tight", pad_inches = 0.1)
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "lig_counts.png"), format = "png", dpi = 200, bbox_inches = "tight", pad_inches = 0.1)
    
    # -----  Structures
    # ax = axs[1]
    
    fig, ax = plt.subplots()
    
    # Get frequent ligand data (SMILES)
    ligand_info = pd.read_csv(INPUT_FILES["ligands_misc_info"], sep = ";")
    
    ligand_smiles = []
    
    for i in frequent_ligands_idxs:
        ligand_smiles += [
            ligand_info.loc[ligand_info["name"] == sorted_lig_names[i], "smiles"].iloc[0]
        ]
        
    ligands_as_mol = [Chem.MolFromSmiles(smiles) for smiles in ligand_smiles]
    
    # Plot
    img = Draw.MolsToGridImage(
        ligands_as_mol,
        molsPerRow = 5,
        legends = [sorted_lig_names[i] for i in frequent_ligands_idxs],
        useSVG = True
    ) 
    
    img0 = Draw.MolsToGridImage(
        ligands_as_mol,
        molsPerRow = 5,
        legends = [sorted_lig_names[i] for i in frequent_ligands_idxs],
        useSVG = False
    )
    
    ax.imshow(img0)
    ax.axis("off")
    
    # Add axis title
    # ax.title.set_text(dset.capitalize())
    
    # Save image
    cairosvg.svg2png(bytestring = img, write_to = os.path.join(OUTPUT_FILES["figures"], "freq_ligands.png"), dpi = 200)
    cairosvg.svg2pdf(bytestring = img, write_to = os.path.join(OUTPUT_FILES["figures"], "freq_ligands.pdf"))
    
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "freq_ligs_structure.pdf"), format = "pdf", bbox_inches = "tight")
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "freq_ligs_structure.png"), format = "png", dpi = 200, bbox_inches = "tight")
        
    plt.show()
    
# %% Main statement

if __name__ == "__main__":
    ligand_count, ligand_count_no_rep = retrieve_ligand_counts()
    plot_ligands(ligand_count, ligand_count_no_rep)
    