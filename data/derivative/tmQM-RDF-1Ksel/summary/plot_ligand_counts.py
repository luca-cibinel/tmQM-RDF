# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header
import sys
sys.path.append(ROOT_DIR)

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Draw
from collections import defaultdict
from matplotlib.gridspec import GridSpec

import matplotlib.transforms as pltran
import matplotlib.patches as pltpatch
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cairosvg
import re

DATASET_VERSION = "v2025dev"

TMQMG = {
        "v2025dev": "v74.637k"
    }

TMQMGL= {
        "v2025dev": "v60k"
    }

INPUT_FILES = {
        "pubchem": os.path.join(ROOT_DIR, "data", "raw", "pubChem", "data", "PubChemElements_all.csv"),
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
        "figures": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "summary", DATASET_VERSION, "figures", "%s")
    }

# Read PubChem information
#pubchem = pd.read_csv("../data/preprocessing/misc/PubChemElements_all.csv", sep = ",")

# Read selection info
# sel_info = pd.read_csv("../data/extract_selection_from_tmQM-RDF/tmc_info_for_selection/detail.csv", index_col = 0)
# sel_df = pd.read_csv("../data/extract_selection_from_tmQM-RDF/1k_selection.csv")

# %% Utility functions

def retrieve_ligand_counts(tm_mode, n_freq_ligands = 100):
    """
    Counts the ligands that appear in different scenarios (i.e., different datasets and different metal centres)
    
    Arguments:
        - tm_mode: whether to process early TMs (mode = 'earlyTM') or the late ones (mode = 'lateTM')
        - n_freq_ligands: how many of the most frequent ligands should be considered
        
    Returns:
        - two structures made of three levels of nested dictionaries indexed by: 
            dataset (tmQM-RDF, train, validation, test), metal centre (including the option 'all'), ligand
          the first structure counts all the occurrences of a ligand, whereas the second only counts the number of distinct 
          TMCs in which each ligand appears
    """
    
    # Validate tm_mode
    if tm_mode == "earlyTM":
        centres = ["Cr", "Mo", "W"]
    elif tm_mode == "lateTM":
        centres = ["Pd", "Ni", "Pt"]
    else:
        raise Exception(f"tm_mode {tm_mode} not recognised!")
    
    # Retrieve 1k selection
    sel_df = pd.read_csv(INPUT_FILES["selection"] % tm_mode)
    
    # Retrieve TMC compositions
    sel_info = pd.read_csv(INPUT_FILES["details"], index_col = 0)
    
    # Initialise working variables
    dsets = ["tmQM-RDF", "train", "validation", "test"]

    ligand_count = {dset: {c: defaultdict(int) for c in ["all"] + centres} for dset in dsets}
    ligand_count_no_rep = {dset: {c: defaultdict(int) for c in ["all"] + centres} for dset in dsets}

    # Count linands within the partitions of the selection
    for dset in dsets[1:]:
        print(f"Preprocessing {dset} set...")
        
        # Read selection
        sel = list(sel_df.loc[sel_df["partition"] == dset, "TMC"])
            
        # Read each TMC in the selection and count the ligands
        for tmc in sel:
            g = sel_info.loc[tmc, :]
            
            ligands = g.loc["ligands"].split(" ")
            
            for lig in ligands:
                ligand_count[dset]["all"][lig] += 1
                ligand_count[dset][g.loc["core"]][lig] += 1
                
            for lig in set(ligands):
                ligand_count_no_rep[dset]["all"][lig] += 1
                ligand_count_no_rep[dset][g.loc["core"]][lig] += 1

    # Count ligands within tmQM-RDF (only requested centres)
    for tmc in tqdm(sel_info.index):
        if sel_info.loc[tmc, "core"] in centres:
            g = sel_info.loc[tmc, :]
            
            ligands = g.loc["ligands"].split(" ")
            
            for lig in ligands:
                ligand_count[dsets[0]]["all"][lig] += 1
                ligand_count[dsets[0]][g.loc["core"]][lig] += 1
                
            for lig in set(ligands):
                ligand_count_no_rep[dsets[0]]["all"][lig] += 1
                ligand_count_no_rep[dsets[0]][g.loc["core"]][lig] += 1

    # Filter ligands counted in tmQM-RDF according to ligand frequency
    freq_ligands = np.argsort(list(ligand_count[dsets[0]]["all"].values()))[-n_freq_ligands:]
    freq_ligands_names = [list(ligand_count[dsets[0]]["all"].keys())[i] for i in freq_ligands]

    ligand_count[dsets[0]] = {el: {l: ligand_count[dsets[0]][el][l] for l in freq_ligands_names} for el in ligand_count[dsets[0]].keys()}
    ligand_count_no_rep[dsets[0]] = {el: {l: ligand_count_no_rep[dsets[0]][el][l] for l in freq_ligands_names} for el in ligand_count_no_rep[dsets[0]].keys()}

    return ligand_count, ligand_count_no_rep

def print_thresholds(ligand_count_no_rep, n_lig_to_highlight = 10, max_threshold = 20):
    """
    Prints a table reporting how many of the most frequent ligands can be found at each threshold.
    
    - Arguments:
        - ligand_count_no_rep: the outputs of retrieve_ligand_counts
        - n_lig_to_highlight: the n most frequent ligands to consider,
        - max_threshold: the maximum threshold to consider (inclusive)
    """
    
    # Retrieve most frequent ligands
    frequent_ligands_idxs = np.argsort(list(ligand_count_no_rep["train"]["all"].values()))[-n_lig_to_highlight:]
    frequent_ligands_names = np.array(list(ligand_count_no_rep["train"]["all"].keys()))[frequent_ligands_idxs]
    frequent_ligands_counts = {
            key: [ligand_count_no_rep["train"][key].get(l, 0) for l in frequent_ligands_names]
            for key in ligand_count_no_rep["train"].keys()
        }
    print(frequent_ligands_counts)
    
    # Print table
    print("\t\t".join(["alpha"] + list(frequent_ligands_counts.keys())))
    for alpha in range(5, max_threshold + 1):
        hits = [str(sum([c >= alpha for c in frequent_ligands_counts[key]])) for key in frequent_ligands_counts.keys()]
        
        print("\t\t".join([str(alpha)] + hits))

# %% Plotting functions

def plot_ligands(tm_mode, ligand_count, ligand_count_no_rep, n_lig_to_highlight = 10):
    """
    Plots the most frequent ligands in the datasets (tmQM-RDF, training, validation, test) by metal centre
    
    - Arguments:
        - tm_mode: whether to process early TMs (mode = 'earlyTM') or the late ones (mode = 'lateTM')
        - ligand_count, ligand_count_no_rep: the outputs of retrieve_ligand_counts
        - n_lig_to_highlight: the n most frequent ligands to highlight in the plots
    """
    
    if not os.path.exists(OUTPUT_FILES["figures"] % tm_mode):
        os.makedirs(OUTPUT_FILES["figures"] % tm_mode)
    
    pubchem = pd.read_csv(INPUT_FILES["pubchem"], sep = ",")
    
    # Initialise working variables
    dsets = ["tmQM-RDF", "train", "validation", "test"]
    
    # Acquire ligand info (size) from tmQMg-L
    ligand_fingerprints = pd.read_csv(INPUT_FILES["ligands_fingerprints"], sep = ";")
    
    # Create subplots
    
    # Process and plot
    for dset in dsets:
        loc_ligand_fingerprints = ligand_fingerprints.loc[ligand_fingerprints["name"].isin(ligand_count[dset]["all"]), ["name", "n_atoms"]]
        loc_ligand_fingerprints = loc_ligand_fingerprints.loc[sorted(loc_ligand_fingerprints.index), :]
        
        ligand_names = np.array(loc_ligand_fingerprints.loc[:, "name"])
        ligand_sizes = np.array(loc_ligand_fingerprints.loc[:, "n_atoms"])
        
        # order ligands by size
        perm = np.argsort(ligand_sizes, stable = True)
        
        sorted_lig_sizes = ligand_sizes[perm]
        sorted_lig_names = ligand_names[perm]

        sorted_lig_count = [ligand_count[dset]["all"][lig] for lig in sorted_lig_names]
        sorted_lig_count_no_rep = [ligand_count_no_rep[dset]["all"][lig] for lig in sorted_lig_names]
        
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
        # Plot ligand counts (all centres)
        # ----------
        main_color = "blue" if dset == "tmQM-RDF" else "orange"
        
        fig, ax = plt.subplots(figsize = (6.2, 4.8), dpi = 200)
        
        ax.bar(
               range(len(sorted_lig_names)), 
               sorted_lig_count, 
               facecolor = "white", 
               edgecolor = main_color,
               hatch = "///",
               linewidth = 0.2,
               label = "Counts"
        )
        
        ax.bar(
               range(len(sorted_lig_names)), 
               sorted_lig_count_no_rep, 
               color = main_color, 
               label = "Counts (multiplicity ignored)"
        )
        
        for i in frequent_ligands_idxs:
            ax.text(i + 1, sorted_lig_count[i], sorted_lig_names[i], size = "xx-small", rotation = 30)
        
        ax.set_ylim([0, max(sorted_lig_count)*1.15])
        
        # Add legend
        if dset == dsets[0]:
            legend_handles = [pltpatch.Patch(color = "black", label = "Counts (multiplicity ignored)")]
            legend_handles += [pltpatch.Patch(edgecolor = "black", facecolor = "white", label = "Counts", hatch = "///")]
            ax.legend(handles = legend_handles, ncols = 2, fontsize = "x-small")
        
        ax.set_xlabel("Ligands (ordered by size)")
        ax.set_ylabel("Counts", labelpad = 10.5 if dset != dsets[0] else 4)
        
        # Set tile (only for display)
        ax.title.set_text(dset.capitalize())
        plt.show()
        ax.title.set_text(None)
    
        fig.savefig(os.path.join(OUTPUT_FILES["figures"] % tm_mode, f"lig_counts_all_{dset}.pdf"), format = "pdf", bbox_inches = "tight")
        fig.savefig(os.path.join(OUTPUT_FILES["figures"] % tm_mode, f"lig_counts_all_{dset}.png"), dpi = 200, bbox_inches = "tight")
    
        # ----------
        # Plot frequent ligands counts and structures (by centre)
        # ----------
        fig = plt.figure(figsize = (6.2*3, 4.8), dpi = 200)
        gs = GridSpec(1, 2, width_ratios = [1, 2])
        
        axs = []
        axs += [fig.add_subplot(gs[0])]
        axs += [fig.add_subplot(gs[1])]
        
        # ----- Counts
        ax = axs[0]
        
        # Filter ligand_count and ligand_count_no_rep for those counts that are above 25
        freq_ligs_count = {
                centre : {
                        sorted_lig_names[i]: ligand_count[dset][centre][sorted_lig_names[i]]
                        for i in frequent_ligands_idxs
                    }
                for centre in ligand_count[dset].keys()
                if centre != "all"
            }
        
        freq_ligs_count_no_rep = {
                centre : {
                        sorted_lig_names[i]: ligand_count_no_rep[dset][centre][sorted_lig_names[i]]
                        for i in frequent_ligands_idxs
                    }
                for centre in ligand_count_no_rep[dset].keys()
                if centre != "all"
            }
        
        # Prepare element based colors
        centre_colors = {centre: pubchem.loc[pubchem["Symbol"] == centre, "CPKHexColor"].iloc[0]
                         for centre in freq_ligs_count.keys()}
        
        # Check that the colors extracted from PubChem are valid HEX colors, using a regex
        #   If a color is not valid, use a default color (pink)
        for centre, color in centre_colors.items():
            if not re.match(r"^(?:[0-9a-fA-F]{2}){3}$", color):
                centre_colors[centre] = "F78BB2"
        
        # x position of bars
        x = np.arange(len(frequent_ligands_idxs))
        
        # graphical parameters
        barwidth = 0.17
        hspace = 0.06
        
        # Add counts
        sorted_counts = sorted(
                freq_ligs_count.items(),
                key = lambda el: pubchem.loc[pubchem["Symbol"] == el[0], "AtomicNumber"].iloc[0])
        
        counter = 0
        for centre, ligs in sorted_counts:
            offset = (barwidth + hspace)*counter
            
            y = list(ligs.values())
            
            ax.bar(
                    x + offset, 
                    y, 
                    barwidth, 
                    facecolor = "white",
                    edgecolor = "#" + centre_colors[centre],
                    hatch = "///"
                )
            
            counter += 1
        
        # Add counts (no rep)
        sorted_counts_no_rep = sorted(
                freq_ligs_count_no_rep.items(),
                key = lambda el: pubchem.loc[pubchem["Symbol"] == el[0], "AtomicNumber"].iloc[0])
        
        counter = 0
        legend_handles = []
        for centre, ligs in sorted_counts_no_rep:
            offset = (barwidth + hspace)*counter
            ax.bar(
                    x + offset, 
                    ligs.values(), 
                    barwidth, 
                    color = "#" + centre_colors[centre],
                    linewidth = 1,
                    label = centre
                )
            
            legend_handles += [pltpatch.Patch(color = "#" + centre_colors[centre], label = centre)]
            
            counter += 1
            
        # Add axis labels
        ax.set_xlabel("Ligands (ordered by size)")
        ax.set_ylabel("Counts", labelpad = 10.5 if dset != dsets[0] else 4)
        
        # Add and rotate x axis tick labels
        ax.set_xticks(x + barwidth + hspace, [sorted_lig_names[i] for i in frequent_ligands_idxs], fontsize = "x-small")
        plt.setp(ax.xaxis.get_majorticklabels(), rotation = 80, rotation_mode = "anchor", ha = "right")
        
        # Shift x axis tick labels by -2 points in x direction and by 3 in y direction
        # Notice: matplotlip uses 75 points per inches
        dx = -2/72.; dy = 3/72. # shifts (in inches)
        offset = pltran.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
        
        for label in ax.xaxis.get_majorticklabels():
            label.set_transform(label.get_transform() + offset)
        
        # Add legend
        if dset == dsets[0]:
            legend_handles += [pltpatch.Patch(color = "black", label = "Counts (multiplicity ignored)")]
            legend_handles += [pltpatch.Patch(edgecolor = "black", facecolor = "white", label = "Counts", hatch = "///")]
            ax.legend(handles = legend_handles, ncols = 2, fontsize = "x-small")
        
        # # Save figure
        # fig.savefig(OUTPUT_FILES["figures"] % (tm_mode, f"lig_counts_centre_{dset}", "pdf"), format = "pdf", bbox_inches = "tight")
        # fig.savefig(OUTPUT_FILES["figures"] % (tm_mode, f"lig_counts_centre_{dset}", "png"), format = "png", dpi = 200, bbox_inches = "tight")
        
        # -----  Structures
        ax = axs[1]
        
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
        cairosvg.svg2png(bytestring = img, write_to = os.path.join(OUTPUT_FILES["figures"] % tm_mode, f"lfreq_ligands_{dset}.png"), dpi = 200)
        cairosvg.svg2pdf(bytestring = img, write_to = os.path.join(OUTPUT_FILES["figures"] % tm_mode, f"lfreq_ligands_{dset}.pdf"))
        
        fig.savefig(os.path.join(OUTPUT_FILES["figures"] % tm_mode, f"freq_ligs_summary_{dset}.pdf"), format = "pdf", bbox_inches = "tight")
        fig.savefig(os.path.join(OUTPUT_FILES["figures"] % tm_mode, f"freq_ligs_summary_{dset}.png"), format = "png", dpi = 200, bbox_inches = "tight")
        
    plt.show()
    
# %% Main statement

if __name__ == "__main__":
    ligand_count, ligand_count_no_rep = retrieve_ligand_counts("lateTM")
    print_thresholds(ligand_count_no_rep)
    plot_ligands("lateTM", ligand_count, ligand_count_no_rep)
    
    ligand_count, ligand_count_no_rep = retrieve_ligand_counts("earlyTM")
    print_thresholds(ligand_count_no_rep)
    plot_ligands("earlyTM", ligand_count, ligand_count_no_rep)
    