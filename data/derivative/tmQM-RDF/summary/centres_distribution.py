"""
This script plots the histogram of centre counts. For simplicity, it relies on the preliminary information computed for tmQM-RDF-1Ksel.
"""
# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header
import matplotlib.transforms as pltran
import matplotlib.pyplot as plt
import pandas as pd
import re

from collections import Counter

DATASET_VERSION = "v2025dev"

INPUT_FILES = {
        "pubchem": os.path.join(ROOT_DIR, "data", "raw", "pubChem", "data", "PubChemElements_all.csv"),
        "details": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", DATASET_VERSION, "details.csv")
    }

OUTPUT_FILES = {
        "figures": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "summary", DATASET_VERSION, "figures")
    }

def plot_centres_distribution():
    """
    Plots the histogram of centres count
    """
    # Read PubChem information
    pubchem = pd.read_csv(INPUT_FILES["pubchem"], sep = ",")
    
    # Read selection info
    sel_info = pd.read_csv(INPUT_FILES["details"], index_col = 0)
    
    # Get centres count
    count = Counter(sel_info.loc[:, "core"])
    count = {
        el: c 
        for el, c in sorted(count.items(), key = lambda ec: pubchem.loc[pubchem["Symbol"] == ec[0], "AtomicNumber"].iloc[0])
    }
    
    # Plot
    fig, ax = plt.subplots(figsize = (6.2, 4.8))
    ax.bar(count.keys(), count.values(), color = "orange")
    
    ax.set_xticks(range(len(count)), count.keys(), fontsize = "x-small")
    plt.setp(ax.xaxis.get_majorticklabels(), rotation = 60, rotation_mode = "anchor", ha = "right")
    
    # Shift x axis tick labels by -2 points in x direction and by 3 in y direction
    # Notice: matplotlip uses 75 points per inches
    dx = -2/72.; dy = 3/72. # shifts (in inches)
    offset = pltran.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
    
    for label in ax.xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)
    
    ax.set_xlabel("Metal centres (ordered by atomic number)")
    ax.set_ylabel("Count")
    
    # Save plot
    if not os.path.exists(OUTPUT_FILES["figures"]):
        os.makedirs(OUTPUT_FILES["figures"])
        
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "centres_distribution.png"), format = "png", dpi = 200)
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "centres_distribution.pdf"), format = "pdf")
    
# %% Main statement
if __name__ == "__main__":
    plot_centres_distribution()