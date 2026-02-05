"""
Creates the scaffolds and the corresponding recosntructions from the TMCs in the test set
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header

import sys
sys.path.append(os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "interface"))

from tqdm import tqdm
from tmQM_RDF_interface import TmQMRDFInterface

import os
import shutil
import numpy as np
import pandas as pd
import scaffold_utils as scut

DEBUG = False
DEBUG_SIZE = 0 # Set value to 0 or lower to disable debugging and start production

INPUT_FILES = {
        "tmQM-RDF": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev", "graphs"),
        "tmQM-RDF-1Ksel": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "data", "%s", "1k_selection.csv"),
        "ligands": os.path.join(ROOT_DIR, "computational", "reconstruction", "intermediate", "ligands", "%s")
    }

OUTPUT_FILES = {
        "scaffolds": os.path.join(ROOT_DIR, "computational", "reconstruction", "intermediate", "scaffold", "%s"),
        "reconstructions": os.path.join(ROOT_DIR, "computational", "reconstruction", "intermediate", "reconstructions", "rdf", "%s"),
        "reconstructions_figures": os.path.join(ROOT_DIR, "computational", "reconstruction", "intermediate", "reconstructions", "figures", "%s")
    }

# %% Utilities
def create_scaffolds_and_reconstructions(tm_mode):
    # Initialise static pointers
    TmQMRDFInterface.path_to_tmQM_RDF = INPUT_FILES["tmQM-RDF"]
    scut.PrunedTMCDataset.path_to_tmQM_RDF = INPUT_FILES["tmQM-RDF"]
    scut.PrunableTMCGraph.path_to_ligand_blocks = INPUT_FILES["ligands"] % tm_mode
    scut.PrunedTMCGraph.path_to_ligand_blocks = INPUT_FILES["ligands"] % tm_mode
    
    # Create dataset
    dataset = scut.PrunedTMCDataset(INPUT_FILES["tmQM-RDF-1Ksel"] % tm_mode, OUTPUT_FILES["scaffolds"] % tm_mode, INPUT_FILES["ligands"] % tm_mode)

    # Prepare directory
    if os.path.exists(OUTPUT_FILES["reconstructions"] % tm_mode):
        shutil.rmtree(OUTPUT_FILES["reconstructions"] % tm_mode)
    
    os.makedirs(OUTPUT_FILES["reconstructions"] % tm_mode)
    
    # Process trial batch (10 graphs)
    loader = dataset.get_loader("test")
    
    for i, datapoint in tqdm(enumerate(loader)):
        #print(f"------------------------------------\n\nProcessing TMC n. {i}\n\n------------------------------------")
        
        local_reconstruction_directory = os.path.join(OUTPUT_FILES["reconstructions"] % tm_mode, f"recon{i}")
        local_graph_directory = os.path.join(local_reconstruction_directory, "graphs")
        local_matches_directory = os.path.join(OUTPUT_FILES["reconstructions"] % tm_mode, f"recon{i}")
        
        os.makedirs(local_reconstruction_directory)
        os.makedirs(local_graph_directory)
        if not os.path.exists(local_matches_directory):
            os.makedirs(local_matches_directory)
        
        datapoint[0].generate_all_reconstructions("train", local_graph_directory)
        
        selection = {
            "TMC": [
                x.replace(".ttl", "") 
                for x in os.listdir(local_graph_directory) 
                if x.endswith(".ttl")
            ]
        }
        selection = pd.DataFrame(selection)
        selection.to_csv(os.path.join(local_reconstruction_directory, "selection.csv"))
        
        with open(os.path.join(local_matches_directory, f"{datapoint[1]}.gt"), "w") as f:
            f.write("")

# %% Main
if __name__ == "__main__":
    create_scaffolds_and_reconstructions("earlyTM")
    create_scaffolds_and_reconstructions("lateTM")