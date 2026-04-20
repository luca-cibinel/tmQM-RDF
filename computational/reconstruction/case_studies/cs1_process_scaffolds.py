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
sys.path.append(os.path.join(ROOT_DIR, "computational", "reconstruction"))

from tqdm import tqdm
from tmQM_RDF_interface import TmQMRDFInterface

import os
import shutil
import numpy as np
import pandas as pd
import scaffold_utils as scut

DEBUG = False
DEBUG_SIZE = 0 # Set value to 0 or lower to disable debugging and start production

CASE_STUDIES = {
        "earlyTM": ["CAKRAH"],
        "lateTM": ["ILUYOD"]
    }

INPUT_FILES = {
        "tmQM-RDF": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev", "graphs"),
        "tmQM-RDF-1Ksel": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "data", "%s", "1k_selection.csv"),
        "ligands": os.path.join(ROOT_DIR, "computational", "reconstruction", "intermediate", "ligands", "%s")
    }

OUTPUT_FILES = {
        "scaffolds": os.path.join(ROOT_DIR, "computational", "reconstruction", "case_studies", "intermediate", "scaffold", "%s"),
        "reconstructions": os.path.join(ROOT_DIR, "computational", "reconstruction", "case_studies", "intermediate", "reconstructions", "rdf", "%s"),
        "reconstructions_figures": os.path.join(ROOT_DIR, "computational", "reconstruction", "case_studies", "intermediate", "reconstructions", "figures", "%s")
    }

# %% Utilities
def create_scaffolds_and_reconstructions():
    """
    For each case study TMC, create a scaffold for each different ligand found in it
    """
    
    # Initialise static pointers
    TmQMRDFInterface.path_to_tmQM_RDF = INPUT_FILES["tmQM-RDF"]
    scut.PrunedTMCDataset.path_to_tmQM_RDF = INPUT_FILES["tmQM-RDF"]
    
    for tm_mode in CASE_STUDIES:
        scut.PrunableTMCGraph.path_to_ligand_blocks = INPUT_FILES["ligands"] % tm_mode
        scut.PrunedTMCGraph.path_to_ligand_blocks = INPUT_FILES["ligands"] % tm_mode
        
        for case_tmc in CASE_STUDIES[tm_mode]:
            if os.path.exists(os.path.join(OUTPUT_FILES["reconstructions"] % tm_mode, case_tmc)):
                print(f"WARNING! Case study {case_tmc} ({tm_mode}) already exists. Skipping...")
                continue
            
            prunable = scut.PrunableTMCGraph(case_tmc)
            
            unique_ligands = []
            
            # prunable.ligands is a dictionary of the form {id : (ligand_instance_URI, ligand_class)}, where id is an integer
            for i, ligand_info in prunable.ligands.items():
                if ligand_info[1] in unique_ligands:
                    continue
                
                unique_ligands += [ligand_info[1]]
                
                pruned = prunable.prune(i)
                
                local_reconstruction_directory = os.path.join(OUTPUT_FILES["reconstructions"] % tm_mode, case_tmc, f"recon_{ligand_info[1]}")
                local_graph_directory = os.path.join(local_reconstruction_directory, "graphs")
                local_matches_directory = local_reconstruction_directory
                
                os.makedirs(local_reconstruction_directory)
                os.makedirs(local_graph_directory)
                if not os.path.exists(local_matches_directory):
                    os.makedirs(local_matches_directory)
                
                pruned.generate_all_reconstructions("train", local_graph_directory)
                
                selection = {
                    "TMC": [
                        x.replace(".ttl", "") 
                        for x in os.listdir(local_graph_directory) 
                        if x.endswith(".ttl")
                    ]
                }
                selection = pd.DataFrame(selection)
                selection.to_csv(os.path.join(local_reconstruction_directory, "selection.csv"))
                
                with open(os.path.join(local_matches_directory, f"{ligand_info[1]}.gt"), "w") as f:
                    f.write("")

# %% Main
if __name__ == "__main__":
    create_scaffolds_and_reconstructions()