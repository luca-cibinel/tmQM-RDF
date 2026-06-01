"""
Compute the scores (i.e. the probabilities assigned by the BN) of the reconstructed TMCs
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header

import sys
sys.path.append(os.path.join(ROOT_DIR, "computational", "bayesian_network"))
sys.path.append(os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "interface"))

import re
import numpy as np
import pandas as pd
#import tmQM_RDF_interface as tmint

from collections import Counter
from openbabel import openbabel, pybel
from tqdm import tqdm

INPUT_FILES = {
        "scores": os.path.join(ROOT_DIR, "computational", "reconstruction", "results", "%s", "scores"),
        "ligands_misc_info": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", "v60k", "ligands_misc_info.csv"),
        "ligands_fingerprints": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", "v60k", "ligands_fingerprints.csv"),
        "ligands": os.path.join(ROOT_DIR, "computational", "reconstruction", "intermediate", "ligands", "%s")
    }

OUTPUT_FILES = {
        "results": os.path.join(ROOT_DIR, "computational", "reconstruction", "results", "%s", "ranks"),
        "reconstructions_figures": os.path.join(ROOT_DIR, "computational", "reconstruction", "results", "%s", "figures", "probable_reconstructions")
    }

# Main component of dataset_tag = {dataset_short_comment}__s_{size_range[0]}_{size_range[1]}
#TM_MODE = "earlyTM"
DATASET_SHORT_ID = "latmod" # Ligand ATtachment MODes

# Range of sizes of interes (extremes included)
SIZE_RANGE = [10, 12]

SETTINGS = {
        "scores_file_name": "scores.txt",
        "filters": ["n", "c", "hd"] # None ; Charge ; Hapticity and Denticity order
    }

# %% Utility functions
def identify_working_paths(dataset_tag, scores_file_name):
    """
    This function identifies all the subdirectories, starting from INPUT_FILES['bn'], that contain a file named as
    dist_matrix_file_name. Each of these subdirectories is marked as a potential working subdirectory.
    The user is then asked to confirm which directories to use.
    
    Arguments:
        - dataset_tag: the dataset_tag  needed to identify the correct intermediate directory
        - scores_file_name: the file name used to identify scores dictionaries and working directories
        
    Returns:
        - The working directories approved by the user
    """
    
    root = INPUT_FILES["scores"] % dataset_tag
    candidates = []
    
    for directory, _, files in os.walk(root):
        if scores_file_name in files:
            candidates += [directory.replace(root, "")[1:]]
            
    prompt = f"Candidate working paths:\n{'\n'.join([f'{i}) {d}' for i, d in enumerate(candidates)])}"
    prompt += "\n\n"
    prompt += " > Type the numbers of the directories to use: "
    
    selection = input(prompt)
    
    selected = []
    for c in selection:
        selected += [candidates[int(c)]]
        
    return selected

def compute_binding_modalities(tm_mode):
    ligands = [x.split(".")[0] for x in os.listdir(os.path.join(INPUT_FILES["ligands"] % tm_mode, "train")) if x.endswith(".nt")]

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

def retrieve_ligand_filtering_info(tm_mode):
    ligands = [x.split(".")[0] for x in os.listdir(os.path.join(INPUT_FILES["ligands"] % tm_mode, "train")) if x.endswith(".nt")]
    
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

def compute_rankings(tm_mode):
    # Assemble dataset tag
    if SIZE_RANGE[1] > SIZE_RANGE[0]:
        size_tag = f"s_{SIZE_RANGE[0]}_{SIZE_RANGE[1]}"
    else:
        size_tag = f"s_{SIZE_RANGE[0]}"
    
    dataset_tag = f"{tm_mode}-{DATASET_SHORT_ID}-{size_tag}"
    
    # Get modalities/filtering info
    ligands_modalities = compute_binding_modalities(tm_mode)
    ligands_filtering_info = retrieve_ligand_filtering_info(tm_mode)
    
    # Get working paths
    working_paths = identify_working_paths(dataset_tag, SETTINGS["scores_file_name"])
    
    recon_sizes = {f: [] for f in SETTINGS["filters"]}
    ranks = {f: {wp: {} for wp in working_paths} for f in SETTINGS["filters"]}
    
    for wp in working_paths:
        print(f"Working with wp {wp}")
        
        with open(os.path.join(INPUT_FILES["scores"] % dataset_tag, wp, SETTINGS["scores_file_name"]), "r") as f:
            scores = eval(f.read())
    
            # Process reconstructions
            for ground_truth, recon_scores in tqdm(scores.items()):
                ligand_gt = ground_truth.split("_")[-1]
                
                for fltr in SETTINGS["filters"]:
                    if ground_truth in recon_scores:
                        viable_ligands = [
                                l 
                                for l, info in ligands_filtering_info.items()
                                if info[fltr] == ligands_filtering_info[ligand_gt][fltr]
                            ]
                        
                        loc_scores = {
                                r: s 
                                for r, s in recon_scores.items()
                                if r.split("_")[-1] in viable_ligands
                            }

                        sscores = sorted(set(loc_scores.values()), reverse = True)
                        
                        # With regards to the notation used in Saha et al. (2022):
                        #  c[sscores[i]] = n_{i + 1} = number of items at rank i
                        #    (Ranks start from 1, indexes start from 0, hence why scores_to_level introduces a +1 shift)
                        #  N[i] = N_i = total number of items retrieved up to rank i (included)
                        c = {
                            s: len(set([
                                ligands_modalities[recon_name.split("__")[1]] for recon_name, x in loc_scores.items() if x == s
                            ])) 
                            for s in sscores
                        }
                        #c = Counter(loc_scores.values())

                        scores_to_level = {
                            s: i + 1 for i, s in enumerate(sscores)
                        }

                        N = {
                            i + 1: sum([c[x] for x in sscores if x >= sscores[i]])
                            for i in range(len(c))
                        }

                        N[0] = 0

                        correct_level = scores_to_level[recon_scores[ground_truth]]

                        tmc_code = ground_truth.split("_")[0]
                        ranks[fltr][wp][tmc_code] = (N[correct_level - 1], c[recon_scores[ground_truth]])
                        
                        if wp == working_paths[0]:
                            recon_sizes[fltr] += [len(viable_ligands)]
    
    for fltr in SETTINGS["filters"]:
        results_folder = os.path.join(OUTPUT_FILES["results"] % dataset_tag, fltr)
        
        if not os.path.exists(results_folder):
            os.makedirs(results_folder)
        
        with open(os.path.join(results_folder, "ranks.txt"), "w") as f:
            f.write(str(ranks[fltr]))
            
        size = np.array(recon_sizes[fltr])
        print(f"Filter: {fltr};\tMean: {np.mean(size)};\tSd: {np.std(size)}")
        
# %% Main statement
if __name__ == "__main__":
    compute_rankings("lateTM")
    compute_rankings("earlyTM")