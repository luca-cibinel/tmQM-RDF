
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
sys.path.append(os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "interface", "v2025dev"))

import numpy as np
import pandas as pd
import tmQM_RDF_interface as tmint

from tqdm import tqdm

INPUT_FILES = {
        "tmQM-RDF": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev", "graphs"),
        "scores": os.path.join(ROOT_DIR, "computational", "reconstruction", "case_studies", "results", "%s", "scores"),
        "ligands_misc_info": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", "v60k", "ligands_misc_info.csv"),
        "ligands_fingerprints": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", "v60k", "ligands_fingerprints.csv"),
        "reconstructions": os.path.join(ROOT_DIR, "computational", "reconstruction", "case_studies", "intermediate", "reconstructions", "rdf", "%s"),
        "ligands": os.path.join(ROOT_DIR, "computational", "reconstruction", "intermediate", "ligands", "%s"),
        "periodic_table": os.path.join(ROOT_DIR, "data", "raw", "pubChem", "data")
    }

OUTPUT_FILES = {
        "results": os.path.join(ROOT_DIR, "computational", "reconstruction", "case_studies", "results", "%s", "figures") # + wp + TMC__gt + filter + [low/high]
    }

CASE_STUDIES = {
        "earlyTM": ["CAKRAH"],
        "lateTM": ["ILUYOD"]
    }

# Main component of dataset_tag = {dataset_short_comment}__s_{size_range[0]}_{size_range[1]}
DATASET_SHORT_ID = "latmod" # Ligand ATtachment MODes

# Range of sizes of interes (extremes included)
SIZE_RANGE = [10, 12]

SETTINGS = {
        "scores_file_name": "scores.txt",
        "filters": ["n", "c", "d", "h"] # None ; Charge ; Denticity order; Hapticity order
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

def retrieve_ligand_filtering_info(tm_mode):
    """
    This function collects, for each ligand in the specified tm_mode (early/late), the information needed to implement the following filters:
    - n: no filter (added for compatibility)
    - c: preserve charge
    - d: preserve denticity order (the candidate ligand must also have 0 hapticity order)
    - h: preserve hapticity order (the candidate ligand must also have 0 denticity order and only 1 haptic binding site)

    Arguments:
        - tm_mode: the tm mode for which information must be collected (either "early" or "late")
    
    Returns:
        - a dictionary of the form {ligand_id: {filter_name: <info for filtering> ...} ...}
        - a dictionary of the form {filter_name: check_filter, ...} where check_filter is a function that takes two arguments: the first is
            the ORIGINAL ligand's filtering info (contained in the first output of retrieve_ligand_filtering_info) and the second is the CANDIDATE
            replacement's filtering info; it returns a boolean (whether the filter is satisfied or not).
    """
    ligands = [x.split(".")[0] for x in os.listdir(os.path.join(INPUT_FILES["ligands"] % tm_mode, "train")) if x.endswith(".nt")]
    
    tmint.TmQMRDFGraph.path_to_tmQM_RDF = INPUT_FILES["tmQM-RDF"]
    for _, case_studies in CASE_STUDIES.items():
        for case_tmc in case_studies:
            g = tmint.TmQMRDFGraph(case_tmc)
            
            ligands += [g.ligands[ln]["class"] for ln in g.ligands if g.ligands[ln]["class"] not in ligands]
    
    misc_info = pd.read_csv(INPUT_FILES["ligands_misc_info"], sep = ";")
    descriptors = pd.read_csv(INPUT_FILES["ligands_fingerprints"], sep = ";")
    
    denticity = lambda idx: sum(len(site) == 1 for site in idx)
    hapticity = lambda idx: sum(len(site) for site in idx if len(site) > 1)
    haptic_sites = lambda idx: sum(len(site) > 1 for site in idx)

    return {
            l : {
                    "n": None,
                    "c": descriptors.loc[descriptors["name"] == l, "charge"].iloc[0],
                    "d": eval(misc_info.loc[misc_info["name"] == l, "smiles_metal_bond_node_idx_groups"].iloc[0]),
                    "h": eval(misc_info.loc[misc_info["name"] == l, "smiles_metal_bond_node_idx_groups"].iloc[0])
                }
            for l in ligands
        }, {
            "n": lambda l0, l_cand: True,
            "c": lambda c0, c_cand: c0 == c_cand,
            "d": lambda idx0, idx_cand: denticity(idx0) == denticity(idx_cand) and hapticity(idx_cand) == 0,
            "h": lambda idx0, idx_cand: hapticity(idx0) == hapticity(idx_cand) and denticity(idx_cand) == 0 and haptic_sites(idx_cand) == 1
        }

def plot_selected_reconstructions():
    for tm_mode in CASE_STUDIES:
        # Assemble dataset tag
        if SIZE_RANGE[1] > SIZE_RANGE[0]:
            size_tag = f"s_{SIZE_RANGE[0]}_{SIZE_RANGE[1]}"
        else:
            size_tag = f"s_{SIZE_RANGE[0]}"
        
        dataset_tag = f"{tm_mode}-{DATASET_SHORT_ID}-{size_tag}"
        
        # Get filtering info
        ligands_filtering_info, check_filter = retrieve_ligand_filtering_info(tm_mode)
        
        # Get working paths
        working_paths = identify_working_paths(dataset_tag, SETTINGS["scores_file_name"])
        
        for wp in working_paths:
            print(f"Working with wp {wp}")
            
            with open(os.path.join(INPUT_FILES["scores"] % dataset_tag, wp, SETTINGS["scores_file_name"]), "r") as f:
                scores = eval(f.read())
        
                # Process reconstructions
                for ground_truth, recon_scores in scores.items():
                    print(f"\tCurrent gt: {ground_truth}")

                    tmc_gt = ground_truth.split("_")[0]
                    ligand_gt = ground_truth.split("_")[-1]
                    
                    for fltr in SETTINGS["filters"]:
                        print(f"\t\tCurrent filter: {fltr}")

                        # Filter scores by ligand properties
                        
                        fltr_gt = ligands_filtering_info[ligand_gt][fltr]

                        viable_ligands = [
                                l 
                                for l, info in ligands_filtering_info.items()
                                if check_filter[fltr](fltr_gt, info[fltr])
                            ]
                    
                        loc_scores = {
                                r: s 
                                for r, s in recon_scores.items()
                                if r.split("_")[-1] in viable_ligands
                            }
                        
                        print(f"\t\t\tNum. of viable ligands: {len(loc_scores)}")

                        if len(loc_scores) == 0:
                            continue

                        # Process reconstructions
                        local_results_root = os.path.join(OUTPUT_FILES["results"] % dataset_tag, wp, ground_truth, fltr)
                        
                        tmint.TmQMRDFGraph.path_to_tmQM_RDF = os.path.join(INPUT_FILES["reconstructions"] % tm_mode, tmc_gt, f"recon_{ligand_gt}", "graphs")
                        
                        for fname, f in [("high_prob", max), ("low_prob", min)]:
                            local_results_dir = os.path.join(local_results_root, fname)
                            print(f"\t\t\t\tMode: {fname}")

                            if not os.path.exists(local_results_dir):
                                os.makedirs(local_results_dir)
                            
                            target_score = f(loc_scores.values())
                            
                            viable_recons = [r for r, s in loc_scores.items() if s == target_score]
                            print(f"\t\t\t\t\tNum of viable recons: {len(viable_recons)}")

                            for r in viable_recons:
                                tmint.TmQMRDFGraph(r).render(
                                        filename = os.path.join(local_results_dir, r), 
                                        layout = "neato"
                                    )
        
# %% Main statement
if __name__ == "__main__":
    tmint.TmQMRDFGraph.path_to_chem_info = INPUT_FILES["periodic_table"]
    plot_selected_reconstructions()
    
    