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

import pickle
import bn_utils # Needed for pickle import of fitted BNs
import numpy as np
import pandas as pd
import tmQM_RDF_interface as tmint

from tqdm import tqdm

INPUT_FILES = {
        "bn": os.path.join(ROOT_DIR, "computational", "bayesian_network", "results", "%s"),
        "substructures": os.path.join(ROOT_DIR, "computational", "pattern_clustering", "results", "%s", "by_metric"),
        "reconstructions": os.path.join(ROOT_DIR, "computational", "reconstruction", "intermediate", "reconstructions", "rdf", "%s"),
        "periodic_table": os.path.join(ROOT_DIR, "data", "raw", "pubChem", "data")
    }

OUTPUT_FILES = {
        "results": os.path.join(ROOT_DIR, "computational", "reconstruction", "results", "%s", "scores"),
        # "reconstructions_figures": os.path.join(ROOT_DIR, "computational", "reconstruction", "results", "%s", "figures", "probable_reconstructions")
    }

# Main component of dataset_tag = {dataset_short_comment}__s_{size_range[0]}_{size_range[1]}
TM_MODE = "earlyTM"
DATASET_SHORT_ID = "latmod" # Ligand ATtachment MODes

# Range of sizes of interes (extremes included)
SIZE_RANGE = [10, 12]

SETTINGS = {
        "fitted_graph_name": "fit__hc100alpha_5.pkl",
        "cluster_info_file_name": "cls_info.csv",
        "set_median_info_file_name": "cls_set_median_graphs.csv",
        #"figures_per_recon": 3
    }

# %% Utility functions
def identify_working_paths(dataset_tag, graph_file_name):
    """
    This function identifies all the subdirectories, starting from INPUT_FILES['bn'], that contain a file named as
    dist_matrix_file_name. Each of these subdirectories is marked as a potential working subdirectory.
    The user is then asked to confirm which directories to use.
    
    Arguments:
        - dataset_tag: the dataset_tag  needed to identify the correct intermediate directory
        - dist_matrix_file_name: the file name used to identify distance matrices and working directories
        
    Returns:
        - The working directories approved by the user
    """
    
    root = INPUT_FILES["bn"] % dataset_tag
    candidates = []
    
    for directory, _, files in os.walk(root):
        if graph_file_name in files:
            candidates += [directory.replace(root, "")[1:]]
            
    prompt = f"Candidate working paths:\n{'\n'.join([f'{i}) {d}' for i, d in enumerate(candidates)])}"
    prompt += "\n\n"
    prompt += " > Type the numbers of the directories to use: "
    
    selection = input(prompt)
    
    selected = []
    for c in selection:
        selected += [candidates[int(c)]]
        
    return selected

def aggregate_raw_matches(dataset_tag, raw_matches, wp):
    # Process raw matches
    matches = pd.read_csv(raw_matches, index_col = 0)
    index = matches.index
    matches = np.sign(np.array(matches))
    
    # Acquire cluster info
    info = pd.read_csv(
        os.path.abspath(os.path.join(INPUT_FILES["substructures"] % dataset_tag, wp, "..", SETTINGS["cluster_info_file_name"])), 
        index_col = 0
    )
    
    n_clusters = np.max(info.loc[:, "cluster"]) + 1
    
    # Determine if set_median info is available
    set_median_info = (SETTINGS["set_median_info_file_name"] in os.listdir(os.path.join(INPUT_FILES["substructures"] % dataset_tag, wp)))
    
    if not set_median_info:
        clusters = {c: np.array(info.loc[info["cluster"] == c, "index_in_df"]) for c in range(n_clusters)}
        aggregated = {str(i): np.apply_along_axis(np.max, 1, matches[:, c], 0) for i, c in clusters.items()}
    else:
        set_median = pd.read_csv(
            os.path.join(INPUT_FILES["substructures"] % dataset_tag, wp, SETTINGS["set_median_info_file_name"]),
            index_col = 0
        )
        
        aggregated = {}
        for i in range(n_clusters):
            set_median_index = set_median.loc[set_median["cluster"] == i, "set_median_graph_index"].iloc[0]
            aggregated[str(i)] = matches[:, set_median_index]
    
    return pd.DataFrame(aggregated, index = index)

def compute_bn_scores():
    # Assemble dataset tag
    if SIZE_RANGE[1] > SIZE_RANGE[0]:
        size_tag = f"s_{SIZE_RANGE[0]}_{SIZE_RANGE[1]}"
    else:
        size_tag = f"s_{SIZE_RANGE[0]}"
    
    dataset_tag = f"{TM_MODE}-{DATASET_SHORT_ID}-{size_tag}"
    
    # Get working paths
    working_paths = identify_working_paths(dataset_tag, SETTINGS["fitted_graph_name"])
    
    for wp in working_paths:
        print(f"Working with wp {wp}")
        
        with open(os.path.join(INPUT_FILES["bn"] % dataset_tag, wp, SETTINGS["fitted_graph_name"]), "rb") as f:
            bn = pickle.load(f)
    
            # Process reconstructions
            recons = [r for r in os.listdir(INPUT_FILES["reconstructions"] % TM_MODE) if r.startswith("recon")]
            
            scores = {}
            
            missing = 0
            for r, recon in tqdm(enumerate(recons)):
                if not recon.startswith("recon"):
                    continue
                
                # Retrieve ground truth
                files = os.listdir(os.path.join(INPUT_FILES["reconstructions"] % TM_MODE, recon))
                
                ground_truth = [x.split(".")[0] for x in files if x.endswith(".gt")][0]
                
                # Process raw matches
                aggregated_matches = aggregate_raw_matches(
                    dataset_tag,
                    os.path.join(INPUT_FILES["reconstructions"] % TM_MODE, recon, "matches.csv"),
                    wp
                )
                
                # Compute scores
                loc_scores = bn.score(aggregated_matches, log = True, unsafe = True)
                
                tmc_code = list(loc_scores.keys())[0].split("_")[0]
                
                scores[f"{tmc_code}__{ground_truth}"] = loc_scores
                
                # Save figures
                # fig_folder = os.path.join(os.path.join(OUTPUT_FILES["reconstructions_figures"] % dataset_tag, wp))
                # tmint.TmQMRDFInterface.path_to_tmQM_RDF = os.path.join(INPUT_FILES["reconstructions"] % TM_MODE, recon, "graphs")
                
                # if not os.path.exists(fig_folder):
                #     os.makedirs(fig_folder)
                
                # for i, k in enumerate(np.array(list(loc_scores.keys()))[perm][-SETTINGS["figures_per_recon"]:]):
                #     tmint.TmQMRDFGraph(k).render(
                #         filename = os.path.join(fig_folder, f"{SETTINGS['figures_per_recon'] - i}__{k.replace('-','_')}"), 
                #         layout = "neato"
                #     )
    
            results_folder = os.path.join(OUTPUT_FILES["results"] % dataset_tag, wp)
            if not os.path.exists(results_folder):
                os.makedirs(results_folder)
            
            with open(os.path.join(results_folder, "scores.txt"), "w") as f:
                f.write(str(scores))
        
# %% Main statement
if __name__ == "__main__":
    tmint.TmQMRDFGraph.path_to_chem_info = INPUT_FILES["periodic_table"]
    compute_bn_scores()