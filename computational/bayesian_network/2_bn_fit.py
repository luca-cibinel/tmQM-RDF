"""
Fit the learned BN to the data
"""
# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header
import pickle
import bn_utils
import pandas as pd

INPUT_FILES = {
        "intermediate": os.path.join(ROOT_DIR, "computational", "bayesian_network", "intermediate", "%s"),
        "substructures": os.path.join(ROOT_DIR, "computational", "pattern_clustering", "results", "%s", "by_metric"),
        "tmQM-RDF-1Ksel": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "data", "v2025dev", "%s", "1k_selection.csv")
    }

OUTPUT_FILES = {
        "results": os.path.join(ROOT_DIR, "computational", "bayesian_network", "results", "%s")
    }

SETTINGS = {
        "graph_file_name": "hc100alpha_5.gml",
        "matches_file_name": "cls_matches.csv",
        "bn_fit_settings": {
            "prior_type": "dirichlet", 
            "pseudo_counts": 10,
            "n_jobs": 8
        },
        "fitted_name": "fit__%s.pkl"
    }

# Main component of dataset_tag = {dataset_short_comment}__s_{size_range[0]}_{size_range[1]}
PATTERN_DATASET_KEYWORD = "earlyTM"
DATASET_SHORT_ID = "latmod" # Ligand ATtachment MODes

# Range of sizes of interes (extremes included)
SIZE_RANGE = [10, 12]

# %% Main functions

def identify_working_paths(dataset_tag, graph_file_name):
    """
    This function identifies all the subdirectories, starting from INPUT_FILES['intermediate'], that contain a file named as
    graph_file_name. Each of these subdirectories is marked as a potential working subdirectory.
    The user is then asked to confirm which directories to use.
    
    Arguments:
        - dataset_tag: the dataset_tag  needed to identify the correct intermediate directory
        - graph_file_name: the file name used to identify distance matrices and working directories
        
    Returns:
        - The working directories approved by the user
    """
    
    root = INPUT_FILES["intermediate"] % dataset_tag
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

def fit_bn_networks():
    # Get dataset selection and train/validation/test partition
    selection = pd.read_csv(INPUT_FILES["tmQM-RDF-1Ksel"] % PATTERN_DATASET_KEYWORD, index_col = 0)
    train_idxs = selection.loc[selection["partition"] == "train", "TMC"]
    
    # Assemble dataset tag
    if SIZE_RANGE[1] > SIZE_RANGE[0]:
        size_tag = f"s_{SIZE_RANGE[0]}_{SIZE_RANGE[1]}"
    else:
        size_tag = f"s_{SIZE_RANGE[0]}"
    
    dataset_tag = f"{PATTERN_DATASET_KEYWORD}-{DATASET_SHORT_ID}-{size_tag}"
    print(f"[] Dataset: {dataset_tag}")
    
    # Identify working paths
    working_paths = identify_working_paths(dataset_tag, SETTINGS["graph_file_name"])
    
    # Process each working path
    for wp in working_paths:
        print(f"Working with wp {wp}")
        
        # Retrieve data and graph
        data = pd.read_csv(os.path.join(INPUT_FILES["substructures"] % dataset_tag, wp, SETTINGS["matches_file_name"]), index_col = 0)
        data = data.loc[train_idxs, :]
        
        graph = os.path.join(INPUT_FILES["intermediate"] % dataset_tag, wp, SETTINGS["graph_file_name"])
        
        # Fit BN
        bn = bn_utils.BayesianNetwork(graph)
        bn.fit(data, **SETTINGS["bn_fit_settings"])
        
        # Save results
        local_results_dir = os.path.join(OUTPUT_FILES["results"] % dataset_tag, wp)
        
        if not os.path.exists(local_results_dir):
            os.makedirs(local_results_dir)
        
        gname = SETTINGS["graph_file_name"].split(".")[0]
        with open(os.path.join(local_results_dir, SETTINGS["fitted_name"] % gname), "wb") as f:
            pickle.dump(bn, f)

# %% Main statement
if __name__ == "__main__":
    fit_bn_networks()