"""
Step 3 in the creation of the "1k" selection.

It creates the 1k selection (1000: train; 300: validation; 300; test) according to the following scheme:
    - First, a "seed" made of the 350 (lateTM) / 1350 (earlyTM) most frequent ligands is extracted
    - Then, a candidate set of ligands is built so that:
        * only TMCs whose metal core is among the three most frequent cores in tmQM-RDF are allowed
        * only TMCs whose ligands are all found in the seed are allowed
    - Finally, 1000/300/300 TMCs are sampled from the candidate set in order to form the 1k selection
        (the relative proportions of the metal cores in each partition of the selection are artificially matched to those of the
         three most frequent cores in tmQM-RDF)
        
The outliers identified by Kneiding et al (2023) are excluded using the outliers.txt file
available in the GitHub repository https://github.com/uiocompcat/tmQMg/
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header
from tqdm import tqdm

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np

INPUT_FILES = {
        "viable": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "viable_tmcs.txt"),
        "centres": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", "v2025dev", "centres.csv"),
        "details": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", "v2025dev", "details.csv"),
        "sizes": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", "v2025dev", "tmc_sizes.txt"),
        "ligands": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", "v2025dev", "%s", "ligands.csv"),
        "outliers": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", "v2025dev", "outliers.txt"),
        "tmQMg_graphs": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg", "v74.637k", "uNatQ_graphs")
    }

OUTPUT_FILES = {
        "selection": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "data", "v2025dev", "%s", "1k_selection.csv")
    }

RANDOM_SEED = 2568576931
RNG = np.random.RandomState(RANDOM_SEED)

SEED_SIZE = {"lateTM": 350, "earlyTM": 1350}

SIZES = [1000, 300, 300]

# %% Utility functions
def read_input_info(tm_mode):
    """
    Read the various information:
        - Metal centre frequency (see extract_info_from_tmc_for_selection.py)
        - Detailed TMC info (centre and ligands) (see extract_info_from_tmc_for_selection.py)
        - Ligands frequency (see compute_ligand_occurrences.py)
        - TMC sizes
    
    Arguments:
        - tm_mode: whether to process early TMs (mode = 'earlyTM') or the late ones (mode = 'lateTM')
    
    Returns:
        - viable: the list of viable TMCs
        - centre_dist: the pd.DataFrame of metal centre counts
        - detailed_info: the pd.DataFrame with detailed TMC info
        - lig_count_info: the pd.DataFrame with the ligand counts
        - size_dist: a dictionary with the size (in n. of atoms) of each TMC
        - outliers: the list of tmQMg outliers
    """

    with open(INPUT_FILES["viable"], "r") as f:
        viables = f.read().split("\n")
    
    centre_dist = pd.read_csv(INPUT_FILES["centres"], index_col = 0)
    
    detailed_info = pd.read_csv(INPUT_FILES["details"], index_col = 0)
    
    lig_count_info = pd.read_csv(INPUT_FILES["ligands"] % tm_mode, index_col = 0)
    
    with open(INPUT_FILES["sizes"], "r") as f:
        size_dist = eval(f.read())
        
    with open(INPUT_FILES["outliers"], "r") as f:
        outliers = f.read().split("\n")
        
    return viables, centre_dist, detailed_info, lig_count_info, size_dist, outliers

def extract_candidate_tmc(lig_seed_size, tm_mode, detailed_info, lig_count_info, viables, outliers):
    """
    Selects from tmQM-RDF the TMCs which are eligible to be part of the 1k selection.
    
    Arguments:
        - lig_seed_size: the size of the ligand seed
        - tm_mode: whether to process early TMs (mode = 'earlyTM') or the late ones (mode = 'lateTM')
        - detailed_info: the pd.DataFrame with detailed TMC info (see read_input_info)
        - lig_count_info: the pd.DataFrame with the ligand counts (see read_input_info)
        - viables: the pd.DataFrame with the ligand counts (see read_input_info)
        - outliers: the list of tmQMg outliers (see read_input_info)
        
    Returns:
        - candidate_tmc: the list of candidate TMCs
    """
    
    # Find the most frequent ligands and store them in a seed
    q = np.argsort(lig_count_info.loc[:, "count"])[-lig_seed_size:]
    seed = np.array(lig_count_info.index)[q]

    # Extract all viable TMCs with compatible cores such that each ligand is found in the seed.
    if tm_mode == "earlyTM":
        centres = ["Cr", "Mo", "W"]
    elif tm_mode == "lateTM":
        centres = ["Pd", "Ni", "Pt"]
    else:
        raise Exception(f"tm_mode {tm_mode} not recognised!")
    
    candidates_tmc = []

    viables = [tmc for tmc in viables if detailed_info.loc[tmc, "core"] in centres and tmc not in outliers]

    for tmc in viables:
        ligs = detailed_info.loc[tmc, "ligands"].split(" ")
        
        if len([ l for l in ligs if l in seed ]) == len(ligs):
            candidates_tmc += [tmc]
            
    return candidates_tmc

def extract_selection(candidates_tmc, tm_mode, sizes, centre_dist, detailed_info, size_dist, random_state):
    """
    Computes the exact selection
    
    Arguments:
        - candidates_tmc: the list of candidate TMCs (see extract_candidate_tmc)
        - tm_mode: whether to process early TMs (mode = 'earlyTM') or the late ones (mode = 'lateTM')
        - sizes: an array with the required sizes of the train/validation/test partitions of the selection
        - centre_dist: the pd.DataFrame of metal centre counts (see read_input_info)
        - detailed_info: the pd.DataFrame with detailed TMC info (see read_input_info)
        - size_dist: a dictionary with the size (in n. of atoms) of each TMC (see read_input_info)
        
    Returns:
        - selection_train: the TMCs in the train partition
        - selection_validation: the TMCs in the validation partition
        - selection_test: the TMCs in the test partition
    """
    
    rng = np.random.RandomState(random_state)
    train_size, validation_size, test_size = sizes
    
    # Determine the number of TMCs to sample per metal centre
    if tm_mode == "earlyTM":
        centres = ["Cr", "Mo", "W"]
    elif tm_mode == "lateTM":
        centres = ["Pt", "Ni", "Pd"]
    else:
        raise Exception(f"tm_mode {tm_mode} not recognised!")
    
    mc_counts = np.array([centre_dist.loc[c, "count"] for c in centres])
    
    sizes_by_mc = np.int64(np.round(np.diag(sizes) @ np.array([mc_counts]*len(mc_counts)) / sum(mc_counts)))
    # sizes_by_mc[partition, element] = # of TMCs with MC = element to put into partition
    # print(sizes_by_mc)
    # Computes the total number of TMC to sample for each metal centre BEFORE the train/validation/test split
    total_sample_size = sizes_by_mc.sum(0)
    
    # For each desired metal core, sample the appropriate number of TMCs
    # with probability inversely proportional to the TMC size
    
    selection_train = []
    selection_validation = []
    selection_test = []
    
    for i in range(len(centres)):
        el = centres[i]
        
        candidates_by_el = [tmc for tmc in candidates_tmc if detailed_info.loc[tmc, "core"] == el]
        
        prob = np.array([1/size_dist[tmc] for tmc in candidates_by_el])
        prob /= np.sum(prob)
        
        selection_by_el = rng.choice(
                candidates_by_el, total_sample_size[i], replace = False, p = prob
            )
        rng.shuffle(selection_by_el)
        
        selection_train += list(selection_by_el[:sizes_by_mc[0,i]])
        selection_validation += list(selection_by_el[sizes_by_mc[0,i]:(sizes_by_mc[0,i] + sizes_by_mc[1,i])])
        selection_test += list(selection_by_el[(sizes_by_mc[0,i] + sizes_by_mc[1,i]):])

    selection = {
            "TMC": selection_train + selection_validation + selection_test,
            "partition": ["train"]*train_size + ["validation"]*validation_size + ["test"]*test_size
        }
    
    # Write the selection to file 
    pd.DataFrame(selection).to_csv(OUTPUT_FILES["selection"] % tm_mode)

    return selection_train, selection_validation, selection_test

def plot_summary(selection_train, selection_validation, selection_test, size_dist):
    """
    Visualizes the distribution of the metal centre counts and TMC sizes inside the selection
    
    Arguments:
        - selection_train: the TMCs in the train partition
        - selection_validation: the TMCs in the validation partition
        - selection_test: the TMCs in the test partition
        - size_dist: a dictionary with the size (in n. of atoms) of each TMC (see read_input_info)
    """

    # Metal centre counts
    fig, (ax_train, ax_validation, ax_test) = plt.subplots(1, 3)

    # Train
    cores_in_sel = {}

    for tmc in selection_train:
        g = nx.read_gml(os.path.join(INPUT_FILES["tmQMg_graphs"], tmc + ".gml"))
        
        met_core = g.graph["meta_data"]["metal_center_element"]

        cores_in_sel[met_core] = cores_in_sel.get(met_core, 0) + 1
        
    ax_train.bar(cores_in_sel.keys(), cores_in_sel.values())

    # Validation
    cores_in_sel = {}

    for tmc in selection_validation:
        g = nx.read_gml(os.path.join(INPUT_FILES["tmQMg_graphs"], tmc + ".gml"))
        
        met_core = g.graph["meta_data"]["metal_center_element"]

        cores_in_sel[met_core] = cores_in_sel.get(met_core, 0) + 1
        
    ax_validation.bar(cores_in_sel.keys(), cores_in_sel.values())

    # Test
    cores_in_sel = {}

    for tmc in selection_test:
        g = nx.read_gml(os.path.join(INPUT_FILES["tmQMg_graphs"], tmc + ".gml"))
        
        met_core = g.graph["meta_data"]["metal_center_element"]

        cores_in_sel[met_core] = cores_in_sel.get(met_core, 0) + 1
        
    ax_test.bar(cores_in_sel.keys(), cores_in_sel.values())

    fig.show()

    # TMC sizes
    fig, (ax_train, ax_validation, ax_test) = plt.subplots(1, 3, sharey = True)

    # Train
    sizes = [size_dist[tmc] for tmc in selection_train]
    ax_train.hist(sizes, density = True, bins = 20)
    ax_train.title.set_text("Train")
    ax_train.set_xlabel("Size")
    ax_train.set_ylabel("Density")

    # Validatioon
    sizes = [size_dist[tmc] for tmc in selection_validation]
    ax_validation.hist(sizes, density = True, bins = 20)
    ax_validation.title.set_text("Validation")
    ax_validation.set_xlabel("Size")

    # Test
    sizes = [size_dist[tmc] for tmc in selection_test]
    ax_test.hist(sizes, density = True, bins = 20)
    ax_test.title.set_text("Test")
    ax_test.set_xlabel("Size")

    fig.show()

# %% Main statement
if __name__ == "__main__":
    print("late TM")
    viables, centre_dist, detailed_info, lig_count_info, size_dist, outliers = read_input_info("lateTM")
    candidates_tmc = extract_candidate_tmc(SEED_SIZE["lateTM"], "lateTM", detailed_info, lig_count_info, viables, outliers)
    selection_train, selection_validation, selection_test = extract_selection(candidates_tmc, "lateTM", SIZES, centre_dist, detailed_info, size_dist, RANDOM_SEED)
    plot_summary(selection_train, selection_validation, selection_test, size_dist)
    
    print("early TM")
    viables, centre_dist, detailed_info, lig_count_info, size_dist, outliers = read_input_info("earlyTM")
    candidates_tmc = extract_candidate_tmc(SEED_SIZE["earlyTM"], "earlyTM", detailed_info, lig_count_info, viables, outliers)
    selection_train, selection_validation, selection_test = extract_selection(candidates_tmc, "earlyTM", SIZES, centre_dist, detailed_info, size_dist, RANDOM_SEED)
    plot_summary(selection_train, selection_validation, selection_test, size_dist)