"""
Compute and process pattern match data based on pattern interestingess
"""
# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header
from collections import defaultdict
from datetime import date
from tqdm import tqdm

import yaml
import shutil
import patterns
import numpy as np
import rdflib as rdf
import rpy2.robjects as robjects

import warnings
np.warnings = warnings

INPUT_FILES = {
        "tmQM-RDF": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev"),
        "patterns": os.path.join(ROOT_DIR, "computational", "pattern_mining", "results", "%s"),
        "tmQM-RDF-1Ksel": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "data", "%s", "1k_selection.csv"),
        "vos_db_dir": os.path.join(ROOT_DIR, "computational", "pattern_clustering", "temp", "vos_db"),
        "R": {
                "tools": os.path.join(ROOT_DIR, "computational", "vos_tools", "vos_interface.R"),
                "domination": os.path.join(ROOT_DIR, "computational", "pattern_clustering", "aux", "compute_pattern_dominations.R"),
                "matches": os.path.join(ROOT_DIR, "computational", "pattern_clustering", "aux", "compute_pattern_matches.R")
            }
    }

OUTPUT_FILES = {
        "results": os.path.join(ROOT_DIR, "computational", "pattern_clustering", "results", "%s")
    }
# path_to_tmQM_RDF = "./../../../../data/tmQM-RDF"

# Main component of dataset_tag = {dataset_short_comment}__s_{size_range[0]}_{size_range[1]}
PATTERN_DATASET_KEYWORD = "lateTM"
DATASET_SHORT_ID = "latmod" # Ligand ATtachment MODes

# Range of sizes of interes (extremes included)
SIZE_RANGE = [10, 12]

# Number of cores to be used for R computations
N_R_CORES = 4

# Directory where processed dataset has to be stored (under {dataset_root_directory}/{dataset_tag}/{version})
# dataset_root_directory = "./../../processed_pattern_data/"

# Directory where raw pattern data is stored
# target_directory = "./../results/ligand_attachment_modes/large/"

# Path to the .csv file containing the desired dataset selection and the train/validation/test split
# tmQM_RDF_selection = "./../../../../data/extract_selection_from_tmQM-RDF/1k_selection.csv"

# Directory where to keep the VOS db files
# vos_db_dir = './temp_/temp_vos_db'

# Auxiliary R scripts
# R_FILES = {
#         "domination": "./vos/compute_pattern_dominations.R",
#         "matches": "./vos/compute_pattern_matches.R"
#     }

# Operational settings
T = True
F = False
DO_EXTRACT_AND_SAVE_PATTERNS = F     # Extract patterns from raw pattern mining data
DO_R_PREPROCESSING_DOMINATION = F    # R preprocessing: compute domination relationship
DO_CHECK_PATTERN_CONTENT = T         # Evaluate pattern content
DO_FILTER_WORKING_DIRECTORY = T      # Filter working directory: remove patterns that do not fit criteria
DO_R_POSTPROCESSING_MATCHES = T      # R postprocessing: compute matches

# %% Main functions
def prepare_working_tree(root, dirs):
    """
    This function prepares the working tree by performing the following operations:
        - Purges the content of all the directories listed in dirs
        - Eliminates all the temp files (filename starts with .) found within root and its subdirectories
    """
    
    # Identify temp files
    temp_to_remove = []
    
    for directory, _, files in os.walk(root):
        for f in files:
            if f.startswith("."):
                temp_to_remove += [os.path.join(directory, f)]
    
    # Confirm user intentions
    temp_file_msg = f"""
    This code will delete the following files, deemed as temporary:
        {'\n'.join(temp_to_remove)}
    """
    
    dirs_msg = f"""
    This code will purge the content of the following directories:
        {'\n'.join(dirs)}
    """
    
    if len(temp_to_remove) == 0:
        temp_file_msg = ""
        
    if len(dirs) == 0:
        dirs_msg = ""
    
    safety_check = input(
        f"""{'-'*20}
            WARNING! 
            {dirs_msg}
            {temp_file_msg}
            > Type 'Yes' to confirm that you understand and you wish to continue: """)
            
    assert safety_check == "Yes"
    
    print("-"*20)
    
    # Purge directories
    for d in dirs:
    
        if os.path.exists(d):
            shutil.rmtree(d)
        
        os.makedirs(d)
        
    # Delete temp files
    for f in temp_to_remove:
        os.remove(f)

def initialise():
    # Assemble dataset tag
    if SIZE_RANGE[1] > SIZE_RANGE[0]:
        size_tag = f"s_{SIZE_RANGE[0]}_{SIZE_RANGE[1]}"
    else:
        size_tag = f"s_{SIZE_RANGE[0]}"
    
    dataset_tag = f"{PATTERN_DATASET_KEYWORD}-{DATASET_SHORT_ID}-{size_tag}"
    print(f"[] Dataset: {dataset_tag}")
    
    # Root directory of current dataset processing
    local_dataset_directory = os.path.join(OUTPUT_FILES["results"] % dataset_tag)
    
    # Directory where processed patterns should be stored
    pattern_storage_directory = os.path.join(local_dataset_directory, "patterns")

    # Directory where clustered patterns should be stored
    cluster_storage_directory = os.path.join(local_dataset_directory, "clusters_figures") # TODO: should this be here?
    
    # Directory where the final results should be stored
    results_directory = os.path.join(local_dataset_directory)
    
    # Prepare R environment
    robjects.globalenv["ROOT.DIR"] = ROOT_DIR
    robjects.globalenv["tmQM.RDF.selection"] = INPUT_FILES["tmQM-RDF-1Ksel"] % PATTERN_DATASET_KEYWORD
    robjects.globalenv["local.dataset.dir"] = local_dataset_directory
    robjects.globalenv["pattern.dir"] = pattern_storage_directory
    robjects.globalenv["results.dir"] = results_directory
    robjects.globalenv["dataset.tag"] = dataset_tag
    robjects.globalenv["vos.db.dir"] = INPUT_FILES["vos_db_dir"]
    robjects.globalenv["path.to.tmQM.RDF"] = INPUT_FILES["tmQM-RDF"]
    robjects.globalenv["vos.interface"] = INPUT_FILES["R"]["tools"]
    robjects.globalenv["N.CORES"] = N_R_CORES
    
    # Prepare working tree
    dirs_to_prepare = []
    
    if not os.path.exists(results_directory):
        dirs_to_prepare += [results_directory]
    if DO_EXTRACT_AND_SAVE_PATTERNS:
        dirs_to_prepare += [pattern_storage_directory]
    if DO_R_PREPROCESSING_DOMINATION or DO_R_POSTPROCESSING_MATCHES:
        dirs_to_prepare += [INPUT_FILES["vos_db_dir"]]
    
    prepare_working_tree(local_dataset_directory, dirs_to_prepare)
    
    # Read pattern batch
    pattern_batch = patterns.PatternBatch(INPUT_FILES["patterns"] % PATTERN_DATASET_KEYWORD, log = True)
    
    # Save config (for reproducibility)
    config = {
            "date": str(date.today()),
            "dataset_short_id": DATASET_SHORT_ID,
            "size_range": SIZE_RANGE,
            "dataset_tag": dataset_tag,
            "children_directories": {
                    "local_dataset_directory": local_dataset_directory,
                    "pattern_storage_directory": pattern_storage_directory,
                    "cluster_storage_directory": results_directory,
                    "results_directory": results_directory
                },
            "INPUT_FILES": INPUT_FILES,
            "dirs_to_prepare": dirs_to_prepare
        }
    
    config_fname = "config_processing.yml"
    
    with open(os.path.join(local_dataset_directory, config_fname), "w") as f:
        yaml.dump(config, f)
        
    return local_dataset_directory, pattern_storage_directory, cluster_storage_directory, results_directory, pattern_batch

def extract_and_save_patterns(pattern_batch, pattern_storage_directory):
    """
    Extract and save all mined patterns to allow for preprocessing via R.
    The patterns mined via the perl script are processed and converted into
    SELECT queries and a RDF graphs (with variables replaced by appropriate URIs; needed
    for dominations computation).
    
    The resulting directory will have the following structure:
        pattern_storage_directory
            |
            |-> queries
            |     |
            |     |-> ...sizes...
            |
            |-> rdf
            |     |
            |     |-> ...sizes...
    
    Arguments:
        - pattern_batch: a patterns.PatternBatch object containing the patterns of interest
        - pattern_storage_directory: the directory in which the patterns have to be stored
    """
    
    sizes_to_process = range(SIZE_RANGE[0], SIZE_RANGE[1] + 1)
    
    for size in sizes_to_process:
    
        print(f"Processing size {size} ({len(pattern_batch[size])} to process)")
    
        """
        Store the patterns as explicit SPARQL queries/RDF graphs
        """
        
        if not os.path.exists(os.path.join(pattern_storage_directory, f"queries/{size}/")):
            os.makedirs(os.path.join(pattern_storage_directory, f"queries/{size}/"))
            
        if not os.path.exists(os.path.join(pattern_storage_directory, f"rdf/{size}/")):
            os.makedirs(os.path.join(pattern_storage_directory, f"rdf/{size}/"))
        
        print("\tSaving candidate patterns...")
        for p_num, p_id in tqdm(enumerate(pattern_batch[size].keys())):
            p = pattern_batch[size][p_id]
            
            # Save as query
            _, query = p._as_sparql_query(return_string = True, assume_named_graphs = True)
            
            with open(os.path.join(pattern_storage_directory, f"queries/{size}/{p_id}.txt"), "w") as f:
                f.write(query)
                
            # Save as RDF
            base_uri = f"resource://integreat/p5/pattern/s{size}_n{p_num}/variable/"
            
            temp_graph = rdf.Graph()
            
            for s, p, o in p._pattern:
                temp_s = s
                if isinstance(s, rdf.Variable):
                    temp_s = rdf.URIRef(base_uri + s)
                
                temp_o = o
                if isinstance(o, rdf.Variable):
                    temp_o = rdf.URIRef(base_uri + o)
                    
                temp_graph.add((temp_s, p, temp_o))
            
            temp_graph.serialize(os.path.join(pattern_storage_directory, f"rdf/{size}", f"{p_id}.nt"), "nt")
        print("\tDone!")

def R_preprocessing_compute_dominations():
    """
    Invokes an R script to preprocess the patterns by computing the domination relationship.
    The domination relations found using this procedure are stored in the 
        {dataset_root_directory}/{dataset_tag}/{date}/dominations.yml 
    file, listing:
        - the sizes considered in the computation (for which dominating patterns are computed)
        - for each pattern, a list:
        - the first entry is the pattern size
        - the remaining entries are the dominating patterns
    """
    robjects.r.source(INPUT_FILES["R"]["domination"])

def check_pattern_content(pattern_batch):
    """
    Computes pattern interestingness based on three rules:
        1. A pattern specifies at least one ligand
        2. A pattern specifies at least one bingind atom
        3. A pattern specifies both of the above, but in a redundant way
        
    Arguments:
        - pattern_batch: a patterns.PatternBatch object containing the patterns of interest
        
    Returns:
        - pattern_info: a defaultdict(lambda: [[], [], []]) in which each pattern size is associated with three list,
            each containing the ids of the patterns that satisfy the corresponding rule
    """
    
    pattern_info = defaultdict(lambda: [[], [], []])
    
    specifies_ligand_query = patterns.Pattern([
            ("?x1", "<resource://integreat/p5/ligand/ligand/isLigand>", "?x2")
        ])
    
    def specifies_ligand(pattern):
        matches = specifies_ligand_query(pattern)
                
        return len( [ m for m in matches if not isinstance(m[rdf.Variable("x2")], rdf.Variable) ] ) > 0
    
    specifies_binding_atom_query = patterns.Pattern([
            ("?y1", "<resource://integreat/p5/ligand/bond/hasBindingAtom>", "?y2"),
            ("?y2", "<resource://integreat/p5/atomic/atom/isAtom>", "?y3")
        ])
    
    def specifies_binding_atom(pattern):
        matches = specifies_binding_atom_query(pattern)
                
        return len( [ m for m in matches if not isinstance(m[rdf.Variable("y3")], rdf.Variable) ] ) > 0
    
    specifies_both_redundant_query = specifies_ligand_query + \
                                            specifies_binding_atom_query + \
                                            patterns.Pattern([
                                                    ("?x1", "<resource://integreat/p5/ligand/structure/bLl>", "?y1")
                                                ])
    def specifies_both_redundant(pattern):
        matches = specifies_both_redundant_query(pattern)
                
        return len( [ 
                m for m in matches 
                if not isinstance(m[rdf.Variable("y3")], rdf.Variable) and 
                    not isinstance(m[rdf.Variable("x2")], rdf.Variable) 
            ] ) > 0
        
    for size, pttns in pattern_batch.items():
        for p_id, p in pttns.items():
            sl = specifies_ligand(p)
            sba = specifies_binding_atom(p)
            sbnr = specifies_both_redundant(p)
            
            if sl:
                pattern_info[size][0] += [p_id]
            if sba:
                pattern_info[size][1] += [p_id]
            if sbnr:
                pattern_info[size][2] += [p_id]
            
        print(f"Size {size}: " + 
              f"{len(pattern_info[size][0])/len(pattern_batch[size]):.2f} [{len(pattern_info[size][0])}] - " +
              f"{len(pattern_info[size][1])/len(pattern_batch[size]):.2f} [{len(pattern_info[size][1])}] - " + 
              f"{len(pattern_info[size][2])/len(pattern_batch[size]):.2f} [{len(pattern_info[size][2])}]")
        
    return pattern_info

def filter_working_directory(pattern_info, pattern_parents, pattern_storage_directory):
    """
    Filters the working directory created by extract_and_save_patterns using the information computed by
    check_pattern_content. The filtered working directory contains all the pattern that specify at least one binding
    atom and have no reduntant information plus all their parents (within the sizes of interest)
    
    Arguments:
        - pattern_info: the output of check_pattern_content
        - pattern_parents: the output of R_preprocessing_compute_dominations
        - pattern_storage_directory: the directory in which the patterns are be stored
    """
    
    sizes_to_process = [int(s) for s in os.listdir(os.path.join(pattern_storage_directory, "queries"))]
    sizes_to_process = sorted(sizes_to_process, reverse = True)
    
    previous_candidates = []
    for size in sizes_to_process:
        candidates = [p_id for p_id in pattern_info[size][1] if not p_id in pattern_info[size][2]]
        
        for q_id in previous_candidates:
            candidates += pattern_parents[q_id][1:]
                
        candidates = set(candidates)
        
        print(f"Filtering patterns of size {size}...")
        pattern_list = [x.split(".")[0] for x in os.listdir(os.path.join(pattern_storage_directory, "queries", f"{size}"))]
        for p_id in tqdm(pattern_list):
            
            if p_id not in candidates:
                os.remove(os.path.join(pattern_storage_directory, "queries", f"{size}", f"{p_id}.txt"))
                os.remove(os.path.join(pattern_storage_directory, "rdf", f"{size}", f"{p_id}.nt"))
        print("\tDone!")
                
        previous_candidates = candidates

def R_postprocessing_compute_matches():
    """
    Invokes an R script to postprocess the patterns by computing the matches onto the training partition of the 1k selection.
    The matches are stored in the file 
        {dataset_root_directory}/{dataset_tag}/{date}/results/matches.csv
    which is an N_graphs x N_patterns dataset containing for each graph and pattern, the number of matches.
    Notice that the script creates the file incrementally with respect to the pattern sizes, i.e., it starts with the columns related to the
    smaller sizes, then moves to the next and so on. Hence, if the sizes of interests are known and the number of patterns per size is also known
    it is possible to reconstruct the size of each pattern associated with each column.
    """
    robjects.r.source(INPUT_FILES["R"]["matches"])

#%% __main__
if __name__ == "__main__":
    local_dataset_directory, pattern_storage_directory, cluster_storage_directory, results_directory, pattern_batch = initialise()
    
    # %%% Extract and save patterns
    if DO_EXTRACT_AND_SAVE_PATTERNS:
        extract_and_save_patterns(pattern_batch, pattern_storage_directory)
    
    # %%% [R] Preprocessing: compute domination relationship
    if DO_R_PREPROCESSING_DOMINATION:
        if SIZE_RANGE[0] < SIZE_RANGE[1]: 
            R_preprocessing_compute_dominations()
            
            # Purge VOS db directory (for safety and memory reasons)
            if os.path.exists(INPUT_FILES["vos_db_dir"]):
                shutil.rmtree(INPUT_FILES["vos_db_dir"])
            
            os.makedirs(INPUT_FILES["vos_db_dir"])
            
            with open(os.path.join(local_dataset_directory, "dominations.yml"), "r") as f:
                pattern_parents = defaultdict(list, yaml.load(f, yaml.CLoader))
        else:
            print("There is only one size of interest. Dominations are irrelevant. Skipping computation.")
            pattern_parents = defaultdict(list)
    else:
        if SIZE_RANGE[0] < SIZE_RANGE[1]:
            with open(os.path.join(local_dataset_directory, "dominations.yml"), "r") as f:
                pattern_parents = defaultdict(list, yaml.load(f, yaml.CLoader))
        else:
            print("There is only one size of interest. Dominations are irrelevant. Skipping computation.")
            pattern_parents = defaultdict(list)
        
    # %%% Check pattern content
    if DO_CHECK_PATTERN_CONTENT:
        pattern_info = check_pattern_content(pattern_batch)
    
    # %%% Remove unnecessary patterns from working directory
    if DO_FILTER_WORKING_DIRECTORY:
        filter_working_directory(pattern_info, pattern_parents, pattern_storage_directory)
    
    # %%% [R] Postprocessing: compute matches
    if DO_R_POSTPROCESSING_MATCHES:
        R_postprocessing_compute_matches()
        
        # Purge VOS db directory (for safety and memory reasons)
        if os.path.exists(INPUT_FILES["vos_db_dir"]):
            shutil.rmtree(INPUT_FILES["vos_db_dir"])
        
        os.makedirs(INPUT_FILES["vos_db_dir"])