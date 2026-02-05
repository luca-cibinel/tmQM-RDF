"""
Compute pattern similarity and cluster similar patterns
"""
# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
from collections import defaultdict
from datetime import date
from tqdm import tqdm

import yaml
import shutil
import pickle
import numpy as np
import pandas as pd
import seaborn as sbn
import matplotlib.pyplot as plt

import warnings
np.warnings = warnings

INPUT_FILES = {
        "intermediate": os.path.join(ROOT_DIR, "computational", "pattern_clustering", "intermediate", "%s"),
        "patterns": os.path.join(ROOT_DIR, "computational", "pattern_mining", "results", "%s"),
        "tmQM-RDF-1Ksel": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "data", "%s", "1k_selection.csv")
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
# Agglomerative clustering settings
AGGLOMERATION_SETTINGS = {
    "dist_matrix_file": "dist_matrix.pkl",
    "linkage": "average",
    "aggregation": ["max", "set_median"],
    "min_n_clusters": 300
}

# Operational settings
T = True
F = False

do_clustering = T                    # Cluster patterns

# %% Utility functions
def identify_working_paths(dataset_tag, dist_matrix_file_name):
    """
    This function identifies all the subdirectories, starting from INPUT_FILES['intermediate'], that contain a file named as
    dist_matrix_file_name. Each of these subdirectories is marked as a potential working subdirectory.
    The user is then asked to confirm which directories to use.
    
    Arguments:
        - dataset_tag: the dataset_tag  needed to identify the correct intermediate directory
        - dist_matrix_file_name: the file name used to identify distance matrices and working directories
        
    Returns:
        - The working directories approved by the user
    """
    
    root = INPUT_FILES["intermediate"] % dataset_tag
    candidates = []
    
    for directory, _, files in os.walk(root):
        if dist_matrix_file_name in files:
            candidates += [directory.replace(root, "")[1:]]
            
    prompt = f"Candidate working paths:\n{'\n'.join([f'{i}) {d}' for i, d in enumerate(candidates)])}"
    prompt += "\n\n"
    prompt += " > Type the numbers of the directories to use: "
    
    selection = input(prompt)
    
    selected = []
    for c in selection:
        selected += [candidates[int(c)]]
        
    return selected

# %% Main functions
def prepare_working_tree(root, dirs):
    """
    This function prepares the working tree by performing the following operations:
        - Purges the content of all the directories listed in dirs
        - Eliminates all the temp files (filename starts with .) found in within root and its subdirectories
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
    
    if len(temp_to_remove) + len(dirs) == 0:
        return
    
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
    cluster_storage_directory = os.path.join(local_dataset_directory, "clusters_figures")
    
    # Directory where the final results should be stored
    results_directory = os.path.join(local_dataset_directory, "by_metric")
    
    # Preliminary operations
    
    # Prepare working tree
    dirs_to_prepare = []
    
    prepare_working_tree(local_dataset_directory, dirs_to_prepare)
    
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
            "dirs_to_prepare": dirs_to_prepare,
            "agglomeration_settings": AGGLOMERATION_SETTINGS
        }
    
    # Determine appropraite config file name
    config_fname = "config_clustering"
    
    with open(os.path.join(local_dataset_directory, config_fname), "w") as f:
        yaml.dump(config, f)

    return local_dataset_directory, pattern_storage_directory, cluster_storage_directory, results_directory, dataset_tag

def cluster_patterns(settings, pattern_parents, dataset_tag, results_directory,
                                  pattern_storage_directory, cluster_storage_directory, tmQM_RDF_selection):
    """
    Clusters patterns based on a similarity metric using the agglomerative clustering algorithm.
    It produces the following files, inside {dataset_root_directory}/{dataset_tag}/{date}/results/{metric}:
        - {aggregation}/cls_matches.csv: an N_graphs x N_clusters dataset containing the aggregated matches of the clusters
        - set_median/cls_set_median_graphs.csv: an N_clusters x 2 dataset containing, for each cluster, the index of the computed
            set median graph (only if aggregation is set_median).
        - cls_info.csv: a dataset containing metadata about each pattern (index in original dataset, id, size, assigned cluster)
        - cls_edges.csv: a two-columns dataset containing the extremes of the "a priori" edges (an "a priori" edge is an edge
            between two clusters such that there exists two patterns, one in each cluster, that are in a domination relationship; the
            direction of the relationship is disregarded).
        - silhouette.csv: a dataframe containing, for each candidate distance_threshold parameter, the computed silhouette score.
        - dist_matrix.pkl: a pickle file containig the distance matrix (1 - similarity).
        - cluster_sizes.png/.pdf: a barplot of cluster size vs cluster id.
        - silhouette.png/.pdf: a plot of the silhouette score profile.
    
    Arguments:
        - settings: the clustring settings
        - pattern_parents: the output of R_preprocessing_compute_dominations (as a defaultdict)
        - dataset_tag: the identifying tag of the dataset
        - results_directory: the path {dataset_root_directory}/{dataset_tag}/{date}/
        - pattern_storage_directory: the directory in which the patterns are stored
        - cluster_storage_directory: the directory where cluster data (images and .nt RDF graphs) are to be stored
        - tmQM_RDF_selection: the path to the .csv file specifying the desired dataset selection and the train/validation/test split
    """
    # Prepare to store used distances (for diagnostic purposes)
    sims = dict()
    
    # Utility step: create a dictionary to retrieve the pattern size based on the position in the DF
    # It exploits the fact that the R script creates the dataset by progresively joining columns associated with increasing pattern size
    sizes = sorted([int(x) for x in os.listdir(os.path.join(pattern_storage_directory, "queries"))])
    lengths = [len(os.listdir(os.path.join(pattern_storage_directory, "queries", str(size)))) for size in sizes]
    
    last_idx = 0
    size_map = dict()
    
    for s, l in zip(sizes, lengths):
        size_map[s] = (last_idx, last_idx + l - 1)
        last_idx += l
    
    idx_to_size = lambda idx: [s for s in size_map if idx >= size_map[s][0] and idx <= size_map[s][1]][0]
            
    # Retrieve working directories
    working_paths = identify_working_paths(dataset_tag, settings["dist_matrix_file"])
    
    # Cluster patterns
    for wp in working_paths:
        print(f"Operating with wp {wp}")
                
        local_results_directory = os.path.join(results_directory, wp)
        
        if not os.path.exists(local_results_directory):
            os.makedirs(local_results_directory)
        
        local_cluster_storage_directory = os.path.join(cluster_storage_directory, wp)
        
        if os.path.exists(local_cluster_storage_directory):
            shutil.rmtree(local_cluster_storage_directory)
            
        os.makedirs(local_cluster_storage_directory)
        
        # Read matches data frame
        df = pd.read_csv(os.path.join(OUTPUT_FILES["results"] % dataset_tag, "matches.csv"), index_col = 0)
        
        # Retrieve distance matrix
        with open(os.path.join(INPUT_FILES["intermediate"] % dataset_tag, wp, settings["dist_matrix_file"]), "rb") as f:
            dist = pickle.load(f)
        
        sims[wp] = 1 - dist[np.triu_indices(dist.shape[0], 1)]
        
        # Identify threshold
        print("\tIdentifying threshold via silhouette coefficient...")
        
        ag_clusters = lambda t: AgglomerativeClustering(
                                    metric = "precomputed", 
                                    linkage = settings["linkage"], 
                                    distance_threshold = t, 
                                    n_clusters = None
                                ).fit_predict(dist)
        
        thr_space = np.linspace(0, 1, 201)
        scores = np.zeros(len(thr_space))
        
        for i, t in enumerate(thr_space):
            cl_loc = ag_clusters(t)
            if len(set(cl_loc)) in [1, dist.shape[0]] or len(set(cl_loc)) < settings["min_n_clusters"]: # score is not defined if n of clusters is 1 or n_samples
                scores[i] = -1
            else:
                scores[i] = silhouette_score(dist, cl_loc, metric = "precomputed")
        
        opt_thr = thr_space[np.argmax(scores)]
        
        fig, ax = plt.subplots()
        ax.plot(thr_space, scores, color = "orange")
        ax.axvline(opt_thr, ls = "--", color = "black")
        ax.set_xlabel("Distance threshold")
        ax.set_ylabel("Silhouette score")
        
        fig.savefig(os.path.join(local_results_directory, "silhouette.png"), format = "png")
        fig.savefig(os.path.join(local_results_directory, "silhouette.pdf"), format = "pdf")
        
        scores_df = pd.DataFrame({"threshold": thr_space, "silhouette": scores})
        scores_df.to_csv(os.path.join(local_results_directory, "silhouette.csv"))
        
        # Cluster patterns
        print(f"\tClustering patterns (optimal threshold = {opt_thr})...")
        
        computed_clusters = ag_clusters(opt_thr)
        
        clusters = [[j for j in range(len(computed_clusters)) if computed_clusters[j] == i] for i in range(max(computed_clusters) + 1)]
    
        # Plot cluster sizes
        cluster_sizes = [len(c) for c in clusters]
        
        fig, ax = plt.subplots()
        ax.bar(range(len(cluster_sizes)), cluster_sizes)
        ax.set_ylabel("Cluster size")
        ax.set_xlabel("Cluster index")
        
        fig.savefig(os.path.join(local_results_directory, "cluster_sizes.png"), format = "png")
        fig.savefig(os.path.join(local_results_directory, "cluster_sizes.pdf"), format = "pdf")
        
        print("\tSaving cluster information...")
        cluster_info = {"index_in_df": [], "id": [], "size": [], "cluster": []}
        
        for i, cluster in enumerate(clusters):
            for idx in cluster:
                cluster_info["index_in_df"] += [idx]
                cluster_info["id"] += [df.columns[idx]]
                cluster_info["size"] += [idx_to_size(idx)]
                cluster_info["cluster"] += [i]
        
        df_cluster_info = pd.DataFrame(cluster_info)
        df_cluster_info.to_csv(os.path.join(local_results_directory, "cls_info.csv"))
    
        print("\tCopying pattern images...")
        id_to_fname = lambda p_id: "pattern-" + p_id.replace("(", "").replace(")", "") + ".png"
        
        for i, cluster in tqdm(enumerate(clusters)):
            for size in sizes:
                os.makedirs(os.path.join(local_cluster_storage_directory, f"cluster-{i}", str(size)), exist_ok = True)
            
            for p_num in cluster:
                fname = id_to_fname(df.columns[p_num])
                p_size = idx_to_size(p_num)
                shutil.copyfile(
                        os.path.join(INPUT_FILES["patterns"] % PATTERN_DATASET_KEYWORD, f"patterns-{p_size}", fname), 
                        os.path.join(local_cluster_storage_directory, f"cluster-{i}", str(p_size), fname)
                    )
            
        print("\tComputing preliminary adjacency matrix...")
        edges = [] # this list contains tuples (i,j) such that there exists a pattern in cluster_i that is a parent of a pattern in cluster_j
        
        for i, pa_cluster in enumerate(clusters):
            pa_cluster_ids = [df.columns[p_num] for p_num in pa_cluster]
            
            for j, ch_cluster in enumerate(clusters):
                if i == j:
                    continue
                
                if (i,j) in edges or (j,i) in edges:
                    continue
                
                ch_cluster_ids = [df.columns[p_num] for p_num in ch_cluster]
                
                for p_id in ch_cluster_ids:
                    if len([q_id for q_id in pattern_parents[p_id] if q_id in pa_cluster_ids]) > 0:
                        edges += [(i,j)]
                        break
                
        df_parents = pd.DataFrame({"v1": [e[0] for e in edges], "v2": [e[1] for e in edges]})
        df_parents.to_csv(os.path.join(local_results_directory, "cls_edges.csv"))
        
        for aggregation_mode in settings["aggregation"]:
            print(f"\tAggregating cluster matches (mode = {aggregation_mode})")
            
            local_results_subdirectory = os.path.join(local_results_directory, f"{aggregation_mode}")
            if not os.path.exists(local_results_subdirectory):
                os.makedirs(local_results_subdirectory)
            
            # Compute clusters matches (disjunction of individual matches)
            X_full = np.sign(np.array(df))
             
            if aggregation_mode == "max":
                clustered_match_data =  {i: np.apply_along_axis(np.max, 1, X_full[:, cluster]) for i, cluster in enumerate(clusters)}
            if aggregation_mode == "set_median":
                print("\t\tComputing and storing median graphs")
                median_graphs = {
                    i: cluster[ np.argmin(dist[np.ix_(cluster, cluster)].sum(1)) ]
                    for i, cluster in enumerate(clusters)
                }
                clustered_match_data =  {i: X_full[:, mg] for i, mg in median_graphs.items()}
                
                # Store median graphs indices
                df_median_graphs = pd.DataFrame({
                    "cluster": list(median_graphs.keys()), 
                    "set_median_graph_index": list(median_graphs.values()),
                    "set_median_graph_pattern_id": [df_cluster_info.loc[df_cluster_info["index_in_df"] == idx, "id"].iloc[0] for idx in median_graphs.values()]
                })
                df_median_graphs.to_csv(os.path.join(local_results_subdirectory, "cls_set_median_graphs.csv"))
                
                # Store median graphs images
                local_cluster_storage_subdirectory = os.path.join(local_cluster_storage_directory, "set_median_graphs")
                if not os.path.exists(local_cluster_storage_subdirectory):
                    os.makedirs(local_cluster_storage_subdirectory)
                    
                for i, _ in tqdm(enumerate(clusters)):
                    
                    fname = id_to_fname(df_median_graphs.loc[df_median_graphs["cluster"] == i, "set_median_graph_pattern_id"].iloc[0])
                    p_size = idx_to_size(df_median_graphs.loc[df_median_graphs["cluster"] == i, "set_median_graph_index"].iloc[0])
                    
                    shutil.copyfile(
                            os.path.join(INPUT_FILES["patterns"] % PATTERN_DATASET_KEYWORD, f"patterns-{p_size}", fname), 
                            os.path.join(local_cluster_storage_subdirectory, f"smg_cluster-{i}.png")
                        )
            
            df_clustered = pd.DataFrame(clustered_match_data, index = df.index)
            df_clustered.to_csv(os.path.join(local_results_subdirectory, "cls_matches.csv"))

    # Diagnostic: pair plot of used distances
    p = sbn.pairplot(
        pd.DataFrame(sims), 
        corner = True,
        diag_kind = "auto", 
        diag_kws={"color": "orange"}, 
        plot_kws = {"color": "orange", "alpha": 0.2},
        grid_kws = {"diag_sharey": False}
    )
    p.set(xlim = (0, 1), ylim = (0, 1))
    plt.subplots_adjust(wspace=0.2, hspace=0.2)
    plt.show()
    
#%% __main__
if __name__ == "__main__":
    local_dataset_directory, pattern_storage_directory, cluster_storage_directory, results_directory, dataset_tag = initialise()
    
    # %%% Cluster patterns based on similarity
    if do_clustering:
        with open(os.path.join(local_dataset_directory, "dominations.yml"), "r") as f:
            pattern_parents = defaultdict(list, yaml.load(f, yaml.CLoader))
            
        cluster_patterns(
            AGGLOMERATION_SETTINGS,
            pattern_parents, 
            dataset_tag,
            results_directory, 
            pattern_storage_directory, 
            cluster_storage_directory, 
            INPUT_FILES["tmQM-RDF-1Ksel"]
        )
        
        