"""
Computes all the pairwise (dis)similarity scores within a given pattern dataset.

NOTE: this code is designed to be run on an external computational server
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header

# Commented imports are imported below conditionally on chosen parameters
import pickle
import numpy as np
# import rdflib as rdf #
import delgsim as sm
import networkx as nx 
# import matplotlib.pyplot as plt #
# import rdflib.extras.external_graph_libs as rdfextra #

from collections import defaultdict
from tqdm import tqdm

INPUT_FILES = {
        "processed_patterns": os.path.join(ROOT_DIR, "computational", "pattern_clustering", "results", "%s")
    }

OUTPUT_FILES = {
        "intermediate": os.path.join(ROOT_DIR, "computational", "pattern_clustering", "intermediate", "%s"),
        "pattern_info": os.path.join("%s", "sizes_and_compressed_patterns.pkl"),
        "weights": os.path.join("%s", "weights.pkl"),
        "dist": os.path.join("%s", "dist_matrix.pkl")
    }

# tmQM_RDF_selection = "./../../../../../data/extract_selection_from_tmQM-RDF/1k_selection.csv"

# processed_data_directory = os.path.join("..", "..", "..", "processed_pattern_data")

PATTERN_DATASET_KEYWORD = "lateTM"
DATASET_SHORT_ID = "latmod" # Ligand ATtachment MODes

# Range of sizes of interes (extremes included)
SIZE_RANGE = [10, 12]

T = True
F = False

DO_REPRESENTATION_COMPUTATION = T
DO_WEIGHT_COMPUTATION = T
USE_LEARNED_WEIGHTS = F
DO_DELGSIM_COMPUTATION = F
DRAW_PLOTS = F

if DO_REPRESENTATION_COMPUTATION:
    import rdflib as rdf
    import rdflib.extras.external_graph_libs as rdfextra
    
if DRAW_PLOTS:
    import matplotlib.pyplot as plt
    
# %% Utility functions

def pattern_to_nx(patt):
    """
    Computes the 'compressed' representation of a pattern as an nx.DiGraph object
    
    Arguments:
        - patt: a rdflib.Graph object representing the pattern
        
    Returns:
        - an nx.DiGraph representing the 'compressed' pattern. The labels are stored under the node/edge attribute 'label'
    """
    g = rdfextra.rdflib_to_networkx_digraph(
            patt,
            False,
            lambda s, p, o: {"label": str(p)} # constructs the edge attribute 'label' by converting the predicate into the label
        )
    
    node_labels = {}
    edge_remove = []
    
    mc_classified = False
    classified_ligands = []
    classified_atoms = []
    
    for e in g.edges(data = "label"):
        predicate = e[2].split("/")[-1]
        
        if predicate == "hasMetalCentre":
            # Domain: TransitionMetalComplex
            # Range: MetalCentre (unless a more specific class is available)
            
            node_labels[e[0]] = "resource://integreat/p5/complex/TMC/TransitionMetalComplex"
            
            if not mc_classified:
                node_labels[e[1]] = "resource://integreat/p5/ligand/centre/MetalCentre"
            
        elif predicate == "hasLigand":
            # Domain: TransitionMetalComplex
            # Range: Ligand (unless a more specific class is available)
            
            node_labels[e[0]] = "resource://integreat/p5/complex/TMC/TransitionMetalComplex"
            
            if e[1] not in classified_ligands:
                node_labels[e[1]] = "resource://integreat/p5/ligand/ligand/Ligand"
                
        elif predicate == "isMetalCentre":
            # Domain: MetalCentre
            # Range: MetalCentreClass
            
            if not str(e[1]).split("/")[-2] == "variable":
                mc_classified = True
                
                node_labels[e[0]] = str(e[1])
        
            edge_remove += [e] 
        
        elif predicate == "bLc":
            # Domain: MetalCentre
            # Range: LigandBond
            
            if not mc_classified:
                node_labels[e[1]] = "resource://integreat/p5/ligand/centre/MetalCentre"
            
            node_labels[e[1]] = "resource://integreat/p5/ligand/bond/LigandBond"
                
        elif predicate == "bLl":
            # Domain: Ligand (unless a more specific class is available)
            # Range: LigandBond
            
            if e[0] not in classified_ligands:
                node_labels[e[0]] = "resource://integreat/p5/ligand/ligand/Ligand"
            
            node_labels[e[1]] = "resource://integreat/p5/ligand/bond/LigandBond"
            
        elif predicate == "isLigand":
            # Domain: Ligand
            # Range: LigandClass
            
            classified_ligands += [e[0]]
            
            node_labels[e[0]] = str(e[1])
            
            edge_remove += [e]
            
        elif predicate == "hasBindingAtom":
            # Domain: LigandBond
            # Range: Atom (unless a more specific class is available)
            
            node_labels[e[0]] = "resource://integreat/p5/ligand/bond/LigandBond"
            
            if e[1] not in classified_atoms:
                node_labels[e[1]] = "resource://integreat/p5/atomic/atom/Atom"
                
        elif predicate == "isAtom":
            # Domain: Atom
            # Range: Element (unless a more specific class is available)
            
            classified_atoms += [e[0]]
            
            node_labels[e[0]] = str(e[1])
            
            edge_remove += [e]
    
    # Add node labels
    g.update(nodes = [(v, {"label": l}) for v, l in node_labels.items()])
    
    # Remove edges of the form is...
    g.remove_edges_from(edge_remove)
    
    # Remove isolated nodes created by edge removal
    g.remove_nodes_from(list(nx.isolates(g)))
    
    # Detach networkx nodes and edges from rdflib
    g_detach = nx.DiGraph()
    
    g_detach.add_nodes_from([(str(v), d) for v, d in g.nodes(data = True)])
    g_detach.add_edges_from([(str(u), str(v), d) for u, v, d in g.edges(data = True)])
    
    return g_detach
          
def wF(x):
    # Remember that
    #  x = (g, node_id, label) if x is a node
    #  x = (g, source_id, target_id, label) if x is an edge
    
    if len(x) == 3:
        if x[2] in [
                    "resource://integreat/p5/atomic/atom/Atom",
                    "resource://integreat/p5/ligand/ligand/Ligand",
                    "resource://integreat/p5/ligand/centre/MetalCentre"
                ]:
            return 0.5
    
    return 1

# %% Main functions

def pattern_dataset_to_nx(patterns_directory, pat_sizes, pat_ids):
    """
    This function converts a pattern dataset into the equivalent nx.DiGraph dataset
    
    Arguments:
        - patterns_directory: the path to directory where patterns are stored (by size) as RDF graphs
        - pat_sizes: a dictionary that associates to each size of interest the number of patterns of that size
        - pat_ids: a list of pattern ids arranged by size, i.e. [...pat of size i - 1...pat of size i...pat of size i + 1...]
        
    Returns:
        - a dictionary that associates to each size of interest a list of patterns compressed into their nx.DiGraph representation
    """
    sizes_of_interest = pat_sizes.keys()
    
    # Networkx graph computation
    buffer1 = 0
    pat_nx = {i: [] for i in sizes_of_interest}
    for i in sizes_of_interest:
        for p_id in pat_ids[buffer1:(buffer1 + pat_sizes[i])]:
            temp1 = rdf.Graph()
            temp1.parse(os.path.abspath(os.path.join(patterns_directory, f"{i}", f"{p_id}.nt")))
            
            pat_nx[i] += [pattern_to_nx(temp1)]
            
        buffer1 += pat_sizes[i]
        
    return pat_nx

def compute_label_weights(pat_nx, labels_hierarchy = dict()):
    """
    Compute /smooth) inverse document frequency weights for each node/edge label:
    
        w(l) = 1 + ln(N_p/(1 + N_p,l))
        
        N_p = total number of patterns
        N_p,l = number of patterns in which l is present
        
    The weights are then renormalised to lie in the range [0, 1]
    """
    count = defaultdict(int)
    edge_labels = []
    
    for size, pat_batch in pat_nx.items():
        for pat in pat_batch:
            labels = []
            
            for v in pat.nodes(data = "label"):
                labels += [v[1]]
                # count[v[1]] += 1
                
            for e in pat.edges(data = "label"):
                labels += [e[2]]
                edge_labels += [e[2]]
                # count[e[2]] += 1
            
            for l in set(labels):
                count[l] += 1
    
    edge_labels = set(edge_labels)
    
    #weights = {l: 1/c for l, c in count.items()}
    n_pat = sum([len(pat_batch) for pat_batch in pat_nx.values()])
    weights = {l: 1 + np.log(n_pat/(1 + c)) for l, c in count.items()}
    
    for root, is_child in labels_hierarchy.items():
        for l in weights.keys():
            if is_child(l) and l != root and l not in edge_labels:
                weights[l] += weights[root]
        
    Z = max(weights.values())
    weights = {l: w/Z for l, w in weights.items()}
    
    return weights

def compute_similarity_matrix(pat_sizes, pat_nx, dataset_tag, wF):
    """
    Compute the similarity matrix
    
    Arguments:
        - pat_sizes: a dictionary that associates to each size of interest the number of patterns of that size
        - pat_nx: a dictionary that associates to each size of interest a list of patterns compressed into their nx.DiGraph representation
        - dataset_tag: a descriptive tag (used to save the resulting matrix)
        - wF: a weighting function
    """
    sims = np.eye(sum(pat_sizes.values()))
    max_size = max(pat_nx.keys())
    
    buffer1 = 0
    for size1 in pat_nx.keys():
        
        for i, pat_i in tqdm(enumerate(pat_nx[size1]), total = len(pat_nx[size1])):
            
            buffer2 = buffer1
            for size2 in range(size1, max_size + 1):
                    
                for j, pat_j in enumerate(pat_nx[size2]):
                    
                    if size1 == size2 and i >= j:
                        continue
                    
                    sim = sm.similarity(pat_i, pat_j, wF = wF)[0]
                    sims[buffer1 + i, buffer2 + j] = sim
                    sims[buffer2 + j, buffer1 + i] = sim
                
                buffer2 += pat_sizes[size2]
            
        buffer1 += pat_sizes[size1]
        
    return np.clip(sims, 0, 1) # Take care of numerical instabilities

# %% __main__
if __name__ == "__main__":
    
    # Compute dataset tag (for dataset identification)
    if SIZE_RANGE[1] > SIZE_RANGE[0]:
        size_tag = f"s_{SIZE_RANGE[0]}_{SIZE_RANGE[1]}"
    else:
        size_tag = f"s_{SIZE_RANGE[0]}"
    
    dataset_tag = f"{PATTERN_DATASET_KEYWORD}-{DATASET_SHORT_ID}-{size_tag}"
    
    # Resolve working directories
    target_directory = os.path.join(INPUT_FILES["processed_patterns"] % dataset_tag)
    
    results_directory = os.path.join(target_directory)
    
    local_output_directory = os.path.join(OUTPUT_FILES["intermediate"] % dataset_tag, "delgsim", "estimated_weights" if USE_LEARNED_WEIGHTS else "fixed_weights")
    
    patterns_directory = os.path.join(target_directory, "patterns", "rdf")
    
    if not os.path.exists(local_output_directory):
        os.makedirs(local_output_directory)
    
    # Parse dataset as nx.DiGraph
    if DO_REPRESENTATION_COMPUTATION:
        # Pattern ids extraction (read column names from matches.csv, notice that first entry is the index column)
        # This guarantees that patterns appear in the same order both in the similarity matrix and in the match-feature matrix
        with open(os.path.join(results_directory, "matches.csv"), "r") as f:
            pat_ids = eval(f.read().split("\n")[0])[1:]
          
        # Pattern sizes extraction
        sizes_of_interest = range(SIZE_RANGE[0], SIZE_RANGE[1] + 1)
        pat_sizes = {i: len(os.listdir(os.path.join(patterns_directory, f"{i}"))) for i in sizes_of_interest}
        
        pat_nx = pattern_dataset_to_nx(patterns_directory, pat_sizes, pat_ids)
        
        with open(OUTPUT_FILES["pattern_info"] % local_output_directory, "wb") as f:
            pickle.dump((sizes_of_interest, pat_sizes, pat_nx), f)
    else:
        with open(OUTPUT_FILES["pattern_info"] % local_output_directory, "rb") as f:
            sizes_of_interest, pat_sizes, pat_nx = pickle.load(f)
    
    # Compute weights
    blank_labels = [
                "resource://integreat/p5/atomic/atom/Atom",
                "resource://integreat/p5/ligand/ligand/Ligand",
                "resource://integreat/p5/ligand/centre/MetalCentre"
            ]
    
    def w_factory():
        return 1
        
    if DO_WEIGHT_COMPUTATION:
        if USE_LEARNED_WEIGHTS:
            # Build hierarchy by comparing URIs:
            #   children of the blank label resource://integreat/p5/*/X/y
            #   will always have the form resource://integreat/p5/*/X/reference/z
            # The default argument bl = bl in lambda forces the expression to evaluate
            # bl at its current value. Otherwise, it would only reference the variable bl
            # whose value, at the end, will be stuck at blank_labels[-1]
            labels_hierarchy = {
                    bl: (lambda cl, bl = bl : (cl.split("/")[-3] == bl.split("/")[-2]))
                    for bl in blank_labels
                }
            
            weights = compute_label_weights(pat_nx, labels_hierarchy)
            
            if DRAW_PLOTS:
                fig, ax = plt.subplots(1, 1)
                
                p = np.argsort(np.array([l.replace("/", "") for l in weights.keys()]))
                x = np.array([l.split("/")[-1] for l in weights.keys()])[p]
                y = np.array(list(weights.values()))[p]
                
                ax.bar(x, y)
                ax.bar(["Atom", "Ligand", "MetalCentre"], [weights[l] for l in [
                            "resource://integreat/p5/atomic/atom/Atom",
                            "resource://integreat/p5/ligand/ligand/Ligand",
                            "resource://integreat/p5/ligand/centre/MetalCentre"
                        ]], color = "orange")
                plt.setp(ax.xaxis.get_majorticklabels(), rotation = 70, rotation_mode = "anchor", ha = "right")
                plt.show()
        else:
            weights = defaultdict(w_factory, {l: 0.5 for l in blank_labels})
        
        with open(OUTPUT_FILES["weights"] % local_output_directory, "wb") as f:
            pickle.dump(weights, f)
    else:
        with open(OUTPUT_FILES["weights"] % local_output_directory, "rb") as f:
            weights = pickle.load(f)
    
    # Compute similarity matrix
    if DO_DELGSIM_COMPUTATION:
        with open(OUTPUT_FILES["pattern_info"] % local_output_directory, "rb") as f:
            sizes_of_interest, pat_sizes, pat_nx = pickle.load(f)
        
        sims = compute_similarity_matrix(pat_sizes, pat_nx, dataset_tag, lambda x: weights[x[-1]])

        with open(OUTPUT_FILES["dist"] % local_output_directory, "wb") as f:
            pickle.dump(1 - sims, f)
        
        if DRAW_PLOTS:
            fig, ax = plt.subplots()
            ax.hist(sims[np.triu_indices(sims.shape[0], 1)], bins = 100)
            ax.set_xlim(0,1)
            ax.set_title(f"s_DELG {'(learned w.)' if DO_WEIGHT_COMPUTATION else '(fixed w.)'}")
        
        # fig, ax = plt.subplots()
        # ax.plot(np.array(cossim)[np.triu_indices(cossim.shape[0], 1)], sims[np.triu_indices(sims.shape[0], 1)], "*")
        # ax.set_title("s_COS vs s_DELG (learned w.)")
        # ax.set_ylim(0,1)
        # ax.set_xlabel("s_COS")
        # ax.set_ylabel("s_DELG")
        
        
