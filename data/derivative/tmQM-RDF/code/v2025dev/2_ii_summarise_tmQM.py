"""
Step 2.ii in the data preprocessing pipeline

The aim is to summarise all of the information in tmQM into .gml files that mirror the structure of tmQMg. The link between the two datasets
is given by the atom IDs.
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header    
import shutil
import pandas as pd
import networkx as nx
from tqdm import tqdm
from tempfile import TemporaryDirectory

INPUT_FILES = {
        "localised_tmcs": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQMg", "localised_tmcs.csv"),
        "tmQM_X_xyz": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQM", "v2024", "tmQM_X%d.xyz"),
        "tmQM_X_bo": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQM", "v2024", "tmQM_X%d.BO"),
        "tmQM_X_q": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQM", "v2024", "tmQM_X.q"),
        "tmQM_y": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQM", "v2024", "tmQM_y.csv"),
        "tmQMg_graphs": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg", "v74.637k", "uNatQ_graphs"),
        "viable": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "step_2_i_viable_tmcs.txt")
	}

OUTPUT_FILES = {
        "tmQM_graphs": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQM", "uNatQ_graphs"),
		"viable": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "viable_tmcs.txt")
    }

# %% Utility functions
def map_tmQM():
    """
    Creates an atlas (pytho dictionary) containing all the information extracted from the files
        - tmQM/tmQM_X*.xyz
        - tmQM/tmQM_X.q
        - tmQM/tmQM_X*.BO
        - tmQM/tmQM_y.csv
    and also the information related to the localisation of TMCs within the .xyz and .BO files (see script locate_in_tmQM.py).

    The structure of the atlas is:
        - localisation: a pandas Data Frame with columns id, xyz, bo representing 
            respectively the CSD code and the number of the .xyz/.BO file
        - xyz_raw: a list of lists, where each of the inner list contains all of the entries of one of the
            three .xyz files. The entries are achieved by splitting by "\n\n"
        - q_raw: a list of the entries of the .q file. The entries are achieved by splitting by "\n\n"
        - bo_raw: a list of lists, where each of the inner list contains all of the entries of one of the
            three .BO files. The entries are achieved by splitting by "\n\n"
        - csv_raw: a pandas Data Frame obtained by reading the file tmQM/tmQM_y.csv (sep = ;)
        
    Returns:
        - a python dictionary representing the atlas described above
    """
    
    atlas = {}
        
    """
    Extract localisation
    """
    atlas["localisation"] = pd.read_csv(INPUT_FILES["localised_tmcs"], index_col=0)

    """
    Extract raw .xyz information
    """
    atlas["xyz_raw"] = []
    for i in range(3):
        with open(INPUT_FILES["tmQM_X_xyz"] % (i + 1), "r") as f:
            atlas["xyz_raw"] += [f.read().split("\n\n")]
        
    """
    Extract raw .q information
    """
    with open(INPUT_FILES["tmQM_X_q"], "r") as f:
        atlas["q_raw"] = f.read().split("\n\n")    
        
    """
    Extract raw .bo information
    """
    atlas["bo_raw"] = []
    for i in range(3):
        with open(INPUT_FILES["tmQM_X_bo"] % (i + 1), "r") as f:
            atlas["bo_raw"] += [f.read().split("\n\n")]
          
    """
    Extract raw .csv information
    """
    atlas["csv_raw"] = pd.read_csv(INPUT_FILES["tmQM_y"], sep = ";")

    return atlas

def listify_xyz(str_raw):
    """
    Converts a text of the form

    "***
    ***
    A       x  y  z
    B       x  y  z
    ..."
    (there are spaces between A and x and between x and y and between y and z)
    with x y z floats
    (format used by tmQM/tmQM_X[n].xyz after splitting by "\n\n")

    into a list of the form
    [('A', [x, y, z]), ('B', [x, y, z]), ...]
    
    Arguments:
        - str_raw: the text to process
        
    Returns:
        - the processed list
    """
    str_list = str_raw.split("\n")[2:]
    
    out = []
    for s in str_list:
        line = [c for c in s.split(" ") if not c == ""]
        out += [(line[0], [float(x) for x in line[1:]])]
        
    return out

def listify_q(str_raw):
    """
    Converts a text of the form

    "***
    A       q
    B       q
    Total charge = y
    ..."
    (format used by tmQM/tmQM_X.q after splitting by "\n\n")

    into a list of the form
    [('A', q), ('B', q), ...]
    and the scalar
    y
    
    Arguments:
        - str_raw: the text to process
        
    Returns:
        - the processed list
        - the scalar y
    """
    str_list = str_raw.split("\n")
    
    total_charge = str_list[-1]
    total_charge = eval(total_charge.replace("Total charge = ", ""))
    
    str_list = str_list[1:-1]
    
    out = []
    for s in str_list:
        line = [c for c in s.split(" ") if not c == ""]
        out += [(line[0], eval(line[1]))]
        
    return out, total_charge

def listify_bo(str_raw):
    """
    Converts a text of the form

    "***
    idA A v   B idB o   C idC o ...
    idB B v   A idA o   D idD o ...
    ..."
    (format used by tmQM/tmQM_X[n].BO after splitting by "\n\n")

    into two lists of the form
    [(idA - 1, 'A', v), (idB - 1, 'B', v), ...]
    [(idA - 1, idB - 1, 'A', 'B', o), (idA - 1, idC - 1, 'A', 'C', o), (idB - 1, idA - 1, 'B', 'A', o), (idB - 1, idD - 1, 'B', 'D', o)]
    
    Arguments:
        - str_raw: the text to process
        
    Returns:
        - the first processed list
        - the second processed list
    """
    
    str_list = str_raw.split("\n")[1:]
    
    out1 = []
    out2 = []
    for s in str_list:
        line = [c for c in s.split(" ") if not c == ""]
        
        id_source = int(line[0]) - 1
        out1 += [(id_source, line[1], eval(line[2]))]
        
        for i in range((len(line) - 3)//3):
            idx_lab = 3 + 3*i
            idx_id_tar = idx_lab + 1
            idx_ord = idx_id_tar + 1
            
            id_tar = int(line[idx_id_tar]) - 1
            
            out2 += [(id_source, id_tar, line[1], line[idx_lab], eval(line[idx_ord]))]
        
    return out1, out2

def collapse_bo_listified(bo_list):
    """
    Collapses the second output of listify_bo

    [(idA - 1, idB - 1, 'A', 'B', o), (idA - 1, idC - 1, 'A', 'C', o), (idB - 1, idA - 1, 'B', 'A', o), (idB - 1, idD - 1, 'B', 'D', o)]

    by removing element wich only differ by a permutation of idA and idB.
    
    bo_list can be seen as a list of directed edges, while the output is that same list but converted into undirected edges.
    
    Arguments:
        - d_list: the list of directed edges
    
    Returns:
        - the corresponding list of undirected edges
    """
    out = []
    
    for e in bo_list:
        duplicates = [f for f in out if (f[1] == e[0]) and (f[0] == e[1])]
        
        if len(duplicates) == 0:
            out += [e]
    
    return out

def merge_listified(xyz_list, q_list, bo_list):
    """
    Merges the output of listify_xyz

    [('A', [x, y, z]), ('B', [x, y, z]), ...],

    the first output of listify_q

    [('A', q), ('B', q), ...]

    and the first output of listify_bo 

    [(idA - 1, 'A', v), (idB - 1, 'B', v), ...]

    into a single list with format:
        
    [(idA - 1, 'A', v, q, [x, y, z]), (idB - 1, 'B', v, q, [x, y, z]), ...]

    based on the order of appearance of the elements

    Arguments:
        - xyz_list: the output of listify_xyz
        - q_list: the output of listify_q
        - bo_list: the output of listify_bo
        
    Returns:
        - the merged list
    """
    
    out = [bo_list[i] + (q_list[i][1], xyz_list[i][1]) for i in range(len(xyz_list))]
    
    return out

def validate_bo_info(tmc_id, bo_info):
    """
    Check whether the information retrieved from tmQM are consistent with those in
    tmQMg (i.e., that the labelled node set inferred from the .BO files coincides
    with the node set found in tmQMg).
    
    Arguments:
        - tmc_id: the CSD code of the TMC to validate
        - bo_info: the output of listify_bo
        
    Returns:
        - True if the sets of labels coincide, false otherwise
    """
    
    g = nx.read_gml(INPUT_FILES["tmQMg_graphs"] + tmc_id + ".gml")
    
    """
    Validate node set (same ids, same labels, same order)
    """
    V = dict(g.nodes(data = "node_label"))
    V_bo = {str(x[0]): x[1] for x in bo_info[0]}
    
    valid = V == V_bo
    
    return valid

def read_from_tmQM(tmc_id, atlas):
    """
    Extract the entry of a given TMC in the .xyz, .BO, .q and .csv files. 

    For .xyz and .bo, the correct file to look into is specified by
    the localisation dataframe (obtained from the script locate_in_tmQM.py).

    Remember :
        
    - The structure of each .xyz file is

    n
    CSD_code = XXYYZZ | q = x | S = y | Stoichiometry = ... | MND = z | yyyy-yyyy CSD
    A x y z
    B x y z

    m
    CSD_code = ZZYYXX | ...
    A x y z
    B x y z
    ...

    - The structure of the .q file is

    CSD_code = XXYYZZ | yyyy-yyyy CSD
    A q
    B q
    ...

    CSD_code = ZZYYXX | yyyy-yyyy CSD
    A q
    B q
    ...

    - The structure of each .BO file is

    CSD_code = XXYYZZ | yyyy-yyyy CSD
    idA A vA   B idB oAB   C idC oAC ...
    idB B vB   A idA oBA   D idD oBD ...

    CSD_code = ZZYYXX | yyyy-yyyy CSD
    idA A vA   B idB oAB   C idC oAC ...
    idB B vB   A idA oBA   D idD oBD ...

    where for each line, the first entry (e.g., idA A v) reports the valency index (v) of the reference atom (A, with id idA)
    and the subsequent entries (B idB oAB, C idC oAC), report the Wiberg bond order of the bond between the reference atom and
    the atom specified in the entry (so oAB is the order of (A, idA) -- (B, idB), oAB = oBA). 
    
    Arguments:
        - tmc_id: the CSD code of the TMC to look up
        - atlas: the atlas produced by map_tmQM
        
    Returns:
        - A dictionary containing all of the necessary tmQM info (empty if compilation fails)
        - A boolean: true if compilation was succesful, false othervise
    """
    
    """
    ---------------------------  .xyz  ----------------------------------
    """
    
    """
    Retrieve the index of the correct .xyz file and retrieve the corresponding entries
    """
    file_idx = atlas["localisation"].loc[atlas["localisation"]["id"] == tmc_id].iloc[0, 1]
    
    entries = atlas["xyz_raw"][file_idx - 1]
        
    """
    Retrieve the correct entry
    """
    selected_entry = None
    for entry in entries:
        if len(entry) > 0:
            """
            Extract the part 'CSD_code = XXYYZZ' 
            """
            temp = entry.split("\n")[1].split(" | ")[0]
            
            """
            Extract the name
            """
            name = temp.replace("CSD_code = ", "")
            
            if name == tmc_id:
                selected_entry = entry
                break
    
    """
    Sanity check
    """
    if selected_entry is None:
        # print("ERROR in read_from_tmQM!!! Unable to recover TMC " + tmc_id + ". It was supposed to be in file tmQM/tmQM_X" + str(file_idx) + ".xyz")
        # exit(1)
        return {}, False
    
    """
    Extract the data: 
    q, S, Stoichiometry, MND (from header, iot has to be converted to python dictionary)
    atom list (from body)
    """
    header = selected_entry.split("\n")[1].split(" | ") # ["q = x", "S = x", "Stoichiometry = zzz", "MND = x", "yyyy-yyyy CSD"]
    
    # header = full_header[1:5] # ["q = x", "S = x", "Stoichiometry = xxx", "MND = x"]
    
    # dict_raw = "{'" + ", '".join(header) + "}" # {'q = x, 'S = x, 'Stoichiometry = xxx, 'MND = x}
    # dict_raw = dict_raw.replace(" =", "':") # {'q': x, 'S': x, 'Stoichiometry': xxx, 'MND': x}
    # dict_raw = dict_raw.replace("'Stoichiometry': ", "'Stoichiometry': '") # {'q': x, 'S': x, 'Stoichiometry': 'xxx, 'MND': x}
    # dict_raw = dict_raw.replace(", 'MND'", "', 'MND'") # {'q': x, 'S': x, 'Stoichiometry': 'xxx', 'MND': x}
    # data = eval(dict_raw)
    
    xyz_info = {}
    xyz_info["q"] = eval(header[1].split(" = ")[1])
    xyz_info["S"] = eval(header[2].split(" = ")[1])
    xyz_info["Stoichiometry"] = header[3].split(" = ")[1]
    xyz_info["MND"] = eval(header[4].split(" = ")[1])
    xyz_info["CSD_years"] = header[5]
    xyz_info["atoms"] = listify_xyz(selected_entry)
    
    """
    ---------------------------  .q  ----------------------------------
    """
    
    entries = atlas["q_raw"]
        
    """
    Retrieve the correct entry
    """
    selected_entry = None
    for entry in entries:
        if len(entry) > 0:
            """
            Extract the part 'CSD_code = XXYYZZ' 
            """
            temp = entry.split("\n")[0].split(" | ")[0]
            
            """
            Extract the name
            """
            name = temp.replace("CSD_code = ", "")
            
            if name == tmc_id:
                selected_entry = entry
                break
    
    """
    Sanity check
    """
    if selected_entry is None:
        # print("ERROR in read_from_tmQM!!! Unable to recover TMC " + tmc_id + ". It was supposed to be in file tmQM/tmQM_X" + str(file_idx) + ".BO")
        # exit(1)
        return {}, False
    
    """
    Extract the information from the entry
    """
    q_info = listify_q(selected_entry)
    
    """
    ---------------------------  .BO  ----------------------------------
    """
    
    """
    Retrieve the index of the correct .BO file and retrieve the corresponding entries
    """
    file_idx = atlas["localisation"].loc[atlas["localisation"]["id"] == tmc_id].iloc[0, 2]
    
    entries = atlas["bo_raw"][file_idx - 1]
        
    """
    Retrieve the correct entry
    """
    selected_entry = None
    for entry in entries:
        if len(entry) > 0:
            """
            Extract the part 'CSD_code = XXYYZZ' 
            """
            temp = entry.split("\n")[0].split(" | ")[0]
            
            """
            Extract the name
            """
            name = temp.replace("CSD_code = ", "")
            
            if name == tmc_id:
                selected_entry = entry
                break
    
    """
    Sanity check
    """
    if selected_entry is None:
        # print("ERROR in read_from_tmQM!!! Unable to recover TMC " + tmc_id + ". It was supposed to be in file tmQM/tmQM_X" + str(file_idx) + ".BO")
        # exit(1)
        return {}, False
    
    """
    Extract and validate information from the entry
    """
    bo_info = listify_bo(selected_entry)
    valid = validate_bo_info(tmc_id, bo_info)
    
    """
    ---------------------------  .csv  ----------------------------------
    """
    csv_row = atlas["csv_raw"].loc[atlas["csv_raw"]["CSD_code"] == tmc_id]
    csv_info = csv_row.to_dict("records")[0]
    
    """
    ---------------------------  MERGE  ----------------------------------
    """
    
    """
    Integrate .xyz, .q, .BO and .csv information
    """
    data = xyz_info
    
    data["total_charge"] = q_info[1]
    data["atoms"] = merge_listified(data["atoms"], q_info[0], bo_info[0])
    data["bonds"] = collapse_bo_listified(bo_info[1])
    
    data.update(csv_info)
    
    return data, valid
    
def build_graphs(atlas, resume = True):
    """
    Builds the .gml graphs that summarise the content of tmQM. It also updates the list of viable TMCs.
    
    Arguments:
        - atlas: the atlas produced by map_tmQM
        - resume: a boolean indicating whether or not the computation should resume from the content of
                intermediate/tmQM/uNatQ_graphs. WARNING! No sanity check is performed on those files!
    """
    """
    Create directory intermediate/tmQM/uNatQ_graphs if it does not exist
    """
    if not os.path.exists(OUTPUT_FILES["tmQM_graphs"]):
        os.makedirs(OUTPUT_FILES["tmQM_graphs"])
    
    """
    Extract list of viable TMCs
    """
    with open(INPUT_FILES["viable"], "r") as f:
         viables = f.read().split("\n")
    
    print(f"Inherited {len(viables)} viable TMCs from step 2.i")
    
    """
    Establish the denomination of the properties (mathces the one used by the tmQMg).
    """
    graph_features_dictionary = {
        "q": "charge",
        "S": "spin",
        "Stoichiometry": "stoichiometry",
        "MND": "metal_node_degree",
        "CSD_years": "CSD_years",
        "Electronic_E": "electronic_energy",
        "Dispersion_E": "dispersion_energy",
        "Dipole_M": "dipole_moment",
        "Metal_q": "metal_node_natural_charge",
        "HL_Gap": "homo_lumo_gap",
        "HOMO_Energy": "homo_energy",
        "LUMO_Energy": "lumo_energy",
        "Polarizability": "polarisability"
    }
    
    # viables = ["WAJJOH"] # DEBUG ONLY!!!!!!!!!!!!!
    # print("DEBUG VERSION ONLY!!!!!!!!!!!!!!!\n"*10)
    
    """
    Convenience option:
        since the process is VERY time consuming, this option ensures that only the TMCs
        which do not already have a tmQM/.gml file will be processed, thus allowing the user to stop
        the script at any time and resuming it in a second moment
    """
    
    if resume:
        print("\nWARNING! Resuming from current state! It is assumed that the files in\n\ttmQM/processed/uNatQ_graphs\nare valid!")
        processed = [x.replace(".gml", "") for x in os.listdir(OUTPUT_FILES["tmQM_graphs"])]
    
        viables = [x for x in viables if x not in processed]
        
    print("\n\n" + str(len(viables)) + " files to process\n\n")
    
    """
    Main cycle
    """
    valid_tmcs = []
    if resume:
        valid_tmcs = [x.replace(".gml", "") for x in processed if x.endswith(".gml")]
    
    print("Processing TMCs...")
    for tmc in tqdm(viables):
        
        """
        Read data from tmQM
        """
        data, valid = read_from_tmQM(tmc, atlas)
        
        """
        If TMC is not valid, move on
        """
        if not valid:
            continue
        
        valid_tmcs += [tmc]
            
        """
        Initialise empty graph
        """
        g = nx.Graph()
        
        """
        Add complex-level features
        """
        complex_features = [f for f in data.keys() if f not in ["atoms", "bonds", "total_charge", "CSD_code"]]
        
        for f in complex_features:
            g.graph[graph_features_dictionary[f]] = data[f]
            
        """
        Add nodes and node-related features
        """
        node_ids = [str(v[0]) for v in data["atoms"]]
        
        g.add_nodes_from(node_ids)
        
        for v in data["atoms"]:
            g.nodes[str(v[0])]["valency_index"] = v[2]
            g.nodes[str(v[0])]["natural_atomic_charge"] = v[3]
            g.nodes[str(v[0])]["node_position"] = v[4]
            
        """
        Add edges and edge-related features
        """
        edges = [(str(e[0]), str(e[1])) for e in data["bonds"]]
        
        g.add_edges_from(edges)
        
        for e in data["bonds"]:
            g.edges[str(e[0]), str(e[1])]["wiberg_bond_order"] = e[4]
            
        """
        Write graph to .gml
        """
        nx.write_gml(g, OUTPUT_FILES["tmQM_graphs"] + tmc + ".gml")
        
    print("Done!")

    """
    Update list of viable TMCs.
    """
    print(f"Final number of valid TMCs: {len(valid_tmcs)}")
    with open(OUTPUT_FILES["viable"], "w") as f:
        f.write("\n".join(valid_tmcs))
        
def compress_graphs():
    """
    Compresses the content of tmQM/uNatQ_graphs into three archives, containing the files from
    A to J, from K to T and from U to Z
    """
    
    with TemporaryDirectory() as dirAJ, TemporaryDirectory() as dirKT, TemporaryDirectory() as dirUZ:
        """
        Sort all of the TMC files in the corresponding directory
        """
        print("Preprocessing files for compression...")
        
        for f in tqdm(os.listdir(OUTPUT_FILES["tmQM_graphs"])):
            if not f.endswith(".gml"):
                continue
            
            initial = f[0]
        
            target = dirAJ if initial < "K" else (dirKT if initial < "U" else dirUZ)
            
            shutil.copy(os.path.join(OUTPUT_FILES["tmQM_graphs"], f), target)
            
        """
        Compress and then delete the three temporary directories
        """
        tmQM_intermediate = os.path.abspath(os.path.join(OUTPUT_FILES["tmQM_graphs"], ".."))
        dirnames = {dirAJ: "_A_J", dirKT: "_K_T", dirUZ: "_U_Z"}
        for d in [dirAJ, dirKT, dirUZ]:
            print(f"Compressing {d}...")
            shutil.make_archive(os.path.join(tmQM_intermediate, "uNatQ_graphs" + dirnames[d]), "zip", d)
        
    print("Compression complete!")
    
#%% Main body
if __name__ == "__main__z":
    atlas = map_tmQM()
    build_graphs(atlas, True)
    compress_graphs()
