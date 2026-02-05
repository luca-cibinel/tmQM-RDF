"""
Step 1.i in the data preprocessing pipeline

From the file tmQMg-L/ligands_xyzs.xyz, for each subgraph extracts the coordinates of the atoms,
then mathces them to the coordinates stored in the .gml files in tmAMg/uNatQ_graphs.
Then creates the file intermediate/tmQMg-L/ligands_atoms_idx.csv file in which each subgraph is matched with a list of indices, which correspond
to the node ids of the atoms whose coordinates were stored in tmQMg-L/ligands_xyzs.xyz

EXAMPLE:
    
    - in tmQMg-L/ligands_xyzs.xyz
    
    2
    XXYYZZ-subgraph-0
    A1 x1  y1  z1
    A2 x2  y2  z2
    
    - in tmQMg/uNatQ_graphs/XXYYZZ.gml
    
    node [
        id "0"
        node_position x1
        node_position y1
        node_position z1    
    ]
    ...
    node [
        id "3"
        node_position x2
        node_position y2
        node_position z2    
    ]
    
    =>
    
    - in intermediate/tmQMg-L/ligands_atoms_idx.csv
    
    XXYYZZ-subgraph-0, 0 3
    
    
NOTE:
    
The format used by tmQMg-L/ligands_xyzs.xyz is

[n atoms]
[XXYYZZ]-subgraph-[n]
A1 x1  y1  z1
A2 x2  y2  z2


[n atoms]
[XXYYZZ]-subgraph-[n]
A1 x1  y1  z1
A2 x2  y2  z2

...

subgraphs are separated by the token "\n\n"
the second line of each subgraph description is the subgraph identifier, with the first six charachters being the TMC id

EXAMPLE OF USAGE:
    import pandas as pd
    
    ddf = pd.read_csv("tmQMg-L/ligands_atoms_idx.csv", sep = ",", header = 0)
    ddf.loc[ddf["subgraph"] == "UKEBAN-subgraph-0"].iloc[0, 1].strip().split(" ")
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header

import networkx as nx
from tqdm import tqdm

INPUT_FILES = {
        "ligands_xyzs": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", "v60k", "ligands_xyzs.xyz"),
        "tmQMg_graphs": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg", "v74.637k", "uNatQ_graphs")
    }

OUTPUT_FILES = {
        "ligands_atoms": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQMg-L", "ligands_atoms_idx.csv")
    }

# %% Utility functions
def listify(str_raw):
    """
    Converts a raw text of the form

    "
    ***
    ***
    A x y z
    B x y z
    ...
    "

    with x y z floats
    (format used by tmQMg-L/ligands_xyzs.xyz after splitting by "\n\n")

    into a list of the form
    [[x, y, z], [x, y, z]]
    (each element is a coordinate set)
    
    Arguments:
        - str_raw: the raw text (as specifeid above)
        
    Returns:
        - the output list described above 
    """
    
    str_list = str_raw.split("\n")[2:]
    
    out = []
    for s in str_list:
        if len(s) > 0:
            out.append([float(x) for x in s.split(" ")[1:]])
        
    return out

def get_subgraphs_coordinates():
    """
    Extracts each subgraph description from tmQMg-L/ligands_xyzs.xyz and store each description under its subgraph id
    
    Returns:
        - a dictionary of the form 
            {
                [XXYYZZ]-subgraph-[n]:
                '***
                 ***
                 A x y z
                 B x y z
                 ...'
            }
    """
    xyzs = {}
    with open(INPUT_FILES["ligands_xyzs"], "r") as fxyz:
        for xyz in fxyz.read().split("\n\n"):
            xyzs[xyz.split("\n")[1]] = xyz
    
    return xyzs


def get_indices_from_coordinates(xyzs):
    """
    Converts each subgraph description into a list of indices, as described in the file header.
    
    Arguments:
        - xyzs: the dictionary produced bu get_subgraphs_coordinates
    
    Returns:
        - a dictionary of the form
            {
                '[XXYYZZ]-subgraph-[n]': 'id0 id1 ...'
            }
    """
    
    l2 = lambda x, y: sum((x[i] - y[i])**2 for i in range(3)) # l2 norm helper function
    
    for subg in tqdm(xyzs.keys()): 
        """
        Extract TMC id and the corresponding .gml graph
        """
        gname = subg.split("-")[0]
        g = nx.read_gml(os.path.join(INPUT_FILES["tmQMg_graphs"], gname + ".gml"))
        
        """
        Extract the list of nodes with the corresponding position (it is a list of the form [(id, [x, y, z]), ...])
        """
        nodes = list(g.nodes(data = "node_position"))
        
        """
        For each coordinate set in the listified subgraph description, look for a match in the list of nodes,
        if match is found, store the id of the match and stop the search (max 1 match per coordinate set).
        Due to numerical rounding, a match between two set of coordinates is defined as being closer than 0.001 in l2 norm
        """
        ids = ""
        for xyz in listify(xyzs[subg]):
            for n in nodes:
                if l2(n[1], xyz) <= 0.001:
                    ids += n[0] + " "
                    break
                
        xyzs[subg] = ids
        
    return xyzs
    
def write_to_csv(xyzs):
    """
    Writes the lists of atom indices associated to each subgraph into a csv file of the form
    
    subgraph, indices
    [XXYYZZ]-subgraph-[n], id0 id1 ...
    [XXYYZZ]-subgraph-[n], id0 id1 ...
    
    Arguments:
        - xyzs: the dictionary produced bu get_indices_from_coordinates
    """
    with open(OUTPUT_FILES["ligands_atoms"], "w") as f:
        text = "subgraph, indices\n"
        for subg in xyzs.keys():
            text += subg + ", " + xyzs[subg] + "\n"
            
        f.write(text)

#%% Main body
if __name__ == "__main__":
    xyzs = get_subgraphs_coordinates()
    xyzs = get_indices_from_coordinates(xyzs)
    write_to_csv(xyzs)

