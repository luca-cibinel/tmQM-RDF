"""
Step 1 in the creation of the 1k selection.

This script extracts from the various tmQM datasets the information used to
select the 1000 complexes that will be included in the 1k selection.
In particular:
    - For each metal centre encountered in the tmQM-RDF dataset, it counts the number of complexes
        which possess that centre. This information is read from the .gml files in tmQMg and 
        it is stored in the file intermediate/centres.csv in the format:
                        count
                El        n
                ...
        
    - For each complex in the tmQM-RDF dataset, it summarises the metal core and the ligands (reduntant 
        information, but it makes it easier to access). This information is read from intermediate/tmQMg-L/uNatQ_graphs 
         and stored in the file intermediate/detail.csv 
        in the format:
                            core    ligands
                XXYYZZ        C       lig1 lig2 ...
                ...
                
NOTE: in both files, the element label/the CSD code (XXYYZZ) is used as index, so when 
        reading them with pandas use index_col = 0
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header
import networkx as nx
import pandas as pd

from tqdm import tqdm

INPUT_FILES = {
        "viable": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "viable_tmcs.txt"),
        "tmQMg_graphs": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg", "v74.637k", "uNatQ_graphs"),
        "tmQMg_L_graphs": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQMg-L", "uNatQ_graphs")
    }

OUTPUT_FILES = {
        "centres": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", "v2025dev", "centres.csv"),
        "details": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", "v2025dev", "details.csv"),
        "sizes": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", "v2025dev", "tmc_sizes.csv")
    }

# %% Utility functions
def extract_tmc_info():
    """
    Extracts info regarding the viable TMCs in tmQM-RDF, namely metal center counts and ligand composition.
    """
    
    """
    Get the TMCs which are viable for tmQM-RDF
    """
    with open(INPUT_FILES["viable"], "r") as f:
        tmc_list = f.read().split("\n")
    
    """
    Main loop: read and summarise the information
    """
    met_cores = {}
    tmc_info = {}
    tmc_sizes = {}
    
    for tmc in tqdm(tmc_list):
        """
        Read the .gml graph in tmQMg
        """
        g = nx.read_gml(os.path.join(INPUT_FILES["tmQMg_graphs"], tmc + ".gml"))
        
        """
        Extract and store the metal core
        """
        met_core = g.graph["meta_data"]["metal_center_element"]
        
        met_cores[met_core] = met_cores.get(met_core, 0) + 1
        
        """
        Extract (from tmQMg-L) and store the ligand information
        """
        tmc_df = pd.read_csv(os.path.join(INPUT_FILES["tmQMg_L_graphs"], tmc + ".csv"), sep = ";")
        ligs = " ".join(list(tmc_df.loc[:, "ligand_id"]))
        
        tmc_info[tmc] = [met_core, ligs]
        
        """
        Extract and store TMC size
        """
        tmc_sizes[tmc] = g.graph["meta_data"]["n_atoms"]
        
    """
    Store the data: metal cores count
    """
    met_cores_df = pd.DataFrame.from_dict(met_cores, orient = "index", columns = ["count"])
    met_cores_df.to_csv(OUTPUT_FILES["centres"])
    
    """
    Store the data: TMC detailed information
    """ 
    tmc_info_df = pd.DataFrame.from_dict(tmc_info, orient = "index", columns = ["core", "ligands"])
    tmc_info_df.to_csv(OUTPUT_FILES["details"])
    
    """
    Store the data: TMC size
    """
    with open(OUTPUT_FILES["sizes"], "w") as f:
        f.write(str(tmc_sizes))

# %% Main statement

if __name__ == "__main__":
    extract_tmc_info()