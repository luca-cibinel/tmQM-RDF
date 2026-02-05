"""
Step 1.2 in the data preprocessing pipeline

Links every TMC in tmQMg to its ligands by combining the data from tmQMg-L/ligands_misc_info.csv and from intermediate/tmQMg-L/ligands_atoms_idx.csv.

For every ligand in tmQMg-L, the file tmQMg-L/ligands_misc_info.csv contains the list of all of its instances, in the form of a dictionary:
    {XXYYZZ-subgraph-n: [...], ...}

The list associated to each subgraph is the list of its coordinating atoms (denticity/hapticity specified by nested lists), where each atom is
identified by its position in the list of atoms specified in tmQMg-L/ligands_xyzs.xyz for the given subgraph. This list has been translated
into a list of tmQMg node ids by 1_i_encode_ligand_subgraphs.py.

This script creates a .csv file for each TMC, containing:
    - the ligand id
    - the subgraph the ligand is instanciated to
    - the list of coordinating atoms for that instance (using the translation into tmQMg node ids)
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header
import shutil
import pandas as pd
from tqdm import tqdm

INPUT_FILES = {
        "tmQMg_graphs": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg", "v74.637k", "uNatQ_graphs"),
        "ligands_misc_info": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", "v60k", "ligands_misc_info.csv"),
        "ligands_atoms": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQMg-L", "ligands_atoms_idx.csv")
    }

OUTPUT_FILES = {
        "checkpoint":  os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQMg-L", "step_1_ii_checkpoint.txt"),
        "tmQMg_L_graphs": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQMg-L", "uNatQ_graphs")
    }

# %% Utility classes
class TqdmCallback(tqdm):
    """
    Helper class: extends tqdm to show a progress bar associated with an iterable object.
    A callback function is called every time that the percentage of completion moves forward of
    a quantity specified by dp.

    NOTE: the callback is called at the END of the iteration.
    """
    
    def __init__(self, iterable, callback, dp = 2):
        self.callback = callback
        
        self.dp = dp
        self.perc = 0
        
        super().__init__(iterable)
    
    def update(self, n=1):
        displayed = super().update(n)
        
        if displayed and (self.n/self.total*100 > self.perc + self.dp):
            self.callback()
            self.perc += self.dp
            
        return displayed

# %% Utility functions
def translate(source, target):
    """
    Takes in input a source with the following nested structure:
        each element is either an integer or a list with the same structure as source
        
    The target is a list of arbitrary elements.

    The output is a list with the same nested structure as the source but with each integer replaced
    by the corresponding element in the target.

    EXAMPLE:
        source = [0, [0, 1], [[2], 3]]
        target = ["a", "b", "c", "d"]
        
        output = ["a", ["a", "b"], [["c"], "d"]]
    """
    
    if not type(source) == type([]):
        return target[source]
    
    return [translate(x, target) for x in source]

def get_baseline_csv_files():
    """
    For each TMC available in tmQMg, prepares a baseline string which represents the content
    of the future .csv file to be later populated with ligand data.
    
    Returns:
        - a dictionary of the form
            {
                XXYYZZ: "ligand_id;subgraph;coordinating_atoms"
            }
          for each CSD code (XXYYZZ) present in tmQMg
    """
    
    """
    Extract the names of the TMC complexes available in tmQMg
    """
    print("Retrieving names...")
    names = []
    for x in os.listdir(INPUT_FILES["tmQMg_graphs"]):
        if x.endswith(".gml"):
            names.append(x.replace(".gml", ""))
    print("Names retrieved!")
    
    """
    For each TMC available in tmQMg create a baseline string
    """
    print("Forging baseline files...")
    uNatQ_graphs = {}
    for n in names:
        uNatQ_graphs[n] = ["ligand_id;subgraph;coordinating_atoms"]
    print("Baseline files forged!")
    
    return uNatQ_graphs

def compile_csv_files(uNatQ_graphs):
    """
    Compiles the .csv files associated to the TMCs present in tmQMg according to the logic described in
    the file header.
    
    Arguments:
        - uNatQ_graphs: the dictionary produced by get_baseline_csv_files
        
    Returns:
        - The same dictionary populated with ligand data
    """
    
    """
    Extract the ligand and subgraph data (coming from tmQMg-L and encode_ligand_subgraphs.py respectively)
    as pandas dataframes
    """
    ligand_data = pd.read_csv(INPUT_FILES["ligands_misc_info"], sep = ";")
    subgraph_data = pd.read_csv(INPUT_FILES["ligands_atoms"], sep = ",", header = 0)
    
    """
    Backup variables: periodically save progress so as to allow the computations
    to be run in different sessions
    """
    total = ligand_data.shape[0]
    checkpoint = 0
    if os.path.exists(OUTPUT_FILES["checkpoint"]):
        with open(OUTPUT_FILES["checkpoint"], "r") as f:
            checkpoint, uNatQ_graphs = eval(f.read())
    
    def backup():
        with open(OUTPUT_FILES["checkpoint"], "w") as f:
            f.write(str((i + 1, uNatQ_graphs)))
    
    """
    Main cycle: for each ligand...
    """
    print(f"Translating coordinating atoms... (resuming from row {checkpoint + 1} of ligands_misc_info.csv)")
    for i in TqdmCallback( range(checkpoint, total), backup ):
        """
        Extract the information about the ligand available in tmQMg-L
        """
        row = ligand_data.iloc[i]
        
        """
        Extract the ligand identifier and the dictionary of the coordinating atoms for each instance of the ligand
        An instance is a subgraph id paired with a list (of lists) of coordinating atoms
        """
        lname = row["name"]
        instances = eval(row["metal_bond_node_idx_groups"])
        
        """
        For each instance of the ligand...
        """
        for subg, atoms in instances.items():
            """
            Extract the name of the TMC from the subgraph id
            """
            tmcname = subg[0:6]
            
            """
            Extract the list of atoms that compose the ligand (indexed by the tmQMg node ids) and convert the ids to integers
            """
            
            structure = [
                int(x) for x in
                subgraph_data.loc[subgraph_data["subgraph"] == subg].iloc[0, 1].strip().split(" ")
            ]
            
            """
            In the (string corresponding to the).csv file associated to the TMC add a line 
            corresponding to the instanced ligand and add the list 
            of coordinating atoms (translated from the index in the .xyz file to the indexes in the tmQMg entries)
            """
            # with open("./tmQMg-L/uNatQ_graphs/" + tmcname + ".csv", "a") as f:
            #     f.write(lname + "; " + subg + "; " + str(translate(atoms, structure)) + "\n")
            uNatQ_graphs[tmcname] += [lname + ";" + subg + ";" + str(translate(atoms, structure))]
                
    print("Coordinating atoms translated!")
    
    return uNatQ_graphs

def write_to_csv(uNatQ_graphs):
    """
    Writes the .csv files previously prepared to actual .csv files
    
    Arguments:
        - uNatQ_graphs: the dictionary produced by compile_csv_files
    """
    
    """
    Create directory tmQMg-L/uNatQ_graphs if it does not exist
    """
    if not os.path.exists(OUTPUT_FILES["tmQMg_L_graphs"]):
        os.makedirs(OUTPUT_FILES["tmQMg_L_graphs"])
    
    """
    Secondary cycle: write the .csv files
    """
    print("Writing .csv files...")
    for tmc, content in tqdm(uNatQ_graphs.items()):
        """
        Write content to file
        """
        with open(OUTPUT_FILES["tmQMg_L_graphs"] + tmc + ".csv", "w") as f:
             f.write("\n".join(content))
    
    print(".csv files created!")

def compress_csv_files():
    """
    Creates a .zip archive with the newly computed .csv files
    """
    
    print("Compressing...")
    shutil.make_archive(OUTPUT_FILES["tmQMg_L_graphs"][:-1], "zip", OUTPUT_FILES["tmQMg_L_graphs"])
    print("Done!")

#%% Main body

if __name__ == "__main__":
    uNatQ_graphs = get_baseline_csv_files()
    uNatQ_graphs = compile_csv_files(uNatQ_graphs)
    write_to_csv(uNatQ_graphs)
    compress_csv_files()
    