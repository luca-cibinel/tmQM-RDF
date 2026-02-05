"""
Step 2 in the creation of the 1k selection.

This script processes the number of occurrences of the ligands encountered in tmQM-RDF, using the
information stored in intermediate/[version]/detail.csv (see extract_info_from_tmc_for_selection.py).
In particular:
    - For each ligand, it counts the number of occurrences and stores the results in intermediate/ligands.csv
      in the format:
                                count
                ligand-id         n
                ...
                
NOTE: the ligand ID/the CSD code is used as index, so when reading them with pandas use index_col = 0
"""
# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header
import pandas as pd

from tqdm import tqdm

INPUT_FILES = {
        "details": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", "v2025dev", "details.csv")
    }

OUTPUT_FILES = {
        "ligands": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", "v2025dev", "%s", "ligands.csv")
    }

# %% Utility functions

def compute_ligand_occurrences(tm_mode):
    """
    Computes the ligand occurrences in the TMCs whose metal centre falls within the list of admissible centres.
    
    - Arguments:
        - tm_mode: whether to process early TMs (mode = 'earlyTM') or the late ones (mode = 'lateTM')
    """
    
    # Retrieve misc infos
    misc_info = pd.read_csv(INPUT_FILES["details"], index_col = 0)
    
    # Filter misc info by metal centres
    if tm_mode == "earlyTM":
        centres = ["Cr", "Mo", "W"]
    elif tm_mode == "lateTM":
        centres = ["Pt", "Ni", "Pd"]
    else:
        raise Exception(f"tm_mode {tm_mode} not recognised!")
        
    misc_info_index_by_centre = [i for i in range(len(misc_info.index)) if misc_info.iloc[i].iloc[0] in centres]

    # Count each ligand met in tmQM-RDF
        
    lig_count = {}
    
    for i in tqdm(misc_info_index_by_centre):
        ligs = misc_info.iloc[i].iloc[1].split(" ")
        
        for l in set(ligs):
            lig_count[l] = lig_count.get(l, 0) + 1
    
    """
    Store the result first as a pandas DF and then as a .csv file
    """
    lig_count_df = pd.DataFrame.from_dict(lig_count, orient = "index", columns = ["count"])
    lig_count_df = lig_count_df.sort_values(by = "count")
    lig_count_df.to_csv(OUTPUT_FILES["ligands"] % tm_mode)
    
# %% Main statement

if __name__ == "__main__":
    compute_ligand_occurrences("lateTM")
    compute_ligand_occurrences("earlyTM")