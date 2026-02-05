"""
Step 2.i in data preprocessing pipeline

Locate each viable TMC (as specified in step_1_viable_tmcs.txt) inside the files tmQM/tmQM_X[1/2/3].[xyz/BO] and produce
a file intermediate/tmQM/localised.csv with structure

id,     xyz, BO
XXYYZZ  n    m

with n, m \in \{1, 2, 3\}.

It also updates the file list of viable TMCs by removing those TMCs for which it has not been possible to find complete information in tmQM.
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
        "viable": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "step_1_viable_tmcs.txt"),
        "tmQM_X_xyz": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQM", "v2024", "tmQM_X%d.xyz"),
        "tmQM_X_bo": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQM", "v2024", "tmQM_X%d.BO")
    }

OUTPUT_FILES = {
        "localised": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQMg", "localised_tmcs.csv"),
        "viable": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "step_2_i_viable_tmcs.txt")
    }

# %% Utility functions
def extract_tmcs():
    """
    Extracts the names of the TMCs referred to in the .xyz/.BO files
    
    Returns:
        - viables: the list of viable TMCs computed in step 1.iii
        - xyz: a list of three nested lists. xyz[i] contains the TMCs referred to in tmQM/tmQM_X[i + 1].xyz
        - bo: a list of three nested lists. bo[i] contains the TMCs referred to in tmQM/tmQM_X[i + 1].bo
    """
    
    """
    Extract list of viable TMCs
    """
    with open(INPUT_FILES["viable"], "r") as f:
        viables = f.read().split("\n")
        
    """
    Read each of the three .xyz files and retreive the names of the TMCs contained in each file
    
    The structure of each file is
    
    n
    CSD_code = XXYYZZ | ...
    ...
    
    m
    CSD_code = ZZYYXX | ...
    ...
    """
    xyz = []
    for i in [1, 2, 3]:
        loc_list = []
        with open(INPUT_FILES["tmQM_X_xyz"] % i, "r") as f:
            """
            Identify the single TMCs
            """
            entries = f.read().split("\n\n")
            
            print(f"Extracting from {INPUT_FILES['tmQM_X_xyz'] % i}...")
            for entry in tqdm(entries):
                if len(entry) > 0:
                    """
                    Extract the part 'CSD_code = XXYYZZ' 
                    """
                    temp = entry.split("\n")[1].split(" | ")[0]
                    
                    """
                    Extract the name
                    """
                    name = temp.replace("CSD_code = ", "")
                    
                    """
                    If TMC is viable, add it to the list
                    """
                    if name in viables:
                        loc_list += [name]
            
        xyz += [loc_list]
        
    """
    Read each of the three .BO files and retreive tha names of the TMCs contained in each file
    
    The structure of each file is
    
    CSD_code = XXYYZZ | ...
    ...
    
    CSD_code = ZZYYXX | ...
    ...
    """
    bo = []
    for i in [1, 2, 3]:
        loc_list = []
        with open(INPUT_FILES["tmQM_X_bo"] % i, "r") as f:
            """
            Identify the single TMCs
            """
            entries = f.read().split("\n\n")
            
            print(f"Extracting from {INPUT_FILES['tmQM_X_bo'] % i}...")
            for entry in tqdm(entries):
                if len(entry) > 0:
                    """
                    Extract the part 'CSD_code = XXYYZZ' 
                    """
                    temp = entry.split("\n")[0].split(" | ")[0]
                    
                    """
                    Extract the name
                    """
                    name = temp.replace("CSD_code = ", "")
                    
                    """
                    If TMC is viable, add it to the list
                    """
                    if name in viables:
                        loc_list += [name]
            
        bo += [loc_list]
    
    return viables, xyz, bo

def locate(viables, xyz, bo):
    """
    Locate viables TMCs in .xyz/.BO files
    
    Arguments:
        - viables: the list of viable TMCs produced in step 1.iii
        - xyz: the list produced by extract_tmcs
        - bo: the list produced by extract_tmcs
        
    Returns:
        - df: a dictionary of the form
            {"id": [], "xyz": [], "bo": []}
            where id is the list of viable TMCs whith available xyz/BO information, xyz and bo are the lists with
            the corresponding indices of the xyz/BO files
        - missing: a list containing the names of the TMCs for which it was not possible to retrieve xyz/bo information
    """
    
    """
    For each viable TMC, locate it inside the .xyz/.BO files (i.e., determine wether they are in file 1, 2 or 3)
    """
    df = {"id": [], "xyz": [], "bo": []}
    missing = []
    print("Locating TMCs...")
    for tmc in tqdm(viables):
        
        xyz_found = False
        bo_found = False
        
        """
        Locate in .xyz
        """
        for i in [1, 2, 3]:
            if tmc in xyz[i - 1]:
                df["xyz"] += [i]
                xyz_found = True
                break
            
        """
        Locate in .BO
        """
        for i in [1, 2, 3]:
            if tmc in bo[i - 1]:
                df["bo"] += [i]
                bo_found = True
                break
        
        """
        Save TMC name
        """
        if xyz_found and bo_found:
            df["id"] += [tmc]
        else:
            missing += [tmc]
    
    print("Done!")
    
    return df, missing

def sanity(viables, df, missing):
    """
    Sanity check: verifies that the list of viable TMCs is partitioned into
    TMCs for which xyz/BO information was found and TMCs for which it was not.
    
    Arguments:
        - viables: the list of viable TMCs produced in step 1.iii
        - df: the dictionary produced by locate
        - missing: the list produced by locate
    """
    print("viables - missing:", len(viables) - len(missing))
    print("df[id]:", len(df["id"]))
    print("df[xyz]", len(df["xyz"]))
    print("df[bo]", len(df["bo"]))

def write(viables, df, missing):
    """
    Stores the localisation of TMCs into the file intermediate/tmQM/localised_tmcs.csv
    Updates the list of viables TMCs by removing those for which xyz/BO information was not found and stores it into the file
    step_2_i_viable_tmcs.txt
    """
    
    """
    Update viable TMCs
    """
    new_viables = [tmc for tmc in viables if tmc not in missing]
    
    with open(OUTPUT_FILES["viable"], "w") as f:
        f.write("\n".join(new_viables))
        
    df = pd.DataFrame(df)
    df.to_csv(OUTPUT_FILES["localised"])

#%% Main body
if __name__ == "__main__":
    viables, xyz, bo = extract_tmcs()
    df, missing = locate(viables, xyz, bo)
    sanity(viables, df, missing)
    write(viables, df, missing)