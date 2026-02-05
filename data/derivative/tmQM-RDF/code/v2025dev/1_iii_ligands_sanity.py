"""
Step 1.iii of the data preprocessing pipeline

Performs a sanity check on the .csv files created by enode_ligand_structure.py.

It checks whether, for each TMC available in tmQMg:
    - there is data about ligands available
    - the recovered ligands cover the entirety of the TMC
    
The results are stored in two files: intermediate/tmQMg-L/no_data.txt, intermediate/tmQMg-L/no_coverage.txt, which contain, 
respectively, the list of the TMCs for which no ligand data is available and the list of TMCs that are not covered by the reconstructed ligands.

A third file, intermediate/step_1_viable_tmcs.txt, is created, containing the names of all the TMCs with full ligand data available.

Finally, the sizes of each processed TMC are stored in intermediate/tmQMg/tmc_sizes.txt, for future convenience
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header
import pandas as pd
import networkx as nx
from tqdm import tqdm

INPUT_FILES = {
		"tmQMg_L_graphs": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQMg-L", "uNatQ_graphs"),
		"tmQMg_graphs": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg", "v74.637k", "uNatQ_graphs"),
		"ligands_atoms": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQMg-L", "ligands_atoms_idx.csv")
	}

OUTPUT_FILES = {
		"no_data": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQMg-L", "no_data.txt"),
		"sizes": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQMg", "tmc_sizes.txt"),
		"no_coverage": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQMg-L", "no_coverage.txt"),
		"viable": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "step_1_viable_tmcs.txt")
	}

def check_availability():
    """
    For each TMC in tmQMg, check that a non-empty .csv file is available in intermediate/tmQMg-L/uNatQ_graphs.
    An empty .csv file only contains the header and if such a file is found, the corresponding TMC
    is stored in a list of complexes with missing ligand information.
    This list is then saved to a file.
    """
    
    no_data = []
    
    print("Checking availability of ligand data...")
    for x in tqdm(os.listdir(INPUT_FILES["tmQMg_L_graphs"])):
        with open(os.path.join(INPUT_FILES["tmQMg_L_graphs"], x), "r") as f:
            text = f.readlines()
            
            if(len(text) == 1):
                no_data += [x.replace(".csv", "")]
    print("Done!")
    
    """
    Write list of ligands with no data to file
    """
    with open(OUTPUT_FILES["no_data"], "w") as f:
        f.write("\n".join(no_data))

def check_coverage():
    """
    Verify that the ligand information extracted during step 1.ii entirely describes the related TMCs. For a given TMC,
    we say that we have complete coverage when for all atoms, except for the metal centre, it is possible to identify the ligand
    it belongs to.
    """
    
    """
    Check if TMC sizes have already been extracted from tmQMg (in a previous run).
    If they are available, read from the file intermediate/tmQMg/tmc_sizes.txt (contains a python dictionary of the
    form {XXYYZZ: size, ...}).
    If not, for each TMC in tmQMg, read the corresponding .gml graph and extract the size
    """
    sizes_available = os.path.isfile(OUTPUT_FILES["sizes"]) 
    if not sizes_available:
        print("Retrieving TMC sizes...")
        sizes = {}
        for x in tqdm(os.listdir(INPUT_FILES["tmQMg_graphs"])):
            if x.endswith(".gml"):
                g = nx.read_gml(os.path.join(INPUT_FILES["tmQMg_graphs"], x))
                sizes[x.replace(".gml", "")] = g.graph["meta_data"]["n_atoms"]
        print("Done!")
    else:
        with open(OUTPUT_FILES["sizes"], "r") as f:
            sizes = eval(f.read())

    """
    Extract the data about subgraph structure (created by 1_i_encode_ligand_subgraphs.py)
    """
    subgraph_data = pd.read_csv(INPUT_FILES["ligands_atoms"], sep = ",", header = 0)
    
    """
    For each TMC for which a size is available...
    """
    covered = {}
    print("Computing ligand coverage...")
    for tmc in tqdm(sizes.keys()):
    
        """
        Read the .csv file in tmQMg-L/uNatQ_graphs associated to the TMC
        """
        tmcdf = pd.read_csv(os.path.join(INPUT_FILES["tmQMg_L_graphs"], tmc + ".csv"), sep = ";")
        
        """
        For each ligand associated to the TMC...
        """
        n_atoms = 0
        for j in range(tmcdf.shape[0]):
            row = tmcdf.iloc[j]
            
            """
            Retrieve the number of atoms that are 'covered' by one ligand
            """
            n_atoms += len(subgraph_data.loc[subgraph_data["subgraph"] == row["subgraph"]].iloc[0, 1].strip().split(" "))
        
        covered[tmc] = n_atoms
    print("Done!")
    
    """
    For each TMC if the number of atoms covered by at least one ligand (plus the metal core, not accounted for previously),
    does not equal the number of atoms in the TMC, save the TMC in the list of TMCs with incomplete coverage
    """
    not_covered = []
    for tmc, size in covered.items():
        if not covered[tmc] + 1 == sizes[tmc]:
            not_covered += [tmc]
            
    """
    Write list of TMCs with incomplete coverage to file
    """
    with open(OUTPUT_FILES["no_coverage"], "w") as f:
        f.write("\n".join(not_covered))

def verify_overlap():
    """
    Check overlap between the list of TMC with no_data and TMc with no coverage (expected: 100%)
    """
    
    with open(OUTPUT_FILES["no_data"], "r") as f:
        no_data = f.read().split("\n")
        
    with open(OUTPUT_FILES["no_coverage"], "r") as f:
        no_coverage = f.read().split("\n")
        
    print(100*len([x for x in no_data if x in no_coverage])/len(no_data), "(Expected: 100)")

def write_viable_list():
    """
    Create the list of viable TMCs (all of those TMCs for which complete ligand information is available)
    """
        
    with open(OUTPUT_FILES["no_coverage"], "r") as f:
        no_coverage = f.read().split("\n")
    
    names = [x.replace(".csv", "") for x in os.listdir(INPUT_FILES["tmQMg_L_graphs"]) if x.endswith(".csv") and x.replace(".csv", "") not in no_coverage]
    
    with open(OUTPUT_FILES["viable"], "w") as f:
        f.write("\n".join(names))
        
#%% Main body
if __name__ == "__main__":
    check_availability()
    check_coverage()
    verify_overlap()
    write_viable_list()