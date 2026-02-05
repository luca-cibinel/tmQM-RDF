"""
Extracts the ligands from the tmQM-RDF 1K selection
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header

from tqdm import tqdm
from pathlib import Path

import rdflib as rdf
import pandas as pd

INPUT_FILES = {
        "tmQM-RDF": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev"),
        "tmQM-RDF-1Ksel": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "data", "%s", "1k_selection.csv")
    }

OUTPUT_FILES = {
        "ligands": os.path.join(ROOT_DIR, "computational", "reconstruction", "intermediate", "ligands", "%s")
    }

# path_to_tmQM_RDF = "../../../../data/tmQM-RDF"
 
# tmQM_RDF_selection = "../../../../data/extract_selection_from_tmQM-RDF/1k_selection.csv"

# ligands_storage_directory = "./ligand_blocks"

# %% Helper functions

def dfs_ligand(tmc, csd_code, ligand_blacklist = []):
    """
    This function navigates a TMC-RDF graph depth first and extracts the ligands.
    """
    
    clean = lambda s, p, o: (rdf.term.URIRef(t.replace(csd_code, "CANDIDATE")) for t in (s, p, o))
    
    ligands = {}
    
    propagation_whitelist = [rdf.term.URIRef(p) for p in [ # Define edges on which propagation is allowed
            "resource://integreat/p5/complex/TMC/hasLigand",
            "resource://integreat/p5/ligand/ligand/isLigand",
            "resource://integreat/p5/ligand/structure/bLl",
            "resource://integreat/p5/ligand/ligand/hasAtom",
            "resource://integreat/p5/atomic/atom/isAtom",
            "resource://integreat/p5/atomic/structure/b",
            "resource://integreat/p5/ligand/bond/hasBindingAtom"
        ]]
    
    closure_whitelist = [rdf.term.URIRef(p) for p in [ # Define incoming edges which are allowed to be included in ligand
            "resource://integreat/p5/complex/TMC/hasLigand",
            "resource://integreat/p5/ligand/structure/bLc",
            "resource://integreat/p5/atomic/structure/b"
        ]]
    
    root = rdf.term.URIRef(f"resource://integreat/p5/complex/TMC/{csd_code}") # Define starting point
    
    to_visit = [root]
    visited = []
    
    key = None
    while len(to_visit) > 0:
        visit = to_visit.pop()
        visited += [visit]
        
        if visit in ligands: # Starting visit of a ligand: prepare to store triples in correct graph
            key = visit
        
        for s, p, o in tmc.triples((None, None, visit)): # Add incoming edges
                                    # (ensures that nodes about the metal centre are added,
                                    # needed for identifiability of binding nodes)
            if p in closure_whitelist:
                if key is not None:
                    ligands[key].add((clean(s, p, o)))
        
        for s, p, o in tmc.triples((visit, None, None)): # Add outgoing edges and propagate
        
            if p == propagation_whitelist[0]: # hasLigand
                
                if o.split("_")[-2] in ligand_blacklist:
                    continue
            
                ligand_blacklist += [o.split("_")[-2]]
                ligands[o] = rdf.Graph() # Prepare to store new ligand
            
            if p in propagation_whitelist:
                if o not in visited:
                    to_visit += [o]
                if key is not None:
                    ligands[key].add((clean(s, p, o)))
    
    for key in ligands: # Complete closure by adding the metal centre path
        tmc_uri, _, mc_ll_uri = next(iter(
            tmc.triples((root, rdf.term.URIRef("resource://integreat/p5/complex/TMC/hasMetalCentre"), None))
            ))
        
        # print(tmc_uri, mc_ll_uri)
        
        _, _, mc_al_uri = next(iter(
            tmc.triples((mc_ll_uri, rdf.term.URIRef("resource://integreat/p5/ligand/centre/hasAtom"), None))
            ))
        
        temp = ligands[key].serialize(format = "nt")
        
        temp = temp.replace(str(tmc_uri).replace(csd_code, "CANDIDATE"), "@TMC@")
        temp = temp.replace(str(mc_ll_uri).replace(csd_code, "CANDIDATE"), "$MC_LL$")
        temp = temp.replace(str(mc_al_uri).replace(csd_code, "CANDIDATE"), "!MC_AL!")
    
        ligands[key] = temp
    
    return ligands

# %% Main function

def extract_ligands_from_selection(selection, path_to_tmQM_RDF, ligands_storage_directory):
    print("Extracting ligands...")
    
    partitions = set(selection.loc[:,"partition"])
    ligs = {part: [] for part in partitions}
    
    for part in partitions:
        print(f"\tProcessing partition: {part}")
        
        local_storage_directory = os.path.join(ligands_storage_directory, part)
        if not os.path.exists(local_storage_directory):
            os.makedirs(local_storage_directory)
        
        tmcs = selection.loc[selection["partition"] == part, "TMC"]
        
        bl = []
        
        progress = tqdm(tmcs)
        for tmc in progress:
            g = rdf.Graph()
            g.parse(Path(os.path.join(path_to_tmQM_RDF, "graphs", f"{tmc}.ttl")))
            
            ligs = dfs_ligand(g, tmc, bl)
            
            for key, lig in ligs.items():
                with open(os.path.join(local_storage_directory, f"{key.split('_')[-2]}.nt"), "w") as f:
                    f.write(lig)
            
            progress.set_description(f"Blacklist: {len(bl)}")
# %% Main

if __name__ == "__main__":
    # earlyTM
    print("Extracting from earlyTM selection")
    selection = pd.read_csv(INPUT_FILES["tmQM-RDF-1Ksel"] % "earlyTM", index_col = 0)

    if not os.path.exists(OUTPUT_FILES["ligands"] % "earlyTM"):
        os.makedirs(OUTPUT_FILES["ligands"] % "earlyTM")
        
    extract_ligands_from_selection(selection, INPUT_FILES["tmQM-RDF"], OUTPUT_FILES["ligands"] % "earlyTM")
    
    # lateTM
    print("Extracting from lateTM selection")
    selection = pd.read_csv(INPUT_FILES["tmQM-RDF-1Ksel"] % "lateTM", index_col = 0)

    if not os.path.exists(OUTPUT_FILES["ligands"] % "lateTM"):
        os.makedirs(OUTPUT_FILES["ligands"] % "lateTM")
        
    extract_ligands_from_selection(selection, INPUT_FILES["tmQM-RDF"], OUTPUT_FILES["ligands"] % "lateTM")