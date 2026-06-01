# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header
import sys
sys.path.append(os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "interface", "v1.0"))

import os
import rdflib
import numpy as np
import tmQM_RDF_interface as tmint

from tqdm import tqdm

DATASET_VERSION = "v1.0.1"

OUTPUT_FILES = {
        "num_triples": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "summary", DATASET_VERSION, "num_triples.txt")
    }

#%% Main statement
if __name__ == "__main__":
    instance = tmint.TmQMRDF(os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", DATASET_VERSION))
    
    lens = {
        "lens_tmc": [],
        "lens_lig": [],
        "lens_cntr": [],
        "lens_el": []
    }

    lens_tmc = []
    lens_lig = []
    lens_cntr = []
    lens_el = []
    
    # Assertions
    for tmc in tqdm(os.listdir(os.path.join(instance.path, "assertions", "TMCs"))):
        if not tmc.endswith(".ttl"):
            continue
        
        g = instance.tmc(tmc.split(".")[0])
        lens["lens_tmc"] += [len(g.rdf)]
        instance.clear()

    for lig in tqdm(os.listdir(os.path.join(instance.path, "assertions", "ligands"))):
        if not lig.endswith(".ttl"):
            continue
        
        g = instance.ligand(lig.split(".")[0])
        lens["lens_lig"] += [len(g.rdf)]
        instance.clear()

    for c in tqdm(os.listdir(os.path.join(instance.path, "assertions", "centres"))):
        if not c.endswith(".ttl"):
            continue
        
        g = instance.centre(c.split(".")[0].split("_")[1])
        lens["lens_cntr"] += [len(g.rdf)]
        instance.clear()

    for el in tqdm(os.listdir(os.path.join(instance.path, "assertions", "elements"))):
        if not el.endswith(".ttl"):
            continue
        
        g = instance.element(os.path.join("..", "elements", el.split(".")[0]))
        lens["lens_el"] += [len(g.rdf)]
        instance.clear()
    
    lens["lens_a"] = sum(lens.values(), [])
    
    # Terminology
    lens["lens_t"] = []
    for prefix, _, files in os.walk(os.path.join(instance.path, "terminology")):
        for tbox in files:
            if not tbox.endswith(".ttl"):
                continue
            
            g = rdflib.Graph()
            g.parse(os.path.join(prefix, tbox))
            lens["lens_t"] += [len(g)]
    
    # Log
    log = []
    for name, n_trip in lens.items():
        line = f"{name}: avg = {np.mean(n_trip)}, sum = {np.sum(n_trip)}"
        log += [line]
        print(line)
    
    with open(OUTPUT_FILES["num_triples"], "w") as f:
        f.write("\n".join(log))
