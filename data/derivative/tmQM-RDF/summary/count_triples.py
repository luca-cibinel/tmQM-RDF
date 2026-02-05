# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header
import sys
sys.path.append(os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "interface"))

import os
import rdflib
import numpy as np
import tmQM_RDF_interface as tmint

from tqdm import tqdm

DATASET_VERSION = "v2025dev"

tmint.TmQMRDFGraph.path_to_tmQM_RDF = os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", DATASET_VERSION, "graphs")

OUTPUT_FILES = {
        "num_triples": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "summary", DATASET_VERSION, "num_triples.txt")
    }

#%% Main statement
if __name__ == "__main__":
    lens = []
    
    for tmc in tqdm(os.listdir(tmint.TmQMRDFGraph.path_to_tmQM_RDF)):
        if not tmc.endswith(".ttl"):
            continue
        
        g = tmint.TmQMRDFGraph(tmc.split(".")[0])
        lens += [len(g.rdf)]
    
    lens_t = []
    for tbox in tqdm(os.listdir(os.path.join(tmint.TmQMRDFGraph.path_to_tmQM_RDF, ".."))):
        if not tbox.endswith(".ttl"):
            continue
        
        g = rdflib.Graph()
        g.parse(os.path.join(tmint.TmQMRDFGraph.path_to_tmQM_RDF, "..", tbox))
        lens_t += [len(g)]
        
    print(f"lens: avg = {np.mean(lens)}, sum = {np.sum(lens)}")
    print(f"lens_t: avg = {np.mean(lens_t)}, sum = {np.sum(lens_t)}")
    print(f"total: avg = {np.mean(lens + lens_t)}, sum = {np.sum(lens + lens_t)}")
    
    with open(OUTPUT_FILES["num_triples"], "w") as f:
        f.write(f"""
TMCs: avg = {np.mean(lens)}, sum = {np.sum(lens)}
TBox + Elements + Ligands: avg = {np.mean(lens_t)}, sum = {np.sum(lens_t)}
Total: avg = {np.mean(lens + lens_t)}, sum = {np.sum(lens + lens_t)}
        """)
