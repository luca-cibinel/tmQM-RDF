"""
This is a utility script designed to decompress all the archives that contain the intermediate tmQM graph data into their original
configuration, namely

tmQM/
    uNatQ_graphs/
        ...
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header
from tempfile import TemporaryDirectory

import zipfile

DATA_DIR = os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "tmQM")

# %% Utility functions

def extract_uNatQ_graphs(data_dir):
    """
    Extracts the tmQM intermediate graphs from [tmQM_intermediate_dir]/uNatQ_graphs_[letter_from]_[letter_to].zip 
    into [tmQM_intermediate_dir]/uNatQ_graphs
    
    - Arguments:
        - data_dir: the directory where the data ought to be stored
    """
    
    graphs_dir = os.path.join(data_dir, "uNatQ_graphs")
    
    if not os.path.exists(graphs_dir):
        os.makedirs(graphs_dir)
    
    for zf in os.listdir(data_dir):
        if not zf.endswith(".zip"):
            continue
        
        print(f"Decompressing {zf}...")
        
        with zipfile.ZipFile(os.path.join(data_dir, zf), "r") as zipf:
            zipf.extractall(graphs_dir)
        
# %% Main statement
if __name__ == "__main__":
    extract_uNatQ_graphs(DATA_DIR)
          