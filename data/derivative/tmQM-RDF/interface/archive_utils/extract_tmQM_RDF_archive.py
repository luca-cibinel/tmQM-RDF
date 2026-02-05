"""
This is a utility script designed to decompress all the archives that contain the tmQM-RDF data into their original
configuration, namely

data/
    graphs/
        ...
    tmQM-RDFS-*.ttl
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header
from tempfile import TemporaryDirectory

import zipfile

DATA_DIR = os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev")

# %% Utility functions

def extract_rdfs_headers(data_dir):
    """
    Extracts the RDFS headers from [data_dir]/archive/tmQM-RDFS.zip into [data_dir]
    
    - Arguments:
        - data_dir: the directory where the data ought to be stored
    """
    
    print("Decompressing RDFS headers...")
    
    archive_dir = os.path.join(data_dir, "archive")
    
    with zipfile.ZipFile(os.path.join(archive_dir, "tmQM-RDFS.zip"), "r") as zipf:
        zipf.extractall(data_dir)

def extract_rdf_graphs(data_dir):
    """
    Extracts the RDF graphs from [data_dir]/archive/graphs/[letter].zip into [data_dir]/graphs
    
    - Arguments:
        - data_dir: the directory where the data ought to be stored
    """
    
    archive_graphs_dir = os.path.join(data_dir, "archive", "graphs")
    graphs_dir = os.path.join(data_dir, "graphs")
    
    if not os.path.exists(graphs_dir):
        os.makedirs(graphs_dir)
    
    for zf in os.listdir(archive_graphs_dir):
        if not zf.endswith(".zip"):
            continue
        
        print(f"Decompressing {zf}...")
        
        with zipfile.ZipFile(os.path.join(archive_graphs_dir, zf), "r") as zipf:
            zipf.extractall(graphs_dir)
        
# %% Main statement
if __name__ == "__main__":
    extract_rdfs_headers(DATA_DIR)
    extract_rdf_graphs(DATA_DIR)
          