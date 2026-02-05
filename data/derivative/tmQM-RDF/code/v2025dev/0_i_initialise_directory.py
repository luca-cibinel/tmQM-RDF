"""
Step 0 of the data preprocessing pipeline

Creates all the necessary subfolders
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% HEADER
DATASET_VERSION = "v2025dev"

# %% Utility functions

def create_intermediate():
    """
    Creates the subtree
    
    tmQM-RDF/
        intermediate/
            tmQM/
            tmQMg/
            tmQMg-L/
    """
    
    for f in ["tmQM", "tmQMg", "tmQMg-L"]:
        os.makedirs(os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", DATASET_VERSION, f))
        
def create_data():
    """
    Creates the subtree
    
    tmQM-RDF/
        data/
    """
    
    os.makedirs(os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", DATASET_VERSION))
    
# %% Main statement
if __name__ == "__main__":
    create_intermediate()
    create_data()