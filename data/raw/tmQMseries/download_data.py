"""
Downloads the raw tmQM series.
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header
from tempfile import TemporaryDirectory
import urllib.request as url
import zipfile
import shutil
import gzip
import git

OUTPUT_FILES = {
        "tmQM": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQM"),
        "tmQMg": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg"),
        "tmQMg-L": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L")
    }

URLS = {
        "tmQM": {
                "v2024": "https://github.com/uiocompcat/tmQM/"
            },
        "tmQMg": {
                "v74.637k": "https://data.archive.sigma2.no/dataset/cc354a73-7398-487f-83f8-4166caa8cc09/download/nird/home/hanneskn/tmQMg/uNatQ_graphs.zip",
                "v74.548k": "https://data.archive.sigma2.no/dataset/4f94f626-b18c-458c-946f-d7052bb05982/download/u-NatQ_graphs.zip"
            },
        "tmQMg-L": {
                "v60k": "https://github.com/hkneiding/tmQMg-L/archive/refs/tags/v60k.zip",
                "v74k": "https://github.com/uiocompcat/tmQMg-L/archive/refs/tags/v74k.zip"
            }
    }

# %% Utility functions
def download_tmQM(xtmQM_dir = "./data/tmQM/", xtmQM_repo = "https://github.com/uiocompcat/tmQM/"):
    """
    Downloads the tmQM dataset from its Github repository into a temporary directory and
    then saves the necessary files into the specified directory. The temporary directory is
    automatically deleted at the end.
    
    Arguments:
        - tmQM_dir: the directory to which the files should be saved
        - tmQM_repo: the url of the Github repository
    """
    
    if not os.path.exists(OUTPUT_FILES["tmQM"]):
        os.makedirs(OUTPUT_FILES["tmQM"])
    
    for v, tmQM_url in URLS["tmQM"].items():
        print(f"v = {v}")
        
        tmQM_dir = os.path.join(OUTPUT_FILES["tmQM"], v)
        if not os.path.exists(tmQM_dir):
            os.makedirs(tmQM_dir)
        
        with TemporaryDirectory() as temp_dir:
            print("\n\nCloning tmQM...")
            git.Repo.clone_from(tmQM_url, temp_dir)
            
            for fname in os.listdir(os.path.join(temp_dir, "tmQM")):
                print(f"Saving {fname}...")
                
                shutil.move(
                    os.path.join(temp_dir, "tmQM", fname), 
                    os.path.join(tmQM_dir, fname)
                )
                
                if fname.endswith(".gz"):
                    print(f"\tDecompressing {fname}...")
                    with gzip.open(os.path.join(tmQM_dir, fname), 'rb') as f_in:
                        with open(os.path.join(tmQM_dir, fname.replace(".gz", "")), 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                            
                    print(f"\tRemoving {fname}...")
                    os.remove(os.path.join(tmQM_dir, fname))

"""
Stores the data from tmQMg into the ./data/tmQMg/uNatQ_graphs folder.
The dataset is a compressed archive which contains .gml graphs inside a folder named
uNatQ_graphs.
"""

def download_tmQMg(xtmQMg_dir = "./data/tmQMg/", xtmQMg_url = "https://data.archive.sigma2.no/dataset/cc354a73-7398-487f-83f8-4166caa8cc09/download/nird/home/hanneskn/tmQMg/uNatQ_graphs.zip"):
    """
    Downloads the tmQMg dataset as a compressed zip archive and extracts it into the specified directory.
    
    Arguments:
        - tmQMg_dir: the directory to which the files should be saved
        - tmQMg_url: the url of the zip archive
    """
    
    if not os.path.exists(OUTPUT_FILES["tmQMg"]):
        os.makedirs(OUTPUT_FILES["tmQMg"])
    
    for v, tmQMg_url in URLS["tmQMg"].items():
        tmQMg_dir = os.path.join(OUTPUT_FILES["tmQMg"], v)
        
        if not os.path.exists(tmQMg_dir):
            os.makedirs(tmQMg_dir)
        
        print(f"\n\nv = {v}")
        print("Downloading tmQMg zip archive...")
        url.urlretrieve(
                tmQMg_url,
                os.path.join(tmQMg_dir, "uNatQ_graphs.zip")
            )
        
        print("Extracting archive...")
        with zipfile.ZipFile(os.path.join(tmQMg_dir, "uNatQ_graphs.zip")) as zipf:
            zipf.extractall(tmQMg_dir)
        os.remove(os.path.join(tmQMg_dir, "uNatQ_graphs.zip"))

"""
Stores the data from tmQM into the ./data/tmQMg-L/ folder.
"""

def download_tmQMg_L(xtmQMg_L_dir = "./data/tmQMg-L/", xtmQMg_L_url = "https://github.com/hkneiding/tmQMg-L/archive/refs/tags/v60k.zip"):
    """
    Downloads the tmQMg-L dataset from its Github repository into a temporary directory and
    then saves the necessary files into the specified directory. The temporary directory is
    automatically deleted at the end.
    
    For reproducibility reasons, the 60k release is downloaded.
    
    Arguments:
        - tmQMg_L_dir: the directory to which the files should be saved
        - tmQMg_L_repo: the url of the Github repository
    """
    
    if not os.path.exists(OUTPUT_FILES["tmQMg-L"]):
        os.makedirs(OUTPUT_FILES["tmQMg-L"])
    
    for v, tmQMg_L_url in URLS["tmQMg-L"].items():
        print(f"v = {v}")
        
        tmQMg_L_dir = os.path.join(OUTPUT_FILES["tmQMg-L"], v)

        if not os.path.exists(tmQMg_L_dir):
            os.makedirs(tmQMg_L_dir)
        
        with TemporaryDirectory() as temp_dir:
            print("Downloading tmQMg-L zip archive...")
            url.urlretrieve(
                    tmQMg_L_url,
                    os.path.join(temp_dir, "tmQMg-L.zip")
                )
            print("Extracting archive...")
            with zipfile.ZipFile(os.path.join(temp_dir, "tmQMg-L.zip")) as zipf:
                zipf.extractall(temp_dir)
            
            # Recover the name of the extracted folder
            zip_folder_name = [fname for fname in os.listdir(temp_dir) if not fname.endswith(".zip")][0]
            
            for fname in [
                    "ligands_descriptors.csv", 
                    "ligands_fingerprints.csv", 
                    "ligands_misc_info.csv", 
                    os.path.join("xyz", "ligands_xyzs.xyz"), 
                    os.path.join("xyz", "ligand_xyzs.xyz")]:
                
                print(f"Saving {fname}...")
                
                try:
                    shutil.move(
                        os.path.join(temp_dir, zip_folder_name, fname), 
                        os.path.join(tmQMg_L_dir, fname.split("/")[-1])
                    )
                except Exception:
                    print("\tFile not found in this version. Skipping...")
                
                
    
#%% Main body
if __name__ == "__main__":
    download_tmQM()
    download_tmQMg()
    download_tmQMg_L()
    