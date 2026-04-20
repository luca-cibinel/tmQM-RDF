"""
Downloads the outliers.txt file from the tmQMg github repository: https://github.com/uiocompcat/tmQMg/
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Main script

from tempfile import TemporaryDirectory
import shutil
import git
import os

def download_tmQMg_outliers(
        download_dir = os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", "v1.0"), 
        tmQMg_repo = "https://github.com/uiocompcat/tmQMg/"
    ):
    """
    Downloads the outliers.txt dataset from the tmQMg Github repository
    
    Arguments:
        - download_dir: the directory to which the file should be saved
        - tmQMg_repo: the url of the Github repository
    """
    
    with TemporaryDirectory() as temp_dir:
        git.Repo.clone_from(tmQMg_repo, temp_dir)
        
        shutil.move(os.path.join(temp_dir, "scripts", "outliers.txt"), download_dir)
        
# %% Main statement
if __name__ == "__main__":
    download_tmQMg_outliers()