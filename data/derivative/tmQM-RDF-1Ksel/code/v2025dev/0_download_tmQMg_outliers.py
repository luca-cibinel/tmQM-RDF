"""
Downloads the outliers.txt file from the tmQMg github repository: https://github.com/uiocompcat/tmQMg/
"""

from tempfile import TemporaryDirectory
import shutil
import git
import os

def download_tmQMg_outliers(download_dir = "../intermediate", tmQMg_repo = "https://github.com/uiocompcat/tmQMg/"):
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