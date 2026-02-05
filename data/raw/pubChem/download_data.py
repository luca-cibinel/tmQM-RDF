"""
Downloads the pubChem csv version of the periodic table
"""

import urllib.request as url
import os

"""
Downloads the pubChem csv version of the periodic table
""" 
def download_misc(misc_dir = "./data/", pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/periodictable/CSV?response_type=save&response_basename=PubChemElements_all"):
    """
    Downloads the PubChem csv version of the periodic table from the official PubChem website and
    also sets up the ./misc/ directory.
    
    Arguments:
        - misc_dir: the directory to which the files should be saved
        - pubchem_url: the url of the csv version of the periodic table
    """
    
    if not os.path.exists(misc_dir):
        os.makedirs(misc_dir)
       
    url.urlretrieve(
            pubchem_url,
            os.path.join(misc_dir, "PubChemElements_all.csv")
        )
    
# %% Main
if __name__ == "__main__":
    download_misc()