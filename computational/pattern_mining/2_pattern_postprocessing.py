# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
     
# %% Header

from tqdm import tqdm
import gzip
import re

INPUT_FILES = {
        "patterns": os.path.join(ROOT_DIR, "computational", "pattern_mining", "results", "%s"),
        "pattern_file_regexp": "output-size-\d*.dat.gz"
    }

# %% Utility functions
def clean_pattern_file(tm_mode, file):
    """
    Cleans a pattern file (of the form output-size%d.dat.gz) by removing all the unnecessary match details
    
    Arguments:
        - tm_mode: whether to process early TMs (mode = 'earlyTM') or the late ones (mode = 'lateTM')
        - file: the file name
    """
    
    with gzip.open(os.path.join(INPUT_FILES["patterns"] % tm_mode, file), "rt") as f:
        lines = f.readlines()
    
    newlines = []
    for l in tqdm(lines):
        if l.startswith("0"):
            newlines += [l]
            
        if l.startswith("6"):
            newlines += [l]
            #newlines += [l.replace("resource://integreat/p5/", "resource://integreat/tmQM/RDF/")]
            
    with gzip.open(os.path.join(INPUT_FILES["patterns"] % tm_mode, file), "wt") as f:
        f.writelines(newlines)
        
def clean_patterns(tm_mode):
    """
    Cleans all the pattern files (of the form output-size%d.dat.gz) inside the results directory specified by tm_mode.
    
    Arguments:
        - tm_mode: whether to process early TMs (mode = 'earlyTM') or the late ones (mode = 'lateTM')
    """
    
    print(f"\n\n--- tm_mode: {tm_mode}")
    
    for f in os.listdir(os.path.join(INPUT_FILES["patterns"] % tm_mode)):
        m = re.match(INPUT_FILES["pattern_file_regexp"], f)
        
        if m is not None:
            print(f"\n\nCleansing {m.string}...")
            clean_pattern_file(tm_mode, m.string)
            
# %% Main statement
if __name__ == "__main__":
    clean_patterns("earlyTM")
    clean_patterns("lateTM")
