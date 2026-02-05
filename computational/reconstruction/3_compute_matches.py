# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header
import os
import rpy2.robjects as robjects
import rpy2.rinterface_lib.embedded as rinterface

TM_MODE = "lateTM"

VOSH_SCRIPT = os.path.join(ROOT_DIR, "computational", "reconstruction", "aux", "3_compute_matches.R")
SAFETY_FILE = os.path.join(ROOT_DIR, "computational", "reconstruction", "DELETE_TO_ENABLE_STOP")

RECON_FOLDER = os.path.join(ROOT_DIR, "computational", "reconstruction", "intermediate", "reconstructions", "rdf", TM_MODE)
RESULTS_DIRECTORY = os.path.join(ROOT_DIR, "computational", "reconstruction", "intermediate", "reconstructions", "rdf", TM_MODE)

N_CORES = 4

def completed(results_directory):
    
    for d in os.listdir(results_directory):
        if os.path.isdir(os.path.join(results_directory, d)):
            if "matches.csv" not in os.listdir(os.path.join(results_directory, d)):
                return False
    
    return True

# %% Main

if __name__ == "__main__":
    
    with open(SAFETY_FILE, "w") as f:
        f.write("Since the package rpy2 masks the behavior of Ctrl + C, if you want to stop the script 3_compute_matches.py, you have to first delete this file and then interrupt the script.")
    
    interrupt = False
    while not completed(RESULTS_DIRECTORY) and not interrupt:
        try:
            robjects.globalenv["TM.MODE"] = TM_MODE
            robjects.globalenv["RECON.FOLDER"] = RECON_FOLDER
            robjects.globalenv["PATTERN.DATASET.TAG"] = f"{TM_MODE}-latmod-s_10_12"
            robjects.globalenv["ROOT.DIR"] = ROOT_DIR
            robjects.globalenv["N.PROC"] = N_CORES
            robjects.r.source(VOSH_SCRIPT)
        except rinterface.RRuntimeError:
            interrupt = (not os.path.exists(SAFETY_FILE))
            
            print("")
            print("")
            print("_"*20)
            print("RRuntimeError caught!")
            robjects.r("""
                       print("   Stopping cluster...")
                       parallel::stopCluster(cl)
                       """)
            print("Restarting...")
            print("-"*20)
            print("")
            print("")