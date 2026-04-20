# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header
import os
import rpy2.robjects as robjects
import rpy2.rinterface_lib.embedded as rinterface

VOSH_SCRIPT = os.path.join(ROOT_DIR, "computational", "reconstruction", "aux", "3_compute_matches.R")
SAFETY_FILE = os.path.join(ROOT_DIR, "computational", "reconstruction", "case_studies", "DELETE_TO_ENABLE_STOP")

RECON_FOLDER = os.path.join(ROOT_DIR, "computational", "reconstruction", "case_studies", "intermediate", "reconstructions", "rdf")
RESULTS_DIRECTORY = os.path.join(ROOT_DIR, "computational", "reconstruction", "case_studies", "intermediate", "reconstructions", "rdf")
PATTERN_DIRECTORY = os.path.join(ROOT_DIR, "computational", "pattern_clustering", "results", "%s", "patterns")

N_CORES = 8

CASE_STUDIES = {
        "earlyTM": ["CAKRAH"],
        "lateTM": ["ILUYOD"]
    }

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
    
    for tm_mode, case_studies in CASE_STUDIES.items():
        print("*" * 10)
        print("*" * 5)
        print(f"PROCESSING TM MODE {tm_mode}")
        for case_tmc in case_studies:
            print(f"\tProcessing TMC {case_tmc}")
            
            local_recon_directory = os.path.join(RECON_FOLDER, tm_mode, case_tmc)
            local_results_directory = os.path.join(RESULTS_DIRECTORY, tm_mode, case_tmc)
            
            interrupt = (not os.path.exists(SAFETY_FILE))
            if not os.path.exists(PATTERN_DIRECTORY % f"{tm_mode}-latmod-s_10_12"):
                continue
            
            while not completed(local_results_directory) and not interrupt:
                try:
                    robjects.globalenv["TM.MODE"] = tm_mode
                    robjects.globalenv["RECON.FOLDER"] = local_recon_directory
                    robjects.globalenv["PATTERN.DATASET.TAG"] = f"{tm_mode}-latmod-s_10_12"
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