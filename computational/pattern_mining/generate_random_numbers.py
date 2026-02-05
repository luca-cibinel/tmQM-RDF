# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
     
# %% Header
import numpy as np

OUTPUT_FILES = {
        "random_numbers": os.path.join(ROOT_DIR, "computational", "pattern_mining", "intermediate", "random_numbers", "%s")
    }

RANDOM_SEED = 2568576931
SEED_SIZE = 10000

# %% Utility functions
def sample_random_numbers(session_keyword, random_state, seed_size):
    if not os.path.exists(OUTPUT_FILES["random_numbers"] % session_keyword):
        os.makedirs(OUTPUT_FILES["random_numbers"] % session_keyword)
    
    rng = np.random.RandomState(RANDOM_SEED)
    
    numbers = ["%.10f" % rng.rand() for _ in range(seed_size)]
    
    with open(os.path.join(OUTPUT_FILES["random_numbers"] % session_keyword, "random_numbers_1.txt"), "w") as f:
        f.write("\n".join(numbers))
    
    numbers = ["%.10f" % rng.rand() for _ in range(seed_size)]
    
    with open(os.path.join(OUTPUT_FILES["random_numbers"] % session_keyword, "random_numbers_2.txt"), "w") as f:
        f.write("\n".join(numbers))
        
# %% Main statement
if __name__ == "__main__":
    sample_random_numbers("earlyTM", RANDOM_SEED, SEED_SIZE)
    sample_random_numbers("lateTM", RANDOM_SEED*2, SEED_SIZE)