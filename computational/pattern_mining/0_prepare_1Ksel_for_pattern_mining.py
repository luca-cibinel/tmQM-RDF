# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
# %% Header
import sys
sys.path.append(os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "interface"))

from tmQM_RDF_interface import TmQMRDFGraph
from tqdm import tqdm
import pandas as pd

TmQMRDFGraph.path_to_tmQM_RDF = os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev", "graphs")
TmQMRDFGraph.path_to_chem_info = os.path.join(ROOT_DIR, "data", "raw", "pubChem")

INPUT_FILES = {
        "selection": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "data", "%s", "1k_selection.csv")
    }

OUTPUT_FILES = {
        "serialized": os.path.join(ROOT_DIR, "computational", "pattern_mining", "intermediate", "local_data", "%s")
    }

# %% Utility functions
def serialize_tmQM_RDF_1Ksel(tm_mode):
    """
    Serializes the RDF graphs in the 1K selection(s) of tmQM-RDF into a format suitable for pattern mining.
    The serialised graphs only include the following predicates:
        
         rdflib.term.URIRef('http://www.w3.org/1999/02/22-rdf-syntax-ns#type'),
         rdflib.term.URIRef('resource://integreat/p5/atomic/atom/isAtom'),
         rdflib.term.URIRef('resource://integreat/p5/atomic/structure/b'),
         rdflib.term.URIRef('resource://integreat/p5/complex/TMC/hasLigand'),
         rdflib.term.URIRef('resource://integreat/p5/complex/TMC/hasMetalCentre'),
         rdflib.term.URIRef('resource://integreat/p5/ligand/bond/hasBindingAtom'),
         rdflib.term.URIRef('resource://integreat/p5/ligand/centre/hasAtom'),
         rdflib.term.URIRef('resource://integreat/p5/ligand/centre/isMetalCentre'),
         rdflib.term.URIRef('resource://integreat/p5/ligand/ligand/hasAtom'),
         rdflib.term.URIRef('resource://integreat/p5/ligand/ligand/isLigand'),
         rdflib.term.URIRef('resource://integreat/p5/ligand/structure/bLc'),
         rdflib.term.URIRef('resource://integreat/p5/ligand/structure/bLl')
         
    Arguments:
        - tm_mode: whether to process early TMs (mode = 'earlyTM') or the late ones (mode = 'lateTM')
    """   
    
    output_directory = OUTPUT_FILES["serialized"] % tm_mode
    
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    df = pd.read_csv(INPUT_FILES["selection"] % tm_mode)
    selection = [x for x in list(df.loc[df["partition"] == "train", "TMC"])]
    
    for tmcname in tqdm(selection):
        g = TmQMRDFGraph(tmcname)    
        
        g.skeleton().serialize(os.path.join(output_directory, tmcname + ".nt"), format = "ntriples")
        
# %% Main statement

if __name__ == "__main__":
    serialize_tmQM_RDF_1Ksel("earlyTM")
    serialize_tmQM_RDF_1Ksel("lateTM")