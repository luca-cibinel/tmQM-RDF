# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header
import sys
sys.path.append(os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "interface"))

from tqdm import tqdm
from pathlib import Path
from collections.abc import Sequence
from tmQM_RDF_interface import TmQMRDFInterface

import os
import copy
import numpy as np
import pandas as pd
import rdflib as rdf

# %% PrunedTMCDataset (classes)

class PrunedTMCGraph():
    
    path_to_ligand_blocks = None#"./ligand_blocks"
    
    def __init__(self, tmc_graph, tmc_name):
        if type(self).path_to_ligand_blocks is None:
            raise Exception("PrunedTMCGraph.path_to_ligand_bocks not set!")
        
        self.rdf = tmc_graph
        self.tmc_name = tmc_name
        
        self._identify_binding_sites()
        # TODO:
            # detect available valency of metal centre
    
    def _identify_binding_sites(self):
        qres = self.rdf.query("""
            SELECT ?tmc ?mcll ?mcal
            WHERE {
                ?tmc <resource://integreat/p5/complex/TMC/hasMetalCentre> ?mcll .
                ?mcll <resource://integreat/p5/ligand/centre/hasAtom> ?mcal .
            }
        """)
        
        qrow = next(iter(qres))
        
        self.binding_sites = (qrow.tmc, qrow.mcll, qrow.mcal)
        
    def attach(self, ligand_class):
        with open(os.path.join(type(self).path_to_ligand_blocks, f"{ligand_class}.nt"), "r") as f:
            raw_ligand_code = f.read()
            
        raw_ligand_code = raw_ligand_code.replace("@TMC@", str(self.binding_sites[0]))
        raw_ligand_code = raw_ligand_code.replace("$MC_LL$", str(self.binding_sites[1]))
        raw_ligand_code = raw_ligand_code.replace("!MC_AL!", str(self.binding_sites[2]))
        
        ligand = rdf.Graph()
        ligand.parse(data = raw_ligand_code)
        
        extended_tmc = copy.deepcopy(self.rdf)
        return extended_tmc + ligand
    
    def generate_all_reconstructions(self, mode, storage_directory):
        tmc_name = self.binding_sites[0].split("/")[-1]
        ligands = [f.replace(".nt", "") for f in os.listdir(os.path.join(self.path_to_ligand_blocks, mode)) if f.endswith(".nt")]
        
        for ligand in ligands:
            recon = self.attach(os.path.join(mode, ligand))
            
            recon.serialize(os.path.join(storage_directory, f"{tmc_name}__{ligand}.ttl"), format = "ttl")

class PrunableTMCGraph(TmQMRDFInterface):
    
    def __init__(self, tmc_name):
        super().__init__(tmc_name)
        
        self.identify_ligands()
        
    def identify_ligands(self):
        qres = self.rdf.query("""
                SELECT ?lig ?lig_class
                WHERE {
                    ?tmc <resource://integreat/p5/complex/TMC/hasLigand> ?lig .    
                    ?lig <resource://integreat/p5/ligand/ligand/isLigand> ?lig_class .
                }
            """)
        
        self.ligands = {i: (r.lig, r.lig_class.split("_")[-1]) for i, r in enumerate(qres)}
        
    def prune(self, lig_id):
        """
        This method uses a DFS-like algorithm to walk along the subgraph related to
        the specified ligand id and delete it. The modification is carried out on a copy of the rdf
        graph, which is later returned.
        """
        
        tmc = copy.deepcopy(self.skeleton())
        
        predicate_whitelist = [rdf.term.URIRef(p) for p in [ # Define edges on which propagation is allowed
                "resource://integreat/p5/ligand/ligand/isLigand",
                "resource://integreat/p5/ligand/structure/bLl",
                "resource://integreat/p5/ligand/ligand/hasAtom",
                "resource://integreat/p5/atomic/atom/isAtom",
                "resource://integreat/p5/atomic/structure/b",
                "resource://integreat/p5/ligand/bond/hasBindingAtom"
            ]]
        
        closure_whitelist = [rdf.term.URIRef(p) for p in [ # Define incoming edges which are to be considered part of the ligand
                "resource://integreat/p5/complex/TMC/hasLigand",
                "resource://integreat/p5/ligand/structure/bLc",
                "resource://integreat/p5/atomic/structure/b"
            ]]
        
        root = self.ligands[lig_id][0] # Define starting point
        
        to_visit = [root]
        visited = []
        
        while len(to_visit) > 0:
            visit = to_visit.pop()
            visited += [visit]
            
            for s, p, o in self.rdf.triples((None, None, visit)): # Remove incoming edges 
                                                # (ensures incoming edges coming from outside the ligand are removed)
                if p in closure_whitelist:
                    tmc.remove((s, p, o))           # FIXME: check if this still removes other isLigand/isAtom edges
            
            for s, p, o in self.rdf.triples((visit, None, None)):
                
                if p in predicate_whitelist:
                    if o not in visited:
                        to_visit += [o]
                
                tmc.remove((s, p, o)) # Remove all outgoing edges
        
        # for s, p, o in tmc.triples((None, rdf.term.URIRef("resource://integreat/p5/complex/TMC/hasLigand"), root)):
        #     print(s, p, o)
        
        # # print(len(tmc))
        # # tmc.remove((None, rdf.term.URIRef("resource://integreat/p5/complex/TMC/hasLigand"), root))        
        # # print(len(tmc))
        
        # for s, p, o in tmc.triples((None, rdf.term.URIRef("resource://integreat/p5/complex/TMC/hasLigand"), root)):
        #     print(s, p, o)
        
        return PrunedTMCGraph(tmc, self.tmc_name)

class PrunedDataLoader(Sequence):

    def __init__(self, files, path_to_files):
        self.files = files
        self.path_to_files = path_to_files
    
    def __getitem__(self, i):
        file = os.path.join(self.path_to_files, self.files[i])
        
        tmc_name, ground_truth = file.split("__")
        ground_truth = ground_truth.split(".")[0]
        
        temp = rdf.Graph()
        temp.parse(file)
        
        return PrunedTMCGraph(temp, tmc_name), ground_truth
    
    def __len__(self):
        return len(self.files)
    
    
class PrunedTMCDataset():
    
    path_to_tmQM_RDF = None # "../../../../data/tmQM-RDF/graphs"
    
    def __init__(self, 
                 tmQM_RDF_selection, 
                 path_to_pruned_tmcs, 
                 path_to_ligand_blocks, 
                 partition_to_prune = "test",
                 pruning = "random", 
                 random_state = 183620575292047
                 ):
        if type(self).path_to_tmQM_RDF is None:
            raise Exception("PrunedTMCDataset.path_to_tmQM_RDF not set!")
        
        self.tmQM_RDF_selection = pd.read_csv(tmQM_RDF_selection, index_col = 0)
        
        self.path_to_ligand_blocks = path_to_ligand_blocks
        
        self.path_to_pruned_tmcs = path_to_pruned_tmcs
        
        self.random_state = random_state
        self.rng = np.random.default_rng(self.random_state)
        
        if not os.path.exists(self.path_to_ligand_blocks):
            os.makedirs(self.path_to_ligand_blocks)
            self._extract_ligands_from_selection()
            
        if not os.path.exists(self.path_to_pruned_tmcs):
            os.makedirs(self.path_to_pruned_tmcs)
            if pruning == "random":
                self._prune_at_random(partition_to_prune)
            else:
                raise Exception(f"Pruning modality {pruning} not implemented!")
                
        self.data = {
            part: os.listdir(os.path.join(self.path_to_pruned_tmcs, part)) 
            for part in os.listdir(self.path_to_pruned_tmcs)
            if os.path.isdir(os.path.join(self.path_to_pruned_tmcs, part))
        }
    
    @staticmethod
    def _dfs_ligand(tmc, csd_code, ligand_blacklist = []):
        """
        This function navigates a TMC-RDF graph depth first and extracts the ligands.
        """
        
        clean = lambda s, p, o: (rdf.term.URIRef(t.replace(csd_code, "CANDIDATE")) for t in (s, p, o))
        
        ligands = {}
        
        propagation_whitelist = [rdf.term.URIRef(p) for p in [ # Define edges on which propagation is allowed
                "resource://integreat/p5/complex/TMC/hasLigand",
                "resource://integreat/p5/ligand/ligand/isLigand",
                "resource://integreat/p5/ligand/structure/bLl",
                "resource://integreat/p5/ligand/ligand/hasAtom",
                "resource://integreat/p5/atomic/atom/isAtom",
                "resource://integreat/p5/atomic/structure/b",
                "resource://integreat/p5/ligand/bond/hasBindingAtom"
            ]]
        
        closure_whitelist = [rdf.term.URIRef(p) for p in [ # Define incoming edges which are allowed to be included in ligand
                "resource://integreat/p5/complex/TMC/hasLigand",
                "resource://integreat/p5/ligand/structure/bLc",
                "resource://integreat/p5/atomic/structure/b"
            ]]
        
        root = rdf.term.URIRef(f"resource://integreat/p5/complex/TMC/{csd_code}") # Define starting point
        
        to_visit = [root]
        visited = []
        
        key = None
        while len(to_visit) > 0:
            visit = to_visit.pop()
            visited += [visit]
            
            if visit in ligands: # Starting visit of a ligand: prepare to store triples in correct graph
                key = visit
            
            for s, p, o in tmc.triples((None, None, visit)): # Add incoming edges
                                        # (ensures that nodes about the metal centre are added,
                                        # needed for identifiability of binding nodes)
                if p in closure_whitelist:
                    if key is not None:
                        ligands[key].add((clean(s, p, o)))
            
            for s, p, o in tmc.triples((visit, None, None)): # Add outgoing edges and propagate
            
                if p == propagation_whitelist[0]: # hasLigand
                    
                    if o.split("_")[-2] in ligand_blacklist:
                        continue
                
                    ligand_blacklist += [o.split("_")[-2]]
                    ligands[o] = rdf.Graph() # Prepare to store new ligand
                
                if p in propagation_whitelist:
                    if o not in visited:
                        to_visit += [o]
                    if key is not None:
                        ligands[key].add((clean(s, p, o)))
        
        for key in ligands: # Complete closure by adding the metal centre path
            tmc_uri, _, mc_ll_uri = next(iter(
                tmc.triples((root, rdf.term.URIRef("resource://integreat/p5/complex/TMC/hasMetalCentre"), None))
                ))
            
            # print(tmc_uri, mc_ll_uri)
            
            _, _, mc_al_uri = next(iter(
                tmc.triples((mc_ll_uri, rdf.term.URIRef("resource://integreat/p5/ligand/centre/hasAtom"), None))
                ))
            
            temp = ligands[key].serialize(format = "nt")
            
            temp = temp.replace(str(tmc_uri).replace(csd_code, "CANDIDATE"), "@TMC@")
            temp = temp.replace(str(mc_ll_uri).replace(csd_code, "CANDIDATE"), "$MC_LL$")
            temp = temp.replace(str(mc_al_uri).replace(csd_code, "CANDIDATE"), "!MC_AL!")
        
            ligands[key] = temp
        
        return ligands

    def _extract_ligands_from_selection(self):
        print("Extracting ligands...")
        
        partitions = set(self.tmQM_RDF_selection.loc[:,"partition"])
        ligs = {part: [] for part in partitions}
        
        for part in partitions:
            print(f"\tProcessing partition: {part}")
            
            local_storage_directory = os.path.join(self.path_to_ligand_blocks, part)
            if not os.path.exists(local_storage_directory):
                os.makedirs(local_storage_directory)
            
            tmcs = self.tmQM_RDF_selection.loc[self.tmQM_RDF_selection["partition"] == part, "TMC"]
            
            bl = []
            
            progress = tqdm(tmcs)
            for tmc in progress:
                g = rdf.Graph()
                g.parse(Path(os.path.join(type(self).path_to_tmQM_RDF, f"{tmc}.ttl")))
                
                ligs = type(self)._dfs_ligand(g, tmc, bl)
                
                for key, lig in ligs.items():
                    with open(os.path.join(local_storage_directory, f"{key.split('_')[-2]}.nt"), "w") as f:
                        f.write(lig)
                
                progress.set_description(f"Ligands found: {len(bl)}")
        
    def _prune_at_random(self, partition_to_prune, available_in_training = True):
        print("Pruning TMCs at random...")
        
        ligs = [l.split(".")[0] for l in os.listdir(os.path.join(self.path_to_ligand_blocks, "train"))]
        
        local_storage_directory = os.path.join(self.path_to_pruned_tmcs, partition_to_prune)
        if not os.path.exists(local_storage_directory):
            os.makedirs(local_storage_directory)
            
        tmcs = self.tmQM_RDF_selection.loc[self.tmQM_RDF_selection["partition"] == partition_to_prune, "TMC"]
        
        progress = tqdm(tmcs)
        reconstruction_not_possible = 0
        for tmc in progress:
            ptmc = PrunableTMCGraph(tmc)
            
            prunable_ligands = ptmc.ligands
            if available_in_training:
                prunable_ligands = {
                        i: x for i, x in prunable_ligands.items()
                        if x[1] in ligs
                    }
                
                if len(prunable_ligands) == 0:
                    prunable_ligands = ptmc.ligands
                    reconstruction_not_possible += 1

            index_to_prune = self.rng.choice(list(prunable_ligands.keys()))
            
            pruned = ptmc.prune(index_to_prune)
                
            fname = f"{tmc}__{ptmc.ligands[index_to_prune][1]}.ttl"
            pruned.rdf.serialize(os.path.join(local_storage_directory, fname), format = "ttl")
    
    def get_loader(self, mode):
        return PrunedDataLoader(self.data[mode], os.path.join(self.path_to_pruned_tmcs, mode))