"""
Step 3.i in the data preprocessing pipeline

Summarise all the information in tmQM, tmQMg and tmQMg-L into single .ttl files
"""

# %% Locate root dir
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header    
from collections import defaultdict
from tqdm import tqdm

import rdf_prefixes as pfx
import networkx as nx
import pandas as pd
import numpy as np
import time

INPUT_FILES = {
		"viable": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "viable_tmcs.txt"),
		"ligands_atoms": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQMg-L", "ligands_atoms_idx.csv"),
		"ligands_misc": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", "v60k", "ligands_misc_info.csv"),
		"ligands_fingerprints": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", "v60k", "ligands_fingerprints.csv"),
		"ligands_descriptors": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", "v60k", "ligands_descriptors.csv"),
		"periodic_table": os.path.join(ROOT_DIR, "data", "raw", "pubChem", "data", "PubChemElements_all.csv"),
		"tmQM_series_properties": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "misc", "tmQM_series_property_table.csv"),
        "tmQM_series_property_descriptions": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "misc", "tmQM_series_property_descriptions.csv"), # TODO
		"tmQM_graphs": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQM", "uNatQ_graphs"),
		"tmQMg_graphs": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg", "v74.637k", "uNatQ_graphs"),
		"tmQMg_L_graphs": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v2025dev", "tmQMg-L", "uNatQ_graphs")
	}

OUTPUT_FILES = {
		"tmQM_RDF_graphs": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev", "graphs"),
		"tmQM_RDFS_C": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev", "tmQM-RDFS-C.ttl"),
		"tmQM_RDFS_L": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev", "tmQM-RDFS-L.ttl"),
		"tmQM_RDFS_E": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev", "tmQM-RDFS-E.ttl"),
		"tmQM_RDFS_M": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev", "tmQM-RDFS-M.ttl"),
		"tmQM_RDFS_P": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v2025dev", "tmQM-RDFS-P.ttl")
	}

#%% Utility functions
def read_global_files():
    """
    Reads the global (not specific to a single TMC) files in the datasets tmQM, tmQMg, tmQMg-L.
    
    Returns:
        - viables: a list of viable TMCs, read from intermediate/viable_tmcs.txt
        - ligands_atlas: a pd.DataFrame read from intermediate/tmQMg-L/ligands_atoms_idx.csv
        - ligands_misc: a pd.DataFrame read from tmQMg-L/ligands_misc_info.csv
        - ligands_fingerprints: a pd.DataFrame read from tmQMg-L/ligands_fingerprints.csv
        - ligands_descriptors: a pd.DataFrame read from tmQMg-L/ligands_descriptors.csv
        - tmQM_series_property_table: a defaultdict(list) built from misc/tmQM_series_property_table.csv
        - periodic_table: a pd.DataFrame read from pubChem/PubChemElements_all.csv
    """
    
    """
    Extract list of viable TMCs
    """
    with open(INPUT_FILES["viable"], "r") as f:
        viables = f.read().split("\n")
    
    """
    Extract ligand composition information
    """
    ligands_atlas = pd.read_csv(INPUT_FILES["ligands_atoms"])
    
    """
    Read reference ligand data
    """
    ligands_misc = pd.read_csv(INPUT_FILES["ligands_misc"], sep = ";")
    ligands_fingerprints = pd.read_csv(INPUT_FILES["ligands_fingerprints"], sep = ";")
    ligands_descriptors = pd.read_csv(INPUT_FILES["ligands_descriptors"], sep = ";")
    
    """
    
    Read the tmQM series property table
    
        property_key    property_name   dataset     optimisation    singlepoint
        ...
        ...
    
    and convert it to a (defaultdict(list)) dictionary of the form
    
    {
        (property_key, reference_db): [
                [inp5("optimisation"), inp5ds(optimisation)],
                [inp5("singlepoint"), inp5ds(singlepoint)]
            ]
    }
    """
    tmQM_series_property_table = defaultdict(list)
    
    csv_property_table = pd.read_csv(INPUT_FILES["tmQM_series_properties"], sep = ";")
    
    for i in range(len(csv_property_table)):
        row = csv_property_table.iloc[i, :]
        
        tmQM_series_property_table[(row.iloc[0], row.iloc[2])] = [
                [pfx.nm("optimisation"), pfx.nm(row.iloc[3].replace("-", "_").replace("/", "__"))], 
                [pfx.nm("singlepoint"), pfx.nm(row.iloc[4].replace("-", "_").replace("/", "__"))]
            ]
    
    """
    Read PubChem periodic table
    """
    periodic_table = pd.read_csv(INPUT_FILES["periodic_table"], sep = ",")
    
    return viables, ligands_atlas, ligands_misc, ligands_fingerprints, ligands_descriptors, tmQM_series_property_table, periodic_table

#%% RDF Helper functions

def helper_xmltype(obj):
    """
    Returns the appropriate xmls type for an object.

    See https://www.w3.org/TR/2004/REC-xmlschema-2-20041028/#built-in-datatypes for a complete list of the available types.
    
    Arguments:
        - obj: a python object. Supported types: int, np.int64, float, np.float64, str, bool
    Returns:
        - a string representing the corresponding xml type
    """
    
    pytype = type(obj)
    
    if pytype in [int, np.int64]:
        return "integer"
    if pytype in [float, np.float64]:
        return "float"
    if pytype == str:
        return "string"
    if pytype == bool:
        return "boolean"

def helper_typed_literal(value):
    """
    Decorates a literal (value) with the corresponding xmls data type. See helper_xmltype
    
    Arguments:
        - value: the value to be decorated. Supported types: int, np.int64, float, np.float64, str, bool
    Returns:
        - a string representing the decorated literal
    """
    
    strv = str(value)
    
    if strv == "nan":
        strv = "NaN"
    
    return "\"" + strv + "\"" + "^^" + pfx.xmls(helper_xmltype(value))

def helper_write_blank_node_sentence(predicate, bn_sentences, indent = 1, raw = False):
    """
    Helper function: creates a sentence (with implicit subject) whose object is a blank node with several properties

    i.e.:
        
        predicate [
            bn_predicate_1 bn_object_1 ;
            bn_predicate_2 bn_object_2 ;
            ...
        ] ;

    The parameter indent defines how many \t characters have to be prefixed.

    If raw is False (default), bn_sentences is a list of lists of the format

    [['bn_predicate_1', 'bn_object_1'], ['bn_predicate_2', 'bn_object_2'], ...]

    If raw is True, bn_sentences is a list of fully formed triples (use with care!)

    Arguments:
        - predicate: a string representing the predicate URI
        - bn_sentences: either a list of lists or a lists of triples (see description)
        - indent: how many indents have to be prefixed (default = 1)
        - raw: whether the provided bn_Sentences are raw triples or not (default = False)
    Returns:
        - A list of triples
    """
    triples = ["\t"*indent + predicate + " ["]
    
    if raw:
        for bn_s in bn_sentences:
            triples += ["\t"*(indent + 1) + bn_s]
    else:
        for bn_s in bn_sentences:
            triples += ["\t"*(indent + 1) + " ".join(bn_s) + " ;"]
        
    triples += ["\t"*indent + "] ;"]
    
    return triples

def helper_write_prefixes():
    """
    Helper function: writes the header of the .ttl file with all the prefixes definitions
    
    Returns:
        - The triples defining the prefixes
    """
    triples = []
    
    path = lambda m: m()
    prefixes = sorted([getattr(pfx, m) for m in dir(pfx) if not m.startswith("__")], key = path)
    
    for p in prefixes:
        triples += ["@prefix " + p("") + " <" + p() + "> ."]
        
    return triples

def helper_write_property(property_type, reference_db, value, method = [], indent = 1):
    """
    Helper function: writes a property (assuming that the subject is already implicit):
        inp5:hasProperty [
            inp5:type property_type ;
            inp5:reference reference_db ;
            inp5:value "value"^^xmls:data_type ;
        ] ;
    
    Arguments:
        - property_type: the URI defining the property type
        - reference_db: the URI defining the reference db
        - value: the term defining the value
        - method: either [] or an entry of the tmQM_property_table defaultdict(list)
        - indent: how many indents have to be prefixed (default = 1)
    Returns:
        - a list of triples
    """
    
    return helper_write_blank_node_sentence(pfx.inp5("hasProperty"), [
            [pfx.rdf("type"), property_type],
            [pfx.inp5("reference"), reference_db],
            [pfx.rdf("value"), helper_typed_literal(value)]
        ] + method, indent)

def helper_reference_centre(metal_symbol):
    """
    Helper function that turns a transition metal chemical symbol into a reference centre URI using the format
        lgCr:MetalCentre_[metal_symbol]
        
    Arguments:
        - metal_symbol: the chemical symbol of the metal centre
    Returns:
        - a URI representing the corresponding centre class
    """
    
    return pfx.lgCr("MetalCentre_" + metal_symbol)

def helper_reference_ligand(ligand_id):
    """
    Helper function that turns a ligand id into a reference ligand URI using the format
        lgLr:Ligand_[ligand_id]
        
    Arguments:
        - ligand_id: the id of the ligand
    Returns:
        - a URI representing the corresponding ligand class
    """
    
    return pfx.lgLr("Ligand_" + ligand_id)

#%% RDF Builder functions

#%%% TMC
def build_atomic_level(tmc_id, tmQM_g, tmQMg_g):
    """
    Writes the triples related to the atomic level
    
    Arguments:
        - tmc_id: the CSD code of the TMC to build
        - tmQM_g: the tmQM .gml graph as a networkx.Graph object
        - tmQMg_g: the tmQMg .gml graph as a networkx.Graph object
    
    Returns:
        - the triples
        - an atom_atlas of the form
            {node_id: atom_uri, ...}
        - a reference_element_set containing the chemical elements used in this TMC
        - a property_atlas
            {prefix: [properties ...], ...}
    """
    
    triples = []
    
    triples += ["\n# ----- ATOMIC LEVEL -----"]
    
    property_atlas = defaultdict(list)
    
    """
    Preprocess atoms
    """
    counts = {}
    
    atom_atlas = dict(tmQMg_g.nodes(data = "node_label"))
    reference_element_set = set(atom_atlas.values())
    
    for at_id, at_lab in atom_atlas.items():
        atom_atlas[at_id] = pfx.tmA(tmc_id + "_" + at_lab + "_" + str(counts.get(at_lab, 0)))
        counts[at_lab] = counts.get(at_lab, 0) + 1
        
    """
    Preprocess bonds
    """
    bonds_atlas = []
    
    for e in list(tmQMg_g.edges()):
        s = e[0]
        t = e[1]
        
        uri_s = atom_atlas[s]
        uri_t = atom_atlas[t]
        
        uri_e = tmc_id + "_bond_" + uri_s.replace(pfx.tmA(tmc_id + "_"), "") + "_" + uri_t.replace(pfx.tmA(tmc_id + "_"), "")
        
        bonds_atlas += [(s, t, uri_e)]
        
    """
    Add atoms
    add properties
    and add atom related bond structures
    """
    
    tmQM_atomic_properties = dict(tmQM_g.nodes(data = True))
    tmQMg_atomic_properties = dict(tmQMg_g.nodes(data = True))
    
    for at_id, uri_atom in atom_atlas.items():
        triples += [uri_atom]
        
        """
        Add tmQMg node_label
        """
        triples += ["\t" + pfx.tmA("isAtom") + " " + pfx.tmAr(tmQMg_atomic_properties[at_id]["node_label"]) + " ;"]
        
        """
        Add tmQMg node id
        """
        triples += helper_write_property(pfx.tmAp("node_id"), pfx.dsG("tmQMg"), int(at_id))
        property_atlas[pfx.tmAp] += ["node_id"]
        
        """
        Add tmQMg related properties (no node_label, node_position)
        """
        for p, v in tmQMg_atomic_properties[at_id].items():
            
            if p in ["node_label", "node_position", "feature_atomic_number"]:
                continue
            
            p_name = p.replace("feature_", "").replace("target_tzvp_", "").replace("target_", "")
            triples += helper_write_property(
                pfx.tmAp(p_name),
                pfx.dsG("tmQMg"),
                v
            )
            
            property_atlas[pfx.tmAp] += [p_name]
           
        """
        Add tmQMg node_position
        """
        xyz_triples = []
        
        xyz_triples += [
            pfx.rdf("type") + " " + pfx.tmAp("node_position") + " ;",
            pfx.inp5("reference") + " " + pfx.dsG("tmQMg") + " ;"
        ]
        
        xyz_triples += helper_write_blank_node_sentence(pfx.rdf("value"), 
            [[pfx.rdf("type"), pfx.inp5("CartesianCoordinates3D")]] + [
                [
                    pfx.tmAp("xyz"[i]), 
                    helper_typed_literal(tmQMg_atomic_properties[at_id]["node_position"][i])
                ] for i in range(3)
            ], indent = 0)
        
        triples += helper_write_blank_node_sentence(
            pfx.inp5("hasProperty"), xyz_triples, raw = True
        )
        
        property_atlas[pfx.tmAp] += ["node_position"]
        
        """
        Add tmQM related properties (no node_label, node_position)
        """
        for p, v in tmQM_atomic_properties[at_id].items():
            
            if p in ["node_label", "node_position"]:
                continue
            
            p_name = p.replace("feature_", "").replace("target_tzvp_", "").replace("target_", "")
            triples += helper_write_property(
                pfx.tmAp(p_name),
                pfx.dsC("tmQM"),
                v
            )
            
            property_atlas[pfx.tmAp] += [p_name]
           
        """
        Add tmQM node_position
        """
        xyz_triples = []
        
        xyz_triples += [
            pfx.rdf("type") + " " + pfx.tmAp("node_position") + " ;",
            pfx.inp5("reference") + " " + pfx.dsC("tmQM") + " ;"
        ]
        
        xyz_triples += helper_write_blank_node_sentence(pfx.rdf("value"),
                [[pfx.rdf("type"), pfx.inp5("CartesianCoordinates3D")]] + [[
                    pfx.tmAp("xyz"[i]), 
                    helper_typed_literal(tmQM_atomic_properties[at_id]["node_position"][i])
                ] for i in range(3)
            ], indent = 0)
        
        triples += helper_write_blank_node_sentence(
            pfx.inp5("hasProperty"), xyz_triples, raw = True
        )
        
        property_atlas[pfx.tmAp] += ["node_position"]
        
        """
        Add atom related bond structure
        """
        for e in bonds_atlas:
            if at_id in e:
                triples += ["\t" + pfx.tmS("b") + " " + pfx.tmB(e[2]) + " ;"]
        
        """
        Close 'paragraph'
        """
        triples += ["."]
        
    """
    Add bonds
    and add proeprties
    """
    
    tmQM_bond_properties = list(tmQM_g.edges(data = True))
    tmQMg_bond_properties = list(tmQMg_g.edges(data = True))
    
    for e in bonds_atlas:
        triples += [pfx.tmB(e[2])]
        
        triples += ["\t" + pfx.rdf("type") + " " + pfx.tmB("AtomicBond") + " ;"]
        
        """
        Add tmQMg related properties (except for feature_nbo_type)
        """
        bond_properties = [f[2] for f in tmQMg_bond_properties if (f[0], f[1]) == (e[0], e[1])][0]
        
        for p, v in bond_properties.items():
            
            if p in ["feature_nbo_type"]:
                continue
            
            p_name = p.replace("feature_", "").replace("target_tzvp_", "").replace("target_", "")
            triples += helper_write_property(
                pfx.tmBp(p_name),
                pfx.dsG("tmQMg"),
                v
            )
            
            property_atlas[pfx.tmBp] += [p_name]
            
        """
        Add tmQMg feature_nbo_type
        """
        triples += helper_write_blank_node_sentence(
            pfx.inp5("hasProperty"), [
                pfx.rdf("type") + " " + pfx.tmBp("nbo_type") + " ;",
                pfx.inp5("reference") + " " + pfx.dsG("tmQMg") + " ;",
                pfx.rdf("value") + " " + pfx.tmBr("NBOType" + bond_properties["feature_nbo_type"]) + " ;",
            ], raw = True)
        
        property_atlas[pfx.tmBp] += ["nbo_type"]
        
        """
        Add tmQM related properties (except for feature_nbo_type)
        """
        bond_properties = [f[2] for f in tmQM_bond_properties if (f[0], f[1]) == (e[0], e[1])]
        
        if len(bond_properties) > 0:
            for p, v in bond_properties[0].items():
                
                p_name = p.replace("feature_", "").replace("target_tzvp_", "").replace("target_", "")
                triples += helper_write_property(
                    pfx.tmBp(p_name),
                    pfx.dsC("tmQM"),
                    v
                )
                
                property_atlas[pfx.tmBp] += [p_name]
                
        """
        Close 'paragraph'
        """    
        triples += ["."]
        
    return triples, atom_atlas, reference_element_set, property_atlas

def build_ligand_level(tmc_id, tmc_ligands, ligand_atlas, atom_atlas, metal_center_element):
    """
    Writes the triples related to the ligand level.

    Uses the tmc_id (CSD_code XXYYZZ) to build the URIs

    The ligands are built using tmc_ligands (pandas dataframe intermediate/tmQMg-L/ligands_atoms_idx.csv built by 1_i_encode_ligand_subgraphs.py) 
    with structure
    ligand_id           subgraph                         coordinating_atoms
    ligandId            XXYYZZ-subgraph-n              [[n], [m, k], ...]
    (NOTE: the indices n, m, k, ... get parsed to integers, but the node IDs are strings! 
    Probably a bad design choiche in encode_ligand_structure.py...)

    Builds the connection with the atomic level using an atom_atlas (returned by build_atomic_level) of the form

    {node_id: atom_uri, ...}

    and a ligand_atlas (pandas dataframe intermediate/tmQMg-L/uNatQ_graphs/XXYYZZ.csv built by 1_ii_encode_ligand_structure.py) of the form

    subgraph                    indices
    XXYYZZ_subgraph_n           n m ...
    ...

    Arguments:
        - tmc_id: the CSD code of the TMC to build
        - tmc_ligands: the pandads.Dataframe version of intermediate/tmQMg-L/ligands_atoms_idx.csv (built by 1_i_encode_ligand_subgraphs.py)
        - ligand_atlas: the pandads.Dataframe version of intermediate/tmQMg-L/uNatQ_graphs/XXYYZZ.csv (built by 1_ii_encode_ligand_structure.py)
        - atom_atlas: the atom atlas returned by build_atomic_level
        - metal_center_element: a string representing the chemical element of the metal centre

    Returns:
        - the triples
        - local_ligand_atlas, a dictionary of the form
            {uri_ligand: (ligand_id, ligand_subgraph, ligand_coordinating_atoms, ligand_bonds_uris), ...}
            (ligand_coordinating_atoms is a list of lists of coordinating atoms, ligand_bonds_uris is a list of uris of the bonds,
             corresponding to the atom sublists in ligand_coordinating_atoms)
        - uri_metal_centre, the URI of the metal centre
        - reference_centre, the dict {reference_centre_uri: metal_centre_element}
    """
    
    triples = []
    
    triples += ["\n# ----- LIGAND LEVEL -----"]
    
    """
    Preprocess ligands:
        - recover ligand id and create ligand uri
        - recover the associated subgraph (needed when connecting to the atomic level)
        - create bond uris
    """
    
    local_ligand_atlas = {}
    count = {}
    
    for i in range(tmc_ligands.shape[0]):
        """
        Extract ligand information
        """
        row = tmc_ligands.iloc[i]
        
        ligand_id = row.iloc[0]
        ligand_subgraph = row.iloc[1]
        ligand_coordinating = eval(row.iloc[2])
        
        """
        Create uri
        """
        uri_ligand = pfx.lgL(tmc_id + "_" + ligand_id + "_" + str(count.get(ligand_id, 0)))
        
        """
        Create bond uris
        Remember that the coordinating atoms list is a list of lists:
            [[a1, a2, ...], [b1, b2, ...], ...]
            Atoms that are in the same sublist (e.g. a1, a2, ...) are haptic with repsect to each other
                and belong to the same bond uri
            Groups of atoms that are in different sublists (e.g. a1, a2, ... and b1, b2, ...) are dentic with respect to each other
                and belong to different bond uris
        """
        uris_bonds = []
        for b in ligand_coordinating:
            uris_bonds += [pfx.lgB(tmc_id + "_bondL_" + ligand_id + "_" + str(count.get(ligand_id, 0))) + "_" + str(len(uris_bonds))]
        
        """
        Update local ligand information
        """
        count[ligand_id] = count.get(ligand_id, 0) + 1
        local_ligand_atlas[uri_ligand] = (ligand_id, ligand_subgraph, ligand_coordinating, uris_bonds)
    
    """
    Add metal centre
    add metal centre related bond structures
    and add connection to atomic level
    """
    uri_metal_centre = pfx.lgC(tmc_id + "_metalCentre_" + metal_center_element)
    reference_centre_uri = helper_reference_centre(metal_center_element)
    
    triples += [uri_metal_centre]
    
    triples += ["\t" + pfx.lgC("isMetalCentre") + " " + reference_centre_uri + " ;"]
    
    for lig_data in local_ligand_atlas.values():
        for b in lig_data[3]:
            triples += ["\t" + pfx.lgS("bLc") + " " + b + " ;"]
    
    uri_metal_centre_atom = [
        u for u in atom_atlas.values() if "_" + metal_center_element + "_" in u
    ][0]
    triples += ["\t" + pfx.lgC("hasAtom") + " " + uri_metal_centre_atom + " ;"]
    
    triples += ["."]
    
    """
    Add ligands
    add ligand related bond structures
    and add connections to atomic level (ligand to composing atoms)
    """
    for uri_lig, lig_data in local_ligand_atlas.items():
        triples += [uri_lig]
        
        """
        Add ligand
        """
        triples += ["\t" + pfx.lgL("isLigand") + " " + helper_reference_ligand(lig_data[0]) + " ;"]
        triples += ["\t" + pfx.lgL("subgraphName") + " " + helper_typed_literal(lig_data[1]) + " ;"]
        
        """
        Add ligand related bond structures
        """
        for b in lig_data[3]:
            triples += ["\t" + pfx.lgS("bLl") + " " + b + " ;"]
        
        """
        Add connections to atomic level (ligand to composing atoms)
        """
        subgraph = lig_data[1]
        
        atoms = ligand_atlas.loc[ligand_atlas["subgraph"] == subgraph].iloc[0, 1].strip().split(" ")
        
        for a in atoms:
            triples += ["\t" + pfx.lgL("hasAtom") + " " + atom_atlas[a] + " ;"]
        
        """
        Close the 'paragraph'
        """
        triples += ["."]
        
    """
    Add ligand bonds
    and add connections to atomic level (ligand-center bond to coordinating atoms)
    """
    for uri_lig, lig_data in local_ligand_atlas.items():
        """
        Remember: the list coordinating_atoms is a list of lists as explained in the ligands preprocessing step above 
        """
        coordinating_atoms = lig_data[2]
        uris_bonds = lig_data[3]
        
        for i in range(len(uris_bonds)):
            coord = coordinating_atoms[i]
            uri_bond = uris_bonds[i]
            
            """
            Add bond
            """
            triples += [uri_bond]
            
            triples += ["\t" + pfx.rdf("type") + " " + pfx.lgB("LigandBond") + " ;"]
            
            """
            Add connection to atomic level (ligand-center bond to coordinating atoms)
            """
            for a in coord:
                triples += ["\t" + pfx.lgB("hasBindingAtom") + " " + atom_atlas[str(a)] + " ;"]
        
            """
            Close the 'paragraph'
            """
            triples += ["."]
        
    return triples, local_ligand_atlas, uri_metal_centre, {reference_centre_uri: metal_center_element}

def build_complex_level(uri_tmc, tmQM_g, tmQMg_g, local_ligand_atlas, uri_metal_centre, tmQM_series_property_table):
    """
    Writes the triples related to the complex level.

    Builds the connection to the ligand level using
    - local_ligand_atlas, a dictionary of the form
        {uri_ligand: xxx, ...}
        (the only thing needed is the list of the keys)
    - uri_metal_centre, the URI of the metal centre
    - the tmQM_series_property_table defaultidct(list):
        {
            (property_key, reference_db): [
                    [inp5("optimisation"), inp5ds(optimisation)],
                    [inp5("singlepoint"), inp5ds(singlepoint)]
                ]
        }
        
    Arguments:
        - tmc_id: the CSD code of the TMC to build
        - tmQM_g: the tmQM .gml graph as a networkx.Graph object
        - tmQMg_g: the tmQMg .gml graph as a networkx.Graph object
        - local_ligand_atlas: the local ligand atlas produced by build_ligand_level
        - uri_metal_centre: the URI of the metal centre produced by build_ligand_level
        - tmQM_series_property_table: the deafaultdict(list) representing the tmQM property table
        
    Returns:
        - the triples
        - property_atlas: the dict {prefix: [properties ...], ...}
    """
    property_atlas = defaultdict(list)
    
    triples = []
    
    triples += ["\n# ----- COMPLEX LEVEL -----"]
    
    """
    Subject
    """
    triples += [uri_tmc]
    
    """
    tmQM properties
    """
    for p, v in tmQM_g.graph.items():
        triples += helper_write_property(
            pfx.cmTp(p),
            pfx.dsC("tmQM"),
            v,
            tmQM_series_property_table[(p, "tmQM")]
        )
        
        property_atlas[pfx.cmTp] += [p]
    
    """
    tmQMg properties (no meta_data)
    """
    for p, v in tmQMg_g.graph.items():
        parsed_p = p.replace("feature_", "").replace("target_tzvp_", "").replace("target_", "")
        
        if p == "meta_data":
            continue
        
        triples += helper_write_property(
            pfx.cmTp(parsed_p),
            pfx.dsG("tmQMg"),
            v,
            tmQM_series_property_table[(parsed_p, "tmQMg")]
        )
        
        property_atlas[pfx.cmTp] += [parsed_p]
    
    """
    tmQMg meta_data (simple properties)
    """
    meta_triples = []
    meta_triples += [pfx.rdf("type") + " " + pfx.cmTp("MetaDataDescription") + " ;"]
    for p, v in tmQMg_g.graph["meta_data"].items():
        if p in ["id", "metal_center_element", "element_counts"]:
            continue
        
        meta_triples += helper_write_property(
            pfx.cmTp(p),
            pfx.dsG("tmQMg"),
            v,
            indent = 0
        )
        
        property_atlas[pfx.cmTp] += [p]
        
    #triples += helper_write_blank_node_sentence(cmT("meta_data"), meta_triples, raw = True)
    
    """
    tmQMg meta_data (CSD_code)
    """
    meta_triples += [pfx.cmT("CSD_code") + " " + helper_typed_literal(tmQMg_g.graph["meta_data"]["id"]) + " ;"]
    
    """
    tmQMg meta_data (metal_center_element)
    """
    meta_triples += helper_write_blank_node_sentence(pfx.inp5("hasProperty"), [
            [pfx.rdf("type"), pfx.cmTp("metal_center_element")],
            [pfx.inp5("reference"), pfx.dsG("tmQMg")],
            [pfx.rdf("value"), pfx.tmAr(tmQMg_g.graph["meta_data"]["metal_center_element"])]
        ], indent = 0)
    
    property_atlas[pfx.cmTp] += ["metal_center_element"]
    
    """
    tmQMg meta_data (element_counts)
    """
    count_triples = []
    
    count_triples += [pfx.rdf("type") + " " + pfx.cmTp("element_counts") + " ;"]
    count_triples += [pfx.inp5("reference") + " " + pfx.dsG("tmQMg") + " ;"]
    
    bag_triples = []
    bag_triples += [pfx.rdf("type") + " " + pfx.inp5("CountDescription") + " ;"]
    for i, el_count in enumerate(tmQMg_g.graph["meta_data"]["element_counts"].items()):
        bag_triples += helper_write_blank_node_sentence(pfx.rdf(f"_{i + 1}"), [
                pfx.rdf("type") + " " + pfx.inp5("CountDescriptionItem") + " ;",
                pfx.inp5("countOf") + " " + pfx.tmAr(el_count[0]) + " ;",
                pfx.rdf("value") + " " + helper_typed_literal(el_count[1]) + " ;"
            ], raw = True, indent = 0)
    
    count_triples += helper_write_blank_node_sentence(pfx.rdf("value"), bag_triples, raw = True, indent = 0)
    
    meta_triples += helper_write_blank_node_sentence(pfx.inp5("hasProperty"), count_triples, raw = True, indent = 0)
    
    property_atlas[pfx.cmTp] += ["element_counts"]
    
    """
    Add tmQMg meta_data
    """
    
    triples += helper_write_blank_node_sentence(pfx.cmTp("meta_data"), meta_triples, raw = True)
    
    """
    Build connection to the ligand level
    """
    triples += ["\t" + pfx.cmT("hasMetalCentre") + " " + uri_metal_centre + " ;"]
    
    for uri_lig in local_ligand_atlas.keys():
        triples += ["\t" + pfx.cmT("hasLigand") + " " + uri_lig + " ;"]
    
    """
    Close the "paragraph"
    """
    triples += ["."]
    
    return triples, property_atlas

def build_tmc(tmc_id, ligands_atlas, tmQM_series_property_table):
    """
    Main bulder function for single RDF graphs

    Arguments:
        - tmc_id: the CSD code of the TMC to build
        - ligand_atlas: the pandads.Dataframe version of intermediate/tmQMg-L/uNatQ_graphs/XXYYZZ.csv (built by 1_ii_encode_ligand_structure.py)
        - tmQM_series_property_table: the deafaultdict(list) representing the tmQM property table
        
    Returns:
        - A list of turtle triples (each entry is a line)
        - The local ligand atlas returned by build_ligand_level (needed to build the ligand reference information)
        - reference_centre, the dict {reference_centre_uri: metal_centre_element}
        - reference_element_set, a set containing the chemical elements used in this TMC, third output of build_atomic_level
        - property_atlas: the dict {prefix: [properties ...], ...}
    """
    
    property_atlas = {}
    
    triples = []
    
    """
    Read .gml graphs
    """
    
    tmQM_g = nx.read_gml(os.path.join(INPUT_FILES["tmQM_graphs"], tmc_id + ".gml"))
    tmQMg_g = nx.read_gml(os.path.join(INPUT_FILES["tmQMg_graphs"], tmc_id + ".gml"))
    
    ligands = pd.read_csv(os.path.join(INPUT_FILES["tmQMg_L_graphs"], tmc_id + ".csv"), sep = ";")
    
    """
    --------------- PREFIXES ---------------
    """
    triples += helper_write_prefixes()
    
    """
    --------------- ATOMIC LEVEL ---------------
    """
    atomic_triples, atom_atlas, reference_element_set, a_property_atlas = build_atomic_level(tmc_id, tmQM_g, tmQMg_g)
    property_atlas.update(a_property_atlas)
    
    """
    --------------- LIGAND LEVEL ---------------
    """
    ligand_triples, local_ligand_atlas, uri_metal_centre, reference_centre = build_ligand_level(
        tmc_id, ligands, ligands_atlas, atom_atlas, tmQMg_g.graph["meta_data"]["metal_center_element"]
    )
    
    """
    --------------- COMPLEX LEVEL ---------------
    """
    uri_tmc = pfx.cmT(tmc_id)
    
    complex_triples , c_property_atlas = build_complex_level(uri_tmc, tmQM_g, tmQMg_g, local_ligand_atlas, uri_metal_centre, tmQM_series_property_table)
    property_atlas.update(c_property_atlas)
    
    """
    --------------- COMBINE LEVELS ---------------
    """
    
    triples += complex_triples
    triples += ligand_triples
    # triples += ligand_reference_triples
    triples += atomic_triples
    
    return triples, local_ligand_atlas, reference_centre, reference_element_set, property_atlas

def assemble_rdf(viables, ligands_atlas, tmQM_series_property_table):
    """
    Main builder function for all viable TMCs.
    While building the RDF graphs, it keeps track of the following items:
    - ligand classes, builds the reference_ligand_atlas
        {reference_ligand_uri: ligand_id}
    - metal centres, builds the reference_centres_atlas
        {reference_metal_centre_uri: metal_centre_element}
    - subgraphs, builds the subgraph_atlas
        {'XXYYZZ-subgraph-n': ligand_uri}
        (it maps the subgraph name to the ligand that occupies that subgraph, it will 
         be used by build_ligand_reference_ontology)
    - properties, builds the property_atlas
        {prefix: set([properties ...]), ...}
    
    Arguments:
        - viables: a list of viable TMCs, read from viable_tmcs.txt by read_global_files
        - ligands_atlas: a pd.DataFrame read from intermediate/tmQMg-L/ligands_atoms_idx.csv by read_global_files
        - ligands_misc: a pd.DataFrame read from tmQMg-L/ligands_misc_info.csv by read_global_files
        - ligands_fingerprints: a pd.DataFrame read from tmQMg-L/ligands_fingerprints.csv by read_global_files
        - ligands_descriptors: a pd.DataFrame read from tmQMg-L/ligands_descriptors.csv by read_global_files
        - tmQM_series_property_table: a defaultdict(list) built from tmQM_series_property_table.csv by read_global_files
        
    Returns:
        - reference_ligand_atlas
        - reference_centres_atlas
        - subgraph_atlas
        - reference_element_set, a set containing the chemical elements used in the TMCs
        - property_atlas
    """

    """
    Create directory tmQM-RDF if it does not exist
    """
    if not os.path.exists(OUTPUT_FILES["tmQM_RDF_graphs"]):
        os.makedirs(OUTPUT_FILES["tmQM_RDF_graphs"])

    """
    Initialise reference_ligand_atlas and subgraph_atlas
    (needed to build reference ligands ontology)
    """
    reference_ligand_atlas = {}
    reference_centres_atlas = {}
    subgraph_atlas = {}
    property_atlas = defaultdict(set)

    reference_element_set = set()

    start_time = time.time()

    print("Building RDF graphs...")
    for tmc in tqdm(viables):
        """
        Compute triples
        """
        triples, local_ligand_atlas, reference_centre, local_reference_element_set, tmc_property_atlas = build_tmc(
                                                            tmc,
                                                            ligands_atlas,
                                                            tmQM_series_property_table
                                                          )
        
        """
        Write to file
        """
        with open(os.path.join(OUTPUT_FILES["tmQM_RDF_graphs"], tmc + ".ttl"), "w") as f:
            f.write("\n".join(triples))
        
        """
        Process reference metal centre class
        """
        reference_centres_atlas.update(reference_centre)
        
        """
        Process local ligands into the reference atlas and subgraph atlas
        """
        for lig_uri, lig_data in local_ligand_atlas.items():
            if lig_data[0] not in reference_ligand_atlas.keys():
                reference_ligand_atlas[helper_reference_ligand(lig_data[0])] = lig_data[0]
            
            subgraph_atlas[lig_data[1]] = lig_uri
        
        """
        Process reference elements set
        """
        reference_element_set = reference_element_set.union(local_reference_element_set)
        
        """
        Process property_atlas
        """
        for prefix, p_list in tmc_property_atlas.items():
            property_atlas[prefix] = property_atlas[prefix].union(set(p_list))
        
    print(f"RDF graphs built in {time.time() - start_time} seconds.")
    
    return reference_ligand_atlas, reference_centres_atlas, subgraph_atlas, reference_element_set, property_atlas

#%%% Reference ontology
def build_classes_ontology():
    """
    Builds the ontology that contains the RDFS definition of the highest level classes used in tmQM-RDF.
    
    Returns:
        - the list of triples that constitute the ontology
    """
    
    print("Building classes ontology...")
    
    """
    Hardwired class definitions.
    
    Each class is represented as a list of 3: subClassOf: list, possibly empty; label: str; comment: str;
    """
    cm_classes = {
            pfx.cmT("TransitionMetalComplex"): [
                    [], 
                    "Transition Metal Complex (TMC)",
                    "This class represents Transition Metal Complexes."
                ],
            pfx.cmTp("MetaDataDescription"): [
                    [], 
                    "Description of meta data",
                    "This class represents an object whose purpose is to collect meta data about a TMC."
                ]
        }
    
    lg_classes = {
            # pfx.lg("LigandLevelObject"): [
            #         [], 
            #         "Ligand level object",
            #         "This class represents an object that exists at the ligand level (i.e. ligands and metal centres)."
            #     ],
            pfx.lgL("Ligand"): [
                    [], 
                    "Ligand",
                    "This class represents a ligand (intended as an abstract object that exists at the ligand level, and not as a molecule made of atoms chemically bonded to each other)."
                ],
            pfx.lgC("MetalCentre"): [
                    [], 
                    "Metal centre",
                    "This class represents a metal centre (intended as an abstract object that exists at the ligand level, and not as an atom)."
                ],
            pfx.lgB("LigandBond"): [
                    [], 
                    "Chemical bond (ligand level)",
                    "This class represents a chemical bond between a metal centre and a binding site within a ligand (intended as an abstract bond between two ligand-level objects, and not as a chemical bond between two atoms). Each LigandBond object can specify several atoms that participate in the bond and when it does, these atoms are considered to be haptic with respect to each other."
                ],
            pfx.lgLr("LigandClass"): [
                    [], 
                    "Ligand class",
                    "This class represents a ligand class, intended as a label that can be assigned to a ligand representation in order to determine its identity."
                ],
            pfx.lgLr("StableOccurrenceDescription"): [
                    [], 
                    "Stable occurrence",
                    "This class represents the description of the most stable occurrence of a ligand class, made of a subgraph name and a ligand instance (if available in the dataset)."
                ],
            pfx.lgLr("BindingAtomsSMILES"): [
                    [pfx.rdf("Bag")], 
                    "Binding atoms in SMILES string",
                    "This class represents the description of the binding atoms of a ligand in terms of indices within a SMILES string. It is a collection of BindingSiteSMILES objects."
                ],
            pfx.lgLr("BindingSiteSMILES"): [
                    [pfx.rdf("Bag")], 
                    "Binding site in SMILES string",
                    "This class represents the description of a group of haptic atoms in a ligand in terms of indices within a SMILES string."
                ],
            pfx.lgLrm("StructuralFeature"): [
                    [pfx.inp5("Countable")], 
                    "Structural feature",
                    "This class represents a structural feature of a ligand that can be counted."
                ],
            pfx.lgCr("MetalCentreClass"): [
                    [], 
                    "Metal centre class",
                    "This class represents a metal centre class, intended as a label that can be assigned to a metal centre representation in order to determine its identity."
                ],
            pfx.lgLro("LigandOccurrenceClass"): [
                    [], 
                    "Ligand occurrence type",
                    "This class represents the types of occurrences/representations of a ligand that can be encountered and used for computations (i.e., SMILES string, most stable occurrence, relaxed structure)."
                ]
        }
    
    tm_classes = {
            pfx.tmA("Atom"): [
                    [], 
                    "Atom",
                    "This class represents an atom."
                ],
            pfx.tmAr("Element"): [
                    [pfx.inp5("Countable")], 
                    "Chemical element",
                    "This class represents a chemical element, intended as a label for an atom."
                ],
            pfx.tmB("AtomicBond"): [
                    [], 
                    "Chemical bond (atomic level)",
                    "This class represents a chemical bond between two atoms."
                ],
            pfx.tmB("NBOType"): [
                    [], 
                    "NBO type",
                    "This class represents a possible NBO type, intended as a label that can be assigned to a chemical bond."
                ]
        }
    
    ds_classes = {
            pfx.ds("Dataset"): [
                    [], 
                    "Dataset",
                    "This class represents a dataset."
                ]
        }
    
    nm_classes = {
            pfx.nm("Optimisation"): [
                    [], 
                    "Optimisation",
                    "Optimizes molecular geometries by relaxing energy gradients"
                ],
            pfx.nm("Singlepoint"): [
                    [], 
                    "Singlepoint",
                    "Computes the energy and other electronic structure properties of optimized geometries"
                ] 
        }
    
    inp5_classes = {
            pfx.inp5("PropertyType"): [
                    [], 
                    "PropertyType",
                    "This class represents a possible property type, intended as one of the many types of properties that are reported in the tmQM series."
                ],
            pfx.inp5("NonSimplePropertyValue"): [
                    [], 
                    "Non-elementary Property value",
                    "This class represents an articulated property value that cannot be represented using a single URI/literal."
                ],
            pfx.inp5("CartesianCoordinates3D"): [
                    [pfx.inp5("NonSimplePropertyValue"), pfx.rdf("Seq")], 
                    "3D Cartesian coordinates",
                    "This class represents a set of coordinates in the Cartesian coordinate system of the three-dimensional Euclidean space R^3."
                ],
            pfx.inp5("CountDescription"): [
                    [pfx.inp5("NonSimplePropertyValue"), pfx.rdf("Bag")], 
                    "Count description",
                    "This class represents the description of a count property, intended as a bag of items, each one representing the count of a specific object."
                ],
            pfx.inp5("CountDescriptionItem"): [
                    [], 
                    "Count description item",
                    "This class represents an item in the description of a count property, representing the count of a specific object."
                ],
            pfx.inp5("ObservedProperty"): [
                    [], 
                    "Observed property",
                    "This class represents an observation of a property."
                ],
            pfx.inp5("Countable"): [
                    [], 
                    "Countable entity",
                    "This class represents a countable entity."
                ] 
        }
    
    sections = {
            "Complex level": cm_classes,
            "Ligand level": lg_classes,
            "Atomic level": tm_classes,
            "Datasets": ds_classes,
            "Numerical": nm_classes,
            "Properties": inp5_classes
        }
    
    """
    Start building the ontology
    """
    triples = []
    
    for header, classes in sections.items():
        triples += [f"\n# -- {header} --"]
        
        for class_name, class_body in classes.items():
            triples += [f"\n{class_name}"]
            
            triples += ["\t" + pfx.rdf("type") + " " + pfx.rdfs("Class") + " ;"]
            
            if len(class_body[0]) > 0:
                triples += ["\t" + pfx.rdfs("subClassOf") + " " + ", ".join(class_body[0]) + " ;"]
                
            triples += ["\t" + pfx.rdfs("label") + " " + helper_typed_literal(class_body[1]) + " ;"]
            
            triples += ["\t" + pfx.rdfs("comment") + " " + helper_typed_literal(class_body[2]) + " ;"]
            
            triples += ["."]
            
    return triples

def build_predicates_ontology():
    """
    Builds the ontology that contains the RDFS vocabulary of the predicates used in tmQM-RDF.
    
    Returns:
        - the list of triples that constitute the ontology
    """
    
    print("Building predicates ontology...")
    
    """
    Hardwired predicates definitions.
    
    Each predicate is represented as a list of 6: 
        type: str: list, possibly empty; 
        subPropertyOf: list, possibly empty; 
        domain: str; 
        range: str; 
        label: str; 
        comment: str;
    """
    cm_pred = {
            pfx.cmT("CSD_code"): [
                    [],
                    [],
                    pfx.cmT("TransitionMetalComplex"),
                    pfx.xmls("string"),
                    "CSD code",
                    "Indicates the Cambridge Structural Database code of a TMC."
                ],
            pfx.cmT("hasMetalCentre"): [
                    [],
                    [],
                    pfx.cmT("TransitionMetalComplex"),
                    pfx.lgC("MetalCentre"),
                    "Has metal centre",
                    "Indicates the metal centre (intended as ligand level object) of a TMC."
                ],
            pfx.cmT("hasLigand"): [
                    [],
                    [],
                    pfx.cmT("TransitionMetalComplex"),
                    pfx.lgL("Ligand"),
                    "Has ligand",
                    "Indicates a ligand (intended as ligand level object) possessed by TMC."
                ],
            pfx.cmTp("meta_data"): [
                    [],
                    [],
                    pfx.cmT("TransitionMetalComplex"),
                    pfx.cmTp("MetaDataDescription"),
                    "Meta data",
                    "Indicates the meta data associated with a given TMC."
                ]
        }
    
    lg_pred = {
            pfx.lgC("isMetalCentre"): [
                    [],
                    [pfx.rdf("type")],
                    pfx.lgC("MetalCentre"),
                    pfx.lgCr("MetalCentreClass"),
                    "Is metal centre",
                    "Specifies the class of a metal centre (intended as ligand level object)."
                ],
            pfx.lgC("hasAtom"): [
                    [],
                    [],
                    pfx.lgC("MetalCentre"),
                    pfx.tmA("Atom"),
                    "Has atom",
                    "Indicates the atom possessed by (i.e., that is) the metal centre (intended as ligand level object)	."
                ],
            pfx.lgL("isLigand"): [
                    [],
                    [pfx.rdf("type")],
                    pfx.lgL("Ligand"),
                    pfx.lgLr("LigandClass"),
                    "Is ligand",
                    "Specifies the class of a ligand (intended as ligand level object)."
                ],
            pfx.lgL("subgraphName"): [
                    [],
                    [],
                    pfx.lgL("Ligand"),
                    pfx.xmls("string"),
                    "Subgraph name",
                    "Specifies the name of the subgraph that corresponds to a ligand instance (intended as ligand level object)."
                ],
            pfx.lgL("hasAtom"): [
                    [],
                    [],
                    pfx.lgL("Ligand"),
                    pfx.tmA("Atom"),
                    "Has atom",
                    "Indicates the atom possessed by (i.e., that is part of) a ligand (intended as ligand level object)	."
                ],
            pfx.lgLr("subgraph"): [
                    [],
                    [],
                    pfx.lgLr("StableOccurrenceDescription"),
                    pfx.xmls("string"),
                    "Stable occurrence subgraph",
                    "Indicates the subgraph name of the most stable occurrence of a ligand class."
                ],
            pfx.lgLr("instance"): [
                    [],
                    [],
                    pfx.lgLr("StableOccurrenceDescription"),
                    pfx.lgL("Ligand"),
                    "Stable occurrence instance",
                    "Indicates the instance which is the most stable occurrence of a ligand class."
                ],
            pfx.lgLrp("basedOn"): [
                    [],
                    [],
                    pfx.inp5("ObservedProperty"),
                    pfx.lgLro("LigandOccurrenceClass"),
                    "Based on occurrence",
                    "Specifies the ligand representation on which a property is based on (SMILES, most stable occurrence or relaxed structure)."
                ],
            pfx.lgS("bLc"): [
                    [],
                    [],
                    pfx.lgC("MetalCentre"),
                    pfx.lgB("LigandBond"),
                    "Chemical bond (ligand level)",
                    "Indicates that a ligand level object participates in a chemical bond represented by the object	."
                ],
            pfx.lgS("bLl"): [
                    [],
                    [],
                    pfx.lgL("Ligand"),
                    pfx.lgB("LigandBond"),
                    "Chemical bond (ligand level)",
                    "Indicates that a ligand level object participates in a chemical bond represented by the object	."
                ],
            pfx.lgB("hasBindingAtom"): [
                    [],
                    [],
                    pfx.lgB("LigandBond"),
                    pfx.tmA("Atom"),
                    "Has binding atom",
                    "Indicates an atom that acts as binding atom in the chemical bond between ligand level objects represented by the subject."
                ]
        }
    
    tm_pred = {
            pfx.tmA("isAtom"): [
                    [],
                    [pfx.rdf("type")],
                    pfx.tmA("Atom"),
                    pfx.tmAr("AtomClass"),
                    "Is atom",
                    "Specifies the class (element) of an atom."
                ],
            pfx.tmAp("x"): [
                    [pfx.rdfs("ContainerMembershipProperty")],
                    [pfx.rdf("_1")],
                    pfx.inp5("CartesianCoordinates3D"),
                    pfx.xmls("float"),
                    "X coordinate",
                    "Specifies the cartesian x coordinate of an atom."
                ],
            pfx.tmAp("y"): [
                    [pfx.rdfs("ContainerMembershipProperty")],
                    [pfx.rdf("_2")],
                    pfx.inp5("CartesianCoordinates3D"),
                    pfx.xmls("float"),
                    "Y coordinate",
                    "Specifies the cartesian y coordinate of an atom."
                ],
            pfx.tmAp("z"): [
                    [pfx.rdfs("ContainerMembershipProperty")],
                    [pfx.rdf("_3")],
                    pfx.inp5("CartesianCoordinates3D"),
                    pfx.xmls("float"),
                    "Z coordinate",
                    "Specifies the cartesian z coordinate of an atom."
                ],
            pfx.tmS("b"): [
                    [],
                    [],
                    pfx.tmA("Atom"),
                    pfx.tmB("AtomicBond"),
                    "Atomic bond",
                    "Indicates that an atom participates in a chemical bond represented by the object."
                ]
        }
    
    inp5_pred = {
            pfx.inp5("hasProperty"): [
                    [],
                    [],
                    pfx.rdfs("Resource"),
                    pfx.inp5("ObservedProperty"),
                    "Has property",
                    "Asserts that an object possesses a property."
                ],
            pfx.inp5("countOf"): [
                    [],
                    [],
                    pfx.inp5("CountDescriptionItem"),
                    pfx.inp5("Countable"),
                    "Is count of",
                    "In count properties (e.g. element counts) specifies the object being counted."
                ],
            pfx.inp5("reference"): [
                    [],
                    [],
                    pfx.inp5("Property"),
                    pfx.ds("Dataset"),
                    "Reference dataset",
                    "Indicates that the an atom participates in the chemical bond represented by the object."
                ]
        }
    
    nm_pred = {
            pfx.nm("optimisation"): [
                    [],
                    [],
                    pfx.inp5("ObservedProperty"),
                    pfx.nm("Optimisation"),
                    "Optimisation",
                    "Indicates the optimisation method used"
                ],
            pfx.nm("singlepoint"): [
                    [],
                    [],
                    pfx.inp5("ObservedProperty"),
                    pfx.nm("Singlepoint"),
                    "Singlepoint",
                    "Indicates the method used to compute this property given an optimised geometry"
                ]
        }
    
    sections = {
            "Complex level": cm_pred,
            "Ligand level": lg_pred,
            "Atomic level": tm_pred,
            "Numerical": nm_pred,
            "Properties": inp5_pred
        }
    
    """
    Start building the ontology
    """
    triples = []
    
    for header, predicates in sections.items():
        triples += [f"\n# -- {header} --"]
        
        for pred_name, pred_body in predicates.items():
            triples += [f"\n{pred_name}"]
            
            if len(pred_body[0]) > 0:
                triples += ["\t" + pfx.rdf("type") + " " + ", ".join(pred_body[0]) + " ;"]
                
            if len(pred_body[1]) > 0:
                triples += ["\t" + pfx.rdfs("subPropertyOf") + " " + ", ".join(pred_body[1]) + " ;"]
            
            triples += ["\t" + pfx.rdfs("domain") + " " + pred_body[2] + " ;"]
        
            triples += ["\t" + pfx.rdfs("range") + " " + pred_body[3] + " ;"]
            
            triples += ["\t" + pfx.rdfs("label") + " " + helper_typed_literal(pred_body[4]) + " ;"]
            
            triples += ["\t" + pfx.rdfs("comment") + " " + helper_typed_literal(pred_body[5]) + " ;"]
            
            triples += ["."]
            
    return triples

def build_ligand_reference_ontology(ligands_misc, ligands_fingerprints, ligands_descriptors, reference_ligand_atlas, reference_centres_atlas, subgraph_atlas, tmQM_series_property_table, periodic_table):
    """
    Builds an ontology with all the reference information of the metal centres and ligands used in tmQM-RDF

    It uses the file tmQMg-L/ligands_misc_info.csv to recover the SMILES representation.

    It uses the file ligands_fingerprints.csv:
        - It adds charge, n_atoms, n_metal_bound, n_dentic_bound, n_haptic_bound as simple properties
        - It adds is_alternative_charge as a simple property after conversion to boolean
        - It treates the three following blocks of atom count (atoms, dentic atoms, haptic atoms; 14 entries each) as an atom count property
            like the one used in the complex level (only non zero entries are added)

    It uses the file ligands_descriptors:
        - Stable occurrence name: points to the URI of the subgraph with the most stable occurrence
        - logp: added as simple property
        - n_aliphatic_rings, n_aromatic_rings, n_saturated_rings, n_rotatable_bonds: treated as a count property like the one used in the
            complex level (zero entries are added anyway)
        - The remaining properties are treated this way: each property is derived from a specific representation, hence they are added as properties
            which have a reference db (as usual) but also a reference representation.
            The references are:
            * No prefix: SMILES
            * L*_: most stable occurrence
            * L_free_: relaxed structure

    It uses the reference ligand atlas:
        {reference_ligand_uri: ligand_id}
        
    It uses the reference centres atlas:
        {reference_metal_centre_uri: metal_centre_element}
    to define the metal centre class
    
    It uses the subgraph atlas:
        {'XXYYZZ-subgraph-n': ligand_uri}
        (it maps the subgraph name to the ligand that occupies that subgraph)

    It uses the tmQM_series_property_table defaultidct(list):
        
        {
            (property_key, reference_db): [
                    [inp5("optimisation"), inp5ds(optimisation)],
                    [inp5("singlepoint"), inp5ds(singlepoint)]
                ]
        }
    
    Arguments:
        - ligands_misc: a pd.DataFrame read from tmQMg-L/ligands_misc_info.csv by read_global_files
        - ligands_fingerprints: a pd.DataFrame read from tmQMg-L/ligands_fingerprints.csv by read_global_files
        - ligands_descriptors: a pd.DataFrame read from tmQMg-L/ligands_descriptors.csv by read_global_files
        - reference_ligand_atlas: first output of assemble_rdf
        - reference_centres_atlas: second output of assemble_rdf
        - subgraph_atlas: third output of assemble_rdf
        - tmQM_series_property_table: a defaultdict(list) built from tmQM_series_property_table.csv by read_global_files
        - periodic_table: a pd.DataFrame read from misc/PubChemElements_all.csv by read_global_files
        
    Returns:
        - the list of triples that constitute the ontology
        - the property_atlas {prefix: set([properties]), ...} of the properties encountered during the construction
    """
    
    property_atlas = defaultdict(list)
    
    start_time = time.time()
    print("Building ligand ontology...")
    
    triples = []
    
    """
    Define metal centre classes
    """
    atnum = lambda uri_el: int(periodic_table.loc[periodic_table["Symbol"] == uri_el[1], "AtomicNumber"].iloc[0])
    
    centre_label = "Metal centre: %s (%s)"
    centre_comment = "This class describes the metal centre %s (%s), intended as an abstract object that exists at the ligand level, and not as an atom."
    
    triples += ["\n# ---- CENTRES ----"]
    print("\tDefining metal centres...")
    for uri_centre, centre_el in tqdm(sorted(reference_centres_atlas.items(), key = atnum)):
        triples += [""]
        triples += [uri_centre]
        
        """
        Add superclass
        """
        triples += ["\t" + pfx.rdfs("subClassOf") + " " + pfx.lgC("MetalCentre") + " ;"]
        
        """
        Add type
        """
        triples += ["\t" + pfx.rdf("type") + " " + pfx.lgCr("MetalCentreClass") + " ;"]
        
        """
        Add human-like description and label
        """
        centre_name = periodic_table.loc[periodic_table["Symbol"] == centre_el, "Name"].iloc[0]
        
        triples += ["\t" + pfx.rdfs("label") + " " + helper_typed_literal( centre_label % (centre_el, centre_name) ) + " ;"]
        triples += ["\t" + pfx.rdfs("comment") + " " + helper_typed_literal( centre_comment % (centre_el, centre_name) ) + " ;"]
        triples += "."
        
    """
    Define ligand classes
    """
    ligsize = lambda uri_id: int(ligands_fingerprints.loc[ligands_fingerprints["name"] == uri_id[1], "n_atoms"].iloc[0])
    
    ligand_label = "Ligand %s"
    ligand_comment = "This class describes the ligand %s, with SMILES string %s. The ligand is intended as an abstract object that exists at the ligand level, and not as a molecule made of atoms chemically bonded with each other."
    
    triples += ["\n# ---- LIGANDS ----"]
    print("\tDefining ligands...")
    for uri_lig, lig_id in tqdm(sorted(reference_ligand_atlas.items(), key = ligsize)):
        
        triples += [f"\n## -- {lig_id} --"]
        
        triples += [uri_lig]
        
        """
        Define ligand class, as subclass of lgL:Ligand, and type
        """
        triples += ["\t" + pfx.rdfs("subClassOf") + " " + pfx.lgL("Ligand") + " ;"]
        triples += ["\t" + pfx.rdf("type") + " " + pfx.lgLr("LigandClass") + " ;"]
        
        smiles = ligands_misc.loc[ligands_misc["name"] == lig_id, "smiles"].iloc[0]
        smiles = smiles.replace("\\", "\\\\")
        
        triples += ["\t" + pfx.rdfs("label") + " " + helper_typed_literal( ligand_label % (lig_id) ) + " ;"]
        triples += ["\t" + pfx.rdfs("comment") + " " + helper_typed_literal( ligand_comment % (lig_id, smiles) ) + " ;"]
        
        """
        Add SMILES
        """
        triples += helper_write_property(pfx.lgLrp("SMILES"), pfx.dsL("tmQMg-L"), smiles)
        property_atlas[pfx.lgLrp] += ["SMILES"]
        
        """
        Add binding atoms (as indices in SMILES string)
        The property has one "value" edge for each group of haptic atoms, pointing to a blank node.
        This blank node has a "value" edge for each atom inside the group.
        
        binding_atoms is a list of lists: each element of the list is a group of haptic atoms
        represented as a list of atoms (indices in the smiles string)
        """
        binding_atoms = eval(ligands_misc.loc[ligands_misc["name"] == lig_id, "smiles_metal_bond_node_idx_groups"].iloc[0])
        
        haptic_group_triples = [pfx.rdf("type") + " " + pfx.lgLr("BindingAtomsSMILES") + " ;"]
        for i, group in enumerate(binding_atoms):
            haptic_group_triples += helper_write_blank_node_sentence(
                    pfx.rdf(f"_{i + 1}"),
                    [[pfx.rdf("type"), pfx.lgLr("BindingSiteSMILES")]] +
                    [[pfx.rdf(f"_{j + 1}"), helper_typed_literal(atom)] for j, atom in enumerate(group)],
                    indent = 0
                )
        
        value_triples = []
        value_triples += [pfx.rdf("type") + " " + pfx.lgLrp("smiles_metal_bond_node_idx_groups") + " ;"]
        value_triples += [pfx.inp5("reference") + " " + pfx.dsL("tmQMg-L") + " ;"]
        value_triples += helper_write_blank_node_sentence(pfx.rdf("value"), haptic_group_triples, raw = True, indent = 0)
        
        triples += helper_write_blank_node_sentence(pfx.inp5("hasProperty"), value_triples, raw = True)
        
        property_atlas[pfx.lgLrp] += ["smiles_metal_bond_node_idx_groups"]
        
        """
        Add ligand fingerprints:
            - basic properties: charge, n_atoms, n_metal_bound, n_dentic_bound, n_haptic_bound
        """
        
        for p in ['charge', 'n_atoms', 'n_metal_bound', 'n_dentic_bound', 'n_haptic_bound']:
            v = ligands_fingerprints.loc[ligands_fingerprints["name"] == lig_id, p].iloc[0]
            
            triples += helper_write_property(pfx.lgLrp(p), pfx.dsL("tmQMg-L"), v)
            property_atlas[pfx.lgLrp] += [p]
            
        """
        Add ligand fingerprints:
            - is_alternative_charge (turn to bool)
        """
        is_alt_charge = bool(ligands_fingerprints.loc[ligands_fingerprints["name"] == lig_id, "is_alternative_charge"].iloc[0])
        triples += helper_write_property(pfx.lgLrp("is_alternative_charge"), pfx.dsL("tmQMg-L"), is_alt_charge)
        property_atlas[pfx.lgLrp] += ["is_alternative_charge"]
        
        """
        Add ligand fingerprints:
            - atom count
            - dentic atom count
            - haptic atom count
        """
        atoms_to_count = ligands_fingerprints.columns[7:(7+14)]
        
        for count_type in ["", "dentic_", "haptic_"]:
            """
            Recover the corresponding count entries from the .csv file
            """
            atoms_count = dict(
                ligands_fingerprints.loc[
                    ligands_fingerprints["name"] == lig_id, 
                    [count_type + el for el in atoms_to_count]
                ].iloc[0])
            
            """
            If there are no atoms to report, skip
            """
            if sum(atoms_count.values()) == 0:
                continue
            
            """
            Start adding the count property
            """
            count_triples = []
            
            count_triples += [pfx.rdf("type") + " " + pfx.lgLrp(count_type + "element_counts") + " ;"]
            count_triples += [pfx.inp5("reference") + " " + pfx.dsL("tmQMg-L") + " ;"]
            
            """
            Add all counts
            """
            bag_triples = []
            bag_triples += [pfx.rdf("type") + " " + pfx.inp5("CountDescription") + " ;"]
            for i, el_count in enumerate(atoms_count.items()):
                bag_triples += helper_write_blank_node_sentence(pfx.rdf(f"_{i + 1}"), [
                        pfx.rdf("type") + " " + pfx.inp5("CountDescriptionItem") + " ;",
                        pfx.inp5("countOf") + " " + pfx.tmAr(el_count[0].replace(count_type, "")) + " ;",
                        pfx.rdf("value") + " " + helper_typed_literal(el_count[1]) + " ;"
                    ], raw = True, indent = 0)
            
            count_triples += helper_write_blank_node_sentence(pfx.rdf("value"), bag_triples, raw = True, indent = 0)
            
            """
            Add the count property
            """
            triples += helper_write_blank_node_sentence(pfx.inp5("hasProperty"), count_triples, raw = True)
            property_atlas[pfx.lgLrp] += [count_type + "element_counts"]
            
        """
        Add ligand descriptors:
            - stable occurence
        """
        stable_occurrence_subgraph = ligands_descriptors.loc[ligands_descriptors["name"] == lig_id, "stable_occurrence_name"]
        stable_occurrence_subgraph = stable_occurrence_subgraph.iloc[0]
        
        stable_occurrence_uri = subgraph_atlas.get(stable_occurrence_subgraph, None)
        
        stable_occurrence_triples = []
        stable_occurrence_triples += [pfx.rdf("type") + " " + pfx.lgLrp("stable_occurrence") + " ;"]
        stable_occurrence_triples += [pfx.inp5("reference") + " " + pfx.dsL("tmQMg-L")  + " ;"]
        
        instance_triple = [[pfx.lgLr("instance"), stable_occurrence_uri]] if stable_occurrence_uri is not None else []
        stable_occurrence_triples += helper_write_blank_node_sentence(
                                        pfx.rdf("value"), 
                                        [
                                            [pfx.rdf("type"), pfx.lgLr("StableOccurrenceDescription")],
                                            [pfx.lgLr("subgraph"), helper_typed_literal(stable_occurrence_subgraph)]
                                        ] + instance_triple,
                                        indent = 0)
        
        triples += helper_write_blank_node_sentence(
                        pfx.inp5("hasProperty"), 
                        stable_occurrence_triples,
                        raw = True
                   )
        
        property_atlas[pfx.lgLrp] += ["stable_occurrence"]
                        
        """
        Add lignad descriptors:
            - logp
        """
        v = ligands_descriptors.loc[ligands_descriptors["name"] == lig_id, "logp"].iloc[0]
        triples += helper_write_blank_node_sentence(pfx.inp5("hasProperty"), [
                        [pfx.rdf("type"), pfx.lgLrp("logp")],
                        [pfx.inp5("reference"), pfx.dsL("tmQMg-L")],
                        [pfx.lgLrp("basedOn"), pfx.lgLro("SMILES")],
                        [pfx.rdf("value"), helper_typed_literal(v)]
                    ])
        
        property_atlas[pfx.lgLrp] += ["logp"]
        
        """
        Add ligand descriptors:
            - n_aliphatic_rings, n_aromatic_rings, n_saturated_rings, n_rotatable_bonds
        """
        to_class_name = lambda s: "".join(w.capitalize() for w in s[:-1].replace("n_", "").split("_"))
        
        struct_count = dict(
            ligands_descriptors.loc[
                ligands_descriptors["name"] == lig_id, 
                ['n_aliphatic_rings', 'n_aromatic_rings', 'n_saturated_rings', 'n_rotatable_bonds']
            ].iloc[0])
        
        count_triples = []
        
        count_triples += [pfx.rdf("type") + " " + pfx.lgLrp("structure_counts") + " ;"]
        count_triples += [pfx.inp5("reference") + " " + pfx.dsL("tmQMg-L") + " ;"]
        
        bag_triples = []
        bag_triples += [pfx.rdf("type") + " " + pfx.inp5("CountDescription") + " ;"]
        for i, str_count in enumerate(struct_count.items()):
            bag_triples += helper_write_blank_node_sentence(pfx.rdf(f"_{i + 1}"), [
                    pfx.rdf("type") + " " + pfx.inp5("CountDescriptionItem") + " ;",
                    pfx.inp5("countOf") + " " + pfx.lgLrm(to_class_name(str_count[0])) + " ;",
                    pfx.rdf("value") + " " + helper_typed_literal(str_count[1]) + " ;"
                ], raw = True, indent = 0)
        
        count_triples += helper_write_blank_node_sentence(pfx.rdf("value"), bag_triples, raw = True, indent = 0)
        
        triples += helper_write_blank_node_sentence(pfx.inp5("hasProperty"), count_triples, raw = True)
        property_atlas[pfx.lgLrp] += ["structure_counts"]
        
        """
        Add ligands_descriptors:
            - SMILES BASED PROPERTIES: 'exact_cone_angle', 'buried_volume', 'solid_angle', 'solid_cone_angle',
            'G_parameter', 'sasa_area_stable', 'sasa_volume_stable', 'sasa_area_free', 'sasa_volume_free'
        """
        
        smiles_properties = [
            'exact_cone_angle', 'buried_volume', 'solid_angle', 'solid_cone_angle',
            'G_parameter', 'sasa_area_stable', 'sasa_volume_stable',
            'sasa_area_free', 'sasa_volume_free'
        ]
        
        for p in smiles_properties:
            v = ligands_descriptors.loc[ligands_descriptors["name"] == lig_id, p].iloc[0]
            
            triples += helper_write_blank_node_sentence(pfx.inp5("hasProperty"), [
                            [pfx.rdf("type"), pfx.lgLrp(p)],
                            [pfx.inp5("reference"), pfx.dsL("tmQMg-L")],
                            [pfx.lgLrp("basedOn"), pfx.lgLro("SMILES")],
                            [pfx.rdf("value"), helper_typed_literal(v)]
                        ])
            property_atlas[pfx.lgLrp] += [p]
            
        """
        Add ligands_descriptors:
            - most stable occurrence(L*)/relaxed(L_free) structure BASED PROPERTIES
        """
        
        ref_occurences = {
            "L*-": pfx.lgLro("most_stable"),
            "L_free-": pfx.lgLro("relaxed_structure")
        }
        
        for ref_occ, uri_occ in ref_occurences.items():
            properties = [x.replace(ref_occ, "") for x in ligands_descriptors.columns if x.startswith(ref_occ)]
            
            for p in properties:
                parsed_p = p.replace("/", "_over_")
                
                v = ligands_descriptors.loc[ligands_descriptors["name"] == lig_id, ref_occ + p].iloc[0]
                
                triples += helper_write_blank_node_sentence(pfx.inp5("hasProperty"), [
                                [pfx.rdf("type"), pfx.lgLrp(parsed_p)],
                                [pfx.inp5("reference"), pfx.dsL("tmQMg-L")],
                                [pfx.lgLrp("basedOn"), uri_occ],
                                [pfx.rdf("value"), helper_typed_literal(v)]
                            ] + tmQM_series_property_table[(parsed_p, "tmQMg-L")])
                
                property_atlas[pfx.lgLrp] += [parsed_p]
                
        """
        Close the 'paragraph'
        """
        triples += ["."]
    
    """
    Parse the property_atlas
    """
    property_atlas = {prefix: set(p_list) for prefix, p_list in property_atlas.items()}
    # """
    # Writes the triples to file
    # """
    # with open(OUTPUT_FILES["tmQM_RDF_L"], "w") as f:
    #     f.write("\n".join(triples))

    print(f"Ligand ontology built in {time.time() - start_time} seconds.")
    
    return triples, property_atlas

def build_atom_reference_ontology(reference_element_set, periodic_table):
    """
    Builds an ontology with all the reference information of the elements used in tmQM-RDF
    
    Arguments:
        - reference_element_set: the set of elements returned by assemble_rdf
        - periodic_table: a pd.DataFrame read from misc/PubChemElements_all.csv by read_global_files
    
    Returns:
        - the list of triples that constitute the ontology
        - the property_atlas {prefix: set([properties]), ...} of the properties encountered during the construction
    """
    
    property_atlas = defaultdict(list)
    
    start_time = time.time()
    print("Building elements ontology...")
    
    triples = []
    
    """
    Build element classes
    """
    atnum = lambda smbl: int(periodic_table.loc[periodic_table["Symbol"] == smbl, "AtomicNumber"].iloc[0])
    
    atom_label = "Chemical element: %s (%s)"
    atom_description = "This class represents the chemical element %s (%s), intended as abstract representation of an atom with atomic number %s."
    for el in sorted(reference_element_set, key = atnum):
        triples += [""]
        triples += [pfx.tmAr(el)]
        
        """
        Add superclass
        """
        triples += ["\t" + pfx.rdfs("subClassOf") + " " + pfx.tmA("Atom") + " ;"]
        
        """
        Add type
        """
        triples += ["\t" + pfx.rdf("type") + " " + pfx.tmAr("Element") + " ;"]
        
        """
        Add human-like description
        """
        el_name = periodic_table.loc[periodic_table["Symbol"] == el, "Name"].iloc[0]
        el_n = periodic_table.loc[periodic_table["Symbol"] == el, "AtomicNumber"].iloc[0]
        
        triples += ["\t" + pfx.rdfs("label") + " " + helper_typed_literal( atom_label % (el, el_name) ) + " ;"]
        triples += ["\t" + pfx.rdfs("comment") + " " + helper_typed_literal( atom_description % (el, el_name, el_n) ) + " ;"]
        
        """
        Add atomic number
        """
        triples += helper_write_blank_node_sentence(pfx.inp5("hasProperty"), [
                                                        [pfx.rdf("type"), pfx.tmArp("atomic_number")],
                                                        [pfx.rdf("value"), helper_typed_literal(atnum(el))]
                                                    ])
        property_atlas[pfx.tmArp] += ["atomic_number"]
        
        """
        Close 'paragraph'
        """
        triples += ["."]

    print(f"Elements ontology built in {time.time() - start_time} seconds.")
    
    """
    Parse the property_atlas
    """
    property_atlas = {prefix: set(p_list) for prefix, p_list in property_atlas.items()}
    
    return triples, property_atlas

def build_misc_reference_ontology():
    """
    Builds an ontology with all the reference information of various miscellanea (NBO types, countable structural features of ligands,
    ligand occurrence types)
    
    Returns:
        - the list of triples that constitute the ontology
    """
    
    opt_obj = {
            pfx.nm("GFN2_xTB"): [
                    [pfx.nm("Optimisation")],
                    "GFN2-xTB",
                    "A semi-empirical method in which the solution to the Schrödinger equation is approached through parametrization"
                ],
            pfx.nm("PBE_D3BJ__def2_SVP"): [
                    [pfx.nm("Optimisation")],
                    "PBE-D3BJ/def2-SVP",
                    "A semi-empirical method in which the solution to the Schrödinger equation is approached through parametrization"
                ]
            }
    
    sgl_obj = {
            pfx.nm("GFN2_xTB"): [
                    [pfx.nm("Singlepoint")],
                    "GFN2-xTB",
                    "A semi-empirical method in which the solution to the Schrödinger equation is approached through parametrization"
                ],
            pfx.nm("PBE0_D3BJ__def2_SVP"): [
                    [pfx.nm("Singlepoint")],
                    "PBE0-D3BJ/def2-SVP",
                    "Hybrid GGA DFT functional with a double-z quality basis set"
                ],
            pfx.nm("PBE0_D3BJ__def2_TZVP"): [
                    [pfx.nm("Singlepoint")],
                    "PBE0-D3BJ/def2-TZVP",
                    "Hybrid GGA DFT functional with a triple-z quality basis set"
                ],
            pfx.nm("PBE_D3BJ__def2_SVP"): [
                    [pfx.nm("Singlepoint")],
                    "PBE-D3BJ/def2-SVP",
                    "GGA DFT functional with a double-z quality basis set"
                ],
            pfx.nm("TPSSh_D3BJ__def2_SVP"): [
                    [pfx.nm("Singlepoint")],
                    "TPSSh-D3BJ/def2-SVP",
                    "Hybrid meta-GGA DFT functional with a triple-z quality basis set"
                ]
        }
    
    nbo_obj = {
            pfx.tmBr("NBOTypeBD"): [
                    [pfx.tmB("NBOType")],
                    "NBO type: BD",
                    "Localized bonding orbital in the framework of NBO theory"
                ],
            pfx.tmBr("NBOTypeNone"): [
                    [pfx.tmB("NBOType")],
                    "NBO type: None",
                    ""
                ]
        }
    
    str_obj = {
            pfx.lgLrm("AliphaticRing"): [
                    [pfx.lgLr("StructuralFeature")], 
                    "Aliphatic ring",
                    "This class represents a structural feature class identifiable as an aliphatic ring."
                ],
            pfx.lgLrm("AromaticRing"): [
                    [pfx.lgLr("StructuralFeature")], 
                    "Aromatic ring",
                    "This class represents a structural feature class identifiable as an aromatic ring."
                ],
            pfx.lgLrm("SaturatedRing"): [
                    [pfx.lgLr("StructuralFeature")], 
                    "Saturated ring",
                    "This class represents a structural feature class identifiable as a saturated ring."
                ],
            pfx.lgLrm("RotatableBond"): [
                    [pfx.lgLr("StructuralFeature")], 
                    "Rotatable bond",
                    "This class represents a structural feature class identifiable as a rotatable bond."
                ]
        }
    
    occ_obj = {
            pfx.lgLrm("SMILES"): [
                    [pfx.lgLro("LigandOccurrenceClass")], 
                    "Occurrence: SMILES",
                    "This class represents a ligand occurrence identified as the SMILES string."
                ],
            pfx.lgLrm("most_stable"): [
                    [pfx.lgLro("LigandOccurrenceClass")], 
                    "Occurrence: most stable",
                    "This class represents a ligand occurrence identified as the most stable occurrence."
                ],
            pfx.lgLrm("relaxes_structure"): [
                    [pfx.lgLro("LigandOccurrenceClass")], 
                    "Occurrence: SMILES",
                    "This class represents a ligand occurrence identified as the relaxed structure."
                ]
        }
    
    sections = {
            "Numerical (optimisation)": opt_obj,
            "Numerical (singlepoint)": sgl_obj,
            "NBO Types": nbo_obj,
            "Structural Features": str_obj,
            "Ligand Occurrences": occ_obj
        }
    
    """
    Start building the ontology
    """
    triples = []
    
    for header, objects in sections.items():
        triples += [f"\n# -- {header} --"]
        
        for obj_name, obj_body in objects.items():
            triples += [f"\n{obj_name}"]
            
            triples += ["\t" + pfx.rdf("type") + " " + ", ".join(obj_body[0]) + " ;"]
                
            triples += ["\t" + pfx.rdfs("label") + " " + helper_typed_literal(obj_body[1]) + " ;"]
            
            triples += ["\t" + pfx.rdfs("comment") + " " + helper_typed_literal(obj_body[2]) + " ;"]
            
            triples += ["."]
            
    return triples

def build_property_ontology(tmc_property_atlas, lig_property_atlas, atom_property_atlas):
    """
    Builds the ontology relative to properties.
    
    Argruments:
        - tmc_property_atlas: fifth output of assemble_rdf
        - lig_property_atlas: second output of build_ligand_reference_ontology
        - atom_property_atlas: second output of build_atom_reference_ontology
        
    Returns:
        - the triples
    """
    
    descriptions = pd.read_csv(INPUT_FILES["tmQM_series_property_descriptions"])
    
    print("Building properties ontology...")
    
    property_atlas = defaultdict(set)
    
    for prefix, p_set in tmc_property_atlas.items():
        property_atlas[prefix] = property_atlas[prefix].union(p_set)
        
    for prefix, p_set in lig_property_atlas.items():
        property_atlas[prefix] = property_atlas[prefix].union(p_set)
        
    for prefix, p_set in atom_property_atlas.items():
        property_atlas[prefix] = property_atlas[prefix].union(p_set)

    triples = []
    
    for prefix, p_set in property_atlas.items():
        for p in p_set:
            p_name = p.split("/")[-1]
            lab = descriptions.loc[descriptions["PROPERTY NAME"] == p_name, ["LABEL"]].iloc[0, 0]
            desc = descriptions.loc[descriptions["PROPERTY NAME"] == p_name, ["COMMENT"]].iloc[0, 0]
            
            triples += [""]
            triples += [prefix(p)]
            
            triples += ["\t" + pfx.rdfs("subClassOf") + " " + pfx.inp5("ObservedProperty") + " ;"]
            
            triples += ["\t" + pfx.rdf("type") + " " + pfx.inp5("PropertyType") + " ;"]
                
            triples += ["\t" + pfx.rdfs("label") + " " + helper_typed_literal(lab) + " ;"]
            
            triples += ["\t" + pfx.rdfs("comment") + " " + helper_typed_literal(desc) + " ;"]
            
            triples += ["."]
            
    return triples

def assemble_concepts():
    """
    This function assembles the concepts (classes and predicates) defined by build_classes_ontology
    and build_predicates_ontology into a single .ttl file.
    """
    
    triples = []
    
    triples += helper_write_prefixes()
    
    triples += [""]
    triples += ["# --------------------"]
    triples += ["# ------ CLASSES"]
    triples += ["# --------------------"]
    triples += build_classes_ontology()
    
    triples += [""]
    triples += ["# --------------------"]
    triples += ["# ------ PREDICATES"]
    triples += ["# --------------------"]
    triples += build_predicates_ontology()
    
    """
    Writes the triples to file
    """
    with open(OUTPUT_FILES["tmQM_RDFS_C"], "w") as f:
        f.write("\n".join(triples))

def assemble_ligands(ligands_misc, ligands_fingerprints, ligands_descriptors, reference_ligand_atlas, reference_centres_atlas, subgraph_atlas, tmQM_series_property_table, periodic_table):
    """
    This function assembles the ligands defined by build_ligand_reference_ontology into a single .ttl file.
    
    Arguments:
        - ligands_misc: a pd.DataFrame read from tmQMg-L/ligands_misc_info.csv by read_global_files
        - ligands_fingerprints: a pd.DataFrame read from tmQMg-L/ligands_fingerprints.csv by read_global_files
        - ligands_descriptors: a pd.DataFrame read from tmQMg-L/ligands_descriptors.csv by read_global_files
        - reference_ligand_atlas: first output of assemble_rdf
        - reference_centres_atlas: second output of assemble_rdf
        - subgraph_atlas: third output of assemble_rdf
        - tmQM_series_property_table: a defaultdict(list) built from tmQM_series_property_table.csv by read_global_files
        - periodic_table: a pd.DataFrame read from misc/PubChemElements_all.csv by read_global_files
        
    Returns:
        - the property atlas built by build_ligand_reference_ontology
    """
    
    lig_triples, lig_property_atlas =  build_ligand_reference_ontology(
                                            ligands_misc, 
                                            ligands_fingerprints, 
                                            ligands_descriptors, 
                                            reference_ligand_atlas, 
                                            reference_centres_atlas, 
                                            subgraph_atlas, 
                                            tmQM_series_property_table, 
                                            periodic_table
                                        )
    
    triples = []
    
    triples += helper_write_prefixes()

    triples += [""]
    triples += ["# --------------------"]
    triples += ["# ------ LIGANDS"]
    triples += ["# --------------------"]
    triples += lig_triples

    """
    Writes the triples to file
    """
    with open(OUTPUT_FILES["tmQM_RDFS_L"], "w") as f:
        f.write("\n".join(triples))
        
    return lig_property_atlas

def assemble_elements(reference_element_set, periodic_table):
    """
    This function assembles the ligands defined by build_atom_reference_ontology into a single .ttl file.
    
    Arguments:
        - reference_element_set: fourth output of assemble_rdf
        - periodic_table: a pd.DataFrame read from misc/PubChemElements_all.csv by read_global_files
        
    Returns:
        - the property atlas built by build_atom_reference_ontology
    """
    
    atom_triples, atom_property_atlas = build_atom_reference_ontology(reference_element_set, periodic_table)
    
    triples = []
    
    triples += helper_write_prefixes()
    
    triples += [""]
    triples += ["# --------------------"]
    triples += ["# ------ ELEMENTS"]
    triples += ["# --------------------"]
    triples += atom_triples
    
    """
    Writes the triples to file
    """
    with open(OUTPUT_FILES["tmQM_RDFS_E"], "w") as f:
        f.write("\n".join(triples))
        
    return atom_property_atlas

def assemble_properties(tmc_property_atlas, lig_property_atlas, atom_property_atlas):
    """
    This function assembles the properties defined by build_property_ontology into a single .ttl file.
    
    Arguments:
        - tmc_property_atlas: the property atlas returned by assemble_rdf
        - lig_property_atlas: the property atlas returned by assmeble_ligands
        - atom_property_atlas: the property atlas returned by assemble atoms
    """
    triples = []
    
    triples += helper_write_prefixes()

    triples += [""]
    triples += ["# --------------------"]
    triples += ["# ------ PROPERTIES"]
    triples += ["# --------------------"]
    triples += build_property_ontology(tmc_property_atlas, lig_property_atlas, atom_property_atlas)

    """
    Writes the triples to file
    """
    with open(OUTPUT_FILES["tmQM_RDFS_P"], "w") as f:
        f.write("\n".join(triples))

def assemble_misc():
    """
    This function assembles the misc entities defined by build_misc_reference_ontology into a single .ttl file.
    """
    
    triples = []
    
    triples += helper_write_prefixes()
    
    triples += [""]
    triples += ["# --------------------"]
    triples += ["# ------ MISC"]
    triples += ["# --------------------"]
    triples += build_misc_reference_ontology()
    
    """
    Writes the triples to file
    """
    with open(OUTPUT_FILES["tmQM_RDFS_M"], "w") as f:
        f.write("\n".join(triples))
        
#%% Main body
if __name__ == "__main__":
    viables, ligands_atlas, ligands_misc, ligands_fingerprints, ligands_descriptors, tmQM_series_property_table, periodic_table = read_global_files()
    
    # viables = ["WAJJOH", "GOKTAD"] #TODO: DEBUG ONLY!!!!!!!!!!!!!
    # viables += list(pd.read_csv("../extract_selection_from_tmQM-RDF/1k_selection.csv").loc[:, "TMC"])
    # print("DEBUG VERSION ONLY!!!!!!!!!!!!!!!\n"*10)
    
    reference_ligand_atlas, reference_centres_atlas, subgraph_atlas, reference_element_set, tmc_property_atlas = assemble_rdf(viables, ligands_atlas, tmQM_series_property_table)
    
    assemble_concepts()
    lig_property_atlas = assemble_ligands(ligands_misc, ligands_fingerprints, ligands_descriptors, reference_ligand_atlas, reference_centres_atlas, subgraph_atlas, tmQM_series_property_table, periodic_table)
    atom_property_atlas = assemble_elements(reference_element_set, periodic_table)
    assemble_properties(tmc_property_atlas, lig_property_atlas, atom_property_atlas)
    assemble_misc()