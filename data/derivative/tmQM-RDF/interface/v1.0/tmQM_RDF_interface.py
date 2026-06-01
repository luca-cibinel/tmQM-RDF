#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script defines a class that converts the RDF representation of TMCs employed in tmQM-RDF
into an networkx-encoded version of one of the three representations defined below (run this script for a visual representation).

For simplicity, it relies on the tmQMg dataset and it expands it using information from tmQMg-L (linked via tmQM-RDF).
"""

from pathlib import Path

import graphviz_syntax_wrapper as gsw
import urllib.request as url
import multiprocessing
import rdflib as rdf
import pandas as pd
import numpy as np
import collections
import tempfile
import graphviz
import re
import os

import warnings
np.warnings = warnings

class TmQMRDF(collections.UserDict):
    """
    A user-friendly interface to the tmQM-RDF knowledge graph.
    Table of contents:
        1. Accessing tmQM-RDF subgraphs
            1.1 Merging subgraphs
        2. Subgraph attributes
        3. Import settings

    1. Accessing tmQM-RDF subgraphs
    ---
    This class is a dictionary-like object that internally stores selected subgraphs of tmQM-RDF.
    Namely, a subgraph can be accessed via the key (<category>, <code>), where:
        - category can be either "TMC", "ligand", "centre" or "element"
        - code is the corresponding code of the desired object (either the CSD code, the tmQMg-L ligand id, or the chemical symbol for both centres and elements)
    Initially, not all subgraphs are available. The desired subgraphs have to be "exposed" using the self.expose method. For quick access
    to single TMCs, ligand species or elements, the methods self.tmc, self.ligand or self.element can also be used.

    Each subgraph (the value returned by self[<category>, <code>]) exposes its rdf triples via the .rdf attribute, which is
    an rdflib.Graph object.
    
        1.1 Merging subgraphs
        See rdflib.Graph for more details on how to merge subgraphs.
        The method self.as_knowledge_graph() can be used to merge all the exposed subgraphs into a single graph.

    2. Subgraph attributes
    ---
    Additional category-specific attributes may be available. Specifically, TMCs expose the following attributes:
        - atoms: a list of pairs (atom_name, atom_element), where
            - atom_name is the name of the object representing the atom, intended
                as the last element on the URI path (e.g. if the URI is https://www.integreat.no/research/rdf/tmqm-rdf-dataset/#/atomic/atom/XXYYZZ_El1,
                the name is XXYYZZ_El1)
            - atom_element is the chemical symbol of the element of the atom, again, extracted as the last element
                of the URI path of the element object
        - atomic_bonds: a list of pairs (atom1_name, atom2_name), where atom[n]_name is the name of the atom at the
            n-th end of the bond (n = 1,2). The atoms are ordered according to the IDs assigned in the
            tmQMg dataset, with id(atom1) < id(atom2).
        - ligands: a dictionary of the form {ligand_name: {'class': ligand_id, 'components': [...]}, ...}
            where:
                - ligand_name is the name of the ligand object (intended as the last element
                    on the URI path)
                - ligand_id is the id of the reference ligand
                - components is a list of the names of the atoms that compose this instance of the ligand
        - metal_centre: a dictionary of the form {'class': centre_element, 'components': [centre_atom]}
            where:
                - centre_element is the chemical element of the metal centre
                - centre_atom is the name of the atom representing the centre at the atomic level
        - ligand_bonds: a dictionary of the form {ligand_id: bonds, }, where ligand_id is the id of the ligand the
                bond refers to, and bonds is a list of lists (in the tmQMg-L sense) of the binding atoms.
    
    3. Import settings
    ---
    The actual import behaviour behaviour of the class is governed by the method self.set_import_setting.
    By default, the following settings are enabled:
        - path_to_chem_info = None,
        - pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/periodictable/CSV?response_type=save&response_basename=PubChemElements_all",
        - request_tmc_atoms = True, 
        - request_tmc_atomic_bonds = True,
        - request_tmc_ligands = True,
        - request_tmc_ligand_bonds = True,
        - auto_expose_tmc_ligands = False
    """

    def __init__(self, path):
        """
        Main interface to tmQM-RDF.

        - Arguments:
            - path: the path to the root directory of tmQM-RDF
        """
        super().__init__()

        self.path = path
        self.import_settings = dict()
        self.set_import_settings(
            None, 
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/periodictable/CSV?response_type=save&response_basename=PubChemElements_all",
            True, 
            True, 
            True, 
            True,
            False, 
            force = False
        )

    def set_import_settings(
                self,
                path_to_chem_info = None,
                pubchem_url = None,
                request_tmc_atoms = None, 
                request_tmc_atomic_bonds = None,
                request_tmc_ligands = None,
                request_tmc_ligand_bonds = None,
                auto_expose_tmc_ligands = None,
                force = False
            ):
        """
        Sets the import settings. Only settings that are not None are saved. If a setting is already set and force = False,
        the corresponding argument is always discarded.

        - Arguments:
            - path_to_chem_info: a string pointing to the .csv version of the pubchem periodic table (if not present and pubchem_url is
                available, it will be downloaded to this location)
            - pubchem_url: the url from which the .csv version of the pubchem periodic table should be downloaded
            - request_tmc_atoms: when importing a TMC, should its composing atoms be extracted in a user-friendly format?
            - request_tmc_atomic_bonds: when importing a TMC, should it atomic bonds be extracted in a user-friendly format?
            - request_tmc_ligands: when importing a TMC, should its composing ligands be extracted in a user-friendly format?
            - request_tmc_ligand_bonds: when importing a TMC, should its ligand-level bonds be extracted in a user-friendly format?
            - auto_expose_tmc_ligands: when importing a TMC via self.expose, should its composing ligands be exposed as well? If True, it
                overrides request_tmc_ligands and the import will always behave as if request_tmc_ligands = True.
            - force: should the old settings be overwritten?
        """
        for key, value in locals().items():
            if key == "force":
                continue

            if value is not None and (force or key not in self.import_settings):
                self.import_settings[key] = value

    def _read(self, args):
        return args[1](self.path, args[0], self.import_settings)

    def expose(self, tmcs = [], ligands = [], centres = [], elements = [], n_cores = 1):
        """
        Extract one or more subgraphs corresponding to TMCs/ligand species/elements and makes them 
        accessible via self.__getitem__

        - Arguments:
            - tmcs: iterable containing CSD codes of TMCs to extract
            - ligands: iterable containing tmQMg-L ids of ligand species to extract
            - centres: iterable containing the chemical symbols of the metal centres to extract
            - elements: iterable containing the chemical symbols of elements to extract
            - n_cores: how many processes should be used to extract the TMCs. Default: 1
        """

        for objects, Category in zip([tmcs, ligands, centres, elements], [_TMC, _Ligand, _Centre, _Element]):
            with multiprocessing.Pool(processes = n_cores) as pool:
                parsed_objects = [(obj, Category) for obj in objects if (Category.name, obj) not in self.data]

                if len(parsed_objects) == 0:
                    continue

                if n_cores > 1:
                    # pool = multiprocessing.Pool(processes = n_cores)
                    parse = pool.imap_unordered
                else:
                    parse = map

                for obj in parse(self._read, parsed_objects):
                    super().__setitem__((Category.name, obj.public_code), obj)

                    if self.import_settings["auto_expose_tmc_ligands"] and Category.name == "TMC":
                        local_ligands = [lig_info["class"] for lig_info in obj.ligands.values() if lig_info["class"] not in ligands]
                        ligands += local_ligands

    def tmc(self, csd_code):
        """
        Wrapper for
            self.expose(tmcs = [csd_code])
            self[("TMC", csd_code)]
        """
        self.expose(tmcs = [csd_code])

        return self["TMC", csd_code]

    def ligand(self, tmqmgl_code):
        """
        Wrapper for
            self.expose(ligands = [tmqmgl_code])

            return self["ligand", tmqmgl_code]
        """
        self.expose(ligands = [tmqmgl_code])

        return self["ligand", tmqmgl_code]

    def centre(self, pubchem_code):
        """
        Wrapper for
            self.expose(centres = [pubchem_code])

            return self["centre", pubchem_code]
        """
        self.expose(centres = [pubchem_code])

        return self["centre", pubchem_code]

    def element(self, pubchem_code):
        """
        Wrapper for
            self.expose(tmcs = [csd_code])
            self[("element", csd_code)]
        """
        self.expose(elements = [pubchem_code])

        return self["element", pubchem_code]

    def as_knowledge_graph(self):
        """
        Returns a single rdflib.Graph RDF graph given by the union of all the exposed subgraphs
        """

        return sum([g.rdf for g in self.values()], rdf.Graph())

class _TmQMRDFSubgraph:
    name = None

    def __init__(self, path, category, code):
        self._rdf_file = Path(os.path.join(path, "assertions", category, f"{code}.ttl")).absolute()
        self.rdf = rdf.Graph()
        self.rdf.parse(self._rdf_file)
        self.code = code
        self.public_code = code

    def query(self, query_object):
        """
        Wrapper for self.rdf.query.
        
        See rdflib.query.
        """
        return self.rdf.query(query_object)

class _TMC(_TmQMRDFSubgraph):
    name = "TMC"
    
    def __init__(self, path, tmc_name, import_settings):
        super().__init__(path, "TMCs", tmc_name)

        # Retrieve information such as CSD code, atoms, bonds and ligands/metal centre
        self.tmc_name = tmc_name
        self.CSD_code = tmc_name.replace("-", "_") # Sanitation is needed in case the interface is needed for constructs not coming from the CSD
        self.path_to_chem_info = import_settings.get("path_to_chem_info", None)
        if import_settings["request_tmc_atoms"]:
            self.atoms = self._get_atoms()
        if import_settings["request_tmc_atomic_bonds"]:
            self.atomic_bonds = self._get_chemical_bonds()
        if import_settings["request_tmc_ligand_bonds"]:
            self.ligand_bonds = self._get_ligand_bonds()
        if import_settings["request_tmc_ligands"] or type(self).import_settings["auto_expose_tmc_ligands"]:
            self.ligands, self.metal_centre = self._get_ligands_components()
    
    def _get_atoms(self):
        """
        (Private method)
        
        Runs a SPARQL query to extract the atoms of the TMC from the RDF graph
        
        Returns:
            A list of pairs (atom_name, atom_element), where
                - atom_name is the name of the object representing the atom, intended
                    as the last element on the URI path (e.g. if the URI is https://www.integreat.no/research/rdf/tmqm-rdf-dataset/#/atomic/atom/XXYYZZ_El1,
                    the name is XXYYZZ_El1)
                - atom_element is the chemical symbol of the element of the atom, again, extracted as the last element
                    of the URI path of the element object
        """
        
        retrieve_atoms = """
            SELECT DISTINCT ?atom ?element
            WHERE {
                ?atom <https://www.integreat.no/research/rdf/tmqm-rdf-dataset/#/atomic/atom/isAtom> ?element .
            }
            """
            
        qres = self.rdf.query(retrieve_atoms)
        
        return [(str(row.atom).split('/')[-1], str(row.element).split('/')[-1]) for row in qres]
    
    def _get_chemical_bonds(self):
        """
        (Private method)
        
        Runs a SPARQL query to extract the chemical bonds within the TMC from the RDF graph.
        
        Returns:
            A list of pairs (atom1_name, atom2_name), where atom[n]_name is the name of the atom at the
                n-th end of the bond (n = 1,2). The atoms are ordered according to the IDs assigned in the
                tmQMg dataset, with id(atom1) < id(atom2).
        """
        
        retrieve_bonds = """
            SELECT DISTINCT ?first ?second
            WHERE {
                   ?first <https://www.integreat.no/research/rdf/tmqm-rdf-dataset/#/atomic/structure/b> ?bnd .
                   ?second <https://www.integreat.no/research/rdf/tmqm-rdf-dataset/#/atomic/structure/b> ?bnd .
                   
                   FILTER (?first != ?second)
            }
            """
        
        qres = self.rdf.query(retrieve_bonds)
        
        all_bonds = [(str(row.first).split('/')[-1], str(row.second).split('/')[-1]) for row in qres]
        
        bonds = []
        # Remove duplicates
        for bond in all_bonds:
            if (bond[1], bond[0]) not in bonds:
                bonds += [bond]
        
        return bonds
    
    def _get_ligands_components(self):
        """
        (Private method)
        
        Runs a SPARQL query to extract the ligands (and the metal centre) of the TMC from the RDF graph.
        Notice that the metal centre is structurally represented as a ligand for simplicity.
        
        Returns:
            A pair of two dictionaries:
                - {ligand_name: {'class': ligand_id, 'components': [...]}, ...}
                    where:
                        - ligand_name is the name of the ligand object (intended as the last element
                          on the URI path)
                        - ligand_id is the id of the reference ligand
                        - components is a list of the names of the atoms that compose this instance of the ligand
                - {'class': centre_element, 'components': [centre_atom]}
                    where:
                        - centre_element is the chemical element of the metal centre
                        - centre_atom is the name of the atom representing the centre at the atomic level
        """
        
        retrieve_ligands = """
            SELECT DISTINCT ?ligand ?atom ?ligand_type
            WHERE {
                ?ligand <https://www.integreat.no/research/rdf/tmqm-rdf-dataset/#/ligand/ligand/hasAtom> ?atom .   
                ?ligand <https://www.integreat.no/research/rdf/tmqm-rdf-dataset/#/ligand/ligand/isLigand> ?ligand_type .
            }
            """
            
        qres = self.rdf.query(retrieve_ligands)
        
        out_lig = {}
        
        for row in qres:
            ligand = str(row.ligand).split("/")[-1]
            atom = str(row.atom).split("/")[-1]
            ligand_type = str(row.ligand_type).split("/")[-1].split("_")[-1]
            
            current_comp_list = out_lig.get(ligand, {"components": []})["components"]
            
            new_comp_list = current_comp_list + [atom]
            
            out_lig[ligand] = {"class": ligand_type, "components": new_comp_list}
        
        retrieve_centre = """
            SELECT DISTINCT ?metal_centre ?atom ?metal_centre_type
            WHERE {
                ?metal_centre <https://www.integreat.no/research/rdf/tmqm-rdf-dataset/#/ligand/centre/hasAtom> ?atom .   
                ?metal_centre <https://www.integreat.no/research/rdf/tmqm-rdf-dataset/#/ligand/centre/isMetalCentre> ?metal_centre_type .
            }
            """
            
        qres = self.rdf.query(retrieve_centre)
        
        for row in qres:
            out_centre = {
                "class": str(row.metal_centre_type).split("/")[-1].split("_")[-1],
                "components": [str(row.atom).split("/")[-1]]
            }
        
        return out_lig, out_centre
    
    def _get_ligand_bonds(self):
        """
        (Private method)
        
        Runs a SPARQL query to extract the ligand-level bonds within the TMC from the RDF graph.
        
        Returns:
            A dictionary of the form {ligand_id: bonds, }, where ligand_id is the id of the ligand the
                bond refers to, and bonds is a list of lists (in the tmQMg-L sense) of the binding atoms.
        """
        
        retrieve_bonds = """
            SELECT DISTINCT ?ligand ?bnd ?atom
            WHERE {
                   ?ligand lgS:bLl ?bnd .
                   ?bnd lgB:hasBindingAtom ?atom .
            }
            """
        
        qres = self.rdf.query(retrieve_bonds)
        
        all_bonds = [(str(row.ligand).split('/')[-1], str(row.bnd).split('/')[-1], str(row.atom).split('/')[-1]) for row in qres]
        
        bonds = dict()
        # Remove duplicates
        for ligand, bond, atom in all_bonds:
            bonds.setdefault(ligand, {})
            bonds[ligand].setdefault(bond, [])
            bonds[ligand][bond] += [atom]
        
        bonds = {ligand: list(bonds.values()) for ligand, bonds in bonds.items()}
        
        return bonds
    
    def skeleton(self):
        """
        This function computes the "skeleton" of a TMC (i.e. the RDF graph obtained
        from the corresponding tmQM-RDF entry via a depth-first search, rooted at the TMC node, allowed to
        move only via URIs) as an auxiliary RDF graph.
        
        Returns:
            - An rdflib.Graph object
        """
        
        # Extract TMC node
        temp = self.query("""
            SELECT ?tmc
            WHERE {
                ?tmc <https://www.integreat.no/research/rdf/tmqm-rdf-dataset/#/complex/TMC/hasMetalCentre> ?c .
            }""")
        source = next(iter(temp)).tmc
        
        # Perform DFS
        skel = rdf.Graph()
        
        to_parse = [source]
        visit = []
        
        while len(to_parse) > 0:
            subj = to_parse.pop()
            visit += [subj]
            
            for s, p, o in self.rdf.triples((subj, None, None)):
                if isinstance(o, rdf.URIRef):
                    skel.add((s, p, o))
                    
                    if o not in visit:
                        to_parse += [o]
                        
        return skel

    def as_graphviz(self, layout = "neato"):
        """
        Encodes the TMC as a graph described via the DOT language.
        In the resulting graph, each ligand (and the metal centre) is represented as a cluster and hence highlighted
        via a box named after the ligand id, so as to emphasize the structure of the ligand level. The entire graph
        is a cluster itself, so it is also highlighted by a box, representing the complex level. The atoms are labelled
        and colored according to their chemcial element, using the information extracted from the PubChem periodic table data
        downloaded from https://pubchem.ncbi.nlm.nih.gov/periodic-table/#view=list .
        
        Parameters:
            layout: the desired layout engine, one of 'dot' and 'neato'.
            
        Returns:
            A graphiv.Source object in which the TMC has been encoded
        """
        
        if self.path_to_chem_info is None:
            raise Exception("TMC import_setting 'path_to_chem_info' not set!")
        
        # Preprocess node statements
        #
        # For each atom, a node statement is created. This statement can then be accessed via the node name
        
        node_statements = {}
        
        for node, _ in self.atoms:
            node_statements[node] = gsw.NodeStatement(node)
        
        # Preprocess edge statements
        #
        # For each edge, an edge statement is created. This statement can be accessed via the pair (atom1_name, atom2_name).
        # See self._get_chemical_bonds for more information about the atom names.
        
        edge_statements = {}
        
        for at1, at2 in self.atomic_bonds:
            edge_statements[at1, at2] = gsw.EdgeStatement((at1, at2))
        
        
        # Define the graph object
        #
        # Create a graph and initialize the graph-level aesthetic properties (font and rankdir)
        #
        # Then add the graph-level attributes for nodes and edges
        
        G = gsw.Graph(
                name = self.CSD_code,
                fontname = "Helvetica,Arial,sans-serif",
                rankdir = "LR",
                compound = "true",
                layout = layout
            )
        
        G.add_statements(gsw.NodeStatement(
                fontname = "Helvetica,Arial,sans-serif",
                penwidth = 1.0, 
                color = "black"
                # width = 0.5,
                # height = 0.5
            ))
        
        G.add_statements(gsw.EdgeStatement(
                fontname = "Helvetica,Arial,sans-serif",
                len = 0.6
            ))
            
        
        # Create elemental subgraphs
        #
        # For each element encountered in the TMC, create a subgraph which collects all the atoms which possess
        # that element, and then for each subgraph define the attributes that apply a unique style to every atom
        # of the given element
        
        # Read PubChem information
        if not os.path.exists(os.path.join(self.path_to_chem_info, "PubChemElements_all.csv")):
            
            if not os.path.exists(self.path_to_chem_info):
                os.makedirs(self.path_to_chem_info)
                
            url.urlretrieve(
                    type(self).pubchem_url,
                    os.path.join(self.path_to_chem_info, "PubChemElements_all.csv")
                )
        
        pubchem = pd.read_csv(os.path.join(self.path_to_chem_info, "PubChemElements_all.csv"), sep = ",")
        
        # Extract the set of the elements found in the TMC
        elements = set([el for _, el in self.atoms])
        
        # Then, for each element...
        subgraphs = []
        for el in elements:
            # Extract all the node statements which pertain to the element in question (they will be
            #   added to the subgraph to assert the belonging of an atom)
            statements_by_el = [node_statements[node] for node, el0 in self.atoms if el0 == el]
            
            # Extract the PubChem information relative to the given element
            el_name = pubchem.loc[pubchem["Symbol"] == el, "Name"].iloc[0]
            el_col = pubchem.loc[pubchem["Symbol"] == el, "CPKHexColor"].iloc[0]
            
            # Check that the color extracted from PubChem is a valid HEX color, using a regex
            #   If the color is not valid, use a default color (pink)
            if not re.match(r"^(?:[0-9a-fA-F]{2}){3}$", el_col):
                el_col = "f78bb2"
            
            # Set up a font color: white for carbon (because of the dark background), black for all the others
            font_col = "black"
            if el == "C":
                font_col = "white"
            
            # Define the node statement that will apply the aesthetic attributes to all nodes in the subgraph
            node_size = "0.3" if el == "H" else "0.5"
            fontsize = "7" if el == "H" else "14"
            header = gsw.NodeStatement(
                style = "filled",
                fillcolor = f"#{el_col}",
                label = el,
                fontcolor = font_col,
                width = node_size,
                height = node_size,
                fontsie = fontsize,
                fixedsize = True
            )
            
            # Create the subgraph and save it
            subgraphs += [
                    gsw.Subgraph(
                            header,
                            *statements_by_el,
                            name = el_name
                        )
                ]
        
        # add all the subgraphs to the graph object
        G.add_statements(*subgraphs)
        
        # Create metal centre cluster (structurally speaking, same as a ligand)
        
        centre_name = self.metal_centre["components"][0]
        centre_cluster_name = "metal_centre_" + self.metal_centre["class"]
        
        node_statements[centre_name].attributes["peripheries"] = 3
        ligands_clusters = [
                gsw.Cluster(
                        node_statements[centre_name],
                        name = centre_cluster_name,
                        label = f"Metal centre: {self.metal_centre['class']}",
                        color = "blue",
                        fontcolor = "blue"
                    )
            ]
        
        
        # Identify bonds from centre to ligands and atoms bonded to the centre
        #   It will be needed also when creating the ligand clusters
        
        centre_lig_bonds = []
        centre_bonded_atoms = []
        
        # For each bond, isolate those that involve the metal centre
        #   Also, identify if the metal centre is the first (tail) or the second (head)
        #   element of the bond and specify ltail/lhead accordingly
        for edge, edge_statement in edge_statements.items():
            if centre_name in edge:
                edge_statement.attributes[
                        "ltail" if edge[0] == centre_name else "lhead"
                    ] = "cluster_" + centre_cluster_name
                edge_statement.attributes["len"] = 3
                edge_statement.attributes["style"] = "dashed"
                
                centre_lig_bonds += [edge_statement]
                
                centre_bonded_atoms += [edge[1] if edge[0] == centre_name else edge[0]]
        
        # For each ligand - centre bond, highlight ligand binding atoms
        for at in centre_bonded_atoms:
            node_statements[at].attributes["peripheries"] = 2
        
        # Create ligand clusters
        #
        # For each ligand, define a cluster that will encompass all the atoms and all the chenical bonds
        # that belong to the ligand. A bond is said to belong to the ligand when both its ends belong to the ligand
        
        # See self._get_ligands_components docstring for a description of ligand
        #   For each ligand...
        for ligand_name, ligand in self.ligands.items():
            # Extract all the node statements associated with the ligand
            #   except for those bonded to the metal centre (those will be added separately)
            nodes_in_ligand = [node_statements[at] for at in ligand["components"] if at not in centre_bonded_atoms]
            
            # Identify the atoms bonded to the centre and put them in their own subgraph
            nil_bonded_to_centre = [node_statements[at] for at in ligand["components"] if at in centre_bonded_atoms]
            
            btc_subgraph = gsw.Subgraph(
                    *nil_bonded_to_centre,
                    rank = "same"
                )
            
            # Extract all the edge statements that belong to the ligand
            edges_in_ligand = []
            
            for edge, edge_statement in edge_statements.items():
                if edge[0] in ligand["components"] and edge[1] in ligand["components"]:
                    edges_in_ligand += [edge_statement]
            
            # Create the cluster using the nodes and the edges and label it according
            #   to the ligand id
            ligands_clusters += [
                    gsw.Cluster(
                            *nodes_in_ligand,
                            *edges_in_ligand,
                            btc_subgraph,
                            name = ligand_name.replace("-", "__"),
                            label = ligand["class"],
                            color = "blue",
                            fontcolor = "blue"
                        )
                ]
    
        # Create TMC cluster using all the ligands clusters and the bonds to the metal centre
        
        tmc_cluster = gsw.Cluster(
                *centre_lig_bonds,
                *ligands_clusters,
                name = self.CSD_code,
                label = self.CSD_code,
                color = "orange",
                fontcolor = "orange"
            )
        
        
        # Add TMC cluster to the graph object
        
        G.add_statements(tmc_cluster)
        
        return G.assemble()["statement"]
    
    def view(self, format = "png", filename = None, layout = "neato"):
        """
        Produces an image of the TMC rendered via the graphviz module and visualises it.
        
        Parameters:
            - format: the desired output format for the resulting graphviz object (pdf, png, svg, ...)
            - filename: the name of the file (without the extension) to which the output should be saved (optional)
            - layout: the desired layout engine, one of 'dot' and 'neato'.
        """
        src = self.as_graphviz(layout)
        src = graphviz.Source(src, format = format)
        
        if filename is None:
            src.view(tempfile.mktemp(".gv"), cleanup = True)
        else:
            src.view(filename, cleanup = True)
            
    def render(self, format = "png", filename = None, layout = "neato"):
        """
        Produces an image of the TMC rendered via the graphviz module
        
        Parameters:
            - format: the desired output format for the resulting graphviz object (pdf, png, svg, ...)
            - filename: the name of the file (without the extension) to which the output should be saved (optional)
            - layout: the desired layout engine, one of 'dot' and 'neato'.
        """
        src = self.as_graphviz(layout)
        src = graphviz.Source(src, format = format)
        
        if filename is None:
            src.render(tempfile.mktemp(".gv"), cleanup = True)
        else:
            src.render(filename, cleanup = True)

class _Ligand(_TmQMRDFSubgraph):
    name = "ligand"
    
    def __init__(self, path, tmqmgl_code, import_settings):
        super().__init__(path, "ligands", tmqmgl_code)

        self.tmqmgl_code = tmqmgl_code

class _Centre(_TmQMRDFSubgraph):
    name = "centre"
    
    def __init__(self, path, pubchem_code, import_settings):
        super().__init__(path, "centres", "MetalCentre_" + pubchem_code)

        self.pubchem_code = pubchem_code
        self.public_code = pubchem_code

class _Element(_TmQMRDFSubgraph):
    name = "element"
    
    def __init__(self, path, pubchem_code, import_settings):
        super().__init__(path, "elements", pubchem_code)

        self.pubchem_code = pubchem_code

# %% Main
if __name__ == "__main__":
    __spec__ = None

    import sys
    bin_path = os.path.abspath(os.path.join(sys.executable, ".."))
    if bin_path not in os.environ["PATH"].split(":"):
        print(f"[ WARNING: the path {bin_path} could not be found inside os.environ['PATH']. It will be added now. ]")
        os.environ["PATH"] = f"{bin_path}:{os.environ['PATH']}"

    ROOT_DIR = os.path.abspath(".")
    while ".prj_root" not in os.listdir(ROOT_DIR):
        ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))
    
    instance = TmQMRDF(os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "data", "v1.0"))
    instance.set_import_settings(
        path_to_chem_info = os.path.join(ROOT_DIR, "data", "raw", "pubChem", "data"),
        auto_expose_tmc_ligands = True,
        force = True
    )

    g = instance.tmc("ABEVAH")
    
    