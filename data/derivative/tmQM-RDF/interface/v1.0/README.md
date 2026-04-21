The file [tmQM_RDF_interface.py](/data/derivative/tmQM-RDF/interface/tmQM_RDF_interface.py) allows easy access to the information contained in tmQM-RDF about a specific TMC. This is achieved through an instance the class `TmQMRDF`, a dictionary-like object that allows to request and access specific subgraphs of the tmQM-RDF KG. 

> [!NOTE]
> Additional helper classes, subclasses of `_TmQMRDFSubgraph`, are introduced to handle different types of subgraphs. As a user, you do not need to initialise these kind of objects, as this is handled directly by `TmQMRDF`. All you need to know is which methods and attributes are available to you (see below).

---
# `TmQMRDF` Class:
---
A user-friendly interface to the tmQM-RDF knowledge graph.
Table of contents:
  1. Accessing tmQM-RDF subgraphs
      1. Merging subgraphs
  2. Import settings
  3. TmQMRDF methods
  
## 1. Accessing tmQM-RDF subgraphs
This class is a dictionary-like object that internally stores selected subgraphs of tmQM-RDF.
Namely, a subgraph can be accessed via the key (\<category\>, \<code\>), where:
  - category can be either "TMC", "ligand", "centre" or "element"
- code is the corresponding code of the desired object (either the CSD code, the tmQMg-L ligand id, or the chemical symbol for both centres and elements)
Initially, not all subgraphs are available. The desired subgraphs have to be "exposed" using the self.expose method. For quick access
to single TMCs, ligand species or elements, the methods self.tmc, self.ligand or self.element can also be used.

Each subgraph (the value returned by self[\<category\>, \<code\>]) exposes its rdf triples via the .rdf attribute, which is
an rdflib.Graph object.

  ### i. Merging subgraphs
  See rdflib.Graph for more details on how to merge subgraphs.
  The method self.as_knowledge_graph() can be used to merge all the exposed subgraphs into a single graph.

## 2. Import settings
The actual import behaviour behaviour of the class is governed by the method self.set_import_setting.
By default, the following settings are enabled:
  - path_to_chem_info = None,
  - pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/periodictable/CSV?response_type=save&response_basename=PubChemElements_all",
  - request_tmc_atoms = True, 
  - request_tmc_atomic_bonds = True,
  - request_tmc_ligands = True,
  - request_tmc_ligand_bonds = True,
  - auto_expose_tmc_ligands = False

## 3. TmQMRDF methods
```
__init__(path)
```
Main interface to tmQM-RDF.  

__Arguments__:
  - path: the path to the root directory of tmQM-RDF

```
set_import_settings(
      path_to_chem_info = None,
      pubchem_url = None,
      request_tmc_atoms = None, 
      request_tmc_atomic_bonds = None,
      request_tmc_ligands = None,
      request_tmc_ligand_bonds = None,
      auto_expose_tmc_ligands = None,
      force = False
  )
```
Sets the import settings. Only settings that are not None are saved. If a setting is already set and force = False, the corresponding argument is always discarded.  

__Arguments__:
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

```
expose(tmcs = [], ligands = [], centres = [], elements = [], n_cores = 1)
```
Extract one or more subgraphs corresponding to TMCs/ligand species/elements and makes them accessible via self.\_\_getitem\_\_  

__Arguments__:
  - tmcs: iterable containing CSD codes of TMCs to extract
  - ligands: iterable containing tmQMg-L ids of ligand species to extract
  - centres: iterable containing the chemical symbols of the metal centres to extract
  - elements: iterable containing the chemical symbols of elements to extract
  - n_cores: how many processes should be used to extract the TMCs. Default: 1

```
tmc(csd_code)
```
Wrapper for:  
    self.expose(tmcs = \[csd_code\])  
    self\["TMC", csd_code\]  
__Returns__:
  - an instance of _TmQMRDFSubgraph (see below)

```
ligand(tmqmgl_code)
```
Wrapper for:  
    self.expose(ligands = \[tmqmgl_code\])
    self\["ligand", tmqmgl_code\]  
__Returns__:
  - an instance of _TmQMRDFSubgraph (see below)

```
centre(pubchem_code)
```
Wrapper for:  
    self.expose(centres = \[pubchem_code\])  
    self\["centre", pubchem_code\]  
__Returns__:
  - an instance of _TmQMRDFSubgraph (see below)

```
element(pubchem_code)
```
Wrapper for:  
    self.expose(elements = \[pubchem_code\])  
    self\["element", pubchem_code\]  
__Returns__:
  - an instance of _TmQMRDFSubgraph (see below)

```
as_knowledge_graph()
```
Merges all the currently exposed subgraphs in a single RDF graph.  
__Returns__:
  - an rdflib.graph.Graph object given by the union of all the exposed subgraphs

---
# `_TmQMRDFSubgraph` Class:
---
Each subgraph of the tmQM-RDF KG is an instance of (a subclass of) the class `_TmQMRDFSubgraph`. This class possesses:
 - the attribute `.rdf`: an rdflib.Graph instance representing the given subgraph;
 - the method `.query(query_object)`: a wrapper for `self.rdf.query(query_object)`.

## 1. TMC subgraphs
Subgraphs representing TMCs expose the following additional attributes:
  - atoms: a list of pairs (atom_name, atom_element), where
      - atom_name is the name of the object representing the atom, intended
          as the last element on the URI path (e.g. if the URI is resource://integreat/p5/atomic/atom/XXYYZZ_El1,
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

Moreover, these objects also possess the following methods:
```
skeleton()
```
A convenience method to compute the "skeleton" of a TMC, intended as the rdflib.graph.Graph object obtained by navigating the tmQM-RDF knowledge graph starting from the `cmT:XXYYZZ` TMC node and proceeding in a depth-first search along the direction of the edges, while only going through the URI nodes (i.e., avoiding blank nodes and literals).  
__Returns__:
  - an rdflib.graph.Graph object

```
as_graphviz(layout: str)
```
Creates the source code of a graphical representation of the TMC data endoded in the graph using the graphviz language.  
__Arguments__:
  - layout: the desired graphviz layout
__Returns__:
  - a graphviz.Source object

```
view(format: str, filename: str, layout: str)
```
Shows a graphical representation of the TMC data encoded by the graph (relies on as_graphviz) and saves it to a file.  
__Arguments__:
  - format: the desired graphical format
  - filename: the path to the file to which the image should be saved
  - layout: the desired graphviz layout

```
render(self, format: str, filename: str, layout: str)
```
Same as view, but it only saves the image withouth showing it in a new window. 
