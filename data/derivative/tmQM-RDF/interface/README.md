This directory contains code that can be used to conveniently manage the tmQM-RDF dataset available in the [data](/data/) folder

# tmQM-RDF interface
The file [tmQM_RDF_interface.py](/data/derivative/tmQM-RDF/interface/tmQM_RDF_interface.py) defines two classes to easily access the information contained in tmQM-RDF about a specific TMC.

### Class: TmQMRDFInterface

A baseline class that manages some base features in accessing tmQM-RDF.
Access to the TMCs in tmQM-RDF is granted by specifying the path to the tmQM-RDF dataset's "graphs" subfolder. The path must be assigned to the static variable TmQMRDFInterface.path_to_tmQM_RDF.

**Attributes:**
- tmc_name: string, the CSD code of the TMC
- rdf: rdflib.graph.Graph, the rdf graph that contains the information extracted from tmQM-RDF

**Methods:**
```
__init__(self, tmc_name: str)
```
Class constructor
- Arguments:
  - tmc_name: the CSD code of the TMC you wish to retrieve

<br>

```
query(self, query_object: Union[str, Query])
```
A wrapper for the method [rdflib.graph.Graph.query](https://rdflib.readthedocs.io/en/7.1.1/apidocs/rdflib.html#rdflib.graph.Graph.query)

<br>

```
skeleton(self)
```
A convenience method to compute the "skeleton" of a TMC, intended as the rdflib.graph.Graph object obtained by navigating the tmQM-RDF knowledge graph starting from the `cmT:XXYYZZ` TMC node and proceeding in a depth-first search along the direction of the edges, while only going through the URI nodes (i.e., avoiding blank nodes and literals)
- Returns:
  - an rdflib.graph.Graph object
  
### Class: TmQMRDFGraph

Extends `TmQMRDFInterface`. This is the main class used to access tmQM-RDF in order to extract information about a specific TMC. In addition to access to tmQM-RDF, this class needs access to the PubChem periodic table, whose location has to be specified by assigning the path to the static variable  `TmQMRDFGraph.path_to_chem_info`. If the variable is assigned, but the data is not present at the given location, the class will automatically download it from [here](https://pubchem.ncbi.nlm.nih.gov/rest/pug/periodictable/CSV?response_type=save&response_basename=PubChemElements_all) when needed.

**Attributes:**
- CSD_code: string, the CSD code of the represented TMC
- atoms: list, a list of two-element tuples, where the first element is the specific part of the URI identifying the atom in the knowledge graph and the second is the chemical symbol of the atom
- bonds: list, a list of two-element tuples. Each tuple represent a chemical bond between two atoms and each element of the tuple is the specific part of the URI identifying wither of the atoms participating in the bond
- ligands: dict, a dictionary where the keys are the specific parts of the URIs describing each ligand in the TMC and the values are themselves dictionaries with the following entries:
  - class: string, the id of the ligand (in the tmQMg-L sense)
  - components: list, a list of the specific parts of the URIs of the atoms that are part of the ligand
- metal_centre: dict, a dictionary edscribing the metal centre, with the same "class" and "components" entries as the values of `ligands`

**Methods:**
```
as_graphviz(self, layout: str)
```
Creates the source code of a graphical representation of the TMC data endoded in the graph using the graphviz language.
- Arguments:
  - layout: the desired graphviz layout
- Returns:
  - a graphviz.Source object

<br>

```
view(self, format: str, filename: str, layout: str)
```
Shows a graphical representation of the TMC data encoded by the graph (relies on as_graphviz) and saves it to a file.
- Arguments:
  - format: the desired graphical format
  - filename: the path to the file to which the image should be saved
  - layout: the desired graphviz layout
 
<br>

```
render(self, format: str, filename: str, layout: str)
```
Same as view, but it only saves the image withouth showing it in a new window. 

# Archive utils
As the copy of tmQM-RDF in this repository is a compressed archive, the file [archive_utils/extract_tmQM_RDF_archive.py](/data/derivative/tmQM-RDF/interface/archive_utils/extract_tmQM_RDF_archive.py) decompressess and organises the archive. Be sure to specify the desired version of tmQM-RDF in the script (line 23)
