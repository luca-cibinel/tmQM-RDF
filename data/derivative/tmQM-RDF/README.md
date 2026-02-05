# The tmQM-RDF dataset: a Knowledge Graph Representing Transition Metal Complexes 
This folder contains the [code](/data/derivative/tmQM-RDF/code) used to create the knowledge graph plus an [archive](/data/derivative/tmQM-RDF/data) containing a compressed local copy of tmQM-RDF (see below). 

# Data availability
*This repository only contains a convenience copy of tmQM-RDF, to be referred to only for reproducibility purposes with respect to the publication. See the [main page](/README.md) for download links.*

The tmQM-RDF dataset is stored in the folder [/data/derivative/tmQM-RDF/data](/data/derivative/tmQM-RDF/data), under the subfolder corresponding to the desired version. For reasons of space, this GitHub repository only stores a compressed version of the dataset, in the data/[version]/archive subfolder. A convenience script is provided to correctly unpack the compressed archives into the correct directories: see [/data/derivative/tmQM-RDF/interface/archive_utils/extract_tmQM_RDF_archive.py](/data/derivative/tmQM-RDF/interface/archive_utils/extract_tmQM_RDF_archive.py). In the script, be sure to specify the correct dataset version in line 23.

## Interface
The file [/data/derivative/tmQM-RDF/interface/tmQM_RDF_interface.py](/data/derivative/tmQM-RDF/interface/tmQM_RDF_interface.py) defines two utility classes that can be used to easily access the information contained in tmQM-RDF. See the [documentation](/data/derivative/tmQM-RDF/interface/README.md).

# Data processing pipeline:
Details on the data-processing pipeline can be found [here](/data/derivative/tmQM-RDF/code/README.md).

# References
- Balcells, D. and B. B. Skjelstad (2020). Tmqm dataset-quantum geometries and properties of 86k transition metal complexes. *Journal of Chemical Information and Modeling 60.*
- Kneiding, H., R. Lukin, L. Lang, S. Reine, T. B. Pedersen, R. De Bin, and D. Balcells (2023). Deep learning metal complex properties with natural quantum graphs. *Digital Discovery 2*, 618–633.
- Kneiding, H., A. Nova, and D. Balcells (2024). Directional multiobjective optimization of metal complexes at the billion-system scale. *Nature Computational Science 4*, 263–273.
