# tmQM-RDF Dataset: a Knowledge Graph Representing Transition Metal Complexes
This repository contains the material associated with the publication presenting the tmQM-RDF knowledge graph (KG). The KG contains data about 47,814 transition metal complexes (TMCs), obtained by coherently synthesising the information contained in the tmQM dataset series, composed of the [tmQM](https://github.com/uiocompcat/tmQM) dataset (Balcells and Skjelstad, 2020), the [tmQMg](https://github.com/uiocompcat/tmQMg) dataset (Kneiding et al, 2023) and the [tmQMg-L](https://github.com/hkneiding/tmQMg-L/tree/main?tab=readme-ov-file) dataset (Kneiding et al, 2024).

![tmQM-RDF-1](/tmQM-RDF-1.png)

## Download links (coming soon)
 - **v2025dev**: development version [download] [documentation]

## Release history:
 - September 2025, *v2025dev*: the development version of tmQM-RDF developed in 2025. All the experiments described in the paper are performed on this version of the dataset. The main features of this version are:
   - It relies on the 2024 release of tmQM, the v74k (September 2024) release of tmQMg and the v60k release of tmQMg-L.
   - Non resolvable URIs.

## This repository
In particular, this repository contains the [code used to generate the data](/data) and the [code used to perform the computational experiments](/computational) mentioned in the paper. _Notice that the experiments have been performed using the v2025dev development version of tmQM-RDF._

# References
- Balcells, D. and B. B. Skjelstad (2020). Tmqm dataset-quantum geometries and properties of 86k transition metal complexes. *Journal of Chemical Information and Modeling 60.*
- Kneiding, H., R. Lukin, L. Lang, S. Reine, T. B. Pedersen, R. De Bin, and D. Balcells (2023). Deep learning metal complex properties with natural quantum graphs. *Digital Discovery 2*, 618–633.
- Kneiding, H., A. Nova, and D. Balcells (2024). Directional multiobjective optimization of metal complexes at the billion-system scale. *Nature Computational Science 4*, 263–273.
