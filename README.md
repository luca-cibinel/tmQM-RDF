# tmQM-RDF Dataset: a Knowledge Graph Representing Transition Metal Complexes
![latest version: 1.0](https://img.shields.io/badge/latest-v1.0-green)  
[![Download link](https://img.shields.io/badge/github-DOWNLOAD-blue?logo=github)](https://github.com/luca-cibinel/tmQM-RDF-archive/releases)  
---
This is the companion repository of the publication "[tmQM-RDF Dataset: a Knowledge Graph Representing Transition Metal Complexes](https://arxiv.org/abs/2602.07091)" (Cibinel et al., 2026). The knowledge graph (KG) contains data about 59,143 transition metal complexes (TMCs), obtained by coherently synthesising the information contained in the tmQM dataset series, composed of the [tmQM](https://github.com/uiocompcat/tmQM) dataset (Balcells and Skjelstad, 2020), the [tmQMg](https://github.com/uiocompcat/tmQMg) dataset (Kneiding et al, 2023) and the [tmQMg-L](https://github.com/hkneiding/tmQMg-L/tree/main?tab=readme-ov-file) dataset (Kneiding et al, 2024).

## Release hisitory
 - ![Static Badge](https://img.shields.io/badge/latest-green) April 2026, **v1.0** [[download](https://github.com/luca-cibinel/tmQM-RDF-archive/releases/tag/23042026)]: the first official release of tmQM-RDF. The main features of this version are:
   - it relies on the 2025 release of tmQM, the [v74.548k release of tmQMg](https://data.archive.sigma2.no/dataset/4f94f626-b18c-458c-946f-d7052bb05982/download/u-NatQ_graphs.zip) and the v74k release of tmQMg-L;
   - URIs are resolvable and will redirect to [https://www.integreat.no/research/rdf/tmqm-rdf-dataset](https://www.integreat.no/research/rdf/tmqm-rdf-dataset).
 - ![Static Badge](https://img.shields.io/badge/pre--release-orange) September 2025, **v2025dev** [[download](https://github.com/luca-cibinel/tmQM-RDF-archive/releases/tag/07022025)]: the development version of tmQM-RDF developed in 2025. All the experiments described in the paper are performed on this version of the dataset. The main features of this version are:
   - it relies on the 2024 release of tmQM, the [v74.637k release of tmQMg](https://data.archive.sigma2.no/dataset/cc354a73-7398-487f-83f8-4166caa8cc09/download/nird/home/hanneskn/tmQMg/uNatQ_graphs.zip) and the v60k release of tmQMg-L;
   - non resolvable URIs.

## This repository
In particular, this repository contains the [code used to generate the data](/data) and the [code used to perform the computational experiments](/computational) mentioned in the paper. _Notice that the experiments have been performed using the v2025dev development version of tmQM-RDF._

# Example
We report here an example of the RDF encoding of a TMC. We consider in particular the TMC which in the Cambridge Structural Database is referenced as KCEYPT. Note that, for simplicity, most of the actual file is omitted and only some of the main features are shown. Furthermore, the example is divided in the three fundamental layers used in tmQM-RDF: the _complex_, _ligand_, and _atomic_ levels.

**Complex layer:**
> cmT:KCEYPT  
>	&emsp;  ms:hasProperty \[  
>	&emsp;&emsp;  rdf:type cmTp:charge ;  
>	&emsp;&emsp; 	ms:reference dsC:tmQM ;  
>	&emsp;&emsp; 	rdf:value "-1"^^xmls:integer ;  
>	&emsp;  \] ;  
>	&emsp;  ms:hasProperty \[  
>	&emsp;&emsp; 	rdf:type cmTp:spin ;  
>	&emsp;&emsp; 	ms:reference dsC:tmQM ;  
>	&emsp;&emsp;rdf:value "0"^^xmls:integer ;  
>	&emsp;  \] ;  
>
> &emsp; ...  
> &emsp; ...  
> &emsp; ...
> 
> &emsp;  cmT:hasMetalCentre lgC:KCEYPT_metalCentre_Pt ;  
>	&emsp;  cmT:hasLigand lgL:KCEYPT_ligand0-0_0 ;  
>
> &emsp; ...  
> &emsp; ...  
> &emsp; ...
> .

**Ligand layer:**
> lgC:KCEYPT_metalCentre_Pt  
>	&emsp;  lgC:isMetalCentre lgCr:MetalCentre_Pt ;  
>	&emsp;  lgS:bLc lgB:KCEYPT_bondL_ligand0-0_0_0 ;   
>	&emsp; ...  
>	&emsp;  lgC:hasAtom tmA:KCEYPT_Pt_0 ;
> .  
> lgL:KCEYPT_ligand0-0_0  
>	&emsp;  lgL:isLigand lgLr:Ligand_ligand0-0 ;  
>	&emsp;  lgL:subgraphName "KCEYPT-subgraph-0"^^xmls:string ;  
>	&emsp;  lgS:bLl lgB:KCEYPT_bondL_ligand0-0_0_0 ;  
>	&emsp;  lgL:hasAtom tmA:KCEYPT_Cl_0 ;  
> .
>
> &emsp; ...  
> &emsp; ...  
> &emsp; ...

**Atomic layer:**
> tmA:KCEYPT_Pt_0  
>	&emsp;  tmA:isAtom tmAr:Pt ;  
>	&emsp;  ms:hasProperty \[  
>	&emsp;&emsp; 	rdf:type tmAp:natural_atomic_charge ;  
>	&emsp;&emsp; 	ms:reference dsG:tmQMg ;  
>	&emsp;&emsp; 	rdf:value "0.73094"^^xmls:float ;  
>	&emsp;&emsp; 	nm:singlepoint nm:PBE0_D3BJ__def2_TZVP ;  
>	&emsp;  \] ;  
> .     
>	&emsp; ...  
>	&emsp; ...  
>	&emsp; ...

# References
- Balcells, D. and B. B. Skjelstad (2020). Tmqm dataset-quantum geometries and properties of 86k transition metal complexes. *Journal of Chemical Information and Modeling 60.*
- Cibinel, L., T. Linjordet, J. Pensar, D. Balcells, R. De Bin, and B. Ell (2026). tmQM-RDF Dataset: a Knowledge Graph Representing Transition Metal Complexes. *arXiv preprint*,	arXiv:2602.07091.
- Kneiding, H., R. Lukin, L. Lang, S. Reine, T. B. Pedersen, R. De Bin, and D. Balcells (2023). Deep learning metal complex properties with natural quantum graphs. *Digital Discovery 2*, 618–633.
- Kneiding, H., A. Nova, and D. Balcells (2024). Directional multiobjective optimization of metal complexes at the billion-system scale. *Nature Computational Science 4*, 263–273.
