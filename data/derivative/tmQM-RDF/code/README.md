This document contains information on the data processing pipeline that is involved in the creation of tmQM-RDF.

**Note:**
 - This pipeline relies on the data coming from the tmQM series and other sources. To acquire the data, run the scripts `download_data.py` in the subfolders of [/data/raw/](/data/raw/).
 - The pipeline is virtually identical for all the versions of tmQM-RDF, unless otherwise specified

All of the scripts mentioned below can be found in /data/derivative/tmQM-RDF/code/[version]/. For brevity, the paths mentioned below (and within the scripts) are not the full paths, but follow the following convention:
 - Paths starting with tmQM/, /tmQMg, /tmQMg-L and /misc are to be interpreted as subpaths of /data/raw/tmQMseries](/data/raw/tmQMseries)
 - Paths starting with pubChem/ are to be interpreted as subpaths of [/data/raw/pubChem](/data/raw/pubChem)
 - Paths starting with intermediate/ are to be interpreted as subpaths of /data/derivative/tmQM-RDF/[version]
 - Paths starting with data/ are to be interpreted as subpaths of /data/derivative/tmQM-RDF/data/[version]

0. Initialise directory:
   1. Run the script 0_i_initialise_directory.py (it will create the subdirectories `intermediate/[tmQM/tmQMg/tmQMg-L]` and `data/`)
1. Link together the tmQMg and tmQMg-L datasets:
   1. Run the script 1_i_encode_ligand_subgraphs.py (it will create the file `intermediate/tmQMg-L/ligands_atoms_idx.csv`, available in this repository)
      > From the file `tmQMg-L/ligands_xyzs.xyz`, for each subgraph `XXYYZZ-subgraph-n` it extracts the coordinates of the atoms, then matches them to the coordinates stored in the corresponding file `tmQMg/uNatQ_graphs/XXYYZZ.gml`. This identifies the atoms which compose the subgraph via the indices used in tmQMg. The resulting lists of indices ares stored in the file `intermediate/tmQMg-L/ligands_atoms_idx.csv`.
   2. Run the script 1_ii_encode_ligand_structure.py (it will create the files `intermediate/tmQMg-L/uNatQ_graphs/XXYYZZ.csv`, available in this repository as `intermediate/tmQMg-L/uNatQ_graphs.zip`)
      > Links every TMC in tmQMg to its ligands by combining the data from `tmQMg-L/ligands_misc_info.csv` and from `intermediate/tmQMg-L/ligands_atoms_idx.csv`. The file `tmQMg-L/ligands_misc_info.csv` contains, for each ligand, a list of its instances, each paired with the list of binding atoms (identified via the position in the file `tmQMg-L/ligands_xyzs.xyz`). This scripts translates the binding atoms to the indices used by tmQMg via the file `intermediate/tmQMg-L/ligands_atoms_idx.csv`. The result is the set of files `intermediate/tmQMg-L/uNatQ_graphs/XXYYZZ.csv`, in which, for each TMC, each composing ligand is identified (via ligand ID), linked to the corresponding subgraph and to its appropriately identified binding atoms.
   3. Sanity check: run the script 1_iii_ligands_sanity.py (it will create the files `intermediate/tmQMg-L/no_coverage.txt`, `intermediate/tmQMg-L/no_data.txt` and `intermediate/step_1_viable_tmcs.txt`, all available in this repository)
      > Performs a sanity check on the .csv files `intermediate/tmQMg-L/uNatQ_graphs/XXYYZZ.csv`. In particular, it checks whether, for each TMC available in tmQMg, 1) it has been possible to find some form of ligand data; 2) the available ligand data cover the entire TMC. Three lists are created: `intermediate/tmQMg-L/no_data.txt` (complexes with no available ligand data), `intermediate/tmQMg-L/no_coverage.txt` (complexes with incomplete ligand data) and `intermediate/step_1_viable_tmcs.txt` (complexes with complete ligand data). It is expected that `intermediate/tmQMg-L/no_data.txt` is a subset of `intermediate/tmQMg-L/no_coverage.txt`.
2. Link together the tmQM and tmQMg datasets:
   1. Run the script 2_i_locate_in_tmQM.py (it will create the file `intermediate/tmQM/localised_tmcs.csv` and the file `intermediate/step_2_viable_tmcs.txt`, all available in this repository)
      > Locates each viable TMC (as specified in `intermediate/step_1_viable_tmcs.txt`) inside the files `tmQM/tmQM_X[1/2/3].[xyz/BO]` and produces the file `intermediate/tmQM/localised_tmcs.csv` which links every complex to the number of the .xyz/.BO file. It also updates the file `intermediate/step_1_viable_tmcs.txt` into `intermediate/step_2_viable_tmcs.txt` by removing those TMCs for which it has not been possible to find complete information in tmQM.
   2. Run the script 2_ii_summarise_tmQM.py (it will create the files `intermediate/tmQM/uNatQ_graphs/XXYYZZ.csv` and it will update the file `intermediate/step_2_viable_tmcs.txt` into `intermediate/viable_tmcs.txt`; the content of the folder `intermediate/tmQM/uNatQ_graphs` is available in this repository within the files `tmQM/uNatQ_graphs_*.zip`; they can be conveniently decompressed via extract_tmQM_graphs.py)
      > It summarises all of the information in tmQM into the files `intermediate/tmQM/uNatQ_graphs/XXYYZZ.csv`, which mirror the structure of those in tmQMg. The indexing of the nodes matches that of tmQMg. It also performs progressive sanity checks to ensure, once more, the correctness of the localisation performed at the previous step, and also to ensure that the labelling of the nodes extracted from the files `tmQM/tmQM_X[1/2/3].BO` is consistent with that of tmQMg. Those complexes for which the labelling is inconsistent are removed from `intermediate/viable_tmcs.txt`. The content of the folder `intermediate/tmQM/uNatQ_graphs/` is compressed into three archives of the form `intermediate/tmQM/uNatQ_graphs_[X]_[Y].zip`, which contain each all of the TMC from letter `[X]` to letter `[Y]`.
3. Merge the information into tmQM-RDF:
   1. Run 3_i_rdf_builder.py (it will create the files `data/graphs/XXYYZZ.ttl` and the files `data/tmQM-RDFS_*.ttl`)
      > It summarises all the information in tmQM, tmQMg and tmQMg-L into single `data/graphsXXYYZZ.ttl` files and also creates an RDFS structure via the header files `data/tmQM-RDFS_*.ttl`.
4. Compress tmQM-RDF for GitHub push:
   1. Run 4_i_compress.py (it will create the files `data/archive/graphs/[X].zip` and `data/archive/tmQM-RDFS.zip`)
      > It compresses all the graphs in `data/graphs` into .zip archives according to the first letter of their CSD code. The RDFS header files are compressed into `data/archive/tmQM-RDFS.ttl`.
