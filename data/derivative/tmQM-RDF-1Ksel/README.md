The 1K selection(s) of tmQM-RDF: a selection of 1000 training TMCs (plus 300 validation and 300 test TMCs), available in two versions, for early transition metals (Cr, Mo, W, tag `earlyTM`) and late transition metal (Pd, Ni, Pt, tag `lateTM`).

# Data availability
The two selections are available in the format of .csv files with two columns: TMC (the CSD codes of the TMCs in the selection) and partition (whether the TMC belongs to the training/validation/test set). The actual graphs are to be retrieved from the [tmQM-RDF dataset](/data/derivative/tmQM-RDF).

# Extraction pipeline
All of the scripts mentioned below can be found in [/data/derivative/tmQM-RDF-1Ksel/code](/data/derivative/tmQM-RDF-1Ksel/code). For brevity, the paths mentioned below (and within the scripts) are not the full paths, but follow the following convention:
 - Paths starting with tmQM/, /tmQMg, /tmQMg-L and /misc are to be interpreted as subpaths of [/data/raw/tmQMseries](/data/raw/tmQMseries)
 - Paths starting with pubChem/ are to be interpreted as subpaths of [/data/raw/pubChem](/data/raw/pubChem)
 - Paths starting with intermediate/tmQM/, intermediate/tmQMg/ and intermediate/tmQMg-L are to be interpreted as subpaths of [/data/derivative/tmQM-RDF/](/data/derivative/tmQM-RDF/))
 - Paths starting with intermediate/ are to be interpreted as subpaths of [/data/derivative/tmQM-RDF-1Ksel](/data/derivative/tmQM-RDF-1Ksel/)
 - Paths starting with data/ are to be interpreted as subpaths of [/data/derivative/tmQM-RDF-1Ksel/data/](/data/derivative/tmQM-RDF-1Ksel/data/)
 - The tag [version] can evaluate either to `earlyTM` or `lateTM`, depending on the nature of the desired selection

0. Run the script [0_download_tmQMg_outliers.py](/data/derivative/tmQM-RDF-1Ksel/code/0_download_tmQMg_outliers.py) (it will create the file `intermediate/outliers.txt`)
   > Downloads the file `outliers.txt` from the GitHub repository of tmQMg ([https://github.com/uiocompcat/tmQMg/](https://github.com/uiocompcat/tmQMg/)) and places it in `intermediate/`
1. Run the script [1_extract_info_from_tmc_for_selection.py](/data/derivative/tmQM-RDF-1Ksel/code/1_extract_info_from_tmc_for_selection.py) (it will create the files `intermediate/centres.csv` and `intermediate/details.csv`)
   > For each metal centre encountered in the tmQM-RDF dataset, it counts the number of complexes which possess that cenrte (`intermediate/cores.csv`). For each complex in the tmQM-RDF dataset, it summarises the metal core and the ligands (`intermediate/details.csv`).
2. Run the script [2_compute_ligand_occurrences.py](/data/derivative/tmQM-RDF-1Ksel/code/2_compute_ligand_occurrences.py) (it will create the files `intermediate/[version]/ligands.csv`)
   > For each ligand, it counts the number of occurrences within the complexes having one of the three requested metal centres, ignoring multiple appearances within a single complex (`intermediate/[version]/ligands.csv`).
3. Run the script [3_extract_selection.py](/data/derivative/tmQM-RDF-1Ksel/code/3_extract_selection.py) (it will create the file `data/[version]/1k_selection.txt`)
   > It creates the 1k selection according to the following scheme:
   >    1. First, a "seed" made of the 350 (lateTM)/1350 (earlyTM) most frequent ligands is extracted
   >    2. Then, a candidate set of ligands is built so that:
   >       - only TMCs whose metal centre is among the three requested centres in tmQM-RDF are allowed
   >       - only TMCs whose ligands are all found in the seed are allowed
   >    3. Finally, 1000 + 300 + 300 TMCs are sampled from the candidate set and then partitioned in order to form the train, validation and test partitions of the 1k selection
   >       (the relative proportions of the metal cores in the selection are artificially matched to those of the three most frequent cores in tmQM-RDF, sampling occurs with probability inversely proportional to complex size)
   
