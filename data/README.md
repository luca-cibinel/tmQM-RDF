This directory contains all the code needed to download, manage and process the data used in the creation of tmQM-RDF. The creation of tmQM-RDF requires the avilability of the tmQM dataset series(see the [main page](/README.md)), hence [appropriate code](/data/raw/tmQMseries/download_data.py) is provided in order to download and store all the required material. The same applies to other prerequisite material (see below).

The derivative dataset tmQM-RDF and the two selections (_lateTM_ and _earlyTM_) are created using the code indexed below.

# Directory content
- [**raw**](/data/raw): original datasets retrieved from external sources (needed for further processing/computation)
    - [tmQM series](/data/raw/tmQMseries): the original chemical datasets on which tmQM-RDF is based
    - [pubChem](/data/raw/pubChem): the csv version of the pubChem periodic table
- [**derivative**](/data/derivative): datasets derived from the raw data
    - [tmQM-RDF](/data/derivative/tmQM-RDF): the tmQM-RDF knowledge graph
        - [code](/data/derivative/tmQM-RDF/code): the code used to derive tmQM-RDF from the tmQM dataset series
        - [data](/data/derivative/tmQM-RDF/data): a local (compressed) copy of tmQM-RDF (**NOTE: if you are only interested in the data, it is reccomended that you use the provided [download links](/data/raw/tmQMseries/download_data.py)**)
    - [tmQM-RDF-1Ksel](/data/derivative/tmQM-RDF-1Ksel): two selections containing 1000 training TMCs each (plus 600 validation/test TMCs)
        - [code](/data/derivative/tmQM-RDF-1Ksel/code): the code used to sample the two selections from the general tmQM-RDF TMC population
        - [data](/data/derivative/tmQM-RDF-1Ksel/data): the two selections, provided as a list of TMCs to be extracted from tmQM-RDF using your prefed managment method 
