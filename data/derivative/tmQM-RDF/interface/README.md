This directory contains code that can be used to conveniently manage the tmQM-RDF dataset available in the [data](/data/) folder

# Archive utils
As the copy of tmQM-RDF in this repository is a compressed archive, the file [archive_utils/extract_tmQM_RDF_archive.py](/data/derivative/tmQM-RDF/interface/archive_utils/extract_tmQM_RDF_archive.py) decompressess and organises the archive. The desired version of tmQM-RDF is spectified at line 23 in the script (notice that only the _v2025dev_ version is available through this channel, please refer to the [main page](/README.md)).

# tmQM-RDF interface
Appropriate python scripts are provided to grant easy access to the data in tmQM-RDF. Each version of the dataset has its own specific interface (packed together with the data itself, available from the [main page](/README.md)), so please refer to the correct version for documentation:
  - [v2025dev](/data/derivative/tmQM-RDF/interface/v2025dev)
  - [v1.0](/data/derivative/tmQM-RDF/interface/v1.0)
