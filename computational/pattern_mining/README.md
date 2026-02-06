This directory contains the code and the utilities needed to mine frequent patterns from [the two 1K selections of tmQM-RDF](/data/derivative/tmQM-RDF-1Ksel).

- [0_prepare_1Ksel_for_pattern_mining.py](/computational/pattern_mining/0_prepare_1Ksel_for_pattern_mining.py): this script prepares the RDF graphs in the 1K selections for pattern mining. Pattern mining requires RDF graphs serialised using the ntriples format and no blank nodes (hence skolemisation is performed). Also, only the triples whose predicates are explicitly allowed, as in the reference paper, are retained. Notice that this predicate filtering is performed in two steps: this scripts produces RDF graphs that only contain triples with the following predicates:
    - http://www.w3.org/1999/02/22-rdf-syntax-ns#type,
    - resource://integreat/p5/atomic/atom/isAtom,
    - resource://integreat/p5/atomic/structure/b,
    - resource://integreat/p5/complex/TMC/hasLigand,
    - resource://integreat/p5/complex/TMC/hasMetalCentre,
    - resource://integreat/p5/ligand/bond/hasBindingAtom,
    - resource://integreat/p5/ligand/centre/hasAtom,
    - resource://integreat/p5/ligand/centre/isMetalCentre,
    - resource://integreat/p5/ligand/ligand/hasAtom,
    - resource://integreat/p5/ligand/ligand/isLigand,
    - resource://integreat/p5/ligand/structure/bLc,
    - resource://integreat/p5/ligand/structure/bLl
    - (The predicates http://www.w3.org/1999/02/22-rdf-syntax-ns#type, resource://integreat/p5/ligand/centre/hasAtom and resource://integreat/p5/ligand/ligand/hasAtom are removed in the next step.)
- [1_fm-gpm-transactional-setting_processpool-keepmappingsinfiles_deterministicsampling.pl](/computational/pattern_mining/1_fm-gpm-transactional-setting_processpool-keepmappingsinfiles_deterministicsampling.pl): **Author: [Basil Ell](https://ekvv.uni-bielefeld.de/pers_publ/publ/PersonDetail.jsp?personId=79220264&lang=EN)**; this perl script performs the actual pattern mining. Among the available configs, a predicate blacklist filters out the forbidden predicates not removed in the previous step. An additional auxiliary script is available in support of this script:
    - [generate_config.pl](/computational/pattern_mining/generate_config.pl): this script generates a config file (in the configs/ folder) that can be fed to the pattern mining script as `perl 1_fm-gpm-transactional-setting_processpool-keepmappingsinfiles_deterministicsampling.pl config/config_file`.
- [2_pattern_postprocessing.py](/computational/pattern_mining/2_pattern_postprocessing.py): as the pattern mining script produces very large output files, which not only include the patterns but also their matches against the training data, it is possible to remove these additional matches (if so desired) by running this script.

# Outputs
- Intermediate:
    - local data (i.e., .nt files containing the preprocessed RDF graphs from which patterns are to be mined)
    - random numbers (i.e., files containing sequences of randomly generated numbers used for reproducibility)
- Results:
    - Mined frequent patterns. Each result directory (one per dataset selection) is organised as follows:
        - *.dat.gz: mined patterns files. These files contain the data of the mined patterns and can be easily read and managed using the code in [patterns.py](/computational/pattern_mining/patterns.py)
        - patterns-[i]/: folders containing visual representations of the mined patterns     

# Notes for reproducibility:
- If you wish to replicate _exactly_ the results found in this directory you need to:
    1. Modify the entry `max_pattern_size` in the corresponding .cfg file in [configs/](/computational/pattern_mining/configs) to 8 (for *earlyTM*) or 9 (for *lateTM*). Leave the other settings untouched.
    2. Run the script [1_fm-gpm-transactional-setting_processpool-keepmappingsinfiles_deterministicsampling.pl](/computational/pattern_mining/1_fm-gpm-transactional-setting_processpool-keepmappingsinfiles_deterministicsampling.pl) using the modified config file. This will mine patterns up to size 8/9. The program should use the random numbers file in [intermediate/random_numbers](/computational/pattern_mining/intermediate/random_numbers) and it should *NOT* populate new random numbers files.
    3. Modify again the config file, setting `max_pattern_size` to 12
    4. Run the script [1_fm-gpm-transactional-setting_processpool-keepmappingsinfiles_deterministicsampling.pl](/computational/pattern_mining/1_fm-gpm-transactional-setting_processpool-keepmappingsinfiles_deterministicsampling.pl) again. Once again, it should use the random number files in [intermediate/random_numbers](/computational/pattern_mining/intermediate/random_numbers).
    5. Run the script [2_pattern_postprocessing.py](/computational/pattern_mining/2_pattern_postprocessing.py).
