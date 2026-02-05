This directory contains the code related to the clustering and features computation step.

- [1_process_candidate_patterns.py](/computational/pattern_clustering/1_process_candidate_patterns.py): this script takes care of two fundamental steps, namely pattern filtering and elementary features computation, for one of the available datasets (`PATTERN_DATASET_KEYWORD = "lateTM"` or `PATTERN_DATASET_KEYWORD = "earlyTM"`)
  - As an initial step, frequent patterns (stored in [../pattern_mining/results](/computational/pattern_mining/results)) are copied and stored in the "[results](/computational/pattern_clustering/results)" folder
  - The domination relationship, on which filtering is based, is computed by calling the R script [aux/compute_pattern_dominations.R](/computational/pattern_clustering/aux/compute_pattern_dominations.R)
    - Dominations are stored in the same "results" folder as the patterns, in the file `dominations.yml`. This file is organised as follows:
        - `sizes`: the list of pattern sizes (i.e., number of triples) for which data is vailable inside the file
        - Each entry other than `size` is the id of a pattern $p$, to which a list of values is associated:
          - The first element of the list is the size $|p|$ of the pattern
          - Every other element is the id of a pattern of size $|p| - 1$ which dominates $p$
  - Patterns are filtered as mentioned in the publication, by removing from the results directory the patterns that do not satisfy the filters
  - Finally, the elementary graph features (i.e., pattern matches) are computed by invoking the R script [aux/compute_pattern_matches.R](/computational/pattern_clustering/aux/compute_pattern_matches.R)
- [2_i_compute_cossim.py](/computational/pattern_clustering/2_i_compute_cossim.py): computes the cosine similarity among the filtered patterns; it allows to specify options such as the dataset to process (`PATTERN_DATASET_KEYWORD = "lateTM"` or `PATTERN_DATASET_KEYWORD = "earlyTM"`) and whether to use the $s_{cos;\,p}$ (`FEATURE_TYPE = "proxy"`) or the $s_{cos;\,s}$ (`FEATURE_TYPE = "semantic"`) configuration. The distance matrix (1 - similarity) is computed and stored in the "[intermediate](/computational/pattern_clustering/intermediate)" folder.
- [2_ii_compute_delgsim.py](/computational/pattern_clustering/2_ii_compute_delgsim.py): computes the $s_{DLG}$ similarity among the filtered patterns. The similarity routines are implemented in [delgsim.py](/computational/pattern_clustering/delgsim.py). The distance matrix (1 - similarity) is computed and stored in the "intermediate" folder. As this computation can be particularly time consuming, the settings are slightly more sophisticated:
    - `PATTERN_DATASET_KEYWORD`: string. The dataset to use (`lateTM` or `earlyTM`).
    - `USE_LEARNED_WEIGHTS`: boolean. Should the script use the $s_{DELG;\,n}$ (`False`) or the $s_{DELG;\,l}$ (`True`) configuration? The type of weighting scheme influences where results and intermediate quantities and utility files will be stored. For $s_{DELG;\,n}$ the folder will be called "fixed_weights", for $s_{DELG;\,l}$, "estimated_weights".
    - `DO_REPRESENTATION_COMPUTATION`: boolean. Should the script compute an appropriate machine-readable representation of the patterns (`True`) or should it be loaded from an external file (`False`)? This representation is saved, using the `pickle` module, in the "intermediate" folder. If `False`, the representation will be loaded from the same folder.
    - `DO_WEIGHT_COMPUTATION`: boolean. Should weights (either fixed or estimated) be initialised from scratch (`True`) or should they be read from an external file (`False`)? If `True`, weights are computed and stored, using pickle, in the "intermediate" folder. If `False` the pickle file in the correct "intermediate" folder is loaded.
    - `DO_DELGSIM_COMPUTATION`: boolean. Should the computation of $s_{DELG}$ (in the chosen configuration) be carried out? The computed distance matrix (1 - similarity) will be stored in the "results" folder.
    - `DRAW_PLOTS`: boolean (diagnostic flag). Should a histogram of the computed similarity values be plotted?
- [3_cluster_patterns.py](/computational/pattern_clustering/3_cluster_patterns.py): clusters the patterns on the basis of the computed similarity values. The script will automatically detect the configurations for which similarity values are available and ask for confirmation of which of them to use. The script will then compute the clusters (stored in `cls_info.csv`), the preliminary undirected edges in the initial DAG $B_0$ (`cls_edges.csv`) and finally the aggreggated features (`cls_matches.csv`, inside a folder that corresponds to the aggregation function). All of this is stored inside the "results" folder, inside a "by_metric" folder.

# Outputs
- Intermediate
    - Similarity values (including intermediate values needed for their computation)
- Results
    - Elementary features, clusters, aggregated features, initial BN DAG $B_0$
