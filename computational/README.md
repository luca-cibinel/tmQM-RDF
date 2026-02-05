This directory contains the code used to perform the computational experiments presented in the publication.

# Training phase
1. Mining a set of frequent patterns (see [pattern_mining](/computational/pattern_mining/))
    1. [Acquiring the graph dataset from tmQM-RDF](/computational/pattern_mining/0_prepare_1Ksel_for_pattern_mining.py)
    2. [Graph pattern mining](/computational/pattern_mining/1_fm-gpm-transactional-setting_processpool-keepmappingsinfiles_deterministicsampling.pl)
        - [Pattern post-processing](/computational/pattern_mining/2_pattern_postprocessing.py)
2. Determining a set of families of typical substructures (see [pattern_clustering](/computational/pattern_clustering))
    1. [Pattern filtering](/computational/pattern_clustering/1_process_candidate_patterns.py) 
    2. [Elementary feature computation](/computational/pattern_clustering/1_process_candidate_patterns.py)
    3. [Identification of structural families](/computational/pattern_clustering/3_cluster_patterns.py) *Note: s_{cos} and s_{DLG} must be computed beforehand:*
        - [Computing s_{cos}](/computational/pattern_clustering/2_i_compute_cossim.py)
        - [Computing s_{DLG}](/computational/pattern_clustering/2_ii_compute_delgsim.py)
    4. [Aggregated features computation](/computational/pattern_clustering/3_cluster_patterns.py)
3. Statistical learning (see [bayesian_network](/computational/bayesian_network/))
    1. Statistical learning (i.e., [structure learning](/computational/bayesian_network/1_bn_learn.R) and [parameter estimation](/computational/bayesian_network/2_bn_fit.py))

# Test phase
See [reconstruction](/computational/reconstruction).

1. [Extracting the set of training ligands](/computational/reconstruction/1_extract_ligands.py)
2. [Create molecular scaffolds](/computational/reconstruction/2_process_scaffolds.py)
3. [Create reconstructions](/computational/reconstruction/2_process_scaffolds.py)
3. Assess typicality:
    1. [Compute elementary features](/computational/reconstruction/3_compute_matches.py)
    2. [Compute aggregated features and feature probabilities](/computational/reconstruction/4_compute_scores.py)
    3. [Rank reconstruction](/computational/reconstruction/5_compute_rankings.py)
    4. [Visualise results](/computational/reconstruction/6_plot_results.py)
