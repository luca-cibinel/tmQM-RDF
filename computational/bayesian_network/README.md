This directory contains the code needed to estimate the DAG and the parameters of a BN which models the joint distribution of the aggregated features.

- [1_bn_learn.R](/computational/bayesian_network/1_bn_learn.R): this script estimates the DAG of a BN using the hill climbing algorithm from the `bnlearn` package. The chosen dataset (`TM.MODE <- "lateTM"` or `TM.MODE <- "earlyTM"`) has to be specified. Also, the path to the root of the repository must be set (`ROOT.DIR <- ...`). The learned graphs are stored in the form of .gml graphs in the "[intermediate](/computational/bayesian_network/intermediate)" folder.
- [2_bn_fit.py](/computational/bayesian_network/2_bn_fit.py): this script estimates the parameters of the BN, on the basis of the learned DAGs, using the pgmpy package. The available DAGs are automatically detected (confirmation will be asked). The fitted BNs, in the form of pickle files, are stored in the "[results](/computational/bayesian_network/results)" folder. The number of cores to be used can be specified via the variable `SETTINGS["bn_fit_settings"]["n_jobs"]`.

# Outputs
- Intermediate:
    - DAGs
- Results:
    - Fitted BNs   
