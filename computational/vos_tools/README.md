This directory contains a script that serves as an easy to use interface to the R package "virtuoso" (Boettiger, 2021).

In particular, the script [vos_interface.R](/computational/vos_tools/vos_interface.R) offers utilities that summarise the procedures to initialise and connect to the VOS database, as well as utilities that directly allow for the for the computations of pattern matches against RDF graphs, by internally handling the conversion to traditional SPARQL queries. In addition, parallelisation via the "doParallel" package (Corporation 2022) is allowed.

# References
Boettiger C (2021). _virtuoso: Interface to 'Virtuoso' using 'ODBC'_. R package version 0.1.8, <https://CRAN.R-project.org/package=virtuoso>.
Corporation M, Weston S (2022). _doParallel: Foreach Parallel Adaptor for the 'parallel' Package_. R package version 1.0.17, <https://CRAN.R-project.org/package=doParallel>.

