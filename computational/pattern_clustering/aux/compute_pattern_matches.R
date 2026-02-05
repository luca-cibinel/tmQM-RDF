#' This code uses the virtuoso package as an interface with Virtuoso Open Source
#' in order to compute the number of matches of the patterns identified via
#' ../1_process_candidate_patterns.py on the 1k selection of tmQM-RDF
#' 
#' The patterns are stored in a query of the form
#'  SELECT *
#'  WHERE { ...pattern... }
#' and the matches are computed by turning each pattern into the query
#'  SELECT *
#'  FROM NAMED <...>
#'  ...
#'  FROM NAMED <...>
#'  WHERE { ...pattern... }
#'  
#'  The following variables have to be specified from the main Python file
#'  via rpy2.robjects.globalenv['varname'] = ...
#'    - tmQM.RDF.selection: path to the file that contains the dataset selection
#'    - pattern.dir: path to the directory that contains the processed patterns
#'    - results.dir: path to the directory where the computed matches have to be stored
#'    - vos.db.dir: path to the directory in which the (temporary) VOS database should be created
#'    - path.to.tmQM.RDF: path to the directory where tmQM-RDF is stored

source(vos.interface)

# Extract the given dataset selection
print(getwd())
print(tmQM.RDF.selection)
path.to.tmQM.RDF.selection <- tmQM.RDF.selection
tmQM.RDF.selection <- read.csv(
  path.to.tmQM.RDF.selection,
  row.names = 1
)[, "TMC"]

# Specify the base URI for the named graphs to be imported into the Virtuoso database
gname.uri <- "resource://integreat/tmQM-RDF/graph/"

# Specify the path to the directory containing the query representation of the patterns of interest
working.pattern.dir <- path.join(pattern.dir, "queries")

# Specify the name of the output file
out.file <- path.join(results.dir, "matches.csv")

# -------- START ====

vosh.cluster <- vosh.cluster.init(
  N.CORES,
  vos.db.dir, 
  path.to.tmQM.RDF, 
  tmQM.RDF.selection, 
  path.to.tmQM.RDF.selection, 
  gname.uri, 
  clean.start = T
)

vosh.cluster.import.dataset(vosh.cluster)

# Attempt pattern matching
pattern.sizes <- sort(as.numeric(
  list.files(working.pattern.dir)
))

results <- NULL

for(size in pattern.sizes){
  
  print(paste("Processing patterns of size", size))
  
  patterns <- list.files(path.join(pattern.dir, "queries", size)) # get list of pattern ids
  patterns <- sapply(patterns, function(p_fname){gsub(".txt$", "", p_fname)})
  
  results <- cbind(results, vosh.cluster.match.patterns(vosh.cluster, patterns, path.join(pattern.dir, "queries", size)))
}

print(paste("Writing results to", out.file))
write.csv(results, out.file)

# Kill VOS instance
vosh.cluster.stop(vosh.cluster)
