#' This code uses the virtuoso package as an interface with Virtuoso Open Source
#' in order to compute the number of matches of the patterns identified via
#' ../identify_candidate_patterns.py on the 1k selection of tmQM-RDF
#' 
#' The patterns are stored in a query of the form
#'  SELECT *
#'  WHERE { GRAPH ?graph {...pattern...}}
#'  
#'  The following variables have to be specified from the main Python file
#'  via rpy2.robjects.globalenv['varname'] = ...
#'    - tmQM.RDF.selection: path to the file that contains the dataset selection
#'    - pattern.dir: path to the directory that contains the processed patterns
#'    - results.dir: path to the directory where the computed matches have to be stored
#'    - vos.db.dir: path to the directory in which the (temporary) VOS database should be created
#'    - path.to.tmQM.RDF: path to the directory where tmQM-RDF is stored

path.join <-  function(...){
  sep <- .Platform$file.sep
  pat <- paste0("\\", sep, "$")
  
  paths <- lapply(list(...), function(x){gsub(pat, "", x)})
  
  return(paste0(paths, collapse = sep))
}

DEBUG <- F

if(!DEBUG){
  source(path.join(ROOT.DIR, "computational", "vos_tools", "vos_interface.R"))
  
  # RECON.FOLDER <- path.join(
  #   ROOT.DIR, 
  #   "computational", 
  #   "reconstruction", 
  #   "intermediate", 
  #   "reconstructions", 
  #   "rdf", 
  #   TM.MODE
  # )
  
  AVAILABLE.RECONSTRUCTIONS <- list.dirs(RECON.FOLDER, recursive = F)
  
  PATTERN.DIR <- path.join(
    ROOT.DIR, 
    "computational", 
    "pattern_clustering", 
    "results", 
    PATTERN.DATASET.TAG, 
    "patterns"
  )
  
  # N.PROC <- 8
}else{
  source("/Users/lucaci/Desktop/InteGreat_p5/paper_1/computational/substructures/vos_tools/vos_interface.R")
  WD0 <- "/Users/lucaci/Desktop/InteGreat_p5/paper_1/computational/substructures/bayes/inference/experiments/recon_experiments"
  
  available.reconstructions <- paste0(
    "/Users/lucaci/Desktop/InteGreat_p5/paper_1/computational/substructures/bayes/inference/recon/recon", 
    c(0, 1)
  )
  
  pattern.dir <- "/Users/lucaci/Desktop/InteGreat_p5/paper_1/computational/substructures/processed_pattern_data/latmod_large__s_10_12/2025-08-25/patterns"
  
  if(dir.exists(path.join(WD0, "recon_matches"))){
    unlink(path.join(WD0, "recon_matches"), recursive = T)
  }
  dir.create(path.join(WD0, "recon_matches"))
  dir.create(path.join(WD0, "recon_matches", "recon0"))
  dir.create(path.join(WD0, "recon_matches", "recon1"))
  
  sanity.check.dir <- "/Users/lucaci/Desktop/InteGreat_p5/paper_1/computational/substructures/bayes/inference/recon_matches"
  
  N.proc <- 8
}

VOS.DB.DIR <- path.join(
  ROOT.DIR, 
  "computational", 
  "reconstruction", 
  "temp", 
  "vos_db"
)
VOS.INI.TEMPLATE <- path.join(ROOT.DIR, "computational", "vos_tools", "virtuoso.ini")

processed <- NULL
for(recon in AVAILABLE.RECONSTRUCTIONS){
  tic <- proc.time()
  
  # Prepare results directory
  recon.name <- tail(strsplit(recon, .Platform$file.sep)[[1]], n = 1)
  results.dir <- recon
  
  processed <- c(processed, recon)
  
  # dir.create(results.dir)
  if("matches.csv" %in% list.files(results.dir)){
    next
  }
  
  print(paste("Processing reconstruction:", recon.name))
  
  # Extract the given dataset selection
  recon.selection <- read.csv(
    path.join(recon, "selection.csv"),
    row.names = 1
  )[, "TMC"]
  
  # Specify the base URI for the named graphs to be imported into the Virtuoso database
  gname.uri <- "resource://integreat/tmQM-RDF/graph/"
  
  # Specify the path to the directory containing the query representation of the patterns of interest
  # (specified from main Python file via rpy2.robjects.globalenv['pattern.dir'] = ...)
  # pattern.dir <- "/Users/lucaci/Desktop/InteGreat_p5/pattern mining (repo)/utils/temp_/temp_patterns/queries/"
  working.pattern.dir <- path.join(PATTERN.DIR, "queries")
  
  # Specify the name of the output file
  out.file <- path.join(results.dir, "matches.csv")
  
  # -------- START ====
  
  # Specify the .ini file of the desired database (created with virtuoso::vos_configure)
  vosh.cluster <- vosh.cluster.init(
    N.PROC, 
    VOS.DB.DIR, 
    recon, 
    recon.selection, 
    gname.uri, 
    template = VOS.INI.TEMPLATE, 
    clean.start = T
  )
  cl <- vosh.cluster$cl
  
  print("  vosh.cluster initialised")
  
  # Import the data
  vosh.cluster.import.dataset(vosh.cluster)
  print("  Data imported")
  
  # Attempt pattern matching
  pattern.sizes <- sort(as.numeric(
    list.files(working.pattern.dir)
  ))
  
  results <- NULL
  
  for(size in pattern.sizes){
    
    loc.pattern.dir <- path.join(PATTERN.DIR, "queries", size)
    patterns <- list.files(loc.pattern.dir) # get list of pattern ids
    patterns <- sapply(patterns, function(p_fname){gsub(".txt$", "", p_fname)})
    
    results <- cbind(results, vosh.cluster.match.patterns(vosh.cluster, patterns, loc.pattern.dir))
  }
  
  print("  Matches computed")
  
  # Save output
  write.csv(results, out.file)
  
  if(DEBUG){
    ground.truth <- read.csv(path.join(sanity.check.dir, recon.name, "matches.csv"))
    to.test <- read.csv(out.file)
    
    print(paste("  Sanity check:", all(ground.truth[,-1] == to.test[,-1])))
  }
  
  # Kill VOS instance
  vosh.cluster.stop(vosh.cluster)
  
  # Timing
  
  toc <- proc.time()
  
  print("  Iteration completed in:")
  print(toc - tic)
  
  print("")
  
  print("  ---")
  print(paste("  To process:", length(setdiff(AVAILABLE.RECONSTRUCTIONS, processed)), "/", length(AVAILABLE.RECONSTRUCTIONS)))
  
  print("")
  print("")
}
