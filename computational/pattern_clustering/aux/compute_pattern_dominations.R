#' This code uses the virtuoso package as an interface with Virtuoso Open Source
#' in order to compute the domination relationships among the patterns identified via
#' ../1_process_candidate_patterns.py.
#' 
#' We say that p2 'dominates' p1 if
#' 
#'  |{matches of p2 in G}| = 0 ==> |{matches of p1 in G}| = 0 for each G
#' In practical terms, this is verified by checking whether p2 matches p1 when p1 is
#' represented as an RDF graph (with variables turned into URIs).
#' 
#' For simplicity, domination relationships are only computed between patterns whose size
#' differs by 1.
#' 
#' The patterns are stored in a query of the form
#'  SELECT *
#'  WHERE { GRAPH ?graph {...pattern... } }
#' and the matches are computed by turning each pattern into the query
#'  SELECT *
#'  FROM NAMED <...>
#'  ...
#'  FROM NAMED <...>
#'  WHERE { GRAPH ?graph {...pattern... } }
#'  
#' The domination relations found using this procedure are stored in the {local.dataset.dir}/dominations.yml file, listing:
#' - the sizes considered in the computation (for which dominating patterns are computed)
#' - for each pattern, a list:
#'  - the first entry is the pattern size
#'  - the remaining entries are the dominating patterns

path.join <-  function(...){
  sep <- .Platform$file.sep
  pat <- paste0("\\", sep, "$")
  
  paths <- lapply(list(...), function(x){gsub(pat, "", x)})
  
  return(paste0(paths, collapse = sep))
}

library(virtuoso)
library(yaml)

# Specify the base URI for the named graphs to be imported into the Virtuoso database
gname.uri <- "resource://integreat/pattern/size/"

# Specify the file to which the results will be written+
out.file <- path.join(local.dataset.dir, "dominations.yml")

pattern.sizes <- sort(as.numeric(
    list.files(path.join(pattern.dir, "queries"))
  ))

# -------- START ====

# Create the vos database directory with virtuoso::vos_configure)
virtuoso.ini <- vos_configure(db_dir = vos.db.dir)

# Start the VOS instance
vos_start(virtuoso.ini)

# Connect to the database
vos.connection <- vos_connect()

# Import patterns (and prepare results list)
results <- list()

results[["sizes"]] <- sapply(pattern.sizes[-1], toString)

for(size in pattern.sizes[-1]){
  print(paste("Importing size", size))
  gname.uri.loc <- paste0(gname.uri, size, "/")
  
  patterns <- list.files(path.join(pattern.dir, "rdf", size)) # get list of pattern ids
  patterns <- sapply(patterns, function(p_fname){gsub(".nt$", "", p_fname)})
  
  pb = txtProgressBar(min = 1, max = length(patterns))
  for(i in 1:length(patterns)){
    setTxtProgressBar(pb, i)
    
    p.id <- patterns[i]
    
    vos_import( # Import the file [tmc].ttl into the named graph '<[gname.uri][size]/[p.id]>'
      vos.connection, 
      files = path.join(pattern.dir, "rdf", size, paste0(p.id, ".nt")),
      graph = paste0(gname.uri.loc, p.id)
    )
    
    results[[p.id]] <- c(toString(size)) # prepare results list
  }
  close(pb)
  
  # Sanity check: ensure import was successful
  graphs <- vos_list_graphs(vos.connection)$g
  for(p.id in patterns){
    if(!(paste0(gname.uri.loc, p.id) %in% graphs)){
      stop(paste("Failed import! Size", size, ":", p.id))
    }
  }
}

tot.n.patterns <- sum(sapply(
    pattern.sizes[-length(pattern.sizes)], function(s){length(list.files(path.join(pattern.dir, "rdf", s)))}
  ))

# Compute dominations

progress <- 0 # monitor completion percentage
processed <- 0
batch.times <- NULL # monitor average batch processing time
tic.global <- proc.time() # record starting time
for(size in pattern.sizes[-1]){
  print(paste("Processing size", size))
  gname.uri.loc <- paste0(gname.uri, size, "/")
  
  patterns.d <- list.files(path.join(pattern.dir, "rdf", size - 1)) # get list of (dominating) pattern ids
  patterns.d <- sapply(patterns.d, function(p_fname){gsub(".nt$", "", p_fname)})
  
  patterns <- list.files(path.join(pattern.dir, "rdf", size)) # get list of pattern ids
  patterns <- sapply(patterns, function(p_fname){gsub(".nt$", "", p_fname)})
  
  for(p.num in 1:length(patterns.d)){
    p.id <- patterns.d[p.num] # get pattern id
    
    p.query <- readLines(path.join(pattern.dir, "queries", size - 1, paste0(p.id, ".txt"))) # get query from pattern file
    p.main.query <- paste(paste0(p.query[-1], collapse = "\n")) # get main query (remove "SELECT *")
    
    tic <- proc.time() # monitor current batch processing time
    
    from.named <- paste(
      paste0(
        "FROM NAMED <", 
        paste0(gname.uri.loc, patterns), # from CSD to URI
        ">"
      ), collapse = "\n"
    ) # prepare "FROM NAMED" statements
    
    effective.query <- paste(p.query[1], from.named, p.main.query, sep = "\n") # assemble effective query
    
    res <- vos_query(vos.connection, effective.query) # run query
    
    if(nrow(res) > 0){
      indices <- sapply(res$graph, function(guri){ gsub(gname.uri.loc, "", guri) }) # get indices of matched patterns (from URIs to IDs)
      indices <- unique(indices)
      
      for(q.id in indices){
        results[[q.id]] <- c(results[[q.id]], p.id)
      }
    }
    batch.times <- c(batch.times, unname((proc.time() - tic)["elapsed"])) # monitor current batch processing time
    
    processed <- processed + 1
    
    if(100*processed/tot.n.patterns >= progress){ # output progress and expected completion time
      progress <- progress + 1
      
      elapsed <- unname((proc.time() - tic.global)["elapsed"])
      expected.remaining <- (tot.n.patterns - processed)*mean(batch.times)
      
      print(paste0(progress, "% - ", round(elapsed/60, 3), " || ", round(expected.remaining/60, 3)))
    }
  }
}

write_yaml(results, out.file)

# Kill VOS instance
vos_kill()
