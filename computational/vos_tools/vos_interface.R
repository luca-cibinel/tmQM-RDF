library(virtuoso)
library(doParallel)

# DIRECTORY UTILS ====

#' path.join
#'
#' @param ... A char vector of paths to be joined
#'
#' @return A path that is the correct union of the provided paths
#' 
path.join <-  function(...){
  sep <- .Platform$file.sep
  pat <- paste0("\\", sep, "$")
  
  paths <- lapply(list(...), function(x){gsub(pat, "", x)})
  
  return(paste0(paths, collapse = sep))
}

#' prepare.temp.wd
#'
#' Initiates a temporary subdirectory in the specified root directory, with the possibility of redirecting
#' standard R output to a log.txt file in the new subdirectory via the sink() function.
#'
#' @param root.dir The root directory where the subdirectory should be created. Default: "."
#' @param sink.to.temp Boolean. Should R output be redirected to a file log.txt in the new subdirectory? Default: FALSE
#'
#' @returns The path to the temporary subdirectory
#' 
prepare.temp.wd <- function(root.dir = ".", sink.to.temp = F){
  temp.wd <- tempfile(tmpdir = root.dir)
  dir.create(temp.wd)
  
  if(sink.to.temp){
    file.create(path.join(temp.wd, "log.txt"))
    sink(path.join(temp.wd, "log.txt"), append = T)
  }
  
  return(temp.wd)
}

# DATABASE UTILS ====

#' vosh.init
#'
#' Configures and starts the VOS process and its handler
#'
#' @param db.dir Directory in which the VOS database will be created
#' @param path.to.dataset the path to the directory were the files to be imported are stored
#' @param dataset a vector containing the names of the files to be imported
#' @param gname.uri the base uri of the named graphs into which future data will be imported. Each named graph will be of the form <[gname.uri][file.name]>
#' @param gigs.ram Max. GB of memory VOS has access to. Default: 2
#' @param template Location of a virtuoso.ini template file. Default: virtuoso:::find_virtuoso_ini()
#' @param clean.start Boolean, should db.dir be created anew? Default: FALSE
#' @param connect Boolean, should a connection to VOS be already established? See vosh.connect. Default: FALSE
#'
#' @returns An object of class VOSH
#'
vosh.init <- function(
      db.dir,
      path.to.dataset,
      dataset,
      game.uri,
      gigs.ram = 2, 
      template = virtuoso:::find_virtuoso_ini(),
      clean.start = F,
      connect = F
    ){
  
  # Initialises the db.dir, if requested
  if(clean.start){
    if(dir.exists(db.dir)){
      unlink(db.dir, recursive = T)
    }
    dir.create(db.dir)
  }
  
  # Configures the VOS instance
  virtuoso.ini <- vos_configure(
    dirs_allowed = c(db.dir, path.to.dataset),
    db_dir = db.dir
  )
  
  # Start the VOS instance
  vos_start(virtuoso.ini)
  
  # Connects to the VOS instance (if required)
  vos.connection <- NULL
  if(connect){
    vos.connection <- vos_connect()
  }
  
  # Creates the VOS handler object
  vosh <- list(
    conn = vos.connection,
    db.dir = db.dir,
    dataset = dataset,
    path.to.dataset = path.to.dataset,
    gname.uri = gname.uri
  )
  
  class(vosh) <- "VOSH"
  
  return(vosh)
}

#' vosh.connect
#'
#' Establish a connection to VOS and link the resulting connection object to the provided VOS handler
#'
#' @param vosh An object of class VOSH to link to the established connection
#'
#' @returns The updated VOSH object
#' 
vosh.connect <- function(vosh){
  if(!is.null(vosh$conn)){
    return(vosh)
  }
  
  vos.connection <- vos_connect()
  
  vosh$conn <- vos.connection
  
  return(vosh)
}

#' vosh.import.dataset
#'
#' Imports a dataset of RDF graphs into the VOS database, assigning each graph to a unique
#' named graph.
#' 
#' Notice that the graphs are imported by copying each individual file to a local cache directory
#' (created by the function within the provided db.dir) and by removing it after the import, so that 
#' the cache always contains a single file at a time. The function virtuoso::vos_import is then used to
#' import the entire cache at once.
#' This is done to ensure that when multiple parallel VOS processes try to import the same data,
#' they don't interfere with each other. The normal procedure used by virtuoso::vos_import to
#' import specific files consists in copying the files to the VOS cache directory and then importing
#' the cache, similarly to what is done by this function, but in doing so, virtuoso::vos_import deletes 
#' any copy of the file that was already present in the cache. Since all the VOS processes share the
#' same global cache, this can likely cause the processes to lose access to the files stored in the
#' cache before they can import them.
#'
#' @param vosh An object of class VOSH that handles the current VOS instance
#' @param N.cores.load the number of cores to be dedicated to data loading. Default: 1.
#' @param f.ext The extension of the files to be imported. See vos_import for a list of supported files. Default: 'ttl'
#' 
vosh.import.dataset <- function(vosh, N.cores.load = 1, f.ext = "ttl"){
  # Create a local cache directory
  cache <- tempfile(tmpdir = vosh$db.dir)
  dir.create(cache)
  
  for(i in 1:length(vosh$dataset)){
    x <- vosh$dataset[i]
    f.name <- paste0(x, ".", f.ext)
    
    # Copy file [x].[f.ext] to the local cache
    file.copy(path.join(vosh$path.to.dataset, "graphs", f.name), cache)
    
    vos_import( # Import the local cache into the graph <[vosh$gname.uri][x]>
      vosh$conn, 
      wd = cache,
      glob = "*",
      graph = paste0(vosh$gname.uri, x),
      n_cores = N.cores.load
    )
    
    # Remove file from cache
    unlink(path.join(cache, f.name))
  }
  
  # Delete local cache
  unlink(cache, recursive = T)
}

## Parallel Computing ====

#' vosh.cluster.init
#'
#' Configures a cluster of foked parallel processes and starts the VOSH in each of them. The connection
#' to VOS is automatically established in each process.
#' Each process will operate inside its own temporary directory created inside the provided db.dir.
#'
#' @param N.proc The number of processes to start
#' @param db.dir Root directory in which the subdirectories of each child process will be created
#' @param path.to.dataset the path to the directory were the files to be imported are stored
#' @param dataset a vector containing the names of the files to be imported
#' @param gname.uri the base uri of the named graphs into which future data will be imported. Each named graph will be of the form <[gname.uri][file.name]>
#' @param gigs.ram Max. GB of memory VOS has access to. Default: 2
#' @param template Location of a virtuoso.ini template file. Default: virtuoso:::find_virtuoso_ini()
#' @param clean.start Boolean, should db.dir be created anew? Default: FALSE
#' 
#' @returns An object of class VOSHcluster
#'
vosh.cluster.init <-  function(
    N.proc,
    db.dir,
    path.to.dataset,
    dataset,
    game.uri,
    gigs.ram = 2, 
    template = virtuoso:::find_virtuoso_ini(),
    clean.start = F
){
  
  # Create and register PSOCK cluster
  cl <- makePSOCKcluster(N.proc)
  registerDoParallel(cl)
  
  # Ensure that VOS is running (facilitate safe initialisation of children VOS processes)
  vosh.init(db.dir, path.to.dataset, dataset, gname.uri, gigs.ram, clean.start = clean.start, connect = T)
  
  # Exports required variables to children processes
  clusterExport(cl, c(
      "db.dir",
      "path.to.dataset",
      "dataset",
      "gname.uri",
      "gigs.ram",
      "template",
      "vosh.init",
      "prepare.temp.wd",
      "path.join"
    ),
    envir = environment()
  )
  
  # Initiates VOS in each children process
  clusterEvalQ(cl, {
    # This code will be run inside the environment of the children processes
    library(virtuoso)
    
    # Prepare temporary wd inside db.dir
    temp.db.dir <- prepare.temp.wd(db.dir, TRUE)
    
    print("VOS process started in directory:")
    print(temp.db.dir)
    
    # Initialise VOSH and connect
    vosh <- vosh.init(temp.db.dir, path.to.dataset, dataset, gname.uri, gigs.ram, template, connect = T)
  })
  
  # Return the created cluster
  vosh.cluster <- list(
    cl = cl,
    db.dir = db.dir,
    dataset = dataset,
    path.to.dataset = path.to.dataset,
    gname.uri = gname.uri
  )
  
  class(vosh.cluster) <- "VOSHcluster"
  
  return(vosh.cluster)
}

#' vosh.cluster.import.dataset
#'
#' Imports a dataset of RDF graphs into the VOS database of each of the forked processes 
#' in the provided cluster. See vosh.import.dataset for details on the import procedure.
#'
#' @param cl An object of class VOSHcluster that contains the forked VOS processes
#' @param N.cores.load the number of cores to be dedicated to data loading. Default: 1.
#' @param f.ext The extension of the files to be imported. See vos_import for a list of supported files. Default: 'ttl'
#'
vosh.cluster.import.dataset <- function(vosh.cluster, N.cores.load = 1, f.ext = "ttl"){
  # Exports required variables to children processes
  clusterExport(vosh.cluster$cl, c(
      "N.cores.load",
      "f.ext",
      "vosh.import.dataset"
    ),
    envir = environment()
  )
  
  # Import dataset in each children process
  clusterEvalQ(vosh.cluster$cl, {
    # This code will be run inside the environment of the children processes
    vosh.import.dataset(vosh, N.cores.load, f.ext)
  })
}

#' vosh.cluster.stop
#'
#' Kills the VOS process and stops the parallel cluster.
#'
#' @param vosh.cluster The VOSH object to stop
#'
vosh.cluster.stop <- function(vosh.cluster){
  vos_kill()
  
  stopCluster(vosh.cluster$cl)
}

# PATTERN UTILS ====

#' vosh.match.pattern
#'
#' Matches a graph pattern against a dataset of RDF graphs.
#'
#' @param vosh An object of class VOSH that handler the current VOS instance
#' @param pattern.id The id of the pattern to match
#' @param pattern.dir The path to the directory where the desired pattern is stored
#'
#' @returns A list with two entries:
#'  - pat.id: The id of the matched pattern
#'  - res: a named vector whose entries are named after the graphs in the dataset for which at least one match
#'      is found. The value of each entry is the number of matches against the corresponding graph.
#'
vosh.match.pattern <- function(vosh, pattern.id, pattern.dir){
  p.query <- readLines(path.join(pattern.dir, paste0(pattern.id, ".txt"))) # get query from pattern file
  p.main.query <- paste(paste0(p.query[-1], collapse = "\n")) # get main query (remove "SELECT *")
  
  from.named <- paste(
    paste0(
      "FROM NAMED <", 
      paste0(vosh$gname.uri, vosh$dataset), # from CSD to URI
      ">"
    ), collapse = "\n"
  ) # prepare "FROM NAMED" statements
  
  effective.query <- paste(p.query[1], from.named, p.main.query, sep = "\n") # assemble effective query
  
  res <- vos_query(vosh$conn, effective.query) # run query
  res.out <-  NULL
  
  if(nrow(res) > 0){
    indices <- sapply(res$graph, function(guri){ gsub(vosh$gname.uri, "", guri) }) # get indices of matched graphs (from URIs to CSD codes)
    
    res.out <- table(indices)[unique(indices)]
    names(res.out) <- unique(indices)
  }
  
  list(pat.id = pattern.id, res = res.out)
}

#' vosh.match.pattern
#'
#' Matches several graph patterns against a dataset of RDF graphs.
#'
#' @param vosh An object of class VOSH that handler the current VOS instance
#' @param pattern.ids A charachter vector containing the ids of the patterns to match
#' @param pattern.dir The path to the directory where the desired pattern is stored
#'
#' @returns A length(vosh$dataset) x length(patterns) numeric matrix, whose rows are indexed by the names of
#'  the elements of the dataset and the columns are indexed by pattern ids. Each entry contains the number
#'  of matches of the associated pattern against the associated graph.
#'
vosh.match.patterns <- function(vosh, pattern.ids, pattern.dir){
  results <- matrix(0, length(vosh$dataset), length(pattern.ids))
  rownames(results) <- vosh$dataset
  colnames(results) <- pattern.ids
  
  for(p.id in pattern.ids){
    p.matches <- vosh.match.pattern(vosh, p.id, pattern.dir)
    
    if(!is.null(p.matches$res)){
      results[names(p.matches$res), p.id] <- p.matches$res
    }
  }
  
  return(results)
}

## Parallel Computing ====

#' vosh.cluster.match.pattern
#'
#' Matches several graph patterns against a dataset of RDF graphs.
#'
#' @param vosh An object of class VOSH that handler the current VOS instance
#' @param pattern.ids A charachter vector containing the ids of the patterns to match
#' @param pattern.dir The path to the directory where the desired pattern is stored
#'
#' @returns A length(vosh$dataset) x length(patterns) numeric matrix, whose rows are indexed by the names of
#'  the elements of the dataset and the columns are indexed by pattern ids. Each entry contains the number
#'  of matches of the associated pattern against the associated graph.
#'
vosh.cluster.match.patterns <- function(vosh.cluster, pattern.ids, pattern.dir){
  results <- matrix(0, length(vosh.cluster$dataset), length(pattern.ids))
  rownames(results) <- vosh.cluster$dataset
  colnames(results) <- pattern.ids
  
  # Computes matches (in parallel)
  partial.results.pll <- foreach(
    p.id = pattern.ids, 
    .inorder = FALSE, 
    .packages = "virtuoso",
    .export = "vosh.match.pattern",
    .noexport = "vosh"
  ) %dopar% {
    # This code will be run inside the environment of the children processes
    vosh.match.pattern(vosh, p.id, pattern.dir)
  }
  
  # Assemble partial results from parallel processes into a single local results matrix
  for(i in 1:length(partial.results.pll)){
    partial.res <- partial.results.pll[[i]]
    
    if(!is.null(partial.res$res)){
      results[names(partial.res$res), partial.res$pat.id] <- partial.res$res
    }
  }
  
  return(results)
}
