#' Update ROOT.DIR (lines 33 and 35) with your own project root directory (location of '.prj_root')
#' This script is meant to be launched via 'RScript 1_bn_learn.R [i]' where i = 1 loads the lateTM dataset selection
#' and i = 2 loads the earlyTM selection
#' Also remember to update the number of cores to be used (lines 54 and 56) to adapt to your own use case

set.seed(863269655)
rm(list = ls())

DEBUG <- F
USE.PCALG <- F

path.join <-  function(...){
  sep <- .Platform$file.sep
  pat <- paste0("\\", sep, "$")
  
  paths <- lapply(list(...), function(x){gsub(pat, "", x)})
  
  return(paste0(paths, collapse = sep))
}

extract.path.tail <- function(path, tail = 3){
  split.path <- strsplit(path, .Platform$file.sep)[[1]]
  l <- length(split.path)
  
  path.tail <- NULL
  for(i in tail:1){
    path.tail <- c(path.tail, split.path[l - i + 1])
  }
  
  return(do.call(path.join, as.list(path.tail)))
}

# Header ====

if(DEBUG){
  ROOT.DIR <- "."
}else{
  ROOT.DIR <- "."
}

if(USE.PCALG){
  library(pcalg)
}else{
  library(bnlearn)
}

library(igraph)

library(foreach)
library(parallel)
library(doParallel)
library(doRNG)

if(DEBUG){
  registerDoParallel(4)
}else{
  registerDoParallel(8)
}

registerDoRNG(863269655)

alpha <- 0.05

TM.MODE <- c("lateTM", "earlyTM")[as.integer(commandArgs(trailingOnly=T)[1])]#"lateTM" # One of earlyTM / lateTM 
PATTERN.DATASET.TAG <- paste0(TM.MODE, "-latmod-s_10_12")
CLS.MATCHES.FILE.NAME <- "cls_matches.csv"

print(TM.MODE)

# Files/directories manipulation ====
DATASET.SELECTION <- path.join(ROOT.DIR, "data", "derivative", "tmQM-RDF-1Ksel", "data", TM.MODE, "1k_selection.csv")
PATTERN.DATASET.ROOT <- path.join(
  ROOT.DIR, 
  "computational", 
  "pattern_clustering", 
  "results",
  PATTERN.DATASET.TAG,
  "by_metric"
)
RESULTS.DIR.ROOT <- path.join(
  ROOT.DIR, 
  "computational", 
  "bayesian_network", 
  "intermediate",
  PATTERN.DATASET.TAG
)

dataset.selection <- read.csv(DATASET.SELECTION, header = T, row.names = 1)

work.dirs <- NULL
for(wd in list.dirs(PATTERN.DATASET.ROOT, recursive = T)){
  if(CLS.MATCHES.FILE.NAME %in% list.files(wd)){
    work.dirs <- c(work.dirs, wd)
    
    local.results.dir <- path.join(RESULTS.DIR.ROOT, extract.path.tail(wd))
    if(!dir.exists(local.results.dir)){
      dir.create(local.results.dir, recursive = T)
    }
  }
}

# Utility ====
sample.dag <- function(e.list, maxn = Inf){
  nodes <- unique(as.vector(e.list))
  
  neighborhoods <- sapply(nodes, function(n){sum(e.list == n)}) # get neighborhood sizes
  crowded.nodes <- nodes[neighborhoods > maxn] # Find nodes with neighborhood sizes greater than maxn
  
  while(length(crowded.nodes) > 0){
    x <- crowded.nodes[1]
    ne.x <- (1:nrow(e.list))[apply(e.list, 1, function(e){x %in% e})] # Find edges that involve x
    
    edges.to.remove <- sample(ne.x, length(ne.x) - maxn, replace = F) # Choose edges to remove
    
    e.list <- e.list[-edges.to.remove, ] # Remove edges
    
    neighborhoods <- sapply(nodes, function(n){sum(e.list == n)}) # Recompute crowded nodes
    crowded.nodes <- nodes[neighborhoods > maxn]
    
    # print(neighborhoods[neighborhoods > maxn])
  }
  
  perm <- sample(1:max(nodes))
  
  oriented.e.list <- NULL
  
  for(i in 1:nrow(e.list)){
    if(perm[e.list[i,1]] < perm[e.list[i,2]]){
      oriented.e.list <- rbind(oriented.e.list, c(e.list[i,1], e.list[i,2]))
    }else{
      oriented.e.list <- rbind(oriented.e.list, c(e.list[i,2], e.list[i,1]))
    }
  }
  
  return(oriented.e.list)
}

# Main cycle ====
fitted.results <- foreach(WD = work.dirs) %dopar% { # use .errorhandling="pass" to continue in case of errors
  local.results.dir <- path.join(RESULTS.DIR.ROOT, extract.path.tail(WD))
  
  file.create(path.join(local.results.dir, "log.txt"))
  sink(path.join(local.results.dir, "log.txt"), append = T)

  print(paste("WD =", WD))
  print("Using doRNG")
  
  MAX.P <- floor(log(.Machine$integer.max, 2)) - 3
  print(paste("MAX.P =", MAX.P))
  
  data.full <- read.csv(path.join(WD, "cls_matches.csv"), header = TRUE, row.names = 1)
  if(DEBUG){
    data <- data.full[dataset.selection$TMC[dataset.selection$partition == "train"], 1:10]
  }else{
    data <- data.full[dataset.selection$TMC[dataset.selection$partition == "train"], ]
  }
  colnames(data) <- 1:ncol(data)
  
  print(paste("   Nodes:", ncol(data)))
  nodes <- ncol(data)
  
  suff.stats <- list(dm = as.matrix(data), adaptDF = FALSE)
  
  scaffold <- matrix(0, ncol(data), ncol(data))
  # cls_edges.csv contains node ids in python enumeration (Python starts at 0, R starts at 1)
  scaffold.whitelist <- 1 + as.matrix(read.csv(
    normalizePath(path.join(WD, "..", "cls_edges.csv")), 
    row.names = 1
  ))
  if(DEBUG){
    scaffold.whitelist <- rbind(c(1,2), c(3,4))
  }
  
  scaffold[scaffold.whitelist] <- 1
  scaffold <- sign(scaffold + t(scaffold))
  
  remove.isolated <- F
  
  # Main cycle ====
  print(paste("   Fitting alpha =", alpha))
  
  ## Fit ====
  if(USE.PCALG){
    fit.pc <- pc(
      suff.stats,
      binCItest,
      alpha = alpha,
      p = ncol(data),
      skel.method = "stable.fast",
      fixedEdges = scaffold
    )
  
    cpdag <- udag2pdagSpecial(fit.pc)
    status <- cpdag$status
    fit.pc <- cpdag$pcObj
    
    g <- graph_from_graphnel(fit.pc@graph)
    graph_attr(g.simple, "status") <- status
  }else{
    for(v in colnames(data)){
      data[, v] <- as.factor(data[, v])
    }
    
    scaffold.whitelist <- sample.dag(scaffold.whitelist, maxn = MAX.P)
    scaffold.whitelist <- data.frame("from" = scaffold.whitelist[, 1], "to" = scaffold.whitelist[, 2])
    
    scaffold.dag <- empty.graph(sapply(1:ncol(data), toString))
    arcs(scaffold.dag) <- apply(scaffold.whitelist, c(1, 2), toString)

    print("  All ready for fit")

    fit.hc <- hc(
      data,
      start = scaffold.dag,
      maxp = MAX.P
    )
    
    g <- as.igraph(fit.hc)
  }
  
  #toc <- proc.time() - tic
  #toc <- round(unname(toc["elapsed"]), 4)
  
  print("  Done!")
  
  ## Output ====
  
  ### Graph conversion ====
  
  ### Graphical formatting ====
  print("Formatting graph for output...")
  for(k in 1:ncol(data)){
    # Convert back into Python enumeration
    V(g)[k]$label <- as.character(k - 1)
  }
  
  print("  Nodes formatted")

  for(k in 1:(ncol(data) - 1)){
    for(j in (k+1):ncol(data)){
      E(g)[k %--% j]$color <- ifelse(scaffold[k,j] == 1, "lightgrey", "black")
      E(g)[k %--% j]$scaffold <- scaffold[k,j]
      E(g)[k %--% j]$mult <- length(E(g)[k %--% j])
    }
  }
  
  print("  Edges formatted")

  ### Cleaning ====
  if(remove.isolated){
    isolated <- which(degree(g) == 0)
    g.simple <- delete_vertices(g, isolated)
    simple.labels <- V(g)[-isolated]
  }else{
    g.simple <- g
    simple.labels <- V(g)
  }
  
  ### Plotting ====
  layout <- layout_on_sphere(g.simple)
  print("  Layout computed")
  # plot(
  #   g.simple, 
  #   #layout = layout_with_graphopt(g.simple, spring.length = 4, charge = 0.04),
  #   layout = layout,
  #   labels = simple.labels,
  #   vertex.size = 7,
  #   edge.width = 0.5,
  #   edge.arrow.size = 0.2,
  #   #edge.color = "black",
  #   vertex.label.cex = 0.5,
  #   main = paste("alpha =", a)
  # )
  
  ### Export graph ====
  # local.results.dir <- path.join(RESULTS.DIR.ROOT, extract.path.tail(WD))
  
  for(k in 1:length(V(g.simple))){
    V(g.simple)[k]$posx <- layout[k, 1]
    V(g.simple)[k]$posy <- layout[k, 2]
  }
  
  print("  Layout added to graph")

  if(DEBUG){
    write_graph(g.simple, path.join(local.results.dir, paste0("DEBUG_100alpha_", 100*alpha, ".gml")), format = "gml")
  }else{
    write_graph(
      g.simple,
      path.join(local.results.dir, paste0(ifelse(USE.PCALG, "pc", "hc"),"100alpha_", 100*alpha, ".gml")),
      format = "gml"
    )
  }
  
  print("Done! Graph written to file")

  list(fit = ifelse(USE.PCALG, fit.pc, fit.hc), dataset = extract.path.tail(WD), nodes = nodes, alpha = alpha, wd = WD)
}

# Save workspace ====
save.image(path.join(RESULTS.DIR.ROOT, "infer_bn.RData"))

