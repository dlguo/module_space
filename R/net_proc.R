##################################################
## Project: module_space
## Script purpose: network processing
## Date: Mon Mar 26 12:36:06 2018
## Author: Dali Guo
##################################################

require(igraph)
require(Hmisc)
# require(ppcor)
# require(corpcor)
require(gdata)

# file loading function (convert to correlation matrix)
GetGlasserNets <- function(ts, windowSize, cutoff, rsn7, rsn17, cen, skip) {
  
  # remove the first and last frames
  ts <- ts[, -c(1:skip, (dim(ts)[2]-skip+1):dim(ts)[2])]
  # get the dimension
  m <- dim(ts)[1]
  n <- dim(ts)[2]
  
  # calculate the correlation matrix
  nwindows <- n-windowSize+1
  glist <- list()
  for (i in 1:nwindows){
    corr_list <- rcorr(t(ts[, i:(i+windowSize-1)]), type='pearson')
    corrmat <- corr_list$r
    if (cutoff == FALSE) {
      corrmat[corrmat<0] <- 0
    }
    else {
      corrmat[corrmat>cutoff] <- 1
      corrmat[corrmat<cutoff] <- 0
    }
    g <- graph_from_adjacency_matrix(corrmat, mode = "undirected", weighted = TRUE, diag=FALSE)
    g <- set_vertex_attr(g, "rsn7", value=rsn7)
    g <- set_vertex_attr(g, "rsn17", value=rsn17)
    g <- set_vertex_attr(g, "cen", value=as.list(data.frame(cen)))
    V(g)$color <- V(g)$rsn7
    glist[[i]] <- g
  }
  return(glist)
}

# Convert correlation matrix into network (thresholding)
convertMST <- function(corrmat, rsn7, rsn17, cen, cost){
  if(is.list(corrmat)) {
    n <- dim(corrmat[[1]])[1]
    glist <- list()
    for(i in 1:length(corrmat)) {
      g <- graph_from_adjacency_matrix(corrmat[[i]], mode = "undirected", weighted=TRUE, diag=FALSE)
      g <- set_vertex_attr(g, "rsn7", value=rsn7)
      g <- set_edge_attr(g, "id", value=1:ecount(g))
      mstg <- mst(g, 1/E(g)$weight)
      diffg <- difference(g, mstg)
      od <- order(E(diffg)$weight, decreasing = T)[round(cost*n*(n-1)/2-n+2):ecount(diffg)]
      del_id <- E(diffg)$id[od]
      osg <- delete_edges(g, del_id)
      V(osg)$color <- V(osg)$rsn7
      glist[[i]] <- osg
    }
    return(glist)
  }
  else {
    n <- dim(corrmat)[1]
    g <- graph_from_adjacency_matrix(corrmat, mode = "undirected", weighted=TRUE, diag=FALSE)
    g <- set_vertex_attr(g, "rsn7", value=rsn7)
    g <- set_vertex_attr(g, "rsn17", value=rsn17)
    g <- set_vertex_attr(g, "cen", value=as.list(data.frame(cen)))
    g <- set_edge_attr(g, "id", value=1:ecount(g))
    mstg <- mst(g, 1/E(g)$weight)
    diffg <- difference(g, mstg)
    od <- order(E(diffg)$weight, decreasing = T)[round(cost*n*(n-1)/2-n+2):ecount(diffg)]
    del_id <- E(diffg)$id[od]
    osg <- delete_edges(g, del_id)
    V(osg)$color <- V(osg)$rsn7
    return(osg)
  }
}

convertSimple <- function(corrmat, rsn7, rsn17, cen, cutoff) {
  if(is.list(corrmat)) {
    n <- dim(corrmat[[1]])[1]
    glist <- list()
    for(i in 1:length(corrmat)) {
      corrmat[[i]][corrmat[[i]] < cutoff[i]] <- 0
      corrmat[[i]][corrmat[[i]] > cutoff[i]] <- 1
      g <- graph_from_adjacency_matrix(corrmat[[i]], mode = "undirected", diag=FALSE)
      g <- set_vertex_attr(g, "rsn7", value=rsn7)
      g <- set_vertex_attr(g, "rsn17", value=rsn17)
      g <- set_vertex_attr(g, "cen", value=as.list(data.frame(cen)))
      V(g)$color <- V(g)$rsn7
      glist[[i]] <- g
      # if(i %% 20 == 0) cat('20 complets\n')
    }
    return(glist)
  }
  else {
    n <- dim(corrmat)[1]
    corrmat[corrmat < cutoff[i]] <- 0
    g <- graph_from_adjacency_matrix(corrmat, mode = "undirected", diag=FALSE)
    g <- set_vertex_attr(g, "rsn7", value=rsn7)
    g <- set_vertex_attr(g, "rsn17", value=rsn17)
    g <- set_vertex_attr(g, "cen", value=as.list(data.frame(cen)))
    V(g)$color <- V(g)$rsn7
    return(g)
  }
}

conbineWeight <- function(w1, w2){
  for (i in 1:length(w1)) {
    if (is.na(w1[i])) w1[i] <- w2[i]
  }
  w1
}

getNetChanges <- function(g1, g2){
  rew <- (g1 %m% g2) %u% (g2 %m% g1)
  rew <- set_vertex_attr(rew, "rsn7", value=vertex_attr(rew, 'rsn7_1'))
  rew <- set_vertex_attr(rew, "color", value=vertex_attr(rew, 'color_1'))
  rew <- set_edge_attr(rew, "weight", value=conbineWeight(edge_attr(rew, "weight_1"), edge_attr(rew, "weight_2")))
  rew <- set_edge_attr(rew, "id", value=conbineWeight(edge_attr(rew, "id_1"), edge_attr(rew, "id_2")))
  rew
}

eigen2 <- function(m, cutoff) {
  m[m >= cutoff] <- 1
  m[m < cutoff] <- 0
  tail(eigen(diag(rowSums(m))-m)$values, 2)[1]
}

GetCutoff <- function(corrmat, prec=1e-4) {
  L <- 0
  U <- 1
  remind <- which(rowSums(corrmat)==0)
  if (length(remind)!=0) m <- corrmat[-remind, -remind] else m <- corrmat
  while(U-L>prec) if(eigen2(m, (L+U)/2) < 1e-6) U <- (L+U)/2 else L <- (L+U)/2
  L
}

GenNet1Comp <- function(subj, sess, windowSize) {
  n <- strsplit(sess, "_")[[1]]
  tr <- n[1]
  task <- n[2]
  phase <- n[3]
  data_loc <- paste(data_folder, "/results_SIFT2/", subj, 
                    "/fMRI/", sess, "/", sess, "_glasser_tseries.csv", sep='')
  op <- LoadGlasser(data_loc, parc_loc, windowSize, type = "full")
  cat(sprintf("Scan %s loaded. \n", subj))
  
  corrmat <- op[[1]]
  corrmat <- lapply(corrmat, abs)
  rsn7 <- op[[2]]
  rsn17 <- op[[3]]
  cen <- op[[4]]
  cutoff <- sapply(corrmat, GetCutoff)
  net_series <- convertSimple(corrmat, rsn7, rsn17, cen, cutoff)
  save(net_series, file=paste("../output/raw_net/", subj, "_", sess, ".RData", sep=""))
}

GenNetFixedGSz <- function(subj, sess, parc_file, cutoff, windowSize) {
  n <- strsplit(sess, "_")[[1]]
  tr <- n[1]
  task <- n[2]
  phase <- n[3]
  if(substr(sess, 1, 1) == "r") {
    data_loc <- paste(data_folder, "/results_SIFT2/", subj, 
                      "/fMRI/", sess, "/", sess, "_glasser_GS_bp_z_tseries.csv", sep='')
  } else {
    data_loc <- paste(data_folder, "/results_SIFT2/", subj, 
                      "/fMRI/", sess, "/", sess, "_glasser_GS_z_tseries.csv", sep='')
  }
  op <- LoadGlasser(data_loc, parc_file, windowSize, type = "full")
  cat(sprintf("Scan %s loaded. \n", subj))
  
  corrmat <- op[[1]]
  corrmat <- lapply(corrmat, abs)
  rsn7 <- op[[2]]
  rsn17 <- op[[3]]
  cen <- op[[4]]
  net_series <- convertSimple(corrmat, rsn7, rsn17, cen, rep(cutoff, length(corrmat)))
  save(net_series, file=paste("../output/raw_net_GSz/", subj, "_", sess, ".RData", sep=""))
}

GenNetFixedGS <- function(subj, sess, parc_file, cutoff, windowSize) {
  n <- strsplit(sess, "_")[[1]]
  tr <- n[1]
  task <- n[2]
  phase <- n[3]
  if(substr(sess, 1, 1) == "r") {
    data_loc <- paste(data_folder, "/results_SIFT2/", subj, 
                    "/fMRI/", sess, "/", sess, "_glasser_GS_bp_tseries.csv", sep='')
  } else {
    data_loc <- paste(data_folder, "/results_SIFT2/", subj, 
                    "/fMRI/", sess, "/", sess, "_glasser_GS_tseries.csv", sep='')
  }
  op <- LoadGlasser(data_loc, parc_file, windowSize, type = "full")
  cat(sprintf("Scan %s loaded. \n", subj))
  
  corrmat <- op[[1]]
  corrmat <- lapply(corrmat, abs)
  rsn7 <- op[[2]]
  rsn17 <- op[[3]]
  cen <- op[[4]]
  net_series <- convertSimple(corrmat, rsn7, rsn17, cen, rep(cutoff, length(corrmat)))
  save(net_series, file=paste("../output/raw_net_GS/", subj, "_", sess, ".RData", sep=""))
}

GenNetFixed <- function(subj, sess, parc_file, cutoff, windowSize) {
  n <- strsplit(sess, "_")[[1]]
  tr <- n[1]
  task <- n[2]
  phase <- n[3]
  data_loc <- paste(data_folder, "/results_SIFT2/", subj, 
                    "/fMRI/", sess, "/", sess, "_glasser_tseries.csv", sep='')
  op <- LoadGlasser(data_loc, parc_file, windowSize, type = "full")
  cat(sprintf("Scan %s loaded. \n", subj))
  
  corrmat <- op[[1]]
  corrmat <- lapply(corrmat, abs)
  rsn7 <- op[[2]]
  rsn17 <- op[[3]]
  cen <- op[[4]]
  net_series <- convertSimple(corrmat, rsn7, rsn17, cen, rep(cutoff, length(corrmat)))
  save(net_series, file=paste("../output/raw_net/", subj, "_", sess, ".RData", sep=""))
}

GetNetsGSz <- function(sess, windowSize) {
  for (subj in subj_list) {
    GenNetFixedGSz(subj, sess, parc_file, cutoff, windowSize)
  }
}

GetNetsGS <- function(sess, windowSize) {
  for (subj in subj_list) {
    GenNetFixedGS(subj, sess, parc_file, cutoff, windowSize)
  }
}

GetNets <- function(sess, windowSize) {
  for (subj in subj_list) {
    GenNetFixed(subj, sess, parc_file, cutoff, windowSize)
  }
}

GenNetSeriesGSbpz <- function(sess, subj_list, windowSize, cutoff, rsn7, rsn17, cen) {
  n <- strsplit(sess, "_")[[1]]
  tr <- n[1]
  task <- n[2]
  phase <- n[3]
  out_dir <- paste("../output/nets_", as.character(windowSize), "f_t", as.character(cutoff*100), sep='')
  if (!file.exists(out_dir)) {
    dir.create(out_dir)
  }
  for (subj in subj_list) {
    if(substr(sess, 1, 1) == "r") {
      data_loc <- paste(data_folder, "/results_SIFT2/", subj, 
                        "/fMRI/", sess, "/", sess, "_glasser_GS_bp_z_tseries.csv", sep='')
    } else {
      data_loc <- paste(data_folder, "/results_SIFT2/", subj, 
                        "/fMRI/", sess, "/", sess, "_glasser_GS_z_tseries.csv", sep='')
    }
    ts <- as.matrix(read.csv(data_loc, header = FALSE))
    glist <- GetGlasserNets(ts, windowSize, cutoff, rsn7, rsn17, cen, skip)
    saveRDS(glist, file=paste(out_dir, "/", subj, "_", sess, ".rds", sep=""))
    cat(paste(subj, 'network series is generated and saved.\n'))
  }
}

# Convert between igraph and networks
iG2Net <- function(g){
  if(class(g)=="list") {
    rsn7 <- V(g[[1]])$rsn7
    rsn17 <- V(g[[1]])$rsn17
    cen <- V(g[[1]])$cen
    netlist <- list()
    for(i in 1:length(g)) {
      adj_mat <- as_adjacency_matrix(g[[i]], sparse = F)
      netlist[[i]] <- network(adj_mat, vertex.attr = list(rsn7, rsn17, as.list(data.frame(cen))), vertex.attrnames = list("rsn7", "rsn17", "cen"), directed = F)
    }
    return(netlist)
  }
  else if(class(g) == "igraph") {
    rsn7 <- V(g)$rsn7
    rsn17 <- V(g)$rsn17
    cen <- V(g)$cen
    adj_mat <- as_adjacency_matrix(g, sparse = F)
    return(network(adj_mat, list(rsn7, rsn17, as.list(data.frame(cen))), vertex.attrnames = list("rsn7", "rsn17", "cen"), directed = F))
  }
}

Net2iG <- function(net) {
  if(class(net)=="list") {
    rsn7 <- get.vertex.attribute(net[[1]], 'rsn7')
    netlist <- list()
    for (i in length(net)) {
      adj_mat <- as.matrix.network(net[[i]])
      g <- graph_from_adjacency_matrix(adj_mat, mode='undirected')
      g <- set_vertex_attr(g, "rsn7", value=rsn7)
      netlist[[i]] <- g
    }
    netlist
  }
  if(class(net)=="network") {
    rsn7 <- get.vertex.attribute(net, 'rsn7')
    adj_mat <- as.matrix.network(net)
    g <- graph_from_adjacency_matrix(adj_mat, mode='undirected')
    g <- set_vertex_attr(g, "rsn7", value=rsn7)
    return(g)
  }
}
