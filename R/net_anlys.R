##################################################
## Project: module_space
## Script purpose: generate statistics/plot for networks
## Date: Mon Mar 26 12:30:48 2018
## Author: Dali Guo
##################################################

require(ggplot2)
require(igraph)
require(dtw)
require(parallel)

# Single network statistics
ClusterCoeff <- function(g) transitivity(g, type = "local", isolates = "zero")

GeoDist <- function(g, i)
  sapply(1:vcount(g), function(g, i, j)
    dist(rbind(V(g)$cen[[i]], V(g)$cen[[j]]))[1], g = g, i = i)

AvgGeoDist <- function(g) {
  n <- vcount(g)
  avg_d <- c()
  for (i in 1:n) avg_d <- c(avg_d, mean(GeoDist(g, i)))
  avg_d
}

GetGeoDistMat <- function(g)
  sapply(1:vcount(g), function(x) GeoDist(g, x))

gd <- function(g, i)  geo_mat[, i]

IsRSN7 <- function(g, i) {
  op <- rep(1e-6, igraph::vcount(g))
  op[which(igraph::V(g)$rsn7 == igraph::V(g)$rsn7[i])] <- 1
  op
}
IsRSN17 <- function(g, i) {
  op <- rep(1e-6, vcount(g))
  op[which(V(g)$rsn17 == V(g)$rsn17[i])] <- 1
  op
}
pr <- function(g) igraph::page.rank(g)$vector
cl <- function(g) closeness(g, normalized = T)
LogDeg <- function(g) log(igraph::degree(g)+2)
LogBtw <- function(g) log(betweenness(g)+2)
LocalEff <- function(g) {
  degs <- igraph::degree(g)
  
  eff <- numeric(length(degs))
  nodes <- which(degs > 1)
  
  eff[nodes] <- simplify2array(lapply(nodes, function(x) {
    neighbs <- neighbors(g, v=x)
    g.sub <- induced.subgraph(g, neighbs)
    Nv <- vcount(g.sub)
    
    paths <- shortest.paths(g.sub, weights=NA)
    paths <- paths[upper.tri(paths)]
    2 / Nv / (Nv - 1) * sum(1 / paths[paths != 0])
  })
  )
  eff
}


InvDeg <- function(g) 1/(igraph::degree(g)+0.001)
InvClose <- function(g) 1/cl(g)
InvBtw <- function(g) 1/(betweenness(g)+0.001)
InvClusterCoeff <- function(g) 1/(ClusterCoeff(g)+0.001)
InvGeoDist <- function(g, node) 1/gd(g, node)
InvIsRSN7 <- function(g, node) 1-IsRSN7(g, node)
InvIsRSN17 <- function(g, node) 1-IsRSN17(g, node)
InvLocalEff <- function(g) 1/(LocalEff(g)+0.001)
MFPTfct <- function(adj) {
  # calculate the mean first-passage time (MFPT) for a fully connected graph from the adjacency matrix
  # note: this function is unable to deal with graphs that are not fully connected
  # 
  # if the adjacency matrix contains only a single gene, return a 1x1 matrix containing 0
  if (is.null(dim(adj))) return(matrix(0, 1, 1))
  ngenes <- nrow(adj)
  
  A <- adj / apply(adj, 1, sum)  # A: the transition probability matrix (there is always movement)
  I <- Diagonal(x=rep(as.integer(1), ngenes))
  pi <- as.numeric(rep(1/ngenes, ngenes) %*% solve(I - A + 1 / ngenes)) # pi: the stationary distribution of the transition matrix
  Z <- solve(t(t(I - A) - pi))
  as.matrix(t(t(I - Z) + Z[cbind(1:nrow(Z), 1:nrow(Z))]) %*% (I * (1 / pi))) # M: the mean first passage matrix
}
GraphMFPT <- function(g) MFPTfct(as_adjacency_matrix(g))
GlobalEff <- function(g){
  D <- distances(g)
  mean(1/D[lower.tri(D)])
}


# Network series statistics function
DensityDiff <- function(glist) {
  t_den <- c()
  for(i in 2:length(glist)) t_den <- c(t_den, graph.density(glist[[i]]-glist[[i-1]]))
  t_den
}
DensitySeries <- function(glist) sapply(glist, graph.density)
ModuleDiff <- function(glist) {
  t_mod <- c()
  com <-V(glist[[1]])$rsn7
  for(i in 2:length(glist)) t_mod <- c(t_mod, modularity(glist[[i]]-glist[[i-1]], com))
  t_mod
}
ModuleSeries <- function(glist) sapply(glist, function(g) modularity(g, V(g)$rsn7, weight=E(g)$weight))
TransitivitySeries <- function(glist) sapply(glist, transitivity)
AssortativitySeries <- function(glist) sapply(glist, assortativity_degree)
KStestDiff <- function(g1, g2, f) {
  ksd <- c()
  for(i in 2:length(g1)) ksd <- c(ksd, ks.test(f(g1[[i]]-g1[[i-1]]), f(g2[[i]]-g2[[i-1]]))$statistic[[1]])
  ksd
}
KStestSeries <- function(g1, g2, f) {
  ksd <- c()
  for(i in 1:length(g1)) ksd <- c(ksd, ks.test(f(g1[[i]]), f(g2[[i]]))$statistic[[1]])
  ksd
}
DiffusionSeries <- function(glist) sapply(glist, GraphMFPT)
GlobalEffSeries <- function(glist) sapply(glist, GlobalEff)

# Modular density
GlobalModEdges <- function(mat, comm){
  mat[mat<0] <- 0
  n <- dim(mat)[1]
  num_of_comm <- length(unique(comm))
  mod_den <- array(0, dim=c(num_of_comm, num_of_comm))
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      a <- max(c(comm[i], comm[j]))
      b <- min(c(comm[i], comm[j]))
      mod_den[a, b] <- mod_den[a, b] + mat[i, j]
    }
  }
  mod_den
}
DynGlobalModDen <- function(mats, comm) {
  num_of_comm <- length(unique(comm))
  allT <- dim(mats)[3]
  # Get regularization matrix
  reg_mat <- matrix(nrow=num_of_comm, ncol=num_of_comm)
  for (i in 1:num_of_comm) {
    ne <- length(which(comm==i))
    reg_mat[i,i] <- ne*(ne-1)
  }
  for (j in 1:(num_of_comm-1)) {
    for (i in(j+1):num_of_comm) {
      ni <- length(which(comm==i))
      nj <- length(which(comm==j))
      reg_mat[i,j] <- ni*nj
    }
  }
  mod_dens <- array(dim=c(num_of_comm,num_of_comm,allT))
  for (i in 1:allT) {
    mod_dens[,,i] <- 2*GlobalModEdges(mats[,,i], comm)/reg_mat
  }
  mod_dens
}

# Modular transitivity
VertexModTrans <- function(m, comm){
  m[m<0] <- 0
  mats_cc <- array(NA, dim=c(length(unique(comm)), length(unique(comm)), length(unique(comm))))
  ig <- list()
  num_of_comm <- length(unique(comm))
  for (g in sort(unique(comm))) {
    ig[[g]] <- which(comm == g)
  }
  for (g in sort(unique(comm))) {
    mcc <- array(NA, dim=c(length(unique(comm)), length(unique(comm)), length(ig[[g]])))
    iv <- 1
    for (v in ig[[g]]) {
      # igv <- ig[[g]][ig[[g]] != v]
      # igv <- igv[m[v, igv] != 0]
      # s <- sum(m[v, igv])
      # k <- sum(m[v, igv] != 0)
      # n_pairs <- combn(igv, 2)
      # sw <- 0
      # for (n in 1:dim(n_pairs)[2]) {
      #   if (m[n_pairs[1,n], n_pairs[2,n]] != 0) {
      #     sw <- sw + (m[v, n_pairs[1,n]] + m[v, n_pairs[2,n]])/2
      #   }
      # }
      # mats_cc[g, g, v] <- sw/s/k
      nbr <- which(m[v,]>0)
      nbr <- nbr[nbr != v]
      s <- array(0, dim=c(num_of_comm, num_of_comm))
      k <- array(0, dim=c(num_of_comm, num_of_comm))
      sw <- array(0, dim=c(num_of_comm, num_of_comm))
      for (ind_n1 in 1:(length(nbr)-1)) {
        n1 <- nbr[ind_n1]
        k[comm[n1], ] <- k[comm[n1], ] + 1
        s[comm[n1], ] <- s[comm[n1], ] + m[v, n1]
        for(ind_n2 in (ind_n1+1):length(nbr)){
          n2 <- nbr[ind_n2]
          if (m[n1, n2] > 0) {
            sw[comm[n1], comm[n2]] <- sw[comm[n1], comm[n2]] + (m[v, n1] + m[v, n2])/2
          }
        }
      }
      k[comm[n2], ] <- k[comm[n2], ] + 1
      s[comm[n2], ] <- s[comm[n2], ] + m[v, n2]
      sw[lower.tri(sw)] <- sw[upper.tri(sw)] + sw[lower.tri(sw)]
      sw[upper.tri(sw)] <- NA
      k[lower.tri(k)] <- k[upper.tri(k)] + k[lower.tri(k)]
      k[upper.tri(k)] <- NA
      s[lower.tri(s)] <- s[upper.tri(s)] + s[lower.tri(s)]
      s[upper.tri(s)] <- NA
      mcc[,,iv] <- sw/s/(k-1)
      iv <- iv + 1
    }
    mcc[is.nan(mcc)] <- 0
    mcc[is.infinite(mcc)] <- 0
    mats_cc[g,,] <- rowMeans(mcc, dims = 2)
  }
  mats_cc
}
GlobalModTrans <- function(adj_mat, comm, type=c('arithmetic', 'geometric', 'new')){
  adj_mat[adj_mat<0] <- 0
  num_of_comm <- length(unique(comm))
  n <- dim(adj_mat)[1]
  closed_tri <- array(0, dim=c(num_of_comm, num_of_comm, num_of_comm))
  all_tri <- array(0, dim=c(num_of_comm, num_of_comm, num_of_comm))
  if (type == 'arithmetic') {
    for (x in 1:(n-2)) {
      for (y in (x+1):(n-1)) {
        for (z in (y+1):n) {
          if (adj_mat[x,y]*adj_mat[x,z]*adj_mat[y,z] != 0) {
            closed_tri[comm[x], comm[y], comm[z]] <- closed_tri[comm[x], comm[y], comm[z]] + (adj_mat[x,y] + adj_mat[x,z])/2
            closed_tri[comm[y], comm[x], comm[z]] <- closed_tri[comm[y], comm[x], comm[z]] + (adj_mat[x,y] + adj_mat[y,z])/2
            closed_tri[comm[z], comm[x], comm[y]] <- closed_tri[comm[z], comm[x], comm[y]] + (adj_mat[x,z] + adj_mat[y,z])/2
            all_tri[comm[x], comm[y], comm[z]] <- all_tri[comm[x], comm[y], comm[z]] + (adj_mat[x,y] + adj_mat[x,z])/2
            all_tri[comm[y], comm[x], comm[z]] <- all_tri[comm[y], comm[x], comm[z]] + (adj_mat[x,y] + adj_mat[y,z])/2
            all_tri[comm[z], comm[x], comm[y]] <- all_tri[comm[z], comm[x], comm[y]] + (adj_mat[x,z] + adj_mat[y,z])/2
          }
          else if (adj_mat[x,y]*adj_mat[x,z] + adj_mat[x,y]*adj_mat[y,z] + adj_mat[y,z]*adj_mat[x,z] != 0) {
            all_tri[comm[x], comm[y], comm[z]] <- all_tri[comm[x], comm[y], comm[z]] + (adj_mat[x,y] + adj_mat[x,z])/2
            all_tri[comm[y], comm[x], comm[z]] <- all_tri[comm[y], comm[x], comm[z]] + (adj_mat[x,y] + adj_mat[y,z])/2
            all_tri[comm[z], comm[x], comm[y]] <- all_tri[comm[z], comm[x], comm[y]] + (adj_mat[x,z] + adj_mat[y,z])/2
          }
        }
      }
    }
  }
  else if (type == 'geometric') {
    for (x in 1:(n-2)) {
      for (y in (x+1):(n-1)) {
        for (z in (y+1):n) {
          if (adj_mat[x,y]*adj_mat[x,z]*adj_mat[y,z] != 0) {
            closed_tri[comm[x], comm[y], comm[z]] <- closed_tri[comm[x], comm[y], comm[z]] + sqrt(adj_mat[x,y] * adj_mat[x,z])
            closed_tri[comm[y], comm[x], comm[z]] <- closed_tri[comm[y], comm[x], comm[z]] + sqrt(adj_mat[x,y] * adj_mat[y,z])
            closed_tri[comm[z], comm[x], comm[y]] <- closed_tri[comm[z], comm[x], comm[y]] + sqrt(adj_mat[x,z] * adj_mat[y,z])
            all_tri[comm[x], comm[y], comm[z]] <- all_tri[comm[x], comm[y], comm[z]] + sqrt(adj_mat[x,y] * adj_mat[x,z])
            all_tri[comm[y], comm[x], comm[z]] <- all_tri[comm[y], comm[x], comm[z]] + sqrt(adj_mat[x,y] * adj_mat[y,z])
            all_tri[comm[z], comm[x], comm[y]] <- all_tri[comm[z], comm[x], comm[y]] + sqrt(adj_mat[x,z] * adj_mat[y,z])
          }
          else if (adj_mat[x,y]*adj_mat[x,z] + adj_mat[x,y]*adj_mat[y,z] + adj_mat[y,z]*adj_mat[x,z] != 0) {
            all_tri[comm[x], comm[y], comm[z]] <- all_tri[comm[x], comm[y], comm[z]] + sqrt(adj_mat[x,y] * adj_mat[x,z])
            all_tri[comm[y], comm[x], comm[z]] <- all_tri[comm[y], comm[x], comm[z]] + sqrt(adj_mat[x,y] * adj_mat[y,z])
            all_tri[comm[z], comm[x], comm[y]] <- all_tri[comm[z], comm[x], comm[y]] + sqrt(adj_mat[x,z] * adj_mat[y,z])
          }
        }
      }
    }
  }
  else if (type == 'new') {
    for (x in 1:(n-2)) {
      for (y in (x+1):(n-1)) {
        for (z in (y+1):n) {
          if (adj_mat[x,y]*adj_mat[x,z]*adj_mat[y,z] != 0) {
            closed_tri[comm[x], comm[y], comm[z]] <- closed_tri[comm[x], comm[y], comm[z]] + adj_mat[y,z]/sqrt(adj_mat[x,y] * adj_mat[x,z])
            closed_tri[comm[y], comm[x], comm[z]] <- closed_tri[comm[y], comm[x], comm[z]] + adj_mat[x,z]/sqrt(adj_mat[x,y] * adj_mat[y,z])
            closed_tri[comm[z], comm[x], comm[y]] <- closed_tri[comm[z], comm[x], comm[y]] + adj_mat[x,y]/sqrt(adj_mat[x,z] * adj_mat[y,z])
            all_tri[comm[x], comm[y], comm[z]] <- all_tri[comm[x], comm[y], comm[z]] + 1/sqrt(adj_mat[x,y] * adj_mat[x,z])
            all_tri[comm[y], comm[x], comm[z]] <- all_tri[comm[y], comm[x], comm[z]] + 1/sqrt(adj_mat[x,y] * adj_mat[y,z])
            all_tri[comm[z], comm[x], comm[y]] <- all_tri[comm[z], comm[x], comm[y]] + 1/sqrt(adj_mat[x,z] * adj_mat[y,z])
          }
          else if (adj_mat[y,z] == 0 & adj_mat[x,y] !=0 & adj_mat[x,z] !=0) {
            all_tri[comm[x], comm[y], comm[z]] <- all_tri[comm[x], comm[y], comm[z]] + 1/sqrt(adj_mat[x,y] * adj_mat[x,z])
          }
          else if (adj_mat[x,z] == 0 & adj_mat[x,y] !=0 & adj_mat[y,z] !=0) {
            all_tri[comm[y], comm[x], comm[z]] <- all_tri[comm[y], comm[x], comm[z]] + 1/sqrt(adj_mat[x,y] * adj_mat[y,z])
          }
          else if (adj_mat[x,y] == 0 & adj_mat[x,z] !=0 & adj_mat[y,z] !=0) {
            all_tri[comm[z], comm[x], comm[y]] <- all_tri[comm[z], comm[x], comm[y]] + 1/sqrt(adj_mat[x,z] * adj_mat[y,z])
          }
        }
      }
    }
  }
  else print("wrong type.")
  closed_tri/all_tri
}

# Dtw mat
TriDtwMat <- function(a,b){
  n <- dim(a)[1]
  dtw_mat <- matrix(nrow=n, ncol=n)
  for (j in 1:n) {
    for (i in j:n) {
      dtw_mat[i,j] <- dtw(a[i,j,],b[i,j,],distance.only = T)$normalizedDistance
    }
  }
  dtw_mat
}
GetAllDtwMats <- function(sess, subj_list, comm){
  n_subj <- length(subj_list)
  num_of_comm <- length(unique(rsn7))
  all_dtw_mats <- array(NA, dim=c(num_of_comm, num_of_comm, n_subj, n_subj))
  for(i in 1:n_subj){
    for(j in 1:n_subj){
      si_mats <- readRDS(paste("../output/mod_den_90f/", subj_list[i], "_", sess, "_LR.rds", sep = ""))
      sj_mats <- readRDS(paste("../output/mod_den_90f/", subj_list[j], "_", sess, "_RL.rds", sep = ""))
      dtw_mat <- TriDtwMat(si_mats, sj_mats)
      all_dtw_mats[,,i,j] <- dtw_mat
    }
  }
  all_dtw_mats
}

# Plot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Identifiablity
IdMat <- function(a, b) {
  n <- dim(a)[2]
  cormat <- matrix(nrow=n, ncol=n)
  for (x in 1:n) {
    for (y in 1:n) {
      cormat[x,y] <- cor(a[,x], b[,y])
    }
  }
  Iself <- mean(diag(cormat))
  diag(cormat) <- NA
  Idiff <- mean(cormat, na.rm = T)
  100*(Iself-Idiff)
}

# PCA
pca_mats <- function(mats, ncomp) {
  n <- dim(mats)[2]
  lt <- dim(mats)[3]
  ts_vec <- matrix(nrow = lt, ncol = (n^2-n)/2)
  k <- 1
  for (j in 1:(n-1)) {
    for (i in (j+1):n) {
      ts_vec[,k] <- mats[i,j,]
      k <- k+1
    }
  }
  ts_pca <- prcomp(ts_vec)
  ts_vec_re <- t(t(ts_pca$x[,1:ncomp] %*% t(ts_pca$rotation[,1:ncomp])) + ts_pca$center)
  mat_re <- array(NA, dim=dim(mats))
  k <- 1
  for (j in 1:(n-1)) {
    mat_re[j,j,] <- rep(1,lt)
    for (i in (j+1):n) {
      # ts_vecs[,k] <- ts_mats[i,j,]
      mat_re[i,j,] <- ts_vec_re[,k]
      mat_re[j,i,] <- ts_vec_re[,k]
      k <- k+1
    }
  }
  mat_re[n,n,] <- rep(1,lt)
  total_var <- (ts_pca$sdev)^2
  pca_result <- list()
  pca_result$mats_recon <- mat_re
  pca_result$pca_comp <- ts_pca
  pca_result$prop_var <- total_var/sum(total_var)
  pca_result
}


# network entropy
entropyProb <- function(prob, base=length(prob)) {
  if (length(prob) == 1) return(0)
  else return(-sum(sapply(prob, function(p) if (p==0) 0 else p * log(p, base))))
}

edgeEntropy <- function(adjmat) {
  rList <- list()
  N <- dim(adjmat)[3]
  n <- dim(adjmat)[1]
  rList[['N']] <- N
  probmat <- rowSums(adjmat, dims=2)/N
  inputmat <- cbind(c(probmat), 1-c(probmat))
  entropy_mat <- matrix(apply(inputmat, MARGIN = 1, function(x) entropyProb(x,base=2)), nrow=n)
  rList[['entropy_mat']]<- entropy_mat
  rList[['entropy_mean']] <- mean(entropy_mat)
  rList[['entropy_mean_offdiag']] <- mean(entropy_mat[lower.tri(entropy_mat, diag = F)])
  rList
}

triadEntropy <- function(adjmat) {
  rList <- list()
  N <- dim(adjmat)[3]
  n <- dim(adjmat)[1]
  rList[["N"]] <- N
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores, type = "FORK")
  tcomb <- combn(n, 3)
  z <- array(1:8, dim=c(2,2,2))
  y <- adjmat+1
  entropy_mat <- parSapply(cl, 1:choose(n,3), function(tri)
    entropyProb(table(sapply(1:N, function(x) z[y[tcomb[1,tri], tcomb[2,tri], x],
                                                y[tcomb[1,tri], tcomb[3,tri], x],
                                                y[tcomb[2,tri], tcomb[3,tri], x]]))/N, base=8))
  stopCluster(cl)
  rList[['entropy_mat']] <- entropy_mat
  rList[['entropy_mean']] <- mean(entropy_mat)
  rList
}

distanceEntropy <- function(adjmat, rm.neighbor=F) {
  rList <- list()
  N <- dim(adjmat)[3]
  n <- dim(adjmat)[1]
  rList[['N']] <- N
  distmat <- array(dim=c(n,n,N))
  for (i in 1:N) {
    distmat[,,i] <- distances(graph_from_adjacency_matrix(adjmat[,,i], mode = "undirected", diag = F))
  }
  entropy_mat <- array(0, dim=c(n,n))
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (rm.neighbor) {
        td <- table(distmat[i,j,])
        td <- td[-which(names(td)=='1')]
        dist_freq <- td/sum(td)
      } else {
        dist_freq <- table(distmat[i,j,])/N
      }
      if (length(dist_freq)>1) entropy_mat[i,j] <- entropyProb(dist_freq)
    }
  }
  entropy_mat <- entropy_mat + t(entropy_mat)
  rList[['entropy_mat']]<- entropy_mat
  rList[['entropy_mean']] <- mean(entropy_mat)
  rList[['entropy_mean_offdiag']] <- mean(entropy_mat[lower.tri(entropy_mat, diag = F)])
  rList
}

TriadMutualInfo <- function(adjmat) {
  N <- dim(adjmat)[3]
  n <- dim(adjmat)[1]
  entropyMat <- array(dim=c(n,n,n))
  for (i in 1:n) {
    for (j in (1:n)[-i]) {
      for (k in (1:n)[-c(i,j)]) {
        neis <- sapply(1:N, function(x) adjmat[i,j,x]*adjmat[i,k,x])
        ys <- adjmat[j,k,]
        entropyMat[i,j,k] <- entropyProb(table(neis)/N) + entropyProb(table(ys)/N) - entropyProb(table((neis+2)*(ys-2))/N)
      }
    }
  }
  entropyMat
}
