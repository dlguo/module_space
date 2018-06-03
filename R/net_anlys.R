##################################################
## Project: module_space
## Script purpose: generate statistics/plot for networks
## Date: Mon Mar 26 12:30:48 2018
## Author: Dali Guo
##################################################

require(ggplot2)
require(igraph)

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
ModuleSeries <- function(glist) sapply(glist, function(g) modularity(g, V(g)$rsn7))
TransitivitySeries <- function(glist) sapply(glist, transitivity)
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
