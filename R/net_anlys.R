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
