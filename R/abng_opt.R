require(igraph)

# Action definition

cluster_coe <- function(g) transitivity(g, type="local", isolates="zero")
geo_dist <- function(g, i) sapply(1:vcount(g), FUN=function(g, i, j) dist(rbind(V(g)$cen[[i]], V(g)$cen[[j]]))[1], g=g, i=i)
avg_geo_dist <- function(g) {
  n <- vcount(g)
  nodes <- 1:n
  avg_d <- c()
  for(i in 1:n) avg_d <- c(avg_d, mean(geo_dist(g, i)))
  avg_d
}
get_geo_mat <- function(g) sapply(1:vcount(g), function(x) geo_dist(g, x))
gd <- function(g, i)  geo_mat[, i]
in_rsn7 <- function(g, i) {
  op <- rep(1e-6, vcount(g))
  op[which(V(g)$rsn7 == V(g)$rsn7[i])] <- 1
  op
}
in_rsn17 <- function(g, i) {
  op <- rep(1e-6, vcount(g))
  op[which(V(g)$rsn17 == V(g)$rsn17[i])] <- 1
  op
}
pr <- function(g) igraph::page.rank(g)$vector
cl <- function(g) closeness(g, normalized = T)
log_deg <- function(g) log(igraph::degree(g)+2)
log_bt <- function(g) log(betweenness(g)+2)
local_eff <- function(g) {
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


inv_deg <- function(g) 1/(igraph::degree(g)+0.001)
inv_cl <- function(g) 1/cl(g)
inv_bt <- function(g) 1/(betweenness(g)+0.001)
inv_cl_coe <- function(g) 1/(cluster_coe(g)+0.001)
inv_gd <- function(g, node) 1/gd(g, node)
inv_rsn7 <- function(g, node) 1-in_rsn7(g, node)
inv_rsn17 <- function(g, node) 1-in_rsn17(g, node)
inv_le <- function(g) 1/(local_eff(g)+0.001)

# {
#   n <- vcount(g)
#   geo_mat <- matrix(nrow = n, ncol = n)
#   for(i in 1:n) geo_mat[, i] = sapply(1:n, function(x) geo_dist(g, x))
# }

est_act <- function(g, global_f, local_f, allow_low = F, verbose = F){
  n <- vcount(g[[1]])
  nwindows <- length(g)
  k <- length(global_f)
  l <- length(local_f)
  act_mat <- matrix(nrow = n, ncol = k+l+1)
  LL <- function(par_vec) -sum(log(par_vec %*% t(prob_mat)))  # log-likelihood function. Negative to minimize
  LL_grad <- function(par_vec) {
    grad_vec <- c()
    for(i in 1:length(par_vec)) {
      grad_value <- 0
      for (j in 1:dim(prob_mat)[1]) grad_value = grad_value + prob_mat[j, i]/(par_vec %*% prob_mat[j, ])
      grad_vec <- c(grad_vec, grad_value)
    }
    -grad_vec
  }  # gradient function for BFGS
  for(node in 1:n){
    if(verbose) cat("\nStart to estimate node #", as.character(node), ":\n", sep="")
    if(!allow_low & sum(sapply(g, function(sg) length(ego(sg, 1, node, mindist=1)[[1]]))) < nwindows) {
      if(verbose) message("no enough neighbors.")
      next
    }
    prob_mat <- matrix(nrow=0, ncol = k+l+1)
    for(i in 1:(nwindows-1)) {
      unconn <- setdiff(1:n, ego(g[[i]], 1, node)[[1]])  # potential neighbors
      new_node <- ego(g[[i+1]]-g[[i]], 1, node, mindist=1)[[1]]  # rewired edges
      no_act <- intersect(ego(g[[i]], 1, node, mindist=1)[[1]], ego(g[[i+1]], 1, node, mindist=1)[[1]])  # non-rewired edges
      all_stats <- matrix(nrow=n, ncol=k+l)
      for(j in 1:k) all_stats[, j] <- global_f[[j]](g[[i]])
      for(j in (k+1):(k+l)) all_stats[, j] <- local_f[[j-k]](g[[i]], node)
      prob_mat <- rbind(prob_mat, cbind(matrix(c(prop.table(all_stats[unconn,], 2)[match(new_node, unconn),], rep(0, length(new_node))), ncol=k+l+1)))
      prob_mat <- rbind(prob_mat, matrix(c(rep(0, length(no_act)*(k+l)), rep(1, length(no_act))), ncol=k+l+1))
      prob_mat[is.na(prob_mat)] <- 0
    }
    if(sum(prob_mat[, k+l+1])!=dim(prob_mat)[1]){
      if(verbose) cat(as.character(dim(prob_mat)[1]), "of nodes are changed.")
      res <- constrOptim(theta=rep(1/(k+l+1)-1e-3, k+l+1), f=LL, grad=LL_grad, ui=rbind(-rep(1, k+l+1), diag(k+l+1)), ci=c(-1, rep(0,k+l+1)))
      act_mat[node, ] <- res$par
      if(verbose) cat("The best value is", as.character(res$value))
    }
    else if(verbose) message("This node has no action.")
  }
  act_mat
}

syn_act <- function(g1, nwindows, global_f, local_f, act_mat, verbose = F){
  glist <- list(g1)
  n <- vcount(glist[[1]])
  k <- length(global_f)
  l <- length(local_f)
  for(tw in 2:nwindows){
    if(verbose) cat("Generating network in time", tw, ":\n")
    rem_list <- c()
    add_list <- c()
    all_stats <- matrix(nrow=n, ncol=k+l)
    for(j in 1:k) all_stats[, j] <- global_f[[j]](glist[[tw-1]])
    for(node in 1:n){
      if(verbose) cat("Process node", node, "\n")
      if(all(is.na(act_mat[node,]))) next
      for(j in (k+1):(k+l)) all_stats[, j] <- local_f[[j-k]](glist[[tw-1]], node)
      old_node <- ego(glist[[tw-1]], 1, node, mindist = 1)[[1]]
      unconn <- setdiff(1:n, ego(glist[[tw-1]], 1, node)[[1]])  # potential neighbors
      if(verbose) cat("Nodes to be changed:", length(old_node), "\n")
      act <- sample(k+l+1, length(old_node), prob=act_mat[node,], replace=T)
      for(indold in seq_along(act)) {
        if(act[indold] != k+l+1) {
          rew <- sample(unconn, 1, prob=all_stats[unconn, act[indold]]/sum(all_stats[unconn, act[indold]]), replace=T)
          add_list <- c(add_list, node, rew)
          rem_list <- c(rem_list, paste(as.character(node), as.character(old_node[indold]), sep = "|"))
        }
      }
      # all_added <- unique(matrix(add_list, ncol=2, byrow=T))
      # judge_num <- dim(unique(cbind(all_added[,1]+all_added[,2], all_added[,1]-all_added[,2])))[1]
      # cat(judge_num, "edges added, ", length(rem_list), "edges removed.\n\n")
    }
    glist[[tw]] <- simplify(glist[[tw-1]] + edge(add_list)) - edge(rem_list)
  }
  glist
}

est_act_all <- function(g, global_f, local_f, allow_low = F, verbose = F){
  n <- vcount(g[[1]])
  nwindows <- length(g)
  k <- length(global_f)
  l <- length(local_f)
  act_mat <- matrix(nrow = n, ncol = k+l+1)
  LL <- function(par_vec) -sum(log(par_vec[] %*% t(prob_mat)))  # log-likelihood function. Negative to minimize
  LL_grad <- function(par_vec) {
    grad_vec <- c()
    for(i in 1:length(par_vec)) {
      grad_value <- 0
      for (j in 1:dim(prob_mat)[1]) grad_value = grad_value + prob_mat[j, i]/(par_vec %*% prob_mat[j, ])
      grad_vec <- c(grad_vec, grad_value)
    }
    -grad_vec
  }  # gradient function for BFGS
  
}
