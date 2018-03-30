library(igraph)
library(plot3Drgl)
library(Matrix)
library(ggplot2)
source("./net_proc.R")


MFPTfct <- function(adj) {
  # calculate the mean first-passage time (MFPT) for a fully connected graph from the adjacency matrix
  # note: this function is unable to deal with graphs that are not fully connected
  
  # if the adjacency matrix contains only a single gene, return a 1x1 matrix containing 0
  if (is.null(dim(adj))) return(matrix(0, 1, 1))
  ngenes <- nrow(adj)
  
  A <- adj / apply(adj, 1, sum)  # A: the transition probability matrix (there is always movement)
  I <- Diagonal(x=rep(as.integer(1), ngenes))
  pi <- as.numeric(rep(1/ngenes, ngenes) %*% solve(I - A + 1 / ngenes)) # pi: the stationary distribution of the transition matrix
  Z <- solve(t(t(I - A) - pi))
  as.matrix(t(t(I - Z) + Z[cbind(1:nrow(Z), 1:nrow(Z))]) %*% (I * (1 / pi))) # M: the mean first passage matrix
}

comm <- function(g){
  D <- shortest.paths(g)
  mean(1/D[lower.tri(D)])
}

DiffusionSeries <- function(glist) sapply(glist, GraphMFPT)
CommSeries <- function(glist) sapply(glist, comm)

subj_list <- c("100307", "100408", "101107")
task_list <- c("rfMRI_REST1", "tfMRI_EMOTION", "tfMRI_GAMBLING", "tfMRI_LANGUAGE", "tfMRI_MOTOR", "tfMRI_RELATIONAL", "tfMRI_SOCIAL", "tfMRI_WM")
M <- data.frame()
for(s in subj_list){
  for(task in task_list){
    load(paste(s, task, "LR.RData", sep='_'))
    m <- cbind(DiffusionSeries(net_series), CommSeries(net_series), rep(s, length(net_series)), rep(task, length(net_series)))
    M <- rbind(M, m)
  }
}
colnames(M) <- c("E_diff", "E_comm", "subject", "task")
M$E_diff <- as.numeric(as.character(M$E_diff))
M$E_comm <- as.numeric(as.character(M$E_comm))
ggplot(M, aes(x=E_diff, y=E_comm, shape=subject, color=task))+geom_point()
ggplot(M[M[,3]=="100408",], aes(x=E_diff, y=E_comm, color=task))+geom_point()+geom_segment(aes(xend=c(tail(E_diff, n=-1), NA), yend=c(tail(E_comm, n=-1), NA)),
                                                                                           arrow=arrow(length=unit(0.3,"cm")))
ggplot(M[M[,4]=="rfMRI_REST1",], aes(x=E_diff, y=E_comm, color=subject))+geom_point()+geom_segment(aes(xend=c(tail(E_diff, n=-1), NA), yend=c(tail(E_comm, n=-1), NA)),
                                                                                                     arrow=arrow(length=unit(0.3,"cm")))
