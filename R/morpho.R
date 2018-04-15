library(igraph)
library(plot3Drgl)
library(Matrix)
library(ggplot2)
setwd("~/Dropbox/Projects/module_space/R")
source("./net_proc.R")
source("./net_anlys.R")


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

GraphMFPT <- function(g) MFPTfct(as_adjacency_matrix(g))

comm <- function(g){
  D <- shortest.paths(g)
  mean(1/D[lower.tri(D)])
}

DiffusionSeries <- function(glist) sapply(glist, GraphMFPT)
CommSeries <- function(glist) sapply(glist, comm)

subj_list <- c("100408", "101107", "101309", "101915", "103111", "103414")
task_list <- c("rfMRI_REST1", "tfMRI_EMOTION", "tfMRI_GAMBLING", "tfMRI_LANGUAGE", "tfMRI_MOTOR", "tfMRI_RELATIONAL", "tfMRI_SOCIAL", "tfMRI_WM")
M <- data.frame()
setwd("/home/dali/Dropbox/Projects/17CN/output/raw_net/")
for(s in subj_list){
  for(task in task_list){
    load(paste(s, task, "LR.RData", sep='_'))
    m <- cbind(rep(s, length(net_series)), rep(task, length(net_series)))
    for (rsn in 1:8) {
      sep_net <- lapply(net_series, function(g) induced.subgraph(g, which(V(g)$rsn7==rsn)))
      m <-  cbind(m, DensitySeries(sep_net), TransitivitySeries(sep_net))
    }
    M <- rbind(M, m)
  }
}
colnames(M) <- c("subject", "task", "visual_den", "visual_trans", "motor_den", "motor_trans", "dorsal_den", "dorsal_trans", "ventral_den", "ventral_trans", "limbic_den", "limbic_trans", "fp_den", "fp_trans", "DMN_den", "DMN_trans", "subcortical_den", "subcortical_trans")
for (x in 3:18) {
  M[,x] <- as.numeric(as.character(M[,x]))
}

ggplot(M, aes(x=visual_den, y=DMN_den, shape=subject, color=task))+geom_point()+geom_path() + xlim(0,1) + ylim(0,1) + coord_fixed()
ggplot(M[M[,1]=="103414",], aes(x=fp_den, y=DMN_den, color=task))+geom_point()+geom_path() + xlim(0,1) + ylim(0,1) + coord_fixed()
ggplot(M[M[,2]=="tfMRI_RELATIONAL",], aes(x=visual_den, y=DMN_den, color=subject))+geom_point()+geom_path() + xlim(0,1) + ylim(0,1) + coord_fixed()
ggplot(M[M[,1]=="101915" & M[,2]=="tfMRI_MOTOR",], aes(x=dorsal_den, y=motor_den, color=task))+geom_point()+geom_path() + xlim(0,1) + ylim(0,1) + coord_fixed()

# Plotly
library(plotly)
data <- M[M[,2]=="tfMRI_MOTOR",]
p <- plot_ly(data, x=~visual_den, y=~fp_den, z=rep(1:(dim(data)[1]/6), 6), type = 'scatter3d', mode = 'lines', opacity = 1, color=~subject, line = list(width = 6, reverscale = FALSE))
p

chart_link = api_create(p, filename="motor")
# chart_link
