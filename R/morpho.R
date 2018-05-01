library(igraph)
library(Matrix)
library(ggplot2)
setwd("~/Dropbox/Projects/module_space/R")
source("./net_proc.R")
source("./net_anlys.R")


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

comm <- function(g){
  D <- shortest.paths(g)
  mean(1/D[lower.tri(D)])
}

DiffusionSeries <- function(glist) sapply(glist, GraphMFPT)
CommSeries <- function(glist) sapply(glist, comm)

subj_list <- c('110411', '135932', '136833', '751348')
# task_list <- c("rfMRI_REST1", "tfMRI_EMOTION", "tfMRI_GAMBLING", "tfMRI_LANGUAGE", "tfMRI_MOTOR", "tfMRI_RELATIONAL", "tfMRI_SOCIAL", "tfMRI_WM")
task_list <- c("rfMRI_REST1", "rfMRI_REST2")
M <- data.frame()
setwd("/home/dali/Dropbox/Projects/module_space/output/raw_net/")
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

plot(1:dim(data)[1], M$dorsal_den, type='l')
ggplot(M, aes(x=fp_den, y=dorsal_den, shape=subject, color=task))+geom_path() + xlim(0,1) + ylim(0,1) + coord_fixed()
p1 <- ggplot(M[M[,1]=="135932",], aes(x=fp_den, y=dorsal_den, color=task))+geom_path() + xlim(0,1) + ylim(0,1) + coord_fixed() + ggtitle("Subject 135932")
p2 <- ggplot(M[M[,1]=="110411",], aes(x=fp_den, y=dorsal_den, color=task))+geom_path() + xlim(0,1) + ylim(0,1) + coord_fixed() + ggtitle("Subject 110411")
p3 <- ggplot(M[M[,1]=="136833",], aes(x=fp_den, y=dorsal_den, color=task))+geom_path() + xlim(0,1) + ylim(0,1) + coord_fixed() + ggtitle("Subject 136833")
p4 <- ggplot(M[M[,1]=="751348",], aes(x=fp_den, y=dorsal_den, color=task))+geom_path() + xlim(0,1) + ylim(0,1) + coord_fixed() + ggtitle("Subject 751348")
multiplot(p1, p2, p3, p4, cols = 2)
ggplot(M[M[,2]=="rfMRI_REST1",], aes(x=visual_den, y=DMN_den, color=subject))+geom_path() + xlim(0,1) + ylim(0,1) + coord_fixed()
ggplot(M[M[,1]=="101915" & M[,2]=="tfMRI_MOTOR",], aes(x=dorsal_den, y=motor_den, color=task))+geom_point()+geom_path() + xlim(0,1) + ylim(0,1) + coord_fixed()

# Plotly
library(plotly)
data <- M[M[,2]=="rfMRI_REST1",]
p <- plot_ly(data, x=~fp_den, y=~dorsal_den, z=rep(1:(dim(data)[1]/length(subj_list)), length(subj_list)), type = 'scatter3d', mode = 'lines', opacity = 1, color=~subject, line = list(width = 6, reverscale = FALSE))

data <- M[M[,2]=="rfMRI_REST1" & M[,1]=='135932',]
p <- plot_ly(data, x=~fp_den, y=~dorsal_den, z=1:(dim(data)[1]), type = 'scatter3d', mode = 'lines', opacity = 1, line = list(width = 6, reverscale = FALSE))
p

chart_link = api_create(p, filename="fp_dorsal_one")
# chart_link
