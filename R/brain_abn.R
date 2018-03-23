library(igraph)
library(Hmisc)
library(ppcor)
library(corpcor)
library(gdata)
library(parallel)

source("./net_proc.R")    # load network processing script
source("./abng_opt.R")    # load abn algorithm
set.seed(2)

data_folder <- Sys.getenv("DATA")
windowSize <- 100
cost <- .04
skip <- 10

subj_list <- list.dirs(paste(data_folder, "results_SIFT2", sep = '/'), recursive = F, full.names = F)
parc_loc <- paste(data_folder, "Table_HCP_20170502.xls", sep='/')

get_nets <- function(sess, windowSize) {
  for (subj in subj_list) {
    gen_net(subj, sess, windowSize)
  }
}

abn <- function(subj, sess) {
  # act_mat_list <- list()
  # syn_net_list <- list()
  load(paste("../output/raw_net/", subj, "_", sess, ".RData", sep=""))
  global_f <- c(igraph::degree,
                inv_cl_coe,
                cl,
                betweenness,
                local_eff)
  local_f <- c(in_rsn7)
  act_mat <- est_act(net_series, global_f, local_f, allow_low = T)
  cat(paste(Sys.time(), "estimation done for", subj, sess), sep="\n", file="../output/msg.txt", append=T)
  save(act_mat, file=paste("../output/act_mat/", subj, "_", sess, ".RData", sep=""))
  #syn_net <- syn_act(net_series[[1]], length(net_series), global_f, local_f, act_mat)
  #cat(paste(Sys.time(), "synthesization done for", subj, sess), sep="\n", file="../output/msg.txt", append=T)
  #save(syn_net, file=paste("../output/syn_net/", subj, "_", sess, ".RData", sep=""))
}

#get_nets("rfMRI_REST1_LR", 100)
#get_nets("tfMRI_EMOTION_LR", 20)
#get_nets("tfMRI_GAMBLING_LR", 20)
#get_nets("tfMRI_MOTOR_LR", 20)
#get_nets("tfMRI_RELATIONAL_LR", 20)
#get_nets("tfMRI_SOCIAL_LR", 20)
#get_nets("tfMRI_WM_LR", 20)
#get_nets("tfMRI_LANGUAGE_LR", 20)

args = commandArgs(trailingOnly=TRUE)
abn(args[1], args[2])

# no_cores <- detectCores()
# cl <- makeCluster(no_cores, type="FORK", outfile="debug140.txt")
# subj_list <- as.character(sapply(subj_list, function(x) paste(x, "LR", 1, sep=''), simplify = "array"))
# est_result <- parLapply(cl, subj_list[1:40], main)

# save(est_result, file="../output/est_result.RData")

# stopCluster(cl)

####### test ########
# Single layer
# op <- load_glasser(data_loc, parc_loc, windowSize = 0, type = "full")
# corrmat <- op[[1]]
# corrmat <- abs(corrmat)
# rsn7 <- op[[2]]
# rsn17 <- op[[3]]
# cen <- op[[4]]
# g <- convertMST(corrmat, rsn7, rsn17, cen, cost)
# net <- iGtoNetwork(g)
# plot.igraph(g, vertex.label=NA, vertex.size=2)
# g.ergm.fit <- ergm(rewnet~edges+gwesp(decay=0.75, fixed=TRUE)+nodematch("rsn7"), verbose=T)

# Multi layer
# op <- load_glasser(data_loc, parc_loc, windowSize, type = "full")
# corrmat <- op[[1]]
# corrmat <- lapply(corrmat, abs)
# rsn7 <- op[[2]]
# rsn17 <- op[[3]]
# cen <- op[[4]]
# cutoff <- sapply(corrmat, get_cutoff)
# glist <- convertSimple(corrmat, rsn7, rsn17, cen, cutoff)
# netlist <- iGtoNetwork(glist)
# plot.igraph(glist[[2]], vertex.label=NA, vertex.size=3)
###### loading end ######

# global_f <- c(igraph::degree, inv_deg,
#               cl_coe, inv_cl_coe,
#               cl, inv_cl,
#               betweenness, inv_bt)
# geo_mat <- get_geo_mat(glist[[1]])
# local_f <- c(gd, inv_gd,
#              in_rsn7, inv_rsn7)
# act_mat <- est_act(glist, global_f, local_f, allow_low = T)
# set.seed(2)
# syn_net <- syn_act(glist[[1]], 12, global_f, local_f, act_mat)
# save(act_mat, file="../output/act_mat_1_LR.RData")
# save(syn_net, file="../output/syn_net_1_LR.RData")
# plot.igraph(syn_net[[12]], vertex.label=NA, vertex.size=3)
