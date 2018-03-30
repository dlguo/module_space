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

GetNets <- function(sess, windowSize) {
  for (subj in subj_list) {
    GenNet(subj, sess, windowSize)
  }
}

abn <- function(subj, sess) {
  # act_mat_list <- list()
  # SynthesizeNets_list <- list()
  load(paste("../output/raw_net/", subj, "_", sess, ".RData", sep=""))
  global_f <- c(igraph::degree,
                InvClusterCoeff,
                cl,
                betweenness,
                LocalEff)
  local_f <- c(IsRSN7)
  act_mat <- EstimateActMat(net_series, global_f, local_f, allow_low = T)
  cat(paste(Sys.time(), "estimation done for", subj, sess), sep="\n", file="../output/msg.txt", append=T)
  save(act_mat, file=paste("../output/act_mat/", subj, "_", sess, ".RData", sep=""))
  #SynthesizeNets <- SynthesizeNets(net_series[[1]], length(net_series), global_f, local_f, act_mat)
  #cat(paste(Sys.time(), "synthesization done for", subj, sess), sep="\n", file="../output/msg.txt", append=T)
  #save(SynthesizeNets, file=paste("../output/SynthesizeNets/", subj, "_", sess, ".RData", sep=""))
}

#GetNets("rfMRI_REST1_LR", 100)
#GetNets("tfMRI_EMOTION_LR", 20)
#GetNets("tfMRI_GAMBLING_LR", 20)
#GetNets("tfMRI_MOTOR_LR", 20)
#GetNets("tfMRI_RELATIONAL_LR", 20)
#GetNets("tfMRI_SOCIAL_LR", 20)
#GetNets("tfMRI_WM_LR", 20)
#GetNets("tfMRI_LANGUAGE_LR", 20)

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
# op <- LoadGlasser(data_loc, parc_loc, windowSize = 0, type = "full")
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
# op <- LoadGlasser(data_loc, parc_loc, windowSize, type = "full")
# corrmat <- op[[1]]
# corrmat <- lapply(corrmat, abs)
# rsn7 <- op[[2]]
# rsn17 <- op[[3]]
# cen <- op[[4]]
# cutoff <- sapply(corrmat, GetCutoff)
# glist <- convertSimple(corrmat, rsn7, rsn17, cen, cutoff)
# netlist <- iGtoNetwork(glist)
# plot.igraph(glist[[2]], vertex.label=NA, vertex.size=3)
###### loading end ######

# global_f <- c(igraph::degree, InvDeg,
#               cl_coe, InvClusterCoeff,
#               cl, InvClose,
#               betweenness, InvBtw)
# geo_mat <- GetGeoDistMat(glist[[1]])
# local_f <- c(gd, InvGeoDist,
#              IsRSN7, InvIsRSN7)
# act_mat <- EstimateActMat(glist, global_f, local_f, allow_low = T)
# set.seed(2)
# SynthesizeNets <- SynthesizeNets(glist[[1]], 12, global_f, local_f, act_mat)
# save(act_mat, file="../output/act_mat_1_LR.RData")
# save(SynthesizeNets, file="../output/SynthesizeNets_1_LR.RData")
# plot.igraph(SynthesizeNets[[12]], vertex.label=NA, vertex.size=3)
