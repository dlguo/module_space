library(igraph)
library(Hmisc)
library(gdata)

args <- commandArgs(trailingOnly=TRUE)
source("./net_proc.R")    # load network processing script
source("./net_optim.R")    # load abn algorithm
set.seed(2)

# Capture from args
windowSize <- as.numeric(args[1])
cutoff <- as.numeric(args[2])
subj_list <- args[3]
sess <- args[4]

data_folder <- Sys.getenv("DATA")
cost <- .04
skip <- 5

parc_loc <- paste(data_folder, "Table_HCP_20170502.xls", sep='/')
parc <- read.xls(parc_loc, sheet=3)
rsn7 <- (parc$glasser_RSN.7)[1:360]
rsn17 <- (parc$glasser_RSN.17)[1:360]
cenX <- parc$centroid.X[1:360]
cenY <- parc$centroid.Y[1:360]
cenZ <- parc$centroid.Z[1:360]
cen <- rbind(cenX, cenY, cenZ)


GenNetSeriesGSbpz(sess, subj_list, windowSize, cutoff, rsn7, rsn17, cen)
