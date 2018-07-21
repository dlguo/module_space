library(igraph)
library(gdata)

args <- commandArgs(trailingOnly=TRUE)
source("./net_anlys.R")    # load network processing script

# Capture from args
subj <- args[1]
sess <- args[2]
windowSize <- args[3]

data_folder <- Sys.getenv("DATA")

parc_loc <- paste(data_folder, "Table_HCP_20170502.xls", sep='/')
parc <- read.xls(parc_loc, sheet=3)
rsn7 <- (parc$glasser_RSN.7)[1:360]
rsn17 <- (parc$glasser_RSN.17)[1:360]
cenX <- parc$centroid.X[1:360]
cenY <- parc$centroid.Y[1:360]
cenZ <- parc$centroid.Z[1:360]
cen <- rbind(cenX, cenY, cenZ)

mats <- readRDS(paste('../output/corrmats_', windowSize, 'f/', subj, '_', sess, '.rds', sep=''))
mod_cc <- array(NA, dim=c(length(unique(rsn7)), length(unique(rsn7)), length(unique(rsn7)), dim(mats)[3]))
for (i in 1:dim(mats)[3]) {
  mod_cc[,,,i] <- GlobalModTrans(mats[,,i], rsn7, 'new')
}
saveRDS(mod_cc, file=paste('../output/mod_cc_', windowSize, 'f/', subj, '_', sess, '.rds', sep=''))