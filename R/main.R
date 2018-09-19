library(igraph)
library(Hmisc)
library(gdata)

args <- commandArgs(trailingOnly=TRUE)
source("./net_anlys.R")
source("./net_proc.R")    # load network processing script
source("./net_optim.R")    # load abn algorithm
set.seed(2)

# Assign parsed arguments and load env var
data_folder <- Sys.getenv("DATA")
func <- args[1]
subj <- args[2]
sess <- args[3]
windowSize <- as.numeric(args[4])
if (is.na(as.logical(args[5]))) {
  cutoff <- as.numeric(args[5])
} else {
  cutoff <- as.logical(args[5])
}
cost <- .04
skip <- 5

# Get the community structure
parc_loc <- paste(data_folder, "Table_HCP_20170502.xls", sep='/')
parc <- read.xls(parc_loc, sheet=3)
rsn7 <- (parc$glasser_RSN.7)[1:360]
rsn17 <- (parc$glasser_RSN.17)[1:360]
cenX <- parc$centroid.X[1:360]
cenY <- parc$centroid.Y[1:360]
cenZ <- parc$centroid.Z[1:360]
cen <- rbind(cenX, cenY, cenZ)

# Get LQ
# LQ_loc <- paste(data_folder, "100_Subj_Geo.xlsx", sep='/')
# LQ_file <- read.xls(LQ_loc)
# LQ_subjs <- LQ_file$Subject
# LQ_value <- LQ_file$PMAT24_A_CR


if(func == 'gen') GenCorrMatsGSbpz(sess, subj, windowSize)
if(func == 'mod') {
    mats <- readRDS(paste('../output/corrmats_', windowSize, 'f/', subj, '_', sess, '.rds', sep=''))
    mod_cc <- array(NA, dim=c(length(unique(rsn7)), length(unique(rsn7)), length(unique(rsn7)), dim(mats)[3]))
    for (i in 1:dim(mats)[3]) {
      mod_cc[,,,i] <- GlobalModTrans(mats[,,i], rsn7, 'new')
    }
    saveRDS(mod_cc, file=paste('../output/mod_cc_', windowSize, 'f/', subj, '_', sess, '.rds', sep=''))
}
if(func == 'den') {
    mats <- readRDS(paste('../output/corrmats_', windowSize, 'f/', subj, '_', sess, '.rds', sep=''))
    mod_den <- DynGlobalModDen(mats, rsn7)
    saveRDS(mod_den, file=paste('../output/mod_den_', windowSize, 'f/', subj, '_', sess, '.rds', sep=''))
}
if(func == "dtw") {
    subj_list <- list.files(paste(data_folder, 'results_SIFT2', sep='/'))
    all_dtw_mats <- GetAllDtwMats(sess, subj_list, rsn7)
    saveRDS(all_dtw_mats, file=paste('../output/dtw_mats_', windowSize, 'f/', sess, '.rds', sep=''))
}
if(func == "idm") {
  subj_list <- list.files(paste(data_folder, 'results_SIFT2/', sep='/'))
  ns <- length(subj_list)
  output_folder <- paste(data_folder, '/output/corrmats_', windowSize, 'f/', sep='')
  m1 <- readRDS(paste(output_folder, subj, '_rfMRI_REST1_LR.rds', sep = ''))
  for (i in 1:ns) {
    m2 <- readRDS(paste(output_folder, subj_list[i], '_rfMRI_REST2_LR.rds', sep = ''))
    file.create(paste('../output/idm/', subj, '_', subj_list[i], '.txt', sep=''))
    write(cor(c(m1), c(m2)), file=paste('../output/idm/', subj, '_', subj_list[i], '.txt', sep=''))
  }
}
if (func == "des") {
  mats <- readRDS(paste('../output/corrmats_', windowSize, 'f/', subj, '_', sess, '.rds', sep=''))
  mean_mat <- rowSums(mats, dims=2)
  sd_mat <- apply(mats, c(1,2), sd)
  saveRDS(mean_mat, file=paste('../output/mean_mat_', windowSize, 'f/', subj, '_', sess, '.rds', sep=''))
  saveRDS(sd_mat, file=paste('../output/sd_mat_', windowSize, 'f/', subj, '_', sess, '.rds', sep=''))
}
if (substr(func,1,3) == "pca") {
  ncom <- as.numeric(substr(func, 4, nchar(func)))
  mats <- readRDS(paste('../output/corrmats_', windowSize, 'f/', subj, '_', sess, '.rds', sep=''))
  pca_output <- pca_mats(mats, ncom)
  saveRDS(pca_output, file=paste('../output/pca_', windowSize, '_comp', ncom, '_f/', subj, '_', sess, '.rds', sep=''))
}
