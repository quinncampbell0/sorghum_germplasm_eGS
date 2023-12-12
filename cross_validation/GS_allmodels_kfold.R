###############################################################################
##########        Run all k-fold cross-validations        #####################
###############################################################################

library(rrBLUP)
library(hibayes)
library(dplyr)

source("/home/quinnc3/GenomicData/xval_kfold_functions.R") # load all functions

### Load Genotypic Data ####
# note- coded for running on HPC cluster, need to change paths to run locally
gd1 <- read.table("/home/quinnc3/GenomicData/sorghum_GD_numeric_rrblupformat.txt", head = T)
envdat <- read.table('/home/quinnc3/GenomicData/envdat_master.txt', head = T) # full environmental dataset
trainingset <- read.table("/home/quinnc3/GenomicData/envdat_minicore_genomicprediction.txt", head = T)
gd2 <- gd1[gd1$taxa %in% envdat$gen_id,]
row.names(gd2) <- gd2$taxa # set taxa as rownames
g.in <- gd2[,-1]
g.in <- as.matrix(g.in)

### Load Environmental Data ###

row.names(envdat) <- envdat$gen_id # set gen_id as rownames
row.names(trainingset) <- trainingset$gen_id

y.trainset <- trainingset[,c(2,8:(ncol(trainingset)-1))]
y.in <- envdat[,c(2,8:ncol(envdat))]

################################################################################
## Run functions

#BayesCPi
#xval_k2_BayesCpi <- k.xval.BayesCpi(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k.fold = 2, reps = 50)
#saveRDS(xval_k2_BayesCpi, "xval_BayesCpi_kfold_2.RData")

#xval_k5_BayesCpi <- k.xval.BayesCpi(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k.fold = 5, reps = 50)
#saveRDS(xval_k5_BayesCpi, "xval_BayesCpi_kfold_5.RData")

#xval_k10_BayesCpi <- k.xval.BayesCpi(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k.fold = 10, reps = 50)
#saveRDS(xval_k10_BayesCpi, "xval_BayesCpi_kfold_10.RData")

#BayesLASSO
#xval_k2_BayesLASSO <- k.xval.BayesLASSO(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k.fold = 2, reps = 50)
#saveRDS(xval_k2_BayesLASSO, "xval_BayesLASSO_kfold_2.RData")

#xval_k5_BayesLASSO <- k.xval.BayesLASSO(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k.fold = 5, reps = 50)
#saveRDS(xval_k5_BayesLASSO, "xval_BayesLASSO_kfold_5.RData")

#xval_k10_BayesLASSO <- k.xval.BayesLASSO(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k.fold = 10, reps = 50)
#saveRDS(xval_k10_BayesLASSO, "xval_BayesLASSO_kfold_10.RData")

### RR-BLUP
y.trainset.rr <- trainingset[,c(8:(ncol(trainingset)-1))]
y.in.rr <- envdat[,8:ncol(envdat)]

y.trainset.mat <- as.matrix(y.trainset.rr)
y.in.mat <- as.matrix(y.in.rr)

xval_k2_rrblup <- k.xval(g.in = g.in, y.in = y.in.mat, y.trainset = y.trainset.mat, k.fold = 2, reps = 50)
saveRDS(xval_k2_rrblup, "xval_rrblup_kfold_2.RData")

xval_k5_rrblup <- k.xval(g.in = g.in, y.in = y.in.mat, y.trainset = y.trainset.mat, k.fold = 5, reps = 50)
saveRDS(xval_k5_rrblup, "xval_rrblup_kfold_5.RData")

xval_k10_rrblup <- k.xval(g.in = g.in, y.in = y.in.mat, y.trainset = y.trainset.mat, k.fold = 10, reps = 50)
saveRDS(xval_k10_rrblup, "xval_rrblup_kfold_10.RData")


# Gaussian Kernel
K <- A.mat(g.in)
k_dist <- dist(K) # Calculate Relationship Matrix

xval_k2_GAUSS <- k.xval.GAUSS(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = k_dist, k.fold = 2, reps = 50)
saveRDS(xval_k2_GAUSS, "xval_GAUSS_kfold_2.RData")

xval_k5_GAUSS <- k.xval.GAUSS(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = k_dist, k.fold = 5, reps = 50)
saveRDS(xval_k5_GAUSS, "xval_GAUSS_kfold_5.RData")

xval_k10_GAUSS <- k.xval.GAUSS(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = k_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_GAUSS, "xval_GAUSS_kfold_10.RData")

### Exponential Kernel
K.Exp=Kernel_computation(X=g.in, name="exponential", degree=NULL, nL=NULL)
row.names(K.Exp) <- rownames(g.in)
colnames(K.Exp) <- rownames(g.in)
exp_dist <- dist(K.Exp) # Calculate Relationship Matrix\

xval_k2_EXP <- k.xval.EXP(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = exp_dist, k.fold = 2, reps = 50)
saveRDS(xval_k2_EXP, "xval_EXP_kfold_2.RData")

xval_k5_EXP <- k.xval.EXP(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = exp_dist, k.fold = 5, reps = 50)
saveRDS(xval_k5_EXP, "xval_EXP_kfold_5.RData")

xval_k10_EXP <- k.xval.EXP(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = exp_dist, k.fold = 10, reps = 50)
saveRDS(xval_k10_EXP, "xval_EXP_kfold_10.RData")

write.table(exp_dist, "distance_matrix_exponential_kernel.txt")
write.table(k_dist, "distance_matrix_gaussian_kernel.txt")
