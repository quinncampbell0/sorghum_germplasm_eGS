###############################################################################
##########        Run all fractional cross-validations        #################
###############################################################################

source("/home/quinnc3/GenomicData/xval_frac_functions.R")
library(rrBLUP)
library(hibayes)
library(dplyr)

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

row.names(envdat)<- envdat$gen_id
row.names(trainingset) <- trainingset$gen_id

y.trainset <- trainingset[,c(2,8:(ncol(trainingset)-1))]
y.in <- envdat[,c(2,8:ncol(envdat))]

############################################################################################
# Run functions

#BayesCPi
#xval_BayesCpi_frac50 <- frac.xval.BayesCpi(g.in = g.in, y.in = y.in, y.trainset = y.trainset, frac.train = 0.5, reps = 50)
#saveRDS(xval_BayesCpi_frac50, "xval_BayesCpi_frac50.RData")

#xval_BayesCpi_frac80 <- frac.xval.BayesCpi(g.in = g.in, y.in = y.in, y.trainset = y.trainset, frac.train = 0.8, reps = 50)
#saveRDS(xval_BayesCpi_frac80, "xval_BayesCpi_frac80.RData")

#xval_BayesCpi_frac90 <- frac.xval.BayesCpi(g.in = g.in, y.in = y.in, y.trainset = y.trainset, frac.train = 0.9, reps = 50)
#saveRDS(xval_BayesCpi_frac90, "xval_BayesCpi_frac90.RData")

#BayesLASSO

#xval_BayesLASSO_frac50 <- frac.xval.BayesLASSO(g.in = g.in, y.in = y.in, y.trainset = y.trainset, frac.train = 0.5, reps = 50)
#saveRDS(xval_BayesLASSO_frac50, "xval_BayesLASSO_frac50.RData")

#xval_BayesLASSO_frac80 <- frac.xval.BayesLASSO(g.in = g.in, y.in = y.in, y.trainset = y.trainset, frac.train = 0.8, reps = 50)
#saveRDS(xval_BayesLASSO_frac80, "xval_BayesLASSO_frac80.RData")

#xval_BayesLASSO_frac90 <- frac.xval.BayesLASSO(g.in = g.in, y.in = y.in, y.trainset = y.trainset, frac.train = 0.9, reps = 50)
#saveRDS(xval_BayesLASSO_frac90, "xval_BayesLASSO_frac90.RData")

# RR-BLUP
y.trainset.rr <- trainingset[,c(2,8:(ncol(trainingset)-1))]
y.in.rr <- envdat[,8:ncol(envdat)]

y.trainset.mat <- as.matrix(y.trainset.rr)
y.in.mat <- as.matrix(y.in.rr)

xval_frac_rrblup80 <- frac.xval(g.in = g.in, y.in = y.in.mat, y.trainset = y.trainset.mat, frac.train = 0.8, reps = 50)
saveRDS(xval_frac_rrblup80, "xval_frac_rrblup80.RData")

xval_frac_rrblup90 <- frac.xval(g.in = g.in, y.in = y.in.mat, y.trainset = y.trainset.mat, frac.train = 0.9, reps = 50)
saveRDS(xval_frac_rrblup90, "xval_frac_rrblup90.RData")

# Gaussian kernel
K <- A.mat(g.in)
k_dist <- dist(K)

xval_GAUSS_frac80 <- frac.GAUSS.xval(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = k_dist, frac.train = 0.8, reps = 50)
saveRDS(xval_GAUSS_frac80, "xval_GAUSS_frac80.RData")

xval_GAUSS_frac90 <- frac.GAUSS.xval(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = k_dist, frac.train = 0.9, reps = 50)
saveRDS(xval_GAUSS_frac90, "xval_GAUSS_frac90.RData")

# Exponential Kernel
K.Exp=Kernel_computation(X=g.in, name="exponential", degree=NULL, nL=NULL)
row.names(K.Exp) <- rownames(g.in)
colnames(K.Exp) <- rownames(g.in)
exp_dist<- dist(K.Exp)

xval_EXP_frac80 <- frac.EXP.xval(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = exp_dist, frac.train = 0.80, reps = 50)
saveRDS(xval_EXP_frac80, "xval_EXP_frac80.RData")

xval_EXP_frac90 <- frac.EXP.xval(g.in = g.in, y.in = y.in, y.trainset = y.trainset, k_dist = exp_dist, frac.train = 0.90, reps = 50)
saveRDS(xval_EXP_frac90, "xval_EXP_frac90.RData")



