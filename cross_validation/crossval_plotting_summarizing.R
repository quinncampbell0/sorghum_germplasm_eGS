#############################################################
##### Cross Validation - Data Organization and Plotting #####
#############################################################

library(ggplot2)
library(tidyr)
#############################################################################
#### RR-BLUP - Fractional Cross Validation ####
rrblup_frac50 <- readRDS("xval_frac_rrblup50.RData")
rrblup_frac50 <- rrblup_frac50$xval.result
rrblup_frac50$r.mean <- as.numeric(rrblup_frac50$r.mean)

rrblup_frac80<- readRDS("xval_frac_rrblup80.RData")
rrblup_frac80<- rrblup_frac80$xval.result
rrblup_frac80$r.mean <- as.numeric(rrblup_frac80$r.mean)

rrblup_frac90 <- readRDS("xval_frac_rrblup90.RData")
rrblup_frac90 <- rrblup_frac90$xval.result
rrblup_frac90$r.mean <- as.numeric(rrblup_frac90$r.mean)

#### RR-BLUP - K-Fold Cross-Validation ####

rrblup_kfold2 <- readRDS("xval_rrblup_kfold_2.RData")
rrblup_kfold5 <- readRDS("xval_rrblup_kfold_5.RData")
rrblup_kfold10 <- readRDS("xval_rrblup_kfold_10.RData")

rrblup_kfold2 <- rrblup_kfold2$xval.result
rrblup_kfold5 <- rrblup_kfold5$xval.result
rrblup_kfold10 <- rrblup_kfold10$xval.result

rrblup_kfold2$r.mean <- as.numeric(rrblup_kfold2$r.mean)
rrblup_kfold5$r.mean <- as.numeric(rrblup_kfold5$r.mean)
rrblup_kfold10$r.mean <- as.numeric(rrblup_kfold10$r.mean)

#############################################################################
# Gaussian Kernel  - Fractional

gauss_frac50 <- readRDS('xval_GAUSS_frac50.RData')
gauss_frac50 <- gauss_frac50$xval.result
gauss_frac50$r.mean <- as.numeric(gauss_frac50$r.mean)

gauss_frac80 <- readRDS('xval_GAUSS_frac80.RData')
gauss_frac80 <- gauss_frac80$xval.result
gauss_frac80$r.mean <- as.numeric(gauss_frac80$r.mean)

gauss_frac90 <- readRDS('xval_GAUSS_frac90.RData')
gauss_frac90 <- gauss_frac90$xval.result
gauss_frac90$r.mean <- as.numeric(gauss_frac90$r.mean)
# Gaussian Kernel - K-fold

gauss_kfold_2 <- readRDS('xval_GAUSS_kfold_2.RData')
gauss_kfold_2 <- gauss_kfold_2$xval.result
gauss_kfold_2$r.mean <- as.numeric(gauss_kfold_2$r.mean)

gauss_kfold_5 <- readRDS('xval_GAUSS_kfold_5.RData')
gauss_kfold_5 <- gauss_kfold_5$xval.result
gauss_kfold_5$r.mean <- as.numeric(gauss_kfold_5$r.mean)

gauss_kfold_10 <- readRDS('xval_GAUSS_kfold_10.RData')
gauss_kfold_10 <- gauss_kfold_10$xval.result
gauss_kfold_10$r.mean <- as.numeric(gauss_kfold_10$r.mean)

################################################################################
# Exponential kernel - Fractional

EXP_frac50 <- readRDS('xval_EXP_frac50.RData')
EXP_frac50 <- EXP_frac50$xval.result
EXP_frac50$r.mean <- as.numeric(EXP_frac50$r.mean)

EXP_frac80 <- readRDS('xval_EXP_frac80.RData')
EXP_frac80 <- EXP_frac80$xval.result
EXP_frac80$r.mean <- as.numeric(EXP_frac80$r.mean)

EXP_frac90 <- readRDS('xval_EXP_frac90.RData')
EXP_frac90 <- EXP_frac90$xval.result
EXP_frac90$r.mean <- as.numeric(EXP_frac90$r.mean)

# Exponential Kernel - K-fold

EXP_kfold_2 <- readRDS('xval_EXP_kfold_2.RData')
EXP_kfold_2 <- EXP_kfold_2$xval.result
EXP_kfold_2$r.mean <- as.numeric(EXP_kfold_2$r.mean)

EXP_kfold_5 <- readRDS('xval_EXP_kfold_5.RData')
EXP_kfold_5 <- EXP_kfold_5$xval.result
EXP_kfold_5$r.mean <- as.numeric(EXP_kfold_5$r.mean)

EXP_kfold_10 <- readRDS('xval_EXP_kfold_10.RData')
EXP_kfold_10 <- EXP_kfold_10$xval.result
EXP_kfold_10$r.mean <- as.numeric(EXP_kfold_10$r.mean)

# BayesCpi - fractional

bayescpi_frac50 <- readRDS("xval_BayesCpi_frac50.RData")
bayescpi_frac50 <- bayescpi_frac50$xval.result
bayescpi_frac50$r.mean <- as.numeric(bayescpi_frac50$r.mean)

bayescpi_frac80 <- readRDS("xval_BayesCpi_frac80.RData")
bayescpi_frac80 <- bayescpi_frac80$xval.result
bayescpi_frac80$r.mean <- as.numeric(bayescpi_frac80$r.mean)

bayescpi_frac90 <- readRDS("xval_BayesCpi_frac90.RData")
bayescpi_frac90 <- bayescpi_frac90$xval.result
bayescpi_frac90$r.mean <- as.numeric(bayescpi_frac90$r.mean)

# BayesCpi - K-fold
bayescpi_kfold_2 <- readRDS('xval_BayesCpi_kfold_2.RData')
bayescpi_kfold_2 <- bayescpi_kfold_2$xval.result
bayescpi_kfold_2$r.mean <- as.numeric(bayescpi_kfold_2$r.mean)

bayescpi_kfold_5 <- readRDS('xval_BayesCpi_kfold_5.RData')
bayescpi_kfold_5 <- bayescpi_kfold_5$xval.result
bayescpi_kfold_5$r.mean <- as.numeric(bayescpi_kfold_5$r.mean)

bayescpi_kfold_10 <- readRDS('xval_BayesCpi_kfold_10.RData')
bayescpi_kfold_10 <- bayescpi_kfold_10$xval.result
bayescpi_kfold_10$r.mean <- as.numeric(bayescpi_kfold_10$r.mean)

# BayesLASSO - fractional

#BayesLASSO_frac50 <- readRDS("xval_BayesLASSO_frac50.RData")
#BayesLASSO_frac50 <- BayesLASSO_frac50$xval.result
#BayesLASSO_frac50$r.mean <- as.numeric(BayesLASSO_frac50$r.mean)

#BayesLASSO_frac80 <- readRDS("xval_BayesLASSO_frac80.RData")
#BayesLASSO_frac80 <- BayesLASSO_frac80$xval.result
#BayesLASSO_frac80$r.mean <- as.numeric(BayesLASSO_frac80$r.mean)

#BayesLASSO_frac90 <- readRDS("xval_BayesLASSO_frac90.RData")
#BayesLASSO_frac90 <- BayesLASSO_frac90$xval.result
#BayesLASSO_frac90$r.mean <- as.numeric(BayesLASSO_frac90$r.mean)

# BayesLASSO - K-fold
#BayesLASSO_kfold_2 <- readRDS('xval_BayesLASSO_kfold_2.RData')
#BayesLASSO_kfold_2 <- BayesLASSO_kfold_2$xval.result
#BayesLASSO_kfold_2$r.mean <- as.numeric(BayesLASSO_kfold_2$r.mean)

#BayesLASSO_kfold_5 <- readRDS('xval_BayesLASSO_kfold_5.RData')
#BayesLASSO_kfold_5 <- BayesLASSO_kfold_5$xval.result
#BayesLASSO_kfold_5$r.mean <- as.numeric(BayesLASSO_kfold_5$r.mean)

#BayesLASSO_kfold_10 <- readRDS('xval_BayesLASSO_kfold_10.RData')
#BayesLASSO_kfold_10 <- BayesLASSO_kfold_10$xval.result
#BayesLASSO_kfold_10$r.mean <- as.numeric(BayesLASSO_kfold_10$r.mean)

 ################################################################################
# Organize all dataframes for merging

## Rename model names
rrblup_frac50$model <- "rrBLUP"
rrblup_frac80$model <- "rrBLUP"
rrblup_frac90$model <- "rrBLUP"
rrblup_kfold2$model <- "rrBLUP"
rrblup_kfold5$model <- "rrBLUP"
rrblup_kfold10$model <- "rrBLUP"

gauss_frac50$model <- "Gaussian Kernel"
gauss_frac80$model <- "Gaussian Kernel"
gauss_frac90$model <- "Gaussian Kernel"
gauss_kfold_2$model <- "Gaussian Kernel"
gauss_kfold_5$model <- "Gaussian Kernel"
gauss_kfold_10$model <- "Gaussian Kernel"

EXP_frac50$model <- "Exponential Kernel"
EXP_frac80$model <- "Exponential Kernel"
EXP_frac90$model <- "Exponential Kernel"
EXP_kfold_2$model <- "Exponential Kernel"
EXP_kfold_5$model <- "Exponential Kernel"
EXP_kfold_10$model <- "Exponential Kernel"

bayescpi_frac50$model <- "BayesCpi"
bayescpi_frac80$model <- "BayesCpi"
bayescpi_frac90$model <- "BayesCpi"
bayescpi_kfold_2$model <- "BayesCpi"
bayescpi_kfold_5$model <- "BayesCpi"
bayescpi_kfold_10$model <- "BayesCpi"

#BayesLASSO_frac50$model <- "BayesLASSO"
#BayesLASSO_frac80$model <- "BayesLASSO"
#BayesLASSO_frac90$model <- "BayesLASSO"
#BayesLASSO_kfold_2$model <- "BayesLASSO"
#BayesLASSO_kfold_5$model <- "BayesLASSO"
#BayesLASSO_kfold_10$model <- "BayesLASSO"


## Input xval type
rrblup_frac50$xval <- "Fractional 50"
rrblup_frac80$xval <- "Fractional 80"
rrblup_frac90$xval <- "Fractional 90"
rrblup_kfold2$xval <- "Two-Fold"
rrblup_kfold5$xval <- "Five-Fold"
rrblup_kfold10$xval <- "Ten-Fold"

gauss_frac50$xval <- "Fractional 50"
gauss_frac80$xval <- "Fractional 80"
gauss_frac90$xval <- "Fractional 90"
gauss_kfold_2$xval <- "Two-Fold"
gauss_kfold_5$xval <- "Five-Fold"
gauss_kfold_10$xval <- "Ten-Fold"

EXP_frac50$xval <- "Fractional 50"
EXP_frac80$xval <- "Fractional 80"
EXP_frac90$xval <- "Fractional 90"
EXP_kfold_2$xval <- "Two-Fold"
EXP_kfold_5$xval <- "Five-Fold"
EXP_kfold_10$xval <- "Ten-Fold"

bayescpi_frac50$xval <- "Fractional 50"
bayescpi_frac80$xval <- "Fractional 80"
bayescpi_frac90$xval <- "Fractional 90"
bayescpi_kfold_2$xval <- "Two-Fold"
bayescpi_kfold_5$xval <- "Five-Fold"
bayescpi_kfold_10$xval <- "Ten-Fold"

#BayesLASSO_frac50$xval <- "Fractional 50"
#BayesLASSO_frac80$xval <- "Fractional 80"
#BayesLASSO_frac90$xval <- "Fractional 90"
#BayesLASSO_kfold_2$xval <- "Two-Fold"
#BayesLASSO_kfold_5$xval <- "Five-Fold"
#BayesLASSO_kfold_10$xval <- "Ten-Fold"

model_list <- list(rrblup_frac50, rrblup_frac80, rrblup_frac90, rrblup_kfold2, rrblup_kfold5, rrblup_kfold10, gauss_frac50, gauss_frac80, gauss_frac90, 
                   gauss_kfold_2, gauss_kfold_5, gauss_kfold_10, EXP_frac50, EXP_frac80, EXP_frac90, EXP_kfold_10, EXP_kfold_2, EXP_kfold_5, 
                   bayescpi_frac50, bayescpi_frac80, bayescpi_frac90, bayescpi_kfold_2, bayescpi_kfold_5, bayescpi_kfold_10)
                   #, BayesLASSO_frac50, BayesLASSO_frac80, BayesLASSO_frac90, BayesLASSO_kfold_2, BayesLASSO_kfold_5)

model_list1 <- lapply(model_list, na.omit)

all_models <- do.call("rbind", model_list1)

all_models$r.sd <- as.numeric(all_models$r.sd)

all_models$xval <- factor(all_models$xval, levels = c("Two-Fold", "Five-Fold", "Ten-Fold", "Fractional 50", "Fractional 80", "Fractional 90"))

all_soil <- all_models[all_models$trait %in% c('topsoil_clay','subsoil_clay','topsoil_silt', "subsoil_silt", 'topsoil_sand',
                                                 "subsoil_sand", 'topsoil_cec', 'subsoil_cec', 'topsoil_nitrogen', 'subsoil_nitrogen',
                                                 'topsoil_soc', 'subsoil_soc', 'topsoil_pH','subsoil_pH'),]
all_precip <- all_models[all_models$trait %in% c('prec_Q1', 'prec_Q2', 'prec_Q3', 'prec_Q4', 'Prc.Wrm.Q', 'Prc.Wet.Q',
                                                 'Prc.Wet.M', 'Prc.Seas', 'Prc.Dry.Q', 'Prc.Dry.M', 'Prc.Cld.Q', 'Ann.Prc'),]
all_temp <- all_models[all_models$trait %in% c('tmean_Q1', 'tmean_Q2', 'tmean_Q3', 'tmean_Q4', 'Tmp.Seas','Min.Tmp.Wrm.Q',
                                               'Mean.Tmp.Wrm.Q','Mean.Tmp.Wet.Q','Mean.Tmp.Dry.Q','Mean.Tmp.Cld.Q','Mean.Diurn.Rng',
                                               'Max.Tmp.Wrm.M','Ann.Mean.Tmp'),]



#####################################################################################################################################
##  X axis by Model Type

ggplot(all_soil, aes(y = r.mean, x = model, color = xval)) +
  # geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  # geom_boxplot(position = 'identity') + 
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15),
        axis.text.y = element_text(size = 15)) +
  labs(color = "Cross-Evaluation Method")
  

ggplot(all_temp, aes(y = r.mean, x = model, color = xval)) +
  # geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15),
        axis.text.y = element_text(size = 15))

ggplot(all_precip, aes(y = r.mean, x = model, color = xval)) +
  # geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15),
        axis.text.y = element_text(size = 15))

########################################################################################################################################
## X axis by cross evaluation method

ggplot(all_soil, aes(y = r.mean, x = factor(xval, level = c("Fractional 50", "Fractional 80", "Fractional 90", "Two-Fold", "Five-Fold", "Ten-Fold")), 
                     color = model)) +
  #geom_errorbar(aes(x = factor(xval, level = c("Fractional", "Two-Fold", "Five-Fold", "Ten-Fold")), ymin = r.mean-r.sd, 
  #                  ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15),
        axis.text.y = element_text(size = 15))

ggplot(all_temp, aes(y = r.mean, x = factor(xval, level = c("Fractional 50", "Fractional 80", "Fractional 90", "Two-Fold", "Five-Fold", "Ten-Fold")), 
                     color = model)) +
  #geom_errorbar(aes(x = factor(xval, level = c("Fractional", "Two-Fold", "Five-Fold", "Ten-Fold")), ymin = r.mean-r.sd, 
  #                  ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15),
        axis.text.y = element_text(size = 15))

ggplot(all_precip, aes(y = r.mean, x = factor(xval, level = c("Fractional 50", "Fractional 80", "Fractional 90", "Two-Fold", "Five-Fold", "Ten-Fold")), 
                       color = model)) +
  #geom_errorbar(aes(x = factor(xval, level = c("Fractional", "Two-Fold", "Five-Fold", "Ten-Fold")), ymin = r.mean-r.sd, 
  #                  ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15),
        axis.text.y = element_text(size = 15))


#######################################################################################################################################
## Plot by # of folds

kfold_models <- all_models[all_models$xval %in% c("Two-Fold","Five-Fold","Ten-Fold"),]

kfold_soil <- kfold_models[kfold_models$trait %in% c('topsoil_clay','subsoil_clay','topsoil_silt', "subsoil_silt", 'topsoil_sand',
                                               "subsoil_sand", 'topsoil_cec', 'subsoil_cec', 'topsoil_nitrogen', 'subsoil_nitrogen',
                                               'topsoil_soc', 'subsoil_soc', 'topsoil_pH','subsoil_pH'),]
kfold_precip <- kfold_models[kfold_models$trait %in% c('prec_Q1', 'prec_Q2', 'prec_Q3', 'prec_Q4', 'Prc.Wrm.Q', 'Prc.Wet.Q',
                                                 'Prc.Wet.M', 'Prc.Seas', 'Prc.Dry.Q', 'Prc.Dry.M', 'Prc.Cld.Q', 'Ann.Prc'),]
kfold_temp <- kfold_models[kfold_models$trait %in% c('tmean_Q1', 'tmean_Q2', 'tmean_Q3', 'tmean_Q4', 'Tmp.Seas','Min.Tmp.Wrm.Q',
                                               'Mean.Tmp.Wrm.Q','Mean.Tmp.Wet.Q','Mean.Tmp.Dry.Q','Mean.Tmp.Cld.Q','Mean.Diurn.Rng',
                                               'Max.Tmp.Wrm.M','Ann.Mean.Tmp'),]

ggplot(kfold_temp, aes(y = r.mean, x = factor(xval, level = c("Two-Fold", "Five-Fold", "Ten-Fold")), color = trait)) +
  geom_point(size = 3) +
  ylim(0,1) 

ggplot(kfold_precip, aes(y = r.mean, x = factor(xval, level = c("Two-Fold", "Five-Fold", "Ten-Fold")), color = trait)) +
  geom_point(size = 3) +
  ylim(0,1)

ggplot(kfold_soil, aes(y = r.mean, x = factor(xval, level = c("Two-Fold", "Five-Fold", "Ten-Fold")), color = trait)) +
  geom_point(size = 3) +
  ylim(0,1)

# Facet by number of folds, model type on x axis

ggplot(kfold_temp, aes(y = r.mean, x = model, color = trait)) +
  geom_point(size = 3) +
  ylim(0,1) +
  facet_wrap(~factor(xval, level = c("Two-Fold", "Five-Fold", "Ten-Fold")), scales="free_x", nrow = 1)

######### Plot by Fractional % ############

frac_models <- all_models[all_models$xval %in% c("Fractional 50","Fractional 80","Fractional 90"),]

frac_soil <- frac_models[frac_models$trait %in% c('topsoil_clay','subsoil_clay','topsoil_silt', "subsoil_silt", 'topsoil_sand',
                                                     "subsoil_sand", 'topsoil_cec', 'subsoil_cec', 'topsoil_nitrogen', 'subsoil_nitrogen',
                                                     'topsoil_soc', 'subsoil_soc', 'topsoil_pH','subsoil_pH'),]
frac_precip <- frac_models[frac_models$trait %in% c('prec_Q1', 'prec_Q2', 'prec_Q3', 'prec_Q4', 'Prc.Wrm.Q', 'Prc.Wet.Q',
                                                       'Prc.Wet.M', 'Prc.Seas', 'Prc.Dry.Q', 'Prc.Dry.M', 'Prc.Cld.Q', 'Ann.Prc'),]
frac_temp <- frac_models[frac_models$trait %in% c('tmean_Q1', 'tmean_Q2', 'tmean_Q3', 'tmean_Q4', 'Tmp.Seas','Min.Tmp.Wrm.Q',
                                                     'Mean.Tmp.Wrm.Q','Mean.Tmp.Wet.Q','Mean.Tmp.Dry.Q','Mean.Tmp.Cld.Q','Mean.Diurn.Rng',
                                                     'Max.Tmp.Wrm.M','Ann.Mean.Tmp'),]

#Plot by % 
ggplot(frac_temp, aes(y = r.mean, x = factor(xval, level = c("Fractional 50","Fractional 80","Fractional 90")), color = trait)) +
  geom_point(size = 3) +
  ylim(0,1) 

ggplot(frac_precip, aes(y = r.mean, x = factor(xval, level = c("Fractional 50","Fractional 80","Fractional 90")), color = trait)) +
  geom_point(size = 3) +
  ylim(0,1)

ggplot(frac_soil, aes(y = r.mean, x = factor(xval, level = c("Fractional 50","Fractional 80","Fractional 90")), color = trait)) +
  geom_point(size = 3) +
  ylim(0,1)


###############################################################################################
#Write supplementary table
all_models$xval <- as.character(all_models$xval)
all_models$xval[which(all_models$xval == "Fractional 50")] <- "frac50"
all_models$xval[which(all_models$xval == "Fractional 80")] <- "frac80"
all_models$xval[which(all_models$xval == "Fractional 90")] <- "frac90"
all_models$xval[which(all_models$xval == "Two-Fold")] <- "kfold2"
all_models$xval[which(all_models$xval == "Five-Fold")] <- "kfold5"
all_models$xval[which(all_models$xval == "Ten-Fold")] <- "kfold10"

all_models$model[which(all_models$model == "Gaussian Kernel")] <- "gblup.gaussian"
all_models$model[which(all_models$model == "Exponential Kernel")] <- "gblup.exponential"

all_models <- all_models[,-4]

all_models_long <- pivot_wider(all_models,
                               names_from = c("model", "xval"), 
                               values_from = "r.mean"
)
write.csv(all_models_long, "supplementary_table_xval_prediction.csv")
