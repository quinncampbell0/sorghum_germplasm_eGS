library(ggplot2)
library(tidyr)

setwd('G:/My Drive/Sorghum/R/')
all_models <- read.csv("./cross_validation/xval_results_allmodels.csv")

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

panel_A <- ggplot(all_soil, aes(y = r.mean, x = model, color = xval)) +
  # geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 2, show.legend = F) +
  # geom_boxplot(position = 'identity') + 
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15),
        strip.text.x = element_text(size = 10),
        legend.position="none") 


panel_B <- ggplot(all_temp, aes(y = r.mean, x = model, color = xval)) +
  # geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 2, show.legend = F) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15),
        strip.text.x = element_text(size = 10))

panel_C <- ggplot(all_precip, aes(y = r.mean, x = model, color = xval)) +
  # geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 2, show.legend = F) +
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15),
        strip.text.x = element_text(size = 10),
        legend.position="none")

#############################################################
##### Cross Validation - Data Organization and Plotting #####
#############################################################

library(ggplot2)
library(tidyr)

setwd('./revision_scripts/xval_results/')
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


model_list <- list(rrblup_frac50, rrblup_frac80, rrblup_frac90, rrblup_kfold2, rrblup_kfold5, rrblup_kfold10, gauss_frac50, gauss_frac80, gauss_frac90, 
                   gauss_kfold_2, gauss_kfold_5, gauss_kfold_10, EXP_frac50, EXP_frac80, EXP_frac90, EXP_kfold_10, EXP_kfold_2, EXP_kfold_5 
                   #,bayescpi_frac50, bayescpi_frac80, bayescpi_frac90, bayescpi_kfold_2, bayescpi_kfold_5, bayescpi_kfold_10
                   #,BayesLASSO_frac50, BayesLASSO_frac80, BayesLASSO_frac90, BayesLASSO_kfold_2, BayesLASSO_kfold_5
)

model_list1 <- lapply(model_list, na.omit)

all_models <- do.call("rbind", model_list1)

all_models$r.sd <- as.numeric(all_models$r.sd)

all_models$xval <- factor(all_models$xval, levels = c("Two-Fold", "Five-Fold", "Ten-Fold", "Fractional 50", "Fractional 80", "Fractional 90"))

all_models$trait[which(all_models$trait == 'pet_min')] <- 'PET_Driest.Month'
all_models$trait[which(all_models$trait == 'ann_heat_index')] <- 'HeatIndex_AnnMean'
all_models$trait[which(all_models$trait == 'max_heat_index')] <- 'HeatIndex_MaxTemp'
all_models$trait[which(all_models$trait == 'coldindex.mean')] <- 'ColdIndex_AnnMean'
all_models$trait[which(all_models$trait == 'coldindex.min')] <- 'ColdIndex_Cold.Month'

all_models$trait <- factor(all_models$trait, 
                           levels = c('PET_Driest.Month', 'AridityIndex', 'ColdIndex_AnnMean', 'ColdIndex_Cold.Month', 'HeatIndex_AnnMean', 'HeatIndex_MaxTemp'))

#####################################################################################################################################
##  X axis by Model Type

panel_D <- ggplot(all_models, aes(y = r.mean, x = model, color = xval)) +
  # geom_errorbar(aes(x = model, ymin = r.mean-r.sd, ymax = r.mean + r.sd), width = 0.3, position = position_dodge(0.2)) +
  geom_point(size = 2) +
  # geom_boxplot(position = 'identity') + 
  facet_wrap(vars(trait), scales = "free_x", nrow = 1) + 
  geom_hline(yintercept = 0.5, color="red", size = 1.5, linetype = "longdash") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 15),
        strip.text.x = element_text(size = 12)) +
  labs(color = "Cross-Evaluation Method")


#########################################################
#Minicore rrBLUP with fixed effect SNPs
##########################################

library(rrBLUP)
library(hibayes)
library(dplyr)

#setwd('/home/quinnc3/')

source("./fixed_effects/xval_kfold_functions_updated_2.R") # load all functions

gd1 <- read.table("./GenomicData/sorghum_GD_numeric_rrblupformat.txt", head = T) # download from Figshare link in Github
envdat <- read.table('./fixed_effects/envdat_master.txt', head = T) # full environmental dataset
trainingset <- read.table("./fixed_effects/envdat_master_mini.txt", head = T)

gd2 <- gd1[gd1$taxa %in% envdat$gen_id,]
row.names(gd2) <- gd2$taxa # set taxa as rownames
g.in <- gd2[,-1]
g.in <- as.matrix(g.in)

### Load Environmental Data ###
row.names(envdat) <- envdat$gen_id # set gen_id as rownames
row.names(trainingset) <- trainingset$gen_id

#y.trainset <- trainingset[,c(2,5:(ncol(trainingset)))] # select unique identifier and environmental data only
#y.in <- envdat[,c(1,4:ncol(envdat))]

y.trainset <- subset(trainingset, select = c('gen_id', 'pet_min', 'Tmp.Seas', 'Prc.Seas', 'topsoil_pH'))
y.in <- subset(envdat, select = c('gen_id', 'pet_min', 'Tmp.Seas', 'Prc.Seas', 'topsoil_pH'))

# Exclude "altitude" from y.in and y.trainset
y.in <- y.in[, !(colnames(y.in) %in% c("gen_id"))]
y.trainset <- y.trainset[, !(colnames(y.trainset) %in% c("gen_id"))]

print(colnames(y.in))  # Check the column names in y.in
print(colnames(y.trainset))  # Check the column names in y.trainset

### RR-BLUP
#y.trainset.rr <- trainingset[,c(5:(ncol(trainingset)))]
#y.in.rr <- envdat[,5:ncol(envdat)]

y.trainset.mat <- as.matrix(y.trainset)
y.in.mat <- as.matrix(y.in)

#### SNP sets for fixed effects

# Maturity1 and Tannin1 polymorphism (respectively)
snps_MatTan <- c("S6_40286721", 'S4_61667908')

# Top 10 APPA for Harvest Index plasticity with Precipitation of Warmest Quarter as a prior
snps_Prc.Wrm.Q <- c('S10_4862610','S1_27310288','S5_1793658','S10_54329625','S3_5891497','S1_31317641',
                    'S1_42444158','S10_9477501','S10_54326475','S1_43500042')

# Top 10 APPA for relative net root growth with topsoil pH as prior
snps_pH <- c('S1_58013535','S3_32314736','S7_64137481','S2_5378283','S2_3386704','S9_59524753','S3_72006525',
             'S3_60807535','S2_20187646','S2_20257191')

# Top 10 APPA for panicle weight plasticity with growing season length as priors
#snps_GS <- c('S10_59797088','S10_59797069','S4_61520903','S4_55202995','S3_802046','S5_61212379','S5_61212384',
#             'S4_60073808','S1_1132237','S2_56923866')

snp_list <- list(MatTan = snps_MatTan, Prc.Wrm.Q = snps_Prc.Wrm.Q, topsoil_pH = snps_pH) # , grow.seas = snps_GS
list_names <- names(snp_list)

lapply(list_names, function(name){
  # get list element matching name
  snps_of_interest <- snp_list[[name]]
  
  # Subset the genotype matrix for the fixed SNPs
  fixed.snps <- g.in[, snps_of_interest] # Ensure snps_of_interest are column names in g.in
  
  # Ensure that rownames of fixed.snps are aligned with g.in
  rownames(fixed.snps) <- rownames(g.in)
  
  # Run k-fold cross-validation using fixed.snps
  xval_k10_rrblup <- k.xval(g.in = g.in, 
                            y.in = y.in.mat, 
                            y.trainset = y.trainset.mat, 
                            k.fold = 10, 
                            reps = 50, 
                            fixed.snps = fixed.snps) # Pass fixed.snps matrix here
  file_name <- paste0('./fixed_effects/output/xval_rrblup_10fold_fixed_',name,'.RData')
  
  saveRDS(xval_k10_rrblup, file_name)
})

###################################################
### Plot Fixed Effects Cross Validation Results ###
library(stringr)
mattan <- readRDS('./revision_scripts/xval_rrblup_10fold_fixed_MatTan.RData')$xval.result
wrmq <- readRDS('./revision_scripts/xval_rrblup_10fold_fixed_Prc.Wrm.Q.RData')$xval.result
topph <- readRDS('./revision_scripts/xval_rrblup_10fold_fixed_topsoil_pH.RData')$xval.result

ecophys_full_xval <- readRDS('./revision_scripts/xval_results/xval_rrblup_kfold_10.RData')

ecophys_fx <- ecophys_full_xval$xval.result

ecophys_pet <- ecophys_fx %>%
  filter(trait %in% c("pet_min")) %>% 
  mutate(model = 'No Fixed Effects')

full_xval <- read.csv('./revision_scripts/xval_rrblup_kfold10_allmarkers.csv')

full_xval_tr <- full_xval %>%
  filter(trait %in% c("topsoil_pH", "Tmp.Seas", "Prc.Seas")) %>% 
  mutate(model = 'No Fixed Effects')

no_fix <- rbind(ecophys_pet, full_xval_tr)

mattan$model <- 'Maturity1 and Tannin1'
topph$model <- 'Top 10 APPA, Harvest Index Plasticity + Prec Warmest Quarter'
wrmq$model <- 'Top 10 APPA, Net Root Growth + topsoil pH'

### Plot All together
all <- rbind(mattan, topph, wrmq, no_fix)

all$r.mean <- as.numeric(all$r.mean)
all$r.sd <- as.numeric(all$r.sd)

# allow for wrapping labels
all$modlab = str_wrap(all$model, width = 20)

# Factor for labelling order
all$modlab <- factor(all$modlab, level = c('No Fixed Effects', "Maturity1 and\nTannin1", 
                                           "Top 10 APPA, Harvest\nIndex Plasticity +\nPrec Warmest Quarter",
                                           "Top 10 APPA, Net\nRoot Growth +\ntopsoil pH"))

#Same Panel
panel_E <- ggplot(all, aes(x = modlab, y = r.mean, color = trait, group = trait)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
  geom_line(linetype = 'dotted', alpha = 0.9, linewidth = 0.8,position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), 
                width = 0.5, linewidth = 1.25, 
                position = position_dodge(width = 0.3)) +
  labs(
    y = "r",
    x = "Fixed Effects"
  ) 



###############################################################################################
#setwd('G:/My Drive/Sorghum/R/')

ld_list <- seq(0.1, 0.95, by = 0.05)
prune_set <- format(ld_list, nsmall = 2)

result_list <- lapply(prune_set, function(x){
  file <- readRDS(paste0('./revision_scripts/xval_ldprune/xval_rrblup_LD_',x,'.RData'))
  file$xval.result$LD_prune <- x
  return(file$xval.result)
})

#full_xval <- read.csv('./revision_scripts/xval_rrblup_kfold10_allmarkers.csv')

#full_xval_tr <- full_xval %>%
#  filter(trait %in% c("topsoil_pH", "Tmp.Seas", "Prc.Seas")) %>% 
#  mutate(LD_prune = 0)

ecophys_full_xval <- readRDS('./revision_scripts/xval_ldprune/xval_rrblup_kfold_10.RData')

ecophys_fx <- ecophys_full_xval$xval.result

ecophys_fx_tr <- ecophys_fx %>%
  filter(trait %in% c('pet_min', 'ann_heat_index', 'coldindex.mean')) %>% 
  mutate(LD_prune = 1)

xval_ld <- do.call(rbind, result_list)

#xval_ld <- xval_ld %>%
#  filter(trait %in% c('pet_min'))

xval_ld <- rbind(xval_ld, ecophys_fx_tr)

xval_ld$r.mean <- as.numeric(xval_ld$r.mean)
xval_ld$r.sd <- as.numeric(xval_ld$r.sd)
xval_ld$LD_prune <- as.numeric(xval_ld$LD_prune)

ggplot(xval_ld, aes(x = LD_prune, y = r.mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = r.mean-r.sd, ymax = r.mean+r.sd))+
  facet_wrap(~trait, scales = 'free_y')

ggplot(xval_ld, aes(x = LD_prune, y = r.mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = r.mean-r.sd, ymax = r.mean+r.sd), fill = 'blue', alpha = 0.2)+
  facet_wrap(~trait, scales = 'free_y') +
  theme_minimal()

# Only PET min

xval_ld_pet <- xval_ld %>%
  filter(trait %in% c('pet_min'))

# Add in SNP number

prune_set <- seq(0.1, 0.95, by = 0.05)
prune <- format(prune_set, nsmall = 2)

snpnums <- lapply(prune, function(x){
  
  filename <- paste0('./revision_scripts/xval_ldprune/prune_r2_',x,'.prune.in')
  
  snps <- read.table(filename)
  
  return(nrow(snps))
})

snpnums <- unlist(snpnums)

snpnums1 <- c(snpnums, 404627)

xval_ld_pet$SNPs <- snpnums1

xval_ld_pet$LD_prune <- format(xval_ld_pet$LD_prune , nsmall = 2)

panel_F <- ggplot(xval_ld_pet, aes(x = SNPs, y = r.mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), fill = 'blue', alpha = 0.2) +
  geom_point(color = 'black') + 
  geom_text(aes(label = LD_prune), vjust = -0.5, size = 3) +  # Add labels for each point
  theme_minimal() +
  scale_x_continuous(labels = scales::comma) +  # Use full numbers instead of scientific notation
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Angle x-axis labels at 45 degrees
  ) +
  labs(
    y = "r",
    x = "# of SNPs"
  )


#######################################################################################################################
#######################################################################################################################
library(ggpubr)

plot_list1 = list(panel_A, panel_B, panel_C, panel_D)
plot_list2 = list(panel_E, panel_F)

g1 <- ggarrange(plotlist = plot_list1, ncol = 2, nrow = 2, labels = c('A', "B", 'C', 'D'), font.label = list(size = 18))
g2 <- ggarrange(plotlist = plot_list2, ncol = 1, labels = c('E', 'F'), font.label = list(size = 18))
ggarrange(g1, g2, nrow = 4)

plot_list <- list(panel_A, panel_B, panel_C, panel_D, panel_E, panel_F)

pdf(file = './output_figures/validation_graphs.pdf', width = 15, height = 25)
ggarrange(plotlist = plot_list, ncol = 1, nrow = 6, labels = c('A', "B", 'C', 'D', 'E', 'F'), font.label = list(size = 18))
dev.off()

