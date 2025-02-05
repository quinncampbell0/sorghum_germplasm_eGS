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

mattan <- readRDS('./fixed_effects/output/xval_rrblup_10fold_fixed_MatTan.RData')$xval.result
wrmq <- readRDS('./fixed_effects/output/xval_rrblup_10fold_fixed_Prc.Wrm.Q.RData')$xval.result
topph <- readRDS('./fixed_effects/output/xval_rrblup_10fold_fixed_topsoil_pH.RData')$xval.result

ecophys_full_xval <- readRDS('./fixed_effects/xval_results/xval_rrblup_kfold_10.RData')

ecophys_fx <- ecophys_full_xval$xval.result

ecophys_pet <- ecophys_fx %>%
  filter(trait %in% c("pet_min")) %>% 
  mutate(model = 'No Fixed Effects')

full_xval <- read.csv('./fixed_effects/xval_rrblup_kfold10_allmarkers.csv')

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
ggplot(all, aes(x = modlab, y = r.mean, color = trait, group = trait)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
  geom_line(linetype = 'dotted', alpha = 0.9, linewidth = 0.8,position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), 
                width = 0.5, linewidth = 1.25, 
                position = position_dodge(width = 0.3)) +
  labs(
    y = "r",
    x = "Fixed Effects"
  ) 

# Only MatTan
mattan_p <- all %>% filter(modlab %in% c('No Fixed Effects', "Maturity1 and\nTannin1"))

ggplot(mattan_p, aes(x = modlab, y = r.mean, color = trait, group = trait)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
  geom_line(linetype = 'dotted', alpha = 0.9, linewidth = 0.8,position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), 
                width = 0.5, linewidth = 1.25, 
                position = position_dodge(width = 0.3)) +
  labs(
    y = "r",
    x = "Fixed Effects"
  ) 


# Separate panels by env trait
ggplot(all, aes(x = modlab, y = r.mean, group = trait, color = trait)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
 
  geom_errorbar(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), 
                width = 0.4, linewidth = 1.25, 
                position = position_dodge(width = 0.3)) +
  labs(
    y = "r",
    x = "Fixed Effects"
  ) +
  facet_wrap(~trait, ncol = 1)

### Plot separately and combine with ggarrange

mattan <- rbind(mattan, no_fix)
topph <- rbind(topph, no_fix)
wrmq<- rbind(topph, no_fix)

mattan$r.mean <- as.numeric(mattan$r.mean)
mattan$r.sd <- as.numeric(mattan$r.sd)

p_mattan <- ggplot(mattan, aes(x = model, y = r.mean, color = trait)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = r.mean-r.sd, ymax = r.mean+r.sd), width = 0.25, linewidth = 1.5)

wrmq$r.mean <- as.numeric(wrmq$r.mean)
wrmq$r.sd <- as.numeric(wrmq$r.sd)

p_wrmq <- ggplot(wrmq, aes(x = model, y = r.mean, color = trait)) +
  geom_point() +
  geom_errorbar(aes(ymin = r.mean-r.sd, ymax = r.mean+r.sd))

topph$r.mean <- as.numeric(topph$r.mean)
topph$r.sd <- as.numeric(topph$r.sd)

p_topph <- ggplot(topph, aes(x = model, y = r.mean, color = trait)) +
  geom_point() +
  geom_errorbar(aes(ymin = r.mean-r.sd, ymax = r.mean+r.sd))


ggarrange(p_mattan, p_wrmq, p_topph, ncol = 3)
