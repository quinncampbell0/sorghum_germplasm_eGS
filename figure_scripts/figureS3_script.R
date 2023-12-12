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
