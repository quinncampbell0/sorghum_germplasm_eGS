###############################################################################
## Figure ED10 - correlation plots

# Spearman correlations for GEAVs for all bioclimatic and ecophysiological indices
library(corrplot)
library(ggcorrplot)
library(ggpubr)

ecophys_GEAV <- read.csv('./revision_scripts/rrblup_GEBV_ecophys.csv')
bioclim_GEAV <- read.csv('./input_data/rrblup_GEBV_data.csv')


all_GEAV <- merge(ecophys_GEAV, bioclim_GEAV, by = "X")

cor_matrix <- cor(all_GEAV[,2:ncol(all_GEAV)], method = "spearman")

corrplot(
  cor_matrix,
  method = "color",
  type = "full",
  tl.cex = 1,
  tl.col = "black",
  tl.srt = 45,
  addCoef.col = "black",
  title = "Spearman Correlation Heatmap",
  number.cex = 0.5
)

panel_A <- ggcorrplot(cor_matrix, 
           lab = TRUE,            # Display correlation coefficients
           lab_size = 2,          # Size of the labels
           method = "square",     # Use circular method for displaying
           colors = c("blue", "white", "red"))
##################################################################
# Panel B

envdat <- read.table('./input_data/envdat_master_revision.txt', head = T)

row.names(envdat) <- envdat$pi

envdat <- envdat[,8:ncol(envdat)]

envdat2 <- na.omit(envdat)

envdat_scaled <- scale(envdat2)

colSums(is.na(envdat_scaled))

corr_matrix <- cor(envdat_scaled)[row.names(cor_matrix),colnames(cor_matrix)]

panel_B <- ggcorrplot(corr_matrix, lab = T, lab_size = 2, )

######
# Combine Plots

ggarrange(panel_A, panel_B, ncol = 1, labels = c('A', 'B'), font.label = list(size = 18))
