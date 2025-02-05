# Spearman correlations for GEAVs for all bioclimatic and ecophysiological indices
library(corrplot)
library(ggcorrplot)

ecophys_GEAV <- read.csv('./input_data/rrblup_GEBV_ecophys.csv')
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

ggcorrplot(cor_matrix, 
           type = "lower",        # Display only lower triangle
           lab = TRUE,            # Display correlation coefficients
           lab_size = 2,          # Size of the labels
           method = "square",     # Use circular method for displaying
           colors = c("blue", "white", "red"), # Color gradient
           title = "Spearman Correlation Matrix") # Title
