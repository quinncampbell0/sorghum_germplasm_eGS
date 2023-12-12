
library(ggcorrplot)
library(corrplot)
envdat <- read.csv('G:/My Drive/Sorghum/R/input_data/envdat_master.csv')

row.names(envdat) <- envdat$pi

envdat <- envdat[,8:46]

envdat2 <- na.omit(envdat)

envdat_scaled <- scale(envdat2)

colSums(is.na(envdat_scaled))

corr_matrix <- cor(envdat_scaled)

ggcorrplot(corr_matrix)
