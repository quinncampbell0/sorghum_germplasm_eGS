
library(ggcorrplot)
library(corrplot)
envdat <- read.table('./input_data/envdat_master_revision.txt', head = T)

row.names(envdat) <- envdat$pi

envdat <- envdat[,8:ncol(envdat)]

envdat2 <- na.omit(envdat)

envdat_scaled <- scale(envdat2)

colSums(is.na(envdat_scaled))

corr_matrix <- cor(envdat_scaled)

ggcorrplot(corr_matrix)
