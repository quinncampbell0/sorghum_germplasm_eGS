library(rrBLUP)
library(dplyr)

envdat <- read.table("./input_data/envdat_minicore_genomicprediction.txt", head = T)
gd1 <- read.table("./input_data/sorghum_GD_numeric_rrblupformat.txt", head = T)
#gd1<- read.csv("sorghum_GM_testdata.csv")

row.names(gd1) <- gd1$taxa
gd2 <- gd1[gd1$taxa %in% envdat$gen_id,]
gd3 <- gd2[,-1]
g.in <- as.matrix(gd3)

# RR-BLUP
row.names(envdat) <- envdat$gen_id
y.in.rr <- envdat[,8:ncol(envdat)]
minicore_entries <- which(y.in.rr$is.minicore == T)
y.in.rr <- y.in.rr[minicore_entries,]
y.in.rr <- y.in.rr[,-ncol(y.in.rr)]
y.in.mat <- as.matrix(y.in.rr)

train <- row.names(y.in.mat) # Names of training lines
g.train <- g.in[train,] # Set training genos

traits <- colnames(y.in.mat)

# Calculate marker effects

pred <- setdiff(row.names(g.in), train)
g.pred <- g.in[pred,] # Set prediction genos
y.pred <- y.in.rr[pred,]

marker.list <- list()
gebv_df <- data.frame(matrix(nrow = nrow(envdat), ncol = length(traits)))
#row.names(gebv_df) <- row.names(envdat)
colnames(gebv_df) <- traits

for(t in 1:length(traits)){
  trait <- traits[t]
  y.train <- as.matrix(y.in.mat[train,trait])
  
  ###### Run Model: RR-BLUP ######
  solve.out <- mixed.solve(y = y.train,  Z = g.train, SE = F, return.Hinv = F)
  u.hat <- solve.out$u
  
  # Calculate GEBVs and correlations
  GEBV <- g.pred %*% u.hat
  GEBV_train <- g.train %*% u.hat
  #outputs
  marker.list[[t]] <- u.hat
  gebv_df[,t] <- rbind(GEBV, GEBV_train)
}
row.names(gebv_df)[1:nrow(y.pred)] <- row.names(GEBV)
row.names(gebv_df)[(nrow(y.pred) + 1):(nrow(y.pred) + nrow(y.train))] <- row.names(GEBV_train)

write.csv(gebv_df, 'rrblup_GEBV_data.csv')
saveRDS(marker.list, 'rrblup_markereffects_data.csv')


