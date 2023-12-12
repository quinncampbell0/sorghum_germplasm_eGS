


## Environmental PCA ##
# Using climate and soil variables
library('corrr')
library(ggplot2)
library("FactoMineR")
library(factoextra)
library(ggcorrplot)
library(tidyverse)
library(magrittr)
library(cluster)
library(cluster.datasets)
library(cowplot)
library(NbClust)
library(clValid)
library(ggfortify)
library(clustree)
library(dendextend)
library(corrplot)
library(GGally)
library(ggiraphExtra)
library(knitr)
library(kableExtra)
library(ggExtra)

envdat <- read.csv('./input_data/envdat_master.csv')

# remove soil variables
envdat1 <- envdat[,-c(1, 3:21)]

# assign gen_id as rownames, remove id column for PCA
row.names(envdat1) <- envdat1$gen_id
envdat1 <- envdat1[,-1]

colSums(is.na(envdat1))

# remove accessions with missing data
envdat2 <- envdat1[-which(is.na(envdat1$Ann.Mean.Tmp)),]
colSums(is.na(envdat_scaled))

# normalize data
envdat_scaled <- scale(envdat2)
head(envdat_scaled)

# Plot correlation matrix
corr_matrix <- cor(envdat_scaled)
ggcorrplot(corr_matrix)

envdat_scaled <- envdat_scaled[,-c(18:25)] # attempt to remove highly correlated variables, ie quarterly precip/temp

corr_matrix <- cor(envdat_scaled)
ggcorrplot(corr_matrix)


# Create PCA
#data.pca <- princomp(corr_matrix)
#summary(data.pca)

#data.pca$loadings[, 1:4]

res.pca <- PCA(envdat_scaled,  graph = FALSE)

saveRDS(res.pca, "./cluster_analysis/environmental_pca_output.rds")
# Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

# Visualizations
#fviz_eig(data.pca, addlabels = TRUE)
#fviz_pca_var(data.pca, col.var = "black")

#fviz_cos2(data.pca, choice = "var", axes = 1:2)

#fviz_pca_var(data.pca, col.var = "cos2",
#             gradient.cols = c("black", "orange", "green"),
#             repel = TRUE)

# Extract the results for variables
var <- get_pca_var(res.pca)
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)

# Control variable colors using their contributions to the principle axis
var_contrib_plot <- fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
) + theme_minimal() + ggtitle("Variables - PCA")

saveRDS(var_contrib_plot, "./cluster_analysis/plot_env.pca_variable_contributions.rds")

# Plot PCA

p<-ggplot(data=envdat_pca, aes(x=EV1, y=EV2,group=region,shape=region,color=region))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  scale_color_manual(values = c("black","blue","red", "green", "grey","orange")) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (6.71%)") +
  ylab("PC 2 (4.43%)")
p
p + theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 30))

# Extract the results for variables
var <- get_pca_var(res.pca)
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# Control variable colors using their contributions to the principle axis
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
) + theme_minimal() + ggtitle("Variables - PCA")
######################################################################################
##### Clustering analysis

##### Hierarchical clustering
d <- dist(envdat_scaled, method = "euclidean")

# Hierarchical clustering using Ward's method
res.hc <- hclust(d, method = "ward.D2" )
# Cut tree into 5 groups
grp <- cutree(res.hc, k = 6)
# Visualize
plot(res.hc, cex = 0.6) # plot tree
rect.hclust(res.hc, k = 6, border = 2:5) # add rectangle

# Plot different numbers of k clusters
kmean_calc <- function(df, ...){
  kmeans(df, scaled = ..., nstart = 30)
}
km2 <- kmean_calc(envdat_scaled, 2)
km3 <- kmean_calc(envdat_scaled, 3)
km4 <- kmeans(envdat_scaled, 4)
km5 <- kmeans(envdat_scaled, 5)
km6 <- kmeans(envdat_scaled, 6)
km7 <- kmeans(envdat_scaled, 7)
km8 <- kmeans(envdat_scaled, 8)
km9 <- kmeans(envdat_scaled, 9)
km10 <- kmeans(envdat_scaled, 10)
km11 <- kmeans(envdat_scaled, 11)
p1 <- fviz_cluster(km2, data = envdat_scaled, ellipse.type = "convex", geom = "point") + theme_minimal() + ggtitle("k = 2") 
p2 <- fviz_cluster(km3, data = envdat_scaled, ellipse.type = "convex", geom = "point") + theme_minimal() + ggtitle("k = 3")
p3 <- fviz_cluster(km4, data = envdat_scaled, ellipse.type = "convex", geom = "point") + theme_minimal() + ggtitle("k = 4")
p4 <- fviz_cluster(km5, data = envdat_scaled, ellipse.type = "convex", geom = "point") + theme_minimal() + ggtitle("k = 5")
p5 <- fviz_cluster(km6, data = envdat_scaled, ellipse.type = "convex", geom = "point") + theme_minimal() + ggtitle("k = 6")
p6 <- fviz_cluster(km7, data = envdat_scaled, ellipse.type = "convex", geom = "point") + theme_minimal() + ggtitle("k = 7")
plot_grid(p1, p2, p3, p4, p5, p6, labels = c("k2", "k3", "k4", "k5", "k6", "k7"))

##### Elbow Method - suggests 5 clusters #####
set.seed(31)
# function to compute total within-cluster sum of squares
fviz_nbclust(envdat_scaled, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")

##### Gap Statistic - suggests 23 cluster #####
gap_stat <- clusGap(envdat_scaled, FUN = kmeans, nstart = 30, K.max = 24, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")

###### Silhouette plot - suggests 3 clusters #####
fviz_nbclust(envdat_scaled, kmeans, method = "silhouette", k.max = 24) + theme_minimal() + ggtitle("The Silhouette Plot")

###### Sum of squares method - suggests 7/8 clusters #####
ssc <- data.frame(
  kmeans = c(2,3,4,5,6,7,8),
  within_ss = c(mean(km2$withinss), mean(km3$withinss), mean(km4$withinss), mean(km5$withinss), mean(km6$withinss), mean(km7$withinss), mean(km8$withinss)),
  between_ss = c(km2$betweenss, km3$betweenss, km4$betweenss, km5$betweenss, km6$betweenss, km7$betweenss, km8$betweenss)
)
ssc %<>% gather(., key = "measurement", value = value, -kmeans)
#ssc$value <- log10(ssc$value)
ssc %>% ggplot(., aes(x=kmeans, y=log10(value), fill = measurement)) + geom_bar(stat = "identity", position = "dodge") + ggtitle("Cluster Model Comparison") + xlab("Number of Clusters") + ylab("Log10 Total Sum of Squares") + scale_x_discrete(name = "Number of Clusters", limits = c("0", "2", "3", "4", "5", "6", "7", "8"))

###### NbClust Method - suggests 2 clusters #####


res.nbclust <- NbClust(envdat_scaled, distance = "euclidean",
                       min.nc = 5, max.nc = 9, 
                       method = "complete", index ="all")
factoextra::fviz_nbclust(res.nbclust) + theme_minimal() + ggtitle("NbClust's optimal number of clusters")


##### Clustree method - #####
tmp <- NULL
for (k in 1:11){
  tmp[k] <- kmeans(envdat_scaled, k, nstart = 30)
}
df <- data.frame(tmp)
# add a prefix to the column names
colnames(df) <- seq(1:11)
colnames(df) <- paste0("k",colnames(df))
# get individual PCA
df.pca <- prcomp(df, center = TRUE, scale. = FALSE)
ind.coord <- df.pca$x
ind.coord <- ind.coord[,1:2]
df <- bind_cols(as.data.frame(df), as.data.frame(ind.coord))
clustree(df, prefix = "k")

df_subset <- df %>% select(1:8,12:13)
clustree_overlay(df_subset, prefix = "k", x_value = "PC1", y_value = "PC2")

overlay_list <- clustree_overlay(df_subset, prefix = "k", x_value = "PC1",
                                 y_value = "PC2", plot_sides = TRUE)
overlay_list$x_side
overlay_list$y_side

# Choosing an algorithm
intern <- clValid(envdat_scaled, nClust = 5:9, 
                  clMethods = c("hierarchical","kmeans","pam"), validation = "internal")
# Summary
summary(intern) %>% kable() %>% kable_styling()


##### Extracting cluster features #####
# Compute dissimilarity matrix with euclidean distances
d <- dist(envdat_scaled, method = "euclidean")
# Hierarchical clustering using Ward's method
res.hc <- hclust(d, method = "ward.D2" )
# Cut tree into 6 groups
grp <- cutree(res.hc, k = 6)
# Visualize
plot(res.hc, cex = 0.6) # plot tree
rect.hclust(res.hc, k = 6, border = 2:5) # add rectangle
# export cluster groups 
grp1 <- as.data.frame(grp)
write.csv(grp, "cluster_groups_hclust_k6.csv", row.names = T)


# Execution of k-means with k=6
# final <- kmeans(envdat_scaled, 6, nstart = 30)
final <- hcut(envdat_scaled, k = 6, hc_method = "ward.D2")

pca_envdat_k6 <- fviz_cluster(final, geom = "point") + 
  theme_minimal() + 
  ggtitle("k = 6") +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")

pca_envdat_k6

hclust_ass <- data.frame(gen_id=names(final$cluster), cluster = as.factor(final$cluster))

hclust_ass %>% group_by(cluster) %>% summarise(count = n())

saveRDS(pca_envdat_k6, "./cluster_analysis/plot_env.pca_clusters_k6.rds")

str(grp1)
clust_ass <- as.factor(grp1$grp)

fviz_pca_biplot(res.pca, label = "var", habillage=grp1$grp,
                addEllipses=FALSE,
                ggtheme = theme_minimal())

# descriptive statistics
as.data.frame(envdat_scaled) %>% mutate(Cluster = final$cluster) %>% group_by(Cluster) %>% summarise_all("mean") %>% kable() %>% kable_styling()

# visualize distributions
envdat_df <- as.data.frame(envdat_scaled) %>% rownames_to_column()
cluster_pos <- as.data.frame(final$cluster) %>% rownames_to_column()
colnames(cluster_pos) <- c("rowname", "cluster")
envdat_final <- inner_join(cluster_pos, envdat_df)

ggRadar(envdat_final[-1], aes(group = cluster), rescale = FALSE, legend.position = "none", size = 1, interactive = FALSE, use.label = TRUE) + 
  facet_wrap(~cluster) + 
  scale_y_discrete(breaks = NULL) + # don't show ticks
  theme(axis.text.x = element_text(size = 10)) 

# Correlation plot
envdat_df <- as.data.frame(envdat_scaled)
envdat_df$cluster <- final$cluster
envdat_df$cluster <- as.character(envdat_df$cluster)
ggpairs(envdat_df, c(1:3,6,7,10,11,14,15), mapping = ggplot2::aes(color = cluster, alpha = 0.5), 
        diag = list(continuous = wrap("densityDiag")), 
        lower=list(continuous = wrap("points", alpha=0.9)))

# plot specific graphs from previous matrix with scatterplot
g <- ggplot(envdat_df, aes(x = Ann.Prc, y = Ann.Mean.Tmp, color = cluster)) +
  geom_point() +
  theme(legend.position = "bottom")
ggExtra::ggMarginal(g, type = "histogram", bins = 20, color = "grey", fill = "blue")
b <- ggplot(envdat_df, aes(x = Mean.Tmp.Dry.Q, y = Ann.Prc, color = cluster)) +
  geom_point() +
  theme(legend.position = "bottom")
ggExtra::ggMarginal(b, type = "histogram", bins = 20, color = "grey", fill = "blue")


ggplot(envdat_df, aes(x = cluster, y = Ann.Prc)) + 
  geom_boxplot(aes(fill = cluster))
ggplot(envdat_df, aes(x = cluster, y = Ann.Mean.Tmp)) + 
  geom_boxplot(aes(fill = cluster))
ggplot(envdat_df, aes(x = cluster, y = Max.Tmp.Wrm.M)) + 
  geom_boxplot(aes(fill = cluster))
ggplot(envdat_df, aes(x = cluster, y = Tmp.Seas)) + 
  geom_boxplot(aes(fill = cluster))
ggplot(envdat_df, aes(x = cluster, y = Prc.Seas)) + 
  geom_boxplot(aes(fill = cluster))


# Parallel coordiante plots allow us to put each feature on seperate column and lines connecting each column
ggparcoord(data = envdat_df, columns = 1:17, groupColumn = 18, alphaLines = 0.4, 
           title = "Parallel Coordinate Plot for Climate Data", 
           scale = "globalminmax", showPoints = TRUE) + 
  theme(legend.position = "bottom")
