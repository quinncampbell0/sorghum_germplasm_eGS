
##### Score partitioning (environmental, spatial, botanical, and taxonomic)
library(terra)
library(sf)
library(Polychrome)
library(dplyr)
library(viridis)
library(maps)
library(paletteer)
library(ggthemes)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggpubr)
#################################################################################################
######### Environmental clustering
clust_groups <- read.csv("./input_data/cluster_groups_hclust_k6.csv")
colnames(clust_groups) <- c('pi', "group")

envdat <- read.csv("./input_data/envdat_master.csv")
envdat_sub <- envdat[,c(1,5,6)]

envdat_cluster <- merge(envdat_sub, clust_groups, by = "pi")
envdat_cluster$group <- as.factor(envdat_cluster$group)

######################################
##### Create Map of climate clusters
######################################

cluster_sf <- st_as_sf(envdat_cluster, coords = c("Longitude", "Latitude"), crs = 4326)
target_crs <- '+proj=moll'
cluster_moll <- st_transform(cluster_sf, crs = target_crs)

worldmap <- ne_countries(scale = 'medium', type = 'map_units', returnclass = 'sf')
worldmap_trans <- st_transform(worldmap, crs = target_crs)

disp_win_wgs84 <- st_sfc(st_point(c(-25, -36)), st_point(c(140, 47)),
                         crs = 4326)
disp_win_trans <- st_transform(disp_win_wgs84, crs = target_crs)
disp_win_coord <- st_coordinates(disp_win_trans)

##### Create Map of climate clusters
clim_map <- ggplot() +
  geom_sf(data = worldmap_trans, fill = "white") +
  geom_sf(aes(color = group), cluster_moll, size = 2) +    
  coord_sf(xlim = disp_win_coord[,'X'], ylim = disp_win_coord[,'Y'],
           datum = target_crs, expand = FALSE) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Climatic Clustering Assignments") +
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) 

######################################
### Partition Score by climate cluster
######################################

### Adaptive Capacity score
ac.score <- read.csv('./score_outputs/score_genomicselection_accession.csv')
ac.env.score <- merge(ac.score, envdat_cluster[c(1,4)], by = "pi")

# Combine scores by cluster
ac.clust.score <- ac.env.score %>% group_by(group) %>% 
  summarise(clust.score = mean(score))
ac.clust.score$group <- as.factor(ac.clust.score$group)

head(ac.clust.score)

clim_ac <- ggplot(data = ac.clust.score, aes(x = group, y = clust.score, fill = group)) +
  geom_col() +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "score", x = "cluster", title = "Adaptive Capacity Score") +
  theme(legend.position = "none")

### Resilience Score
r.score <- read.csv('./score_outputs/score_futureclimates_sorghum_accession.csv')
r.env.score <- merge(r.score, envdat_cluster[c(1,4)], by = "pi")

# Combine scores by cluster
r.env.score <- na.omit(r.env.score)
r.clust.score <- r.env.score %>% group_by(group) %>% 
  summarise(clust.score = mean(score))
r.clust.score$group <- as.factor(r.clust.score$group)

head(r.clust.score)

clim_r <- ggplot(data = r.clust.score, aes(x = group, y = clust.score, fill = group)) +
  geom_col() +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "score", x = "cluster", title = "Climate Resilience Score") +
  theme(legend.position = "none")

##### Create combined plot #####

# read pca variable contribution plot
var_contrib_plot <- readRDS("./cluster_analysis/plot_env.pca_variable_contributions.rds")

# read pca for K=6
pca_envdat_k6 <- readRDS("./cluster_analysis/plot_env.pca_clusters_k6.rds")

ggarrange(
  ggarrange(
    clim_map,                # First row with line plot
    # Second row with box and dot plots
    ggarrange(clim_r, clim_ac, nrow = 2, labels = c("B", "C")), 
    ncol = 2, widths = c(3,1),
    labels = "A", legend = "bottom" # Label of the line plot
  ) ,
  ggarrange(pca_envdat_k6, var_contrib_plot, ncol = 2, labels = c("D", "E")),
  nrow = 2, heights = c(1.5,1)
)



########################################################################################################################
### Partition Score by botanical race
######################################

### Adaptive Capacity score
ac.score <- read.csv('./score_outputs/score_genomicselection_accession.csv')
race <- read.csv('./input_data/envdat_master_races.csv')
ac.race.score <- merge(ac.score, race[c(1,7)], by = "pi")

# Combine scores by cluster
ac.race.score <- na.omit(ac.race.score)

ac.racepart.score <- ac.race.score %>% group_by(race) %>% 
  summarise(clust.score = mean(score), count = n()) %>%
  filter(count >= 10)
ac.racepart.score$race <- as.factor(ac.racepart.score$race)

race_count <- ac.race.score %>% group_by(race) %>% 
  count()

head(ac.racepart.score)

### Resilience Score
r.score <- read.csv('./score_outputs/score_futureclimates_sorghum_accession.csv')

r.race.score <- merge(r.score, race[c(1,7)], by = "pi")

# Combine scores by cluster
r.race.score <- na.omit(r.race.score)

r.racepart.score <- r.race.score %>% group_by(race) %>% 
  summarise(clust.score = mean(score), count = n()) %>%
  filter(count >= 10)

r.racepart.score$race <- as.factor(r.racepart.score$race)

race_count <- r.race.score %>% group_by(race) %>% 
  count()

head(ac.racepart.score)

# Color palette
set.seed(12123)
P26 <- createPalette(26, c("#800000", "#228B22", "#0000FF"), range = c(30, 90))

names(P26) <- unique(race$race)

# Plot bar chart - Adaptive Capacity 
ac_race_bar <- ggplot(data = ac.racepart.score, aes(x = race, y = clust.score, fill = race)) +
  geom_col() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust=1), legend.position = "none") +
  scale_fill_manual(values = P26) +
  labs(y = "score")

# Plot bar chart - Resilience
r_race_bar <- ggplot(data = r.racepart.score, aes(x = race, y = clust.score, fill = race)) +
  geom_col() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust=1), legend.position = "none") +
  scale_fill_manual(values = P26) +
  labs(y = "score")

### Plot Map of races ###
race_sf <- st_as_sf(race, coords = c("Longitude", "Latitude"), crs = 4326)
target_crs <- '+proj=moll'
race_moll <- st_transform(race_sf, crs = target_crs)

race_map <- ggplot() +
  geom_sf(data = worldmap_trans, fill = "white") +
  geom_sf(aes(color = race), race_moll, size = 2) +    
  coord_sf(xlim = disp_win_coord[,'X'], ylim = disp_win_coord[,'Y'],
           datum = target_crs, expand = FALSE) +
  labs(title = "Score partitioned by botanical races") +
  theme_bw()+ 
  scale_color_manual(values = P26) + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) 

# Combine all plots
ggarrange(
  race_map,                # First row with line plot
  # Second row with box and dot plots
  ggarrange(r_race_bar, ac_race_bar, nrow = 2, labels = c("B", "C")), 
  ncol = 2, widths = c(3,1),
  labels = "A", legend = "bottom", common.legend = T # Label of the line plot
) 

########################################################################################################################
### Partition Score by spatial cluster
######################################

# Adaptive Capacity score
spat_clust <- read.csv('./cluster_analysis/spatial_clusters.csv')
ac.score <- read.csv('./score_outputs/score_genomicselection_accession.csv')

ac.spat.score <- merge(ac.score, spat_clust[c(1,2)], by = "pi")

envdat <- read.csv('./input_data/envdat_master.csv')
envdat_ll <- subset(envdat, select = c(pi, Latitude, Longitude))

ac.spat.sp <- merge(ac.spat.score, envdat_ll, by = "pi")

hulls <- ac.spat.sp %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  group_by(clust_country) %>%
  summarize(geometry = st_union(geometry), count = n(), score = mean(score)) %>%
  filter(count > 10) %>%
  st_convex_hull()
plot(hulls)

# Transform

target_crs <- '+proj=moll'
spat_moll <- st_transform(hulls, crs = target_crs)

spat_moll$clust_country <- as.factor(spat_moll$clust_country)
  
#### Plot score - spatial clusters ####
spatclust_map <- ggplot() +
  geom_sf(data = worldmap_trans, fill = "white") +
  geom_sf(aes(fill = score), spat_moll, size = 2) +
  coord_sf(xlim = disp_win_coord[,'X'], ylim = disp_win_coord[,'Y'],
           datum = target_crs, expand = FALSE) +
  labs(title = "Score Partitioned by Spatial Clusters") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_fill_gradientn(colours = c("#001242", "#5010ff", "#087d83", "#deee12", "yellow"), name = "Score", limits = c(0,0.45))
  
# Country chloropleth map for comparison
ac.score.country <- read.csv('./score_outputs/score_genomicselection_country.csv')

ac.score.country[ac.score.country$country %in% worldmap_trans$admin,]
merge <- merge(worldmap_trans, ac.score.country, by.x = "admin", by.y = "country", all.x = T)

country_map <- ggplot() +
  geom_sf(data = worldmap_trans, fill = "white") +
  geom_sf(data = merge, aes(fill = gs.score)) +
  coord_sf(xlim = disp_win_coord[,'X'], ylim = disp_win_coord[,'Y'],
           datum = target_crs, expand = FALSE) +
  labs(title = "Score Partitioned by Country") +
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  scale_fill_gradientn(colours = c("#001242", "#5010ff", "#087d83", "#deee12", "yellow"), name = "Score", 
                       na.value = "white", limits = c(0,0.45))

#


##################################
### Combined figure
##################################

ggarrange(country_map, spatclust_map, ncol = 2, labels = c("A", "B"), legend = "bottom", common.legend = T)

# Combine scores by cluster
ac.spat.score <- na.omit(ac.spat.score)

ac.spat.scorepart <- ac.spat.score %>% group_by(clust_country) %>% 
  summarise(clust.score = mean(score), count = n()) %>% filter(count > 10)

ac.spat.scorepart$clust_country <- as.factor(ac.spat.scorepart$clust_country)

ggplot(ac.spat.scorepart, aes(x = clust_country, y = clust.score)) +
  geom_col()

########################################################################################################################
### Genomic Clustering- Score partitioning

gen_clust <- read.csv("./cluster_analysis/cluster_groups_genomic.csv")
envdat <- read.csv("./input_data/envdat_master.csv")
ac.score <- read.csv('./score_outputs/score_genomicselection_accession.csv')
r.score <- read.csv('./score_outputs/score_futureclimates_sorghum_accession.csv')

gen_clust_geo <- merge(gen_clust[,2:3], envdat[,c(1,5,6)], by = "pi")

gen_clust_ac.score <- merge(gen_clust_geo, ac.score, by = "pi")
gen_clust_r.score <- merge(gen_clust_geo, r.score, by = "pi")


####################################
### Create map of genetic clusters
####################################

genclust_sf <- st_as_sf(gen_clust_geo, coords = c("Longitude", "Latitude"), crs = 4326)
target_crs <- '+proj=moll'
genclust_moll <- st_transform(genclust_sf, crs = target_crs)

worldmap <- ne_countries(scale = 'medium', type = 'map_units', returnclass = 'sf')
worldmap_trans <- st_transform(worldmap, crs = target_crs)

disp_win_wgs84 <- st_sfc(st_point(c(-25, -36)), st_point(c(140, 47)),
                         crs = 4326)
disp_win_trans <- st_transform(disp_win_wgs84, crs = target_crs)
disp_win_coord <- st_coordinates(disp_win_trans)

genclust_moll$cluster <- as.factor(genclust_moll$clust)

gen_map <- ggplot() +
  geom_sf(data = worldmap_trans, fill = "white") +
  geom_sf(aes(color = cluster), genclust_moll, size = 2) +    
  coord_sf(xlim = disp_win_coord[,'X'], ylim = disp_win_coord[,'Y'],
           datum = target_crs, expand = FALSE) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Genomic Clustering Assignments") +
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) 

################################
### Combine scores by cluster
################################

ac.score_bygenclust <- gen_clust_ac.score %>% group_by(clust) %>% 
  summarise(clust.score = mean(score))
ac.score_bygenclust$clust <- as.factor(ac.score_bygenclust$clust)

ac.gen <- ggplot(data = ac.score_bygenclust, aes(x = clust, y = clust.score, fill = clust)) +
  geom_col() +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "score", x = "cluster", title = "Adaptive Capacity Score") +
  theme(legend.position = "none")

### Resilience Score

# Combine scores by cluster
gen_clust_r.score <- na.omit(gen_clust_r.score)
r.score_bygenclust <- gen_clust_r.score %>% group_by(clust) %>% 
  summarise(clust.score = mean(score))
r.score_bygenclust$clust <- as.factor(r.score_bygenclust$clust)

head(r.score_bygenclust)

r.gen <- ggplot(data = r.score_bygenclust, aes(x = clust, y = clust.score, fill = clust)) +
  geom_col() +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "score", x = "cluster", title = "Climate Resilience Score") +
  theme(legend.position = "none")

### Combine genomic clustering figures into one

ggarrange(
  gen_map,                # First row with line plot
  # Second row with box and dot plots
  ggarrange(r.gen, ac.gen, nrow = 2, labels = c("B", "C")), 
  ncol = 2, widths = c(3,1),
  labels = "A", legend = "bottom" # Label of the line plot
) 

