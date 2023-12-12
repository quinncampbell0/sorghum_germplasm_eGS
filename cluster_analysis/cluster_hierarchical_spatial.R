

library(sf)
library(terra)
library(geosphere)
library(rnaturalearth)
library(dplyr)

# read data, convert to spatial point data frame
envdat <- read.csv("./input_data/envdat_master.csv")
envdat_sp <- st_as_sf(envdat, coords = c("Longitude", "Latitude"), crs = 4326)

# use the distm function to generate a geodesic distance matrix in meters
mdist <- st_distance(envdat_sp)

# cluster all points using a hierarchical clustering approach
hc <- hclust(as.dist(mdist), method="complete")

# define the distance threshold, in this case 40 m
d=1100000

# define clusters based on a tree "height" cutoff "d" and add them to the SpDataFrame
envdat_sp$clust_buffer <- cutree(hc, h=d)
envdat_sp$clust_buffer <- as.factor(envdat_sp$clust_buffer)
envdat_sp$clust_country <- cutree(hc, k=51)

spat_clust <- data.frame(pi = envdat_sp$pi, clust_country = envdat_sp$clust_country, clust_buffer = envdat_sp$clust_buffer)

write.csv(spat_clust, "./cluster_analysis/spatial_clusters.csv", row.names=F)

# count number of accessions per cluster
count_country <- envdat_sp %>% group_by(clust_country) %>% count()

worldmap <- ne_countries(scale = 'medium', type = 'map_units', returnclass = 'sf')
head(worldmap[c('name', 'continent')])

disp_win_wgs84 <- st_sfc(st_point(c(-27, -36)), st_point(c(180, 48)),
                         crs = 4326)
disp_win_trans <- st_transform(disp_win_wgs84, crs = target_crs)
disp_win_coord <- st_coordinates(disp_win_trans)

target_crs <- '+proj=moll'
worldmap_trans <- st_transform(worldmap, crs = target_crs)
envdat_moll <- st_transform(envdat_sp, crs = target_crs)

#pdf("plot_accession_futureclimate.pdf", width = 14, height=11)
set.seed(935234)
P48 <- createPalette(48, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))

names(P48) <- NULL
# Cluster by distance (cut tree at specific height)
ggplot() +
  geom_sf(data = worldmap_trans, fill = "white") +
  geom_sf(aes(color = clust_buffer), envdat_moll, size = 2) +    
  coord_sf(xlim = disp_win_coord[,'X'], ylim = disp_win_coord[,'Y'],
           datum = target_crs, expand = FALSE) +
  scale_colour_manual(values = P48) +
  #geom_point(aes(x = Longitude, y = Latitude, colour = Genomic.Selection.Score), size = 2) +
  #scale_colour_gradientn(colours = c("#001242", "#5010ff", "#087d83", "#deee12", "#f0fa01"), name = "Score", na.value = "grey") +
  #scale_color_viridis() +
  labs(x = "Longitude", y = "Latitude", title = "Future Climates") +
  theme_bw()


# Cluster by k=51 (number of countries in study)

