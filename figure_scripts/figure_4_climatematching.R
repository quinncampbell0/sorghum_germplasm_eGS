library(ggplot2)
library(maps)
library(network)
library(sna)
library(RColorBrewer)
library(dplyr)
library(GGally)
library(ggnetwork)
library(ggplot2)
library(dplyr)
library(rgeos)
library(sf)
library(ggrepel)
library(ggpubr)

#setwd('G:/My Drive/Sorghum/R/)

geav <- read.csv('./input_data/rrblup_GEBV_data.csv')
envdat <- read.csv('./input_data/envdat_master.csv')
scores <- read.csv('./score_outputs/score_futureclimates_sorghum_accession_updatedregions.csv')

metad <- envdat %>% subset(select = c('gen_id','country', 'Latitude', 'Longitude'))

top10_AnnMeanTmp <- geav %>%
  arrange(desc(Ann.Mean.Tmp)) %>%
  slice_head(n = 10) %>% 
  mutate(rank = 1:10) %>% 
  subset(select = c("X", "rank"))

top10_AnnMeanTmp <- merge(top10_AnnMeanTmp, metad, by.y = 'gen_id', by.x = 'X')

top10_Mean.Tmp.Dry.Q <- geav %>%
  arrange(desc(Mean.Tmp.Dry.Q)) %>%
  slice_head(n = 10) %>% 
  mutate(rank = 1:10) %>% 
  subset(select = c("X", "rank"))

top10_Mean.Tmp.Dry.Q <- merge(top10_Mean.Tmp.Dry.Q, metad, by.y = 'gen_id', by.x = 'X')


top10_Max.Tmp.Wrm.M <- geav %>%
  arrange(desc(Max.Tmp.Wrm.M)) %>%
  slice_head(n = 30) %>% 
  mutate(rank = 1:30) %>% 
  subset(select = c("X", "rank"))

top10_Max.Tmp.Wrm.M <- merge(top10_Max.Tmp.Wrm.M, metad, by.y = 'gen_id', by.x = 'X')


library(dplyr)

#write.csv(top10_AnnMeanTmp, 'G:/My Drive/Sorghum/R/map_climatematch/top10_AnnMeanTmp.csv')
#top10_AnnMeanTmp<-read.csv('G:/My Drive/Sorghum/R/map_climatematch/top10_AnnMeanTmp.csv')[,-1]

#coords <- read.csv('G:/My Drive/Sorghum/R/map_climatematch/coordinates_df.csv', row.names = 1)

set.seed(123)  # For reproducibility

plot_coords <- read.csv('./input_data/coordinates_for_plotting.csv')

network <- data.frame(orig = plot_coords$X, dest = "Kenya")

# convert to network
flights <- network(network, directed = F, loops = F, multiple = F)

# add geographic coordinates
#flights %v% "lat" <- coords[network.vertex.names(flights), "lat"]
#flights %v% "lon" <- coords[network.vertex.names(flights), "lon"]
flights %v% "lat" <- plot_coords$Latitude
flights %v% "lon" <- plot_coords$Longitude
flights %v% 'pi' <- plot_coords$pi

# compute degree centrality
flights %v% "degree" <- degree(flights, gmode = "digraph")
flights %v% 'rank' <- plot_coords$rank
# add random groups
#flights %v% "mygroup" <- sample(letters[1:4], network.size(flights), replace = TRUE)

# Convert world map data to sf object for centroid calculation
world_sf <- st_as_sf(map_data("world"), coords = c("long", "lat"), crs = 4326) %>%
  group_by(region) %>%
  summarize(geometry = st_union(geometry)) %>%
  ungroup()

# Compute centroids of countries in data frame
country_labels1 <- world_sf %>%
  filter(region %in% plot_coords$country) %>%
  mutate(centroid = st_centroid(geometry)) %>%
  st_as_sf() %>%
  mutate(lon = st_coordinates(centroid)[,1],
         lat = st_coordinates(centroid)[,2])

# Manually shift country labels
country_labels1$lat[6] <- 5.57
country_labels1$lon[6] <- 46.5

country_labels1$lat[4] <- 17
country_labels1$lon[4] <- -3

country_labels1$lat[2] <- 6

country_labels1$lat[1] <- 11
country_labels1$lon[1] <- -6


#pdf(width = 15, height = 12, res = 300)
# Plot the map with country labels
g1 <- ggplot() +
  # Base map
  geom_polygon(data = world_sf, aes(x = long, y = lat, group = group, 
                                 fill = ifelse(region == "Kenya", "Kenya", "Other")),  
               color = "white") +
  scale_fill_manual(values = c("Kenya" = "lightblue", "Other" = "gray90")) +  
  coord_quickmap(xlim = c(-15, 50), ylim = c(-10, 25)) +
  
  # Draw edges with jittered starting points
  geom_curve(data = flights, 
             aes(x = lon, y = lat, 
                 xend = 36.8219, yend = -1.2921),  #as.factor(rank)
             curvature = -0.2, ncp = 100, linewidth = 1,
             color = 'black') +
  
  geom_point(data = flights, aes(x = lon, y = lat), color = 'red') +
  #geom_text(data = flights, aes(x = lon, y = lat, label = pi), hjust = 1.1, vjust = 1, size = 3, color = 'darkred') + #, position = position_dodge2(width = 0.5)
  geom_label_repel(data = plot_coords, aes(x = Longitude, y = Latitude, label = pi), size = 3, color = 'darkred', box.padding = 0.2, seed = 123) +
  # Add country labels at polygon centroids
  geom_text(data = country_labels1, aes(x = lon, y = lat, label = region), 
            color = "black", size = 4, fontface = "bold", hjust = 0.5) +
  
  #scale_color_manual(values = green_palette) +  
  labs(fill = "Country", color = "Rank") +
  theme_void() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) 


#########################################################################################################################
# Plot B - Score plot

geodat <- envdat %>% subset(select = c('pi','country', 'Latitude', 'Longitude'))
top10_score <- scores %>%
  arrange(desc(score)) %>%
  slice_head(n = 10) %>% 
  mutate(rank = 1:10) %>% 
  subset(select = c("pi", "rank"))

top10_score <- merge(top10_score, geodat, by = 'pi')


network_score <- data.frame(orig = top10_score$pi, dest = "Kenya")

# convert to network
score_edges <- network(network_score, directed = F, loops = F, multiple = F)

# add geographic coordinates
#flights %v% "lat" <- coords[network.vertex.names(flights), "lat"]
#flights %v% "lon" <- coords[network.vertex.names(flights), "lon"]
score_edges %v% "lat" <- top10_score$Latitude
score_edges %v% "lon" <- top10_score$Longitude
score_edges %v% 'pi' <- top10_score$pi

# compute degree centrality
score_edges %v% "degree" <- degree(score_edges, gmode = "digraph")
score_edges %v% 'rank' <- top10_score$rank
# add random groups
#score_edges %v% "mygroup" <- sample(letters[1:4], network.size(score_edges), replace = TRUE)

# Convert world map data to sf object for centroid calculation
world_sf <- st_as_sf(map_data("world"), coords = c("long", "lat"), crs = 4326) %>%
  group_by(region) %>%
  summarize(geometry = st_union(geometry)) %>%
  ungroup()

# Compute centroids of countries in data frame
country_labels <- world_sf %>%
  filter(region %in% top10_score$country) %>%
  mutate(centroid = st_centroid(geometry)) %>%
  st_as_sf() %>%
  mutate(lon = st_coordinates(centroid)[,1],
         lat = st_coordinates(centroid)[,2])

# Manually shift country labels
country_labels$lat[2] <- 8.2 # Central African Republic
country_labels$lon[2] <- 22.6 # Central African Republic

country_labels$lat[3] <- -29 # Lesotho
country_labels$lon[3] <- 30.2

country_labels$lat[1] <- -3.4 # Burundi
country_labels$lon[1] <- 27.5 # Burundi

country_labels$lat[4] <- 0.7# Uganda
country_labels$lon[4] <- 29.2 # Uganda

country_labels$lon[5] <- 27 # Zambia
country_labels$lat[5] <- -14.25 

#pdf(width = 15, height = 12, res = 300)
# Plot the map with country labels
g2 <- ggplot() +
  # Base map
  geom_polygon(data = world_sf, aes(x = long, y = lat, group = group, 
                                 fill = ifelse(region == "Kenya", "Kenya", "Other")),  
               color = "white") +
  scale_fill_manual(values = c("Kenya" = "lightblue", "Other" = "gray90")) +  
 coord_quickmap(xlim = c(0, 40), ylim = c(-30, 10)) +
  
  # Draw edges with jittered starting points
  geom_curve(data = score_edges, 
             aes(x = lon, y = lat, 
                 xend = 36.8219, yend = -1.2921),  #as.factor(rank)
             curvature = -0.2, ncp = 100, linewidth = 1,
             color = 'black') +
  
  geom_point(data = score_edges, aes(x = lon, y = lat), color = 'red') +
  #geom_text(data = score_edges, aes(x = lon, y = lat, label = pi), size = 3, color = 'darkred') + #, hjust = 1.1, vjust = 1, position = position_dodge2(width = 0.5), ,, position = position_jitter(width = 0.8)
  geom_label_repel(data = top10_score, aes(label = pi, x = Longitude, y = Latitude), size = 3, color = 'darkred', box.padding = 0.4, seed = 123) +
  # Add country labels at polygon centroids
  geom_text(data = country_labels, aes(x = lon, y = lat, label = region), 
            color = "black", size = 4, fontface = "bold", hjust = 0.5) +
  
  #scale_color_manual(values = green_palette) +  
  labs(fill = "Country", color = "Rank") +
  theme_void() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) 

pdf('G:/My Drive/Sorghum/R/map_climatematch/map_figure_climatematch.pdf', width  = 10, height = 12)
ggarrange(g1, g2, ncol = 1, labels = c('A', 'B'), heights = c(1,1.5))
dev.off()



