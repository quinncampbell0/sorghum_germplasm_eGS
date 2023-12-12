# Recreate Figure 2 (maps and bar charts of scores) #

library(viridis)
library(ggplot2)
library(sf)
library(rnaturalearth)

#########################################################################################################################
### Bar Chart - Country/Region/Core Scores for inset graph ###

###### Genomic Selection Graph ######
gs.scores <- read.csv("./score_outputs/score_genomicselection_country.csv")
gs.regionalscores <- read.csv("./score_outputs/genomic_selection_score_region.csv")

gs.scores <- gs.scores %>% filter(num_accessions > 10)

gs.scores$country[which(gs.scores$country == "United Republic of Tanzania")] <- "Tanzania"

gs.acc.score.mini <- read.csv("./score_outputs/score_genomicselection_minicore_accession.csv")

mini_mean_score <- mean(gs.acc.score.mini$score)

ggplot(gs.scores, aes(x=country, y=gs.score, label = country, fill = region)) +
  geom_col(position = "dodge", color = "black") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 25), legend.position = "none", axis.text.x = element_blank(), 
        strip.text.x = element_text(size = 20)) +
  # scale_y_continuous(breaks = seq(0, 0.3, by = 0.05)) +
  labs(y = "", x = "") +
  geom_hline(yintercept = mean(gs.scores$gs.score), color="red", size = 1.5, linetype = "longdash") + # add global mean
  geom_hline(yintercept = mini_mean_score, color="purple", size = 1.5, linetype = "longdash") + # add minicore mean score
  geom_hline(data = gs.regionalscores, aes(yintercept = gs.score), color = "darkblue", size = 1.5, linetype = "longdash") +
  geom_text(data = subset(gs.scores, gs.score > 0.041), position = position_dodge2(width = 0.9, preserve = "single"), 
            angle = 90, hjust = "right", size = 10) + 
  geom_text(data = subset(gs.scores, gs.score < 0.041), position = position_dodge2(width = 0.9, preserve = "single"), 
            angle = 90, hjust = "left", size = 10) +
  facet_grid(cols = vars(region), scale="free_x", space="free_x")


##### Future Climate Graph #####
fc.scores <- read.csv("./score_outputs/score_futureclimates_sorghum_country.csv")
fc.regionalscores <- read.csv("./score_outputs/score_futureclimates_sorghum_region.csv")
# fc_nozeros <- fc.scores[which(fc.scores$score.norm > 0),]
fc.scores$country[which(fc.scores$country == "United Republic of Tanzania")] <- "Tanzania"
fc.regionalscores.core <- na.omit(read.csv("./score_outputs/score_futureclimates_sorghum_core_region.csv"))
fc.acc.score <- read.csv('./score_outputs/score_futureclimates_sorghum_accession.csv')
fc.acc.score <- na.omit(fc.acc.score)

fc.acc.score$region[which(fc.acc.score$region == "East Asia" | fc.acc.score$region == "South Asia")] <- "Asia"

fc.score.country <- fc.acc.score %>% group_by(country) %>% 
  summarise(region = region, count = n(), score = mean(score.acc)) %>%
  filter(!duplicated(country), count > 10)

fc.score.region <- fc.acc.score %>% group_by(region) %>% 
  summarise(score = mean(score.acc))

fc.mini_score <- fc.acc.score[(fc.acc.score$pi %in% gs.acc.score.mini$pi),]
fc.mini_score.mean <- mean(fc.mini_score$score.acc)

ggplot(fc.score.country, aes(x=country, y=score, label = country, fill = region)) +
  geom_col(position = "dodge", color = "black") + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 25), legend.position = "none", axis.text.x = element_blank(), 
        strip.text.x = element_text(size = 20)) +
  # scale_y_continuous(breaks = seq(0, 0.3, by = 0.05)) +
  labs(y = "", x = "") +
  geom_hline(yintercept = mean(fc.acc.score$score.acc), color="red", size = 1.5, linetype = "longdash") + # add global mean
  geom_hline(yintercept = fc.mini_score.mean, color="orange", size = 1.5, linetype = "longdash") +
  geom_hline(data = fc.score.region, aes(yintercept = score), color = "darkblue", size = 1.5, linetype = "longdash") +
  geom_hline(data = fc.regionalscores.core, aes(yintercept = score), color = "green", size = 1.5, linetype = "longdash") +
  geom_text(data = fc.score.country, 
            angle = 90, hjust = "right", size = 10) + 
  #geom_text(data = subset(fc.score.country, score < 0.031), position = position_dodge2(width = 0.9, preserve = "single"), 
  #          angle = 90, hjust = "left", size = 7) +
  facet_grid(cols = vars(region), scale="free_x", space="free_x")

#########################################################################################################################
### Point maps - accession scores plotted on map of study region ###

# Load Data
geodat <- read.csv('./input_data/envdat_master.csv')

geodat <- subset(geodat, select = c("pi", "Latitude", "Longitude"))

acc_score <- read.csv("./score_outputs/score_genomicselection_accession.csv")

colnames(acc_score)[1] <- "pi"

geo_score <- merge(acc_score, geodat, by = "pi")

##### Spatial transformations of base map and zoom

worldmap <- ne_countries(scale = 'medium', type = 'map_units', returnclass = 'sf')
head(worldmap[c('name', 'continent')])

target_crs <- '+proj=moll'
worldmap_trans <- st_transform(worldmap, crs = target_crs)
disp_win_wgs84 <- st_sfc(st_point(c(-27, -36)), st_point(c(180, 48)),
                         crs = 4326)

disp_win_trans <- st_transform(disp_win_wgs84, crs = target_crs)

disp_win_coord <- st_coordinates(disp_win_trans)


# transformed maps - Mollweide proj
####################################

# Genomic Selection Score
score_sf <- st_as_sf(geo_score, coords = c("Longitude", "Latitude"), crs = 4326)
score_moll <- st_transform(score_sf, crs = target_crs) # spatial transformation of point coordinates

ggplot() +
  geom_sf(data = worldmap_trans, fill = "white") +
  geom_sf(aes(color = score), score_moll, size = 2) +    
  coord_sf(xlim = disp_win_coord[,'X'], ylim = disp_win_coord[,'Y'],
           datum = target_crs, expand = FALSE) +
  scale_colour_gradientn(colours = c("#001242", "#5010ff", "#087d83", "#deee12", "#f0fa01"), name = "Score", na.value = "grey") +
  labs(x = "Longitude", y = "Latitude", title = "Climate Adaptive Capacity") +
  theme_bw()

##############################
# Future Climate Score

fc.dat <- read.csv("./score_outputs/score_futureclimates_sorghum_accession.csv")
geo_fc_score <- merge(fc.dat, geodat, by = "pi") # load and prep data

score_fc <- st_as_sf(geo_fc_score, coords = c("Longitude", "Latitude"), crs = 4326)
score_moll_fc <- st_transform(score_fc, crs = target_crs) # spatial transformation of point coordinates

ggplot() +
  geom_sf(data = worldmap_trans, fill = "white") +
  geom_sf(aes(color = score.acc), score_moll_fc, size = 2) +    
  coord_sf(xlim = disp_win_coord[,'X'], ylim = disp_win_coord[,'Y'],
           datum = target_crs, expand = FALSE) +
  #scale_colour_gradientn(colours = c("#001242", "#5010ff", "#087d83", "#deee12", "#f0fa01"), name = "Score", na.value = "grey") +
  scale_color_viridis() +
  labs(x = "Longitude", y = "Latitude", title = "Future Climates") +
  theme_bw()
