# Recreate Figure 2 (maps and bar charts of scores) #

library(viridis)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(dplyr)
library(scico)
library(paletteer)

#########################################################################################################################
### Bar Chart - Country/Region/Core Scores for inset graph ###

###### Genomic Selection Graph ######
fc.score <- read.csv("./score_outputs/score_futureclimates_sorghum_accession_updatedregions.csv")
gs.score <- read.csv('./score_outputs/score_genomicselection_accession_updatedregions.csv')

gs.score.country <- gs.score %>% group_by(country) %>% 
  summarise(region = region, count = n(), score = mean(score)) %>%
  filter(!duplicated(country), count > 10)

gs.score.region <- gs.score.country %>% group_by(region) %>% 
  summarise(count = sum(count), score = mean(score))

gs.acc.score.mini <- read.csv("./score_outputs/score_genomicselection_minicore_accession.csv")
mini_mean_score <- mean(gs.acc.score.mini$score)

gs.score.country$country[which(gs.score.country$country == "United Republic of Tanzania")] <- "Tanzania"

#gs.score.country$region <- factor(gs.score.country$region, levels = c('West Africa', 'Central Africa','South Africa', 'East Africa', 'West Asia', 'South Asia', 'East Asia'))
tiff('./revision_scripts/gs_barplot_revised.tif', width = 1000, height = 500)
ggplot(gs.score.country, aes(x=country, y=score, label = country)) +
  geom_col(position = "dodge", color = "black", fill = '#d1e5f0') + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 25), legend.position = "none", axis.text.x = element_blank(), 
        strip.text.x = element_text(size = 15)) +
  # scale_y_continuous(breaks = seq(0, 0.3, by = 0.05)) +
  labs(y = "", x = "") +
  geom_hline(yintercept = mean(gs.score.country$score), color="#5b8efc", size = 3, linetype = "longdash") + # add global mean
  geom_hline(yintercept = mini_mean_score, color="#de217c", size = 3, linetype = "longdash") + # add minicore mean score
  geom_hline(data = gs.score.region, aes(yintercept = score), color = "#ffae0d", size = 3, linetype = "longdash") +
  geom_text(data = subset(gs.score.country, score > 0.041), position = position_dodge2(width = 0.9, preserve = "single"), 
            angle = 90, hjust = "right", size = 9.5) + 
  geom_text(data = subset(gs.score.country, score < 0.041), position = position_dodge2(width = 0.9, preserve = "single"), 
            angle = 90, hjust = "left", size = 9.5) +
  facet_grid(~factor(region, levels = c('West Africa', 'Central Africa','South Africa', 'East Africa', 'West Asia', 'South Asia', 'East Asia')), 
             scale="free_x", space="free_x", labeller = (region = label_wrap_gen(8)))
dev.off()

######## Future Climates ###########

fc.score.country <- fc.score %>% group_by(country) %>% 
  summarise(region = region, count = n(), score = mean(score)) %>%
  filter(!duplicated(country), count > 10)

fc.score.region <- fc.score.country %>% group_by(region) %>% 
  summarise(count = sum(count), score = mean(score))

fc.core <- read.csv("./score_outputs/score_futureclimates_sorghum_core_accession_updatedregions.csv")
fc.core <- na.omit(fc.core)
fc.core.region <- fc.core %>% group_by(region) %>% summarize(count = n(), score = mean(score))
fc.core.region <- fc.core.region[-1,]

fc.mini_score <- fc.score[(fc.score$pi %in% gs.acc.score.mini$pi),]
fc.mini_score.mean <- mean(fc.mini_score$score)

fc.score.country$country[which(fc.score.country$country == "United Republic of Tanzania")] <- "Tanzania"
tiff('./revision_scripts/fc_score_barplots_Revised.tif', width = 1000, height = 500)
ggplot(fc.score.country, aes(x=country, y=score, label = country)) +
  geom_col(position = "dodge", color = "black", fill = '#d1e5f0') + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 25), legend.position = "none", axis.text.x = element_blank(), 
        strip.text.x = element_text(size = 15)) +
  labs(y = "", x = "") +
  geom_hline(yintercept = mean(na.omit(fc.score$score)), color="#5b8efc", size = 3, linetype = "longdash") + # add global mean
  geom_hline(yintercept = fc.mini_score.mean, color="#de217c", size = 3, linetype = "longdash") +
  geom_hline(data = fc.score.region, aes(yintercept = score), color = "#ffae0d", size = 3, linetype = "longdash") +
  geom_hline(data = fc.core.region, aes(yintercept = score), color = "#735df0", size = 3, linetype = "longdash") +
  geom_text(data = fc.score.country, angle = 90, hjust = "right", size = 9.5) + 
  facet_grid(~factor(region, levels = c('West Africa', 'Central Africa','South Africa', 'East Africa', 'West Asia', 'South Asia', 'East Asia')),
             scale="free_x", space="free_x", labeller = (region = label_wrap_gen(8)))
dev.off()

#########################################################################################################################
### Point maps - accession scores plotted on map of study region ###

# Load Data
geodat <- read.csv('./input_data/envdat_master_updatedregions.csv')

geodat <- subset(geodat, select = c("pi", "Latitude", "Longitude"))

acc_score <- read.csv("./score_outputs/score_genomicselection_accession_updatedregions.csv")

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
  geom_sf(data = worldmap_trans, fill = "#e9e9e9") +
  geom_sf(aes(color = score), score_moll, size = 3, alpha = 1) +    
  coord_sf(xlim = disp_win_coord[,'X'], ylim = disp_win_coord[,'Y'],
           datum = target_crs, expand = FALSE) +
  #scale_colour_gradientn(colours = c("#001242", "#5010ff", "#087d83", "#c7e022", "#fce726"), name = "Score", na.value = "grey") +
  scico::scale_colour_scico(palette = "lajolla") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

##############################
# Future Climate Score

fc.dat <- read.csv("./score_outputs/score_futureclimates_sorghum_accession_updatedregions.csv")
geo_fc_score <- merge(fc.dat, geodat, by = "pi") # load and prep data

score_fc <- st_as_sf(geo_fc_score, coords = c("Longitude", "Latitude"), crs = 4326)
score_moll_fc <- st_transform(score_fc, crs = target_crs) # spatial transformation of point coordinates

ggplot() +
  geom_sf(data = worldmap_trans, fill = "#e9e9e9") +
  geom_sf(aes(color = score), score_moll_fc, size = 3) +    
  coord_sf(xlim = disp_win_coord[,'X'], ylim = disp_win_coord[,'Y'],
           datum = target_crs, expand = FALSE) +
  #scale_colour_gradientn(colours = c("#001242", "#5010ff", "#087d83", "#80d44c", "#fce726"), name = "Score", na.value = "grey") +
  scale_colour_paletteer_c('viridis::inferno')+
  labs(x = "Longitude", y = "Latitude", title = "Future Climates") +
  theme_bw()
