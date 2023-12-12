#### Future Climate Score Calculation ####

library(terra)
library(sf)
library(ggplot2)
setwd("G:/My Drive/Sorghum/R/")

# Load raster files
prc.max <- rast("./geospatial/prc_max.tif")
prc.min <- rast("./geospatial/prc_min.tif")
tmean.max <- rast("./geospatial/tmean_max.tif")
tmean.min <- rast("./geospatial/tmean_min.tif")
bioclim.prc.max <- rast("./geospatial/bioclim_prc_max.tif")
bioclim.prc.min <- rast("./geospatial/bioclim_prc_min.tif")
bioclim.temp.min <- rast("./geospatial/bioclim_temp_min.tif")
bioclim.temp.max <- rast("./geospatial/bioclim_temp_max.tif")

# Load environment-of-origin data
envdat<- read.csv("./input_data/envdat_master.csv")
# envdat <- read.csv("./input_data/envdat_core.csv") # for creating core collection score
regions <- envdat[,c(1,3,4)]

#################################################################################################################################
# Load and prep sorghum cover data
# Import cropland raster
# cropland.og <- rast("Global_cropland_3km_2019.tif")
sorghum.og <- rast("./geospatial/sorghum_HarvestedAreaFraction.tif")

# Resample to align rasters to future climate data

prc.da <- disagg(prc.max, fact = 7)

sorghum.rs <- resample(sorghum.og, prc.da)

sorghum <- aggregate(sorghum.rs, fact = 7)

#### Get polygon of all countries in study
country_region <- as.data.frame(cbind(envdat$country, envdat$region))
colnames(country_region) <- c("country", "region")
country_region <- country_region[!duplicated(country_region$country),]
countries<- st_read(dsn = "./geospatial/world-administrative-boundaries.geojson", crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

#### Further prep
# crop cropland raster
# rasterize country polygon
v <- vect(countries)
r <- rast(v, ncols=2160, nrows=1080)
terra::rasterize(v, sorghum, "name", filename = "borders.tif", overwrite = T)
writeRaster(sorghum, "./geospatial/sorghum.tif", overwrite = TRUE)

#################################################################################################################################
## Variable Subsetting and Data Prep
# Temperature Variables
temp.q <- subset(envdat, select = c(pi, tmean_Q1:tmean_Q4, country, region))
temp.q <- temp.q %>% # get data into same form as TIF files (currently is degrees C * 10)
  mutate_at(c("tmean_Q1", "tmean_Q2", "tmean_Q3", "tmean_Q4"), function(x) {x/10}) 
temp.bioclim <- subset(envdat, 
                       select = c(pi:region, Ann.Mean.Tmp, Tmp.Seas, Mean.Diurn.Rng, Max.Tmp.Wrm.M, Min.Tmp.Cld.M, 
                                  Mean.Tmp.Wet.Q, Mean.Tmp.Dry.Q, Mean.Tmp.Wrm.Q,Mean.Tmp.Cld.Q))
colnames(temp.bioclim)[5:ncol(temp.bioclim)] <- c("bio01", "bio04", "bio02", "bio05", "bio06", "bio08", "bio09", "bio10", "bio11")
temp.bioclim <- temp.bioclim[,c(1, 5:13, 3, 4)]
temp.bioclim <- temp.bioclim %>% # get data into same form as TIF files (currently is degrees C * 10)
  mutate_at(colnames(temp.bioclim)[c(2:10)], function(x) {x/10}) 

# Precipitation Variables
prc.q <- subset(envdat, select = c(pi, prec_Q1:prec_Q4, country, region))
prc.bioclim <- subset(envdat, select = c(pi:region, Ann.Prc, Prc.Seas, Prc.Wet.M, Prc.Dry.M, Prc.Wet.Q, 
                                         Prc.Dry.Q, Prc.Cld.Q, Prc.Wrm.Q))
colnames(prc.bioclim)[5:ncol(prc.bioclim)] <- c("bio12", "bio15", "bio13", "bio14", "bio16", "bio17", "bio19", "bio18")
prc.bioclim <- prc.bioclim[,c(1, 5:12, 3, 4)]

#################################################################################################################################
# Write function to calculate score

future_climate_score <- function(df, rast_max, rast_min){ # enter uncombined dataframe i.e. prc.q or prc.extra
  score.df <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(df)))
  colnames(score.df) <- colnames(df)
  borders <- rast("borders.tif") # raster with country names as cell value
  crop <- rast("sorghum.tif") # raster with percentage of cropland in cell as cell value
  crop[is.na(crop)] <- 0
  for(i in 1: nrow(df)){
    for(j in 2:(ncol(df)-2)){
      country_sub <- borders %in% df[i,which(colnames(df)=="country")]
      country_sub[is.na(country_sub)] <- 0
      cxc <- crop * country_sub # mask everything outside of country by setting to zero
      rast2sum <- rast_max[[colnames(df)[j]]] >= df[i,j] & rast_min[[colnames(df)[j]]] <= df[i,j]
      rast2sum[is.na(rast2sum)] <- 0
      crop_corr <- cxc * rast2sum # correct for percentage cropland in cell
      sum <- global(crop_corr, fun = "sum")
      denom <- global(cxc, fun = "sum")
      score <- sum[1,1]/denom[1,1]
      score.df[i,j] <- score
    } 
  }
  score.df$pi <- df$pi
  score.df$region <- df$region
  score.df$country <- df$country
  return(score.df)
}

#################################################################################################################################
# Run function and write files

prc.out <- future_climate_score(df = prc.q, rast_max = prc.max, rast_min = prc.min)
bioclim.prc.out <- future_climate_score(prc.bioclim, bioclim.prc.max, bioclim.prc.min)
temp.out <- future_climate_score(temp.q, tmean.max, tmean.min)
bioclim.temp.out <- future_climate_score(temp.bioclim, bioclim.temp.max, bioclim.temp.min)

#write.csv(prc.out, "future_climate_prc_out.csv", row.names = FALSE)
#write.csv(temp.out, "future_climate_temp_out.csv", row.names = FALSE)
#write.csv(bioclim.prc.out, "future_climate_prc_bioclim_out.csv", row.names = FALSE)
#write.csv(bioclim.temp.out, "future_climate_temp_bioclim_out.csv", row.names = FALSE)


#################################################################################################################################
# Calculate Country Score

# Precipitation
prc.out <- merge(bioclim.prc.out, prc.out, by = c("pi", "region", "country"))
prc.out$score.prc <- apply(prc.out[,4:(ncol(prc.out))], 1, FUN = "mean")

#write.csv(prc.out, "score_futureclimates_precip_byaccession.csv")

prc.score <- aggregate(prc.out$score.prc, by=list(prc.out$country), FUN=sum)
colnames(prc.score) <- c("country", "score.sum")
prc.score$total <- NA
for(i in 1:length(prc.score$total)){
  count <- length(which(envdat$country==prc.score$country[i]))
  prc.score$total[i] <- count
}
prc.score$prc.score <- prc.score$score.sum / prc.score$total

# Temperature
temp.out <- merge(bioclim.temp.out, temp.out, by = c("pi", "region", "country"))
temp.out$score.temp <- apply(temp.out[,5:(ncol(temp.out))], 1, FUN = "mean")

#write.csv(temp.out, "score_futureclimates_temp_byaccession.csv")

temp.score <- aggregate(temp.out$score.temp, by=list(temp.out$country), FUN=sum)
colnames(temp.score) <- c("country", "score.sum")
temp.score$total <- NA
for(i in 1:length(temp.score$total)){
  count <- length(which(envdat$country==temp.score$country[i]))
  temp.score$total[i] <- count
}
temp.score$temp.score <- temp.score$score.sum / temp.score$total

# Aggregate
country_region <- data.frame(country = envdat$country, region = envdat$region)
country_region <- country_region[!duplicated(country_region$country),]

prc <- data.frame(country = prc.score$country, prc.score = prc.score$prc.score, num_accessions = prc.score$total)
temp <- data.frame(country = temp.score$country, temp.score = temp.score$temp.score)
fc.scores <- merge(prc, temp, by = "country")
fc.scores <- merge(fc.scores, country_region, by = "country")
fc.scores <- fc.scores[,c(1,5,3,2,4)]
fc.scores$fc.score <- apply(fc.scores[,4:5], 1, FUN="sum")

# Remove countries with less than 5 accessions
fc.scores <- fc.scores[which(fc.scores$num_accession > 5),]

# Combine South and East Asia
fc.scores$region[which(fc.scores$region == "East Asia" | fc.scores$region == "South Asia")] <- "Asia"

# Replace NAs with 0s
fc.scores <- fc.scores %>% replace(is.na(.), 0)

write.csv(fc.scores, "./score_outputs/score_futureclimates_country.csv", row.names = F)
# write.csv(fc.scores, "./score_outputs/score_futureclimates_core_country.csv", row.names = F) # write core collection ouput

#######################################################################################################################################################
# Calculate regional score

# Bind score categories into single dataframe
fc.regional <- as.data.frame(cbind(country = temp.out$country, region = temp.out$region, temp_all = temp.out$score.temp, prc_all = prc.out$score.prc))

# Write scores as numeric
fc.regional$temp_all <- as.numeric(fc.regional$temp_all)
fc.regional$prc_all <- as.numeric(fc.regional$prc_all)

# Combine South and East Asia
fc.regional$region[which(fc.regional$region == "East Asia" | fc.regional$region == "South Asia")] <- "Asia"

# Sum all scores by region
fc.regionalscores <- fc.regional %>% group_by(region) %>% 
  summarise(total = n(), # counts number of accessions in each region
            temp_sum = sum(temp_all),
            prc_sum = sum(prc_all),
            .groups = 'drop')

# Regional Score- Normalize by number of accessions
fc.regionalscores$temp_score <- fc.regionalscores$temp_sum/fc.regionalscores$total
fc.regionalscores$prc_score <- fc.regionalscores$prc_sum/fc.regionalscores$total

fc.regionalscores <- as.data.frame(fc.regionalscores)

fc.regionalscores$score <- fc.regionalscores$temp_score + fc.regionalscores$prc_score

write.csv(fc.regionalscores, "./score_outputs/score_futureclimates_region.csv", row.names = F)
# write.csv(fc.regionalscores, "./score_outputs/score_futureclimates_core_region.csv", row.names = F) # write core collection score
#################################################################################################################################
# Output accession score

fc.score.acc <- data.frame(pi = temp.out$pi, country = temp.out$country, region= temp.out$region, prc.score = prc.out$score.prc, temp.score = temp.out$score.temp)

fc.score.acc$score.acc <-  fc.score.acc$temp.score +  fc.score.acc$prc.score

write.csv(fc.score.acc, "./score_outputs/score_futureclimates_accession.csv", row.names = F)
# write.csv(fc.score.acc, "./score_outputs/score_futureclimates_core_accession.csv", row.names = F) # write core collection score