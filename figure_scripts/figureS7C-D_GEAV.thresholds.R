
library(sf)
library(dplyr)
library(terra)
library(ggplot2)
library(reshape2)

### Threshold Analysis

## Bioclim
#get list of rasters (one for each model)
filenames <- list.files("D:/future_climates/bioclimatic/") # Downloaded bioclimatic data for SSP585, 2041-2060 from WorldClim

bioclim.rast.list <- lapply(filenames, function(x){
  file <- paste0("D:/future_climates/bioclimatic/", x)
  raster <- rast(file)
})

bioclim.sds <- sds(bioclim.rast.list)

bioclim.sds.mean <- app(bioclim.sds, mean)

# Define the function to find agreement means for each cell across all CMIP6 models
process_raster_agreement <- function(raster, tolerance = 0.3, agreement_threshold = 0.75) {
  # Function to compute agreement and mean
  process_cell <- function(values) {
    values <- values[!is.na(values)]  # Remove NA values
    if (length(values) == 0) return(NA)
    
    # Calculate the cell-wise mean or median (e.g., use median for robustness)
    cell_stat <- mean(values)  # You can replace with `median(values)`
    
    # Calculate the tolerance range
    lower_bound <- cell_stat * (1 - tolerance)
    upper_bound <- cell_stat * (1 + tolerance)
    
    # Count how many values fall within the tolerance range
    in_range <- values >= lower_bound & values <= upper_bound
    agreement_ratio <- sum(in_range) / length(values)
    
    # If agreement exceeds the threshold, return the mean of agreeing values
    if (agreement_ratio >= agreement_threshold) {
      return(mean(values[in_range]))
    } else {
      return(NA)  # Insufficient agreement
    }
  }
  
  # Apply the function to the raster using app()
  result <- app(raster, process_cell)
  
  return(result)
}

# Apply the agreement function with a 30% tolerance and 75% agreement threshold
agreement_bioclim <- process_raster_agreement(bioclim.sds, tolerance = 0.3, agreement_threshold = 0.75)


### Precipitation
filenames <- list.files("D:/future_climates/precipitation/") # Downloaded precipitation data for SSP585, 2041-2060 from WorldClim

prc.rast.list <- lapply(filenames, function(x){
  file <- paste0("D:/future_climates/precipitation/", x)
  raster <- rast(file)
})

prc.sds <- sds(prc.rast.list)

agreement_prc <- process_raster_agreement(prc.sds, tolerance = 0.3, agreement_threshold = 0.75)

## Temperature Max

filenames <- list.files("D:/future_climates/temperature_max/") # Downloaded max temperature data for SSP585, 2041-2060 from WorldClim

tmp_max.rast.list <- lapply(filenames, function(x){
  file <- paste0("D:/future_climates/temperature_max/", x)
  raster <- rast(file)
})

tmp.max.sds <- sds(tmp_max.rast.list)

agreement_tmp.max <- process_raster_agreement(tmp.max.sds, tolerance = 0.3, agreement_threshold = 0.75)

# Temperature Min
filenames <- list.files("D:/future_climates/temperature_min/") # Downloaded min temperature data for SSP585, 2041-2060 from WorldClim

tmp_min.rast.list <- lapply(filenames, function(x){
  file <- paste0("D:/future_climates/./temperature_min/", x)
  raster <- rast(file)
})

tmp.min.sds <- sds(tmp_min.rast.list)

agreement_tmp.min <- process_raster_agreement(tmp.min.sds, tolerance = 0.3, agreement_threshold = 0.75)

## Extract climate variable for envdat

envdat <- read.csv('./input_data/future_climate/envdat_master_revision.csv')

envdat_sp <- st_as_sf(envdat, coords = c("Longitude", "Latitude"), crs = 4326)

bioclim1 <- extract(agreement_bioclim, envdat_sp)

prc1 <- extract(agreement_prc, envdat_sp)

tmp.min1 <- extract(agreement_tmp.min, envdat_sp)

tmp.max1 <- extract(agreement_tmp.max, envdat_sp)

all_future <- cbind(bioclim1, prc1, tmp.min1, tmp.max1)

future_clim <- all_future[,-which(colnames(all_future)=='ID')]

future_clim$pi <- envdat$pi

# convert temperature values (multiply by 10)

future_clim1 <- future_clim %>%
  mutate(across(c(bio01, bio02, bio04, bio05:bio11), ~ . * 10))

# Rename columns
mapping <- read.csv('./input_data/var_rename.csv')

# Convert the mapping to a named vector
col_mapping <- setNames(mapping$new.name, mapping$original.name)


# Rename columns
names(future_clim1) <- ifelse(names(future_clim1) %in% names(col_mapping), col_mapping[names(future_clim1)], names(future_clim1))


# Threshold calculation function
calculate_thresholds_perc <- function(var) {
  lower_threshold <- quantile(var, probs = 0.3, na.rm = TRUE) # XXth percentile
  upper_threshold <- quantile(var, probs = 0.7, na.rm = TRUE) # XXth percentile
  return(c(lower = lower_threshold, upper = upper_threshold))
}

# Apply the function across all columns
thresholds_perc <- apply(future_clim1[,-ncol(future_clim1)], 2, calculate_thresholds_perc)

# Transpose 
thresholds_perc <- t(thresholds_perc)

colnames(thresholds_perc) <- c('lower','upper')

up_vars <- c('Ann.Mean.Tmp', 'Mean.Diurn.Rng', 'Tmp.Seas', 'Max.Tmp.Wrm.M', 'Mean.Tmp.Wet.Q', 'Mean.Tmp.Dry.Q', 'Mean.Tmp.Wrm.Q', 'Prc.Seas')
down_vars <- c('Min.Tmp.Cld.M', 'Prc.Wet.M', 'Ann.Prc', 'Prc.Wet.Q', 'Prc.Dry.Q', 'Prc.Wet.M', 'Prc.Dry.M', 'Prc.Wrm.Q', 'Prc.Cld.Q')


# separate variables with lower vs. upper thresholds of interest

uppers_perc <- thresholds_perc %>% 
  as.data.frame() %>% 
  rownames_to_column('var') %>% 
  filter(var %in% up_vars) %>% 
  column_to_rownames('var')

downers_perc <- thresholds_perc %>% 
  as.data.frame() %>% 
  rownames_to_column('var') %>% 
  filter(var %in% down_vars) %>% 
  column_to_rownames('var')

# 3. Identify accessions that are over/under this threshold for current climates
#   A. for current climate

up.current.perc <- lapply(row.names(uppers_perc), function(x){
  threshold <- uppers_perc[x, 'upper']
  above<- envdat[which(envdat[,x] > threshold),]
  return(above$pi)
})

down.current.perc <- lapply(row.names(downers_perc), function(x){
  threshold <- downers_perc[x, 'lower']
  above<- envdat[which(envdat[,x] < threshold),]
  return(above$pi)
})

names(up.current.perc) <- row.names(uppers_perc)
names(down.current.perc) <- row.names(downers_perc)

######
# Function to create above/below dataframe
threshold_above <- function(accession_list, variable_list) {
  # Get all unique accessions from the list
  all_accessions <- envdat$pi
  
  # Initialize an empty dataframe
  result_df <- data.frame(pi = all_accessions)
  
  # Loop through each variable in the list
  for (var in names(variable_list)) {
    # Create a column for the variable
    result_df[[var]] <- ifelse(result_df$pi %in% variable_list[[var]], "above", "below")
  }
  
  return(result_df)
}

threshold_below <- function(accession_list, variable_list) {
  # Get all unique accessions from the list
  all_accessions <- envdat$pi
  
  # Initialize an empty dataframe
  result_df <- data.frame(pi = all_accessions)
  
  # Loop through each variable in the list
  for (var in names(variable_list)) {
    # Create a column for the variable
    result_df[[var]] <- ifelse(result_df$pi %in% variable_list[[var]], "below", "above")
  }
  
  return(result_df)
}

# Run function and melt result
thresh_above <- threshold_above(up.current.perc, up.current.perc)
thresh_below <- threshold_below(down.current.perc, down.current.perc)

thresh_above_m <- melt(thresh_above, id.vars = 'pi', value.name = 'thresh_pos', variable.name = 'var')
thresh_below_m <- melt(thresh_below, id.vars = 'pi', value.name = 'thresh_pos', variable.name = 'var')

# Prepare GEAV dataframe for merge
geav <- read.csv('./input_data/rrblup_GEBV_data.csv', row.names = 1)

geav_thresh_up <- geav[,colnames(geav) %in% colnames(thresh_above)]
geav_thresh_down <- geav[,colnames(geav) %in% colnames(thresh_below)]

geav_thresh_up <- geav_thresh_up[envdat$gen_id,]
geav_thresh_down <- geav_thresh_down[envdat$gen_id,]

geav_thresh_up$pi <- envdat$pi
geav_thresh_down$pi <- envdat$pi

geav_thresh_up_m <- melt(geav_thresh_up, id.vars = 'pi', value.name = "GEAV", variable.name = 'var')
geav_thresh_down_m <- melt(geav_thresh_down, id.vars = 'pi', value.name = "GEAV", variable.name = 'var')

# Merge dataframes
thresh_full_up <- merge(thresh_above_m, geav_thresh_up_m, by = c('pi', 'var'))
thresh_full_down <- merge(thresh_below_m, geav_thresh_down_m, by = c('pi', 'var'))

## Plot results

ggplot(thresh_full_up, aes(y = GEAV, x = as.factor(thresh_pos))) +
  geom_boxplot() +
  labs(x = "Position relative to threshold", y = "GEAV") +
  facet_wrap(~var, scales = 'free_y')

ggplot(thresh_full_down, aes(y = GEAV, x = as.factor(thresh_pos))) +
  geom_boxplot() +
  labs(x = "Position relative to threshold", y = "GEAV") +
  facet_wrap(~var, scales = 'free_y')

#################################################################################################
### Function to calculate percent of future environments where an accession's environment of origin
### is more extreme (e.g. higher temperature or lower precipitation)

# Function to calculate percentages
calculate_percentage <- function(envdat, future_clim, vars, direction = "higher") {
  # Select only the columns for the specified variables
  env_subset <- envdat[, vars, drop = FALSE]
  future_subset <- future_clim[, vars, drop = FALSE]
  
  # Calculate percentages
  if (direction == "higher") {
    percentage <- sapply(vars, function(var){
      # Initialize an empty vector to store percentages
      percentages <- numeric(length(envdat[[var]]))
      
      # Loop through each value in envdat
      for (i in 1:length(envdat[[var]])) {
        # Calculate the percentage of values in future_clim greater than envdat[i]
        percentages[i] <- mean(envdat[[var]][i] > future_clim[[var]], na.rm = T) * 100
      }
      
      return(percentages)
    }
    )
  } else if (direction == "lower") {
    percentage <- sapply(vars, function(var){
      # Initialize an empty vector to store percentages
      percentages <- numeric(length(envdat[[var]]))
      
      # Loop through each value in envdat
      for (i in 1:length(envdat[[var]])) {
        # Calculate the percentage of values in future_clim greater than envdat[i]
        percentages[i] <- mean(envdat[[var]][i] < future_clim[[var]], na.rm = T) * 100
      }
      
      return(percentages)
    }
    )
  }
  
  # Create a new dataframe with percentages
  return(as.data.frame(percentage))
}


# Calculate for up_vars (higher percentage)
up_vars_higher <- calculate_percentage(envdat, future_clim1, up_vars, direction = "higher")

# Calculate for down_vars (lower percentage)
down_vars_lower <- calculate_percentage(envdat, future_clim1, down_vars, direction = "lower")

#   B. for GEAV

# read in GEAV file
geav <- read.csv('./input_data/rrblup_GEBV_data.csv', row.names = 1)

### 4. plot GEAV vs. percent of overlap
# calculate percentage of future environments where accession environment is beyond
# envdat, future_clim

# Plot GEAV vs. overlap

geav_lower <- geav[,colnames(geav) %in% colnames(down_vars_lower)]
geav_upper <- geav[,colnames(geav) %in% colnames(up_vars_higher)]

pi_genid <- envdat[,1:2]

geav_lower_ordered <- geav_lower[pi_genid$gen_id,]
geav_upper_ordered <- geav_upper[pi_genid$gen_id,]

geav_merge_lower <- cbind(geav_lower_ordered, pi = envdat$pi)
geav_merge_upper <- cbind(geav_upper_ordered, pi = envdat$pi)

# melt table to bring variable columns to rows
geav_melt_lower <- melt(geav_merge_lower, id.vars = 'pi', value.name = "GEAV", variable.name = 'var')
geav_melt_upper <- melt(geav_merge_upper, id.vars = 'pi', value.name = "GEAV", variable.name = 'var')

# melt overlap dataframe
up_vars_higher$pi <- envdat$pi
down_vars_lower$pi <- envdat$pi

overlap_high <- melt(up_vars_higher, id.vars = 'pi', value.name = "overlap.perc", variable.name = 'var')
overlap_low <- melt(down_vars_lower, id.vars = 'pi', value.name = "overlap.perc", variable.name = 'var')

upper_all <- merge(overlap_high, geav_melt_upper, by = c('pi','var'))
lower_all <- merge(overlap_low, geav_melt_lower, by = c('pi','var'))

### Add a column with top and bottom 10% for plotting color
# Function to classify values into top 10%, bottom 10%, and middle
# Function to classify values into top 10%, bottom 10%, and middle

# Function to classify values into top 10%, bottom 10%, and middle
add_top_bottom_column <- function(df, vars_column) {
  # Get the unique variables from the vars column
  unique_vars <- unique(df$var)
  df$pos <- NA
  # Loop through each variable in the vars column
  for (var in unique_vars) {
    # Subset the data for the current variable
    subset_values <- df[which(df$var == var), "GEAV"]
    
    # Calculate the 10th and 90th percentiles for the subset
    top_threshold <- quantile(subset_values, 0.9, na.rm = TRUE)
    bottom_threshold <- quantile(subset_values, 0.1, na.rm = TRUE)
    
    # Classify the values and update a new column in the original dataframe
    df$pos[which(df$var == var)] <- ifelse(
      subset_values >= top_threshold, "top10",
      ifelse(subset_values <= bottom_threshold, "bottom10", "middle")
    )
  }
  #browser()
  return(df)
}


out_up <- add_top_bottom_column(df = upper_all, "var")
out_down <- add_top_bottom_column(lower_all, "var")

# Set names for labelling

biovar_names <- as_labeller(c(`Ann.Mean.Tmp` = "Annual Mean Temperature", `Max.Tmp.Wrm.M` = "Maximum Temperature, Warmest Month", 
                              `Mean.Tmp.Dry.Q` = 'Mean Temperature, Driest Quarter', `Mean.Diurn.Rng` = "Mean Diurnal Range", 
                              `Tmp.Seas` = 'Temperature Seasonality',`Mean.Tmp.Wrm.M` = 'Mean Temperature, Warmest Month',
                              `Mean.Tmp.Wet.Q` = 'Mean Temperature, Wettest Quarter', `Mean.Tmp.Wrm.Q` = 'Mean Temperature, Warmest Quarter',
                              `Prc.Seas` = 'Precipitation Seasonality'))

ggplot(out_up, aes(x = GEAV, y = overlap.perc, color = pos)) +
  geom_point() +
  ylim(0,100) +
  scale_color_manual(
    labels = c("Bottom 10%","Center","Top 10%"),
    values = c(
      "top10" = "red",         # Set color for "top10"
      "bottom10" = "blue",     # Set color for "bottom10"
      "middle" = "gray"        # Set color for "middle"
    )
  ) +
  labs(x = "GEAV", y = "Percent of Future Environments Exceeding", color = "GEAV Value") +
  facet_wrap(~var, scales = 'free_x', nrow = 2, labeller = biovar_names)

ggplot(out_down, aes(x = GEAV, y = overlap.perc, color = pos)) +
  geom_point() +
  ylim(0,100) +
  scale_color_manual(
    labels = c("Bottom 10%","Center","Top 10%"),
    values = c(
      "top10" = "red",         # Set color for "top10"
      "bottom10" = "blue",     # Set color for "bottom10"
      "middle" = "gray"        # Set color for "middle"
    )
  ) +
  labs(x = "GEAV", y = "Percent of Future Environments Below", color = "GEAV Value") +
  facet_wrap(~var, scales = 'free_x')

#save.image(file='./revision_scripts/future_climate/threshold_workspace.RData')

# Load future climate data
#load(file='./revision_scripts/future_climate/threshold_workspace.RData')

# Future climate agreement values (extracted for each accession) stored in future_clim1
# GEAV values stored in object geav

## Get 70% threshold for Max Temperature of Warmest Month from Kenya sorghum growing areas
bio5 <- raster(agreement_bioclim[[5]])
threshold_bio5 <- quantile(bio5, probs=c(0.75), na.rm = T)*10

bio1 <- raster(agreement_bioclim[[1]])
threshold_bio1 <- quantile(bio1, probs=c(0.75), na.rm = T)*10

bio9 <- raster(agreement_bioclim[[9]])
threshold_bio9 <- quantile(bio9, probs=c(0.75), na.rm = T)*10

#geav <- read.csv('./input_data/rrblup_GEBV_data.csv', row.names = 1)
#envdat <- read.csv('./revision_scripts/future_climate/envdat_master_koppen-geiger.csv')
pi_genid <- envdat[,1:2]

geav_ordered <- geav[pi_genid$gen_id,]

geav_merge <- cbind(geav_ordered, pi = envdat$pi)

geav_bio <- geav_merge[,c(40,15,18,21)]

#colnames(geav_merge1)[2:4] <- c('GEAV_Ann.Mean.Tmp', 'GEAV_Max.Tmp.Wrm.M'

geav_bio_m <- reshape2::melt(geav_bio, measure.vars = 2:4, id.vars = 'pi', value.name = 'GEAV', variable.name = 'var')

#### read in current climate data and create above, below classification

envdat <- read.csv('G:/My Drive/Sorghum/R/revision_scripts/future_climate/envdat_master_koppen-geiger.csv')

envdat_s <- subset(envdat, select = c('pi', 'Ann.Mean.Tmp', 'Max.Tmp.Wrm.M', 'Mean.Tmp.Dry.Q'))

envdat_m <- reshape2::melt(envdat_s, id.vars = 'pi', value.name = 'value', variable.name = 'var')

envdat_m$pos[envdat_m$value >= threshold_bio5 & envdat_m$var == 'Max.Tmp.Wrm.M'] <- 'above'
envdat_m$pos[envdat_m$value < threshold_bio5 & envdat_m$var == 'Max.Tmp.Wrm.M'] <- 'below'

envdat_m$pos[envdat_m$value >= threshold_bio1 & envdat_m$var == 'Ann.Mean.Tmp'] <- 'above'
envdat_m$pos[envdat_m$value < threshold_bio1 & envdat_m$var == 'Ann.Mean.Tmp'] <- 'below'

envdat_m$pos[envdat_m$value >= threshold_bio9 & envdat_m$var == 'Mean.Tmp.Dry.Q'] <- 'above'
envdat_m$pos[envdat_m$value < threshold_bio9 & envdat_m$var == 'Mean.Tmp.Dry.Q'] <- 'below'

pos_envdat <- envdat_m[,c(1,2,4)]

geav_envdat <- merge(geav_bio_m, pos_envdat, by = c('pi', 'var'))
geav_envdat <- na.omit(geav_envdat)
# plot results

var_names <- as_labeller(c(`Ann.Mean.Tmp` = "Annual Mean Temperature", `Max.Tmp.Wrm.M` = "Maximum Temperature, Warmest Month", 
                      `Mean.Tmp.Dry.Q` = 'Mean Temperature, Driest Quarter'))

ggplot(geav_envdat, aes(y = GEAV, x = as.factor(pos))) +
  geom_boxplot() +
  ggplot2::labs(x = "Current climate relative to future climate threshold", y = "GEAV") +
  facet_wrap(~var, scales = 'free_y', labeller = var_names)

