################################################################################
################   RR-BLUP Genomic Selection Score ##############################
################################################################################
library(dplyr)

# Read data
gebv_dat <- read.csv('./input_data/rrblup_GEBV_data.csv', head = T)
markers <- readRDS('./input_data/rrblup_markereffects_data.csv')
envdat <- read.csv('./input_data/envdat_master.csv')


geo <- envdat[,1:4]

colnames(gebv_dat)[1] <- 'gen_id'

gebv1 <- merge(gebv_dat, geo, by = 'gen_id')

gebv_list <- list()

gs_find_topx <- function(gs_dat, top_x){
  gs_dat1 <- subset(gs_dat, select = c(topsoil_pH:prec_Q4))
  apply(gs_dat1, MARGIN = 2, function (x){
    x1 <- gs_dat[order(x, decreasing = TRUE),]
    x2 <- x1[1:top_x,]
    return(x2)
  })
}

test <- gs_find_topx(gs_dat = gebv1, top_x = 100)

test1 <- list()

for(i in 1:length(test)){
  test1[[i]] <- test[[i]][,c(1,41:43,i+1)]
  colnames(test1[[i]])[5] <- 'Value'
}

names(test1) <- names(test)

######## Score Function ########

gs_score_calc <- function(gebv_list, var_names, geodat){
  list_topx <- gebv_list[var_names]
  df <- data.frame(matrix(data = NA, nrow = 0, ncol = ncol(gebv_list[[1]])))
  colnames(df) <- colnames(gebv_list[[1]])
  
  for(i in 1:length(list_topx)){ # combine the list of top 100 accessions per variable into one dataframe
   comb <- cbind(list_topx[[i]], var = var_names[[i]])
   df <- rbind(df, comb)
  }
  
  count.table <- df %>% group_by(country) %>% tally() # count number of accessions from each country
  names(count.table)[2] <- "count"
  
  envdat_vars <- geodat[,var_names]
  envdat_meta <- geodat[,1:4]
  geodat1 <- cbind(envdat_meta, envdat_vars)
  
  score.table <- geodat1 %>% group_by(country) %>% count() # count total number of accessions per country
  names(score.table)[2] <- "total"
  score.table <- merge(score.table, count.table, by = "country", all.x = T)
  score.table$count[which(is.na(score.table$count))] <- 0
  score.table$nvars <- length(var_names)
  score.table$score <- score.table$count/(score.table$total * score.table$nvars) # create normalized score

  # create empty count table for accession-level score
  count.acc <- df %>% group_by(pi) %>% tally()
  names(count.acc)[2] <- "count"
  count.acc$score <- count.acc$count / length(var_names)
  scores.acc <- merge(geodat1, count.acc, all.x = T)
  
  # count accession per variable
  varxacc <- df %>% count(pi, var)
  
  return(list(score.table, scores.acc, varxacc))
}


# Temperature Score
# select variables that came out of cross-validation over 50% correlation
temp_out_vars <- c("Ann.Mean.Tmp", "Tmp.Seas", "Max.Tmp.Wrm.M", "Min.Tmp.Cld.M", 
                     "Mean.Tmp.Dry.Q", "Mean.Tmp.Wrm.Q", "Mean.Tmp.Cld.Q", 
                    "tmean_Q1", "tmean_Q2", "tmean_Q3", "tmean_Q4")

gs_temp <- gs_score_calc(gebv_list = test1, var_names = temp_out_vars, geodat = envdat)

tempscore_country <- gs_temp[[1]]
tempscore_acc <- gs_temp[[2]]

write.csv(tempscore_country, './score_outputs/score_genomicselection_temperature_country.csv')

write.csv(tempscore_acc, './score_outputs/score_genomicselection_temperature_accession.csv')

###### Precipitation Score #######

precip_out_vars <- c('Prc.Dry.Q', 'Prc.Dry.M', 'Prc.Seas', 'prec_Q1', 'prec_Q3', 'prec_Q4')

gs_precip <- gs_score_calc(gebv_list = test1, var_names = precip_out_vars, geodat = envdat)

precipscore_country <- gs_precip[[1]]
precipscore_acc <- gs_precip[[2]]

write.csv(precipscore_country, './score_outputs/score_genomicselection_precipitation_country.csv')

write.csv(precipscore_acc, './score_outputs/score_genomicselection_precipitation_accession.csv')

#### Overall Accession Score Table #####
acc.score <- data.frame(pi = precipscore_acc$pi, country = precipscore_acc$country, region = precipscore_acc$region,
                        score.prc = precipscore_acc$score, score.temp = tempscore_acc$score)

acc.score$score.prc[which(is.na(acc.score$score.prc))] <- 0
acc.score$score.temp[which(is.na(acc.score$score.temp))] <- 0

acc.score$score <- acc.score$score.prc + acc.score$score.temp

write.csv(acc.score, "./score_outputs/score_genomicselection_accession.csv", row.names = F)

###### Overall Country Score ######

gs.scores <- as.data.frame(cbind(country = tempscore_country$country, temp_score = tempscore_country$score, prc_score = precipscore_country$score))

# Write Region into score dataframe
country_region <- data.frame(country = envdat$country, region = envdat$region)
country_region <- country_region[!duplicated(country_region$country),]
gs.scores <- merge(gs.scores, country_region, by = "country")
gs.scores <- gs.scores[,c(1, 4, 2, 3)] # reorder data frame

# Combine East Asia and South Asia
gs.scores$region[which(gs.scores$region == "East Asia" | gs.scores$region == "South Asia")] <- "Asia"

# Format Data (NAs -> 0, write as numeric data)
gs.scores <- gs.scores %>% replace(is.na(.), 0)
gs.scores$temp_score <- as.numeric(gs.scores$temp_score)
gs.scores$prc_score <- as.numeric(gs.scores$prc_score)

# Remove countries with less than 5 accessions
count.table <- envdat %>% group_by(country) %>% tally()
colnames(count.table)[2] <- "num_accessions"
gs.scores <- merge(gs.scores, count.table, by = "country")
gs.scores <- gs.scores[which(gs.scores$num_accession > 5),] 

# Sum across temp/precip scores
gs.scores$gs.score <- gs.scores$temp_score + gs.scores$prc_score

write.csv(gs.scores, "./score_outputs/score_genomicselection_country.csv")


################## Regional Score ####################
# Bind score categories into single dataframe
gs.regional <- as.data.frame(cbind(country = tempscore_country$country, total = tempscore_country$total, temp_count = tempscore_country$count, 
                                   prc_count = precipscore_country$count))

# Write Region into score dataframe
country_region <- data.frame(country = envdat$country, region = envdat$region)
country_region <- country_region[!duplicated(country_region$country),]
gs.regional <- merge(gs.regional, country_region, by = "country")
gs.regional <- gs.regional[,c(1, 5, 2:4)] # reorder data frame

# Combine South and East Asia
gs.regional$region[which(gs.regional$region == "East Asia" | gs.regional$region == "South Asia")] <- "Asia"

# Format as numeric
gs.regional[,3:5] <- apply(gs.regional[,3:5], 2, as.numeric)

# Combine counts by region
gs.regionalscores <- gs.regional %>% group_by(region) %>% 
  summarise(temp_count = sum(temp_count),
            prc_count = sum(prc_count),
            total = sum(total),
            .groups = 'drop')


# Regional Score
gs.regionalscores$temp_score <- gs.regionalscores$temp_count/(gs.regionalscores$total * mean(tempscore_country$nvars))
gs.regionalscores$prc_score <- gs.regionalscores$prc_count/(gs.regionalscores$total * mean(precipscore_country$nvars))

gs.regionalscores <- as.data.frame(gs.regionalscores)

gs.regionalscores$gs.score <- gs.regionalscores$temp_score + gs.regionalscores$prc_score

write.csv(gs.regionalscores, "./score_outputs/genomic_selection_score_region.csv", row.names = F)

#Minicore score
minicore <- read.csv('./input_data/envdat_minicore_genomicprediction.csv')

mini_score <- acc.score[which(minicore$is.minicore),]

write.csv(mini_score, "./score_outputs/score_genomicselection_minicore_accession.csv")

