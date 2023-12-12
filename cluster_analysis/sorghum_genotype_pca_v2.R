#### Re-running PCA on SNP Data using:
# https://rpubs.com/Danya/983464

library(vcfR)
library(vegan)
library(ggplot2)
library(ggpubr)
library(dataPreparation)
library(data.table)
setwd("G:/My Drive/Sorghum/R/")

snps <- vcfR::read.vcfR("./sequence_data/genotypes_sorghum_landraces_all.vcf", convertNA  = TRUE)

snps_num <- vcfR::extract.gt(snps, 
                             element = "GT",
                             IDtoRowNames  = T,
                             as.numeric = T,
                             convertNA = T,
                             return.alleles = F)
rm(snps)
snps_num_t <- t(snps_num)
rm(snps_num)
#snps_num_t <- data.frame(snps_num_t) 

find_NAs <- function(x){
  NAs_TF <- is.na(x)
  i_NA <- which(NAs_TF == TRUE)
  N_NA <- length(i_NA)
  
  cat("Results:",N_NA, "NAs present\n.")
  return(i_NA)
}

# N_rows
# number of rows (individuals)
N_rows <- nrow(snps_num_t)

# N_NA
# vector to hold output (number of NAs)
N_NA   <- rep(x = 0, times = N_rows)

# N_SNPs
# total number of columns (SNPs)
N_SNPs <- ncol(snps_num_t)

# the for() loop
for(i in 1:N_rows){
  
  # for each row, find the location of
  ## NAs with snps_num_t()
  i_NA <- find_NAs(snps_num_t[i,]) 
  
  # then determine how many NAs
  ## with length()
  N_NA_i <- length(i_NA)
  
  # then save the output to 
  ## our storage vector
  N_NA[i] <- N_NA_i
}

# 50% of N_SNPs
cutoff50 <- N_SNPs*0.5

hist(N_NA)            
abline(v = cutoff50, 
       col = 2, 
       lwd = 2, 
       lty = 2)

percent_NA <- N_NA/N_SNPs*100

# Call which() on percent_NA
i_NA_50percent <- which(percent_NA > 50) 

snps_num_t02 <- snps_num_t[-i_NA_50percent, ]

row_names <- row.names(snps_num_t02) # Key



row_names02 <- gsub("sample_","",row_names)

sample_id <- gsub("^([ATCG]*)(_)(.*)",
                  "\\3",
                  row_names02)
pop_id <- gsub("[01-9]*",    
               "",
               sample_id)

table(pop_id)  

invar_omit <- function(x){
  cat("Dataframe of dim",dim(x), "processed...\n")
  sds <- apply(x, 2, sd, na.rm = TRUE)
  i_var0 <- which(sds == 0)
  
  
  cat(length(i_var0),"columns removed\n")
  
  if(length(i_var0) > 0){
    x <- x[, -i_var0]
  }
  
  ## add return()  with x in it
  return(x)                      
}


snps_no_invar <- invar_omit(snps_num_t02) 

snps_noNAs <- snps_no_invar

N_col <- ncol(snps_no_invar)
for(i in 1:N_col){
  
  # get the current column
  column_i <- snps_noNAs[, i]
  
  # get the mean of the current column
  mean_i <- mean(column_i, na.rm = TRUE)
  
  # get the NAs in the current column
  NAs_i <- which(is.na(column_i))
  
  # record the number of NAs
  N_NAs <- length(NAs_i)
  
  # replace the NAs in the current column
  column_i[NAs_i] <- mean_i
  
  # replace the original column with the
  ## updated columns
  snps_noNAs[, i] <- column_i
  
}

#write.csv(snps_noNAs, file = "./sequence_data/SNPs_cleaned.csv",
#          row.names = F)

rm(snps_no_invar)
rm(snps_num_t)
rm(snps_num_t02)
###################################################################################################################
# Part 2: Analyze Data
# https://rstudio-pubs-static.s3.amazonaws.com/981420_4bd07f0f519f4b56bcd235eab87a706b.html
SNPs_cleaned <- snps_noNAs
rm(snps_noNAs)
#SNPs_cleaned <- read.csv(file = "./sequence_data/SNPs_cleaned.csv")

#saveRDS(SNPs_cleaned, "./sequence_data/SNPs_cleaned.rds")
#SNPs_cleaned <- readRDS("./sequence_data/SNPs_cleaned.rds")

SNPs_scaled <- scale(SNPs_cleaned)

# remove zero variance columns
# requires data.table and dataPreparation packages
SNPs_dt <- as.data.table(SNPs_scaled)
cons_col <- which_are_constant(SNPs_dt, verbose=FALSE)
SNPs_scaled1 <- SNPs_scaled[,-cons_col]

############################################
## Run PCA 
############################################
pca_scaled <- prcomp(SNPs_scaled1)
#saveRDS(pca_scaled, "pca_imputed_all_SNPs.rds")

#pca_scaled <- readRDS("pca_imputed_all_SNPs.rds")
screeplot(pca_scaled, 
          ylab  = "Relative importance",
          main = "Scaled PCA Screeplot")

summary_out_scaled <- summary(pca_scaled)

PCA_variation <- function(pca_summary, PCs = 2){
  var_explained <- pca_summary$importance[2,1:PCs]*100
  var_explained <- round(var_explained,1)
  return(var_explained)
}

var_out <- PCA_variation(summary_out_scaled,PCs = 10)

N_columns <- ncol(SNPs_scaled)
barplot(var_out,
        main = "Percent variation Scree plot",
        ylab = "Percent variation explained")
abline(h = 1/N_columns*100, col = 2, lwd = 2)

# biplot(pca_scaled, cex = 0)

pca_scores <- vegan::scores(pca_scaled)

# Load metadata to add country/Score/Region/etc. for plotting

row_names <- read.csv('gen_ids_for_pca.csv')

pca_scores2 <- data.frame(gen_id = row_names,
                          pca_scores)

pca_scores3 <- pca_scores2[,1:11]

envdat <- read.csv('./input_data/envdat_master_races.csv')

geodat <- envdat[,c(1:4,7)]

geo_pca <- merge(envdat, pca_scores3, by = "gen_id")

# by landrace
ggpubr::ggscatter(data = geo_pca,
                  y = "PC2",
                  x = "PC1",
                  color = "race",
                  xlab = "PC1 (2.8% variation)",
                  ylab = "PC2 (1.1% variation)",
                  main = "PC1 and PC2 % Variation")

# by region
ggpubr::ggscatter(data = geo_pca,
                  y = "PC2",
                  x = "PC1",
                  color = "region",
                  xlab = "PC1 (2.8% variation)",
                  ylab = "PC2 (1.1% variation)",
                  main = "PC1 and PC2 % Variation")

# by future climate index
r.score <- read.csv('./score_outputs/score_futureclimates_sorghum_accession.csv')
r_geo_pca <- merge(geo_pca, ac.score, by = c("pi", "region", "country"))

ggpubr::ggscatter(data = r_geo_pca,
                  y = "PC2",
                  x = "PC1",
                  color = "score",
                  xlab = "PC1 (2.8% variation)",
                  ylab = "PC2 (1.1% variation)",
                  main = "PC1 and PC2 % Variation")

# by genomic selection index
ac.score <- read.csv('./score_outputs/score_genomicselection_accession.csv')
ac_geo_pca <- merge(geo_pca, ac.score, by = c("pi", "region", "country"))

ggpubr::ggscatter(data = r_geo_pca,
                  y = "PC2",
                  x = "PC1",
                  color = "score",
                  xlab = "PC1 (2.8% variation)",
                  ylab = "PC2 (1.1% variation)",
                  main = "PC1 and PC2 % Variation")
