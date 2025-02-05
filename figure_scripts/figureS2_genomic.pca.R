library(ggplot2)
library(ggpubr)

# Read PCA output and geographical/botanical race data
pca_scores3 <- read.csv("./cluster_analysis/genomic_pca_scores.csv")
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
r.score <- read.csv('./score_outputs/score_futureclimates_sorghum_accession_updatedregions.csv')
r_geo_pca <- merge(geo_pca, r.score, by = c("pi", "region", "country"))

ggpubr::ggscatter(data = r_geo_pca,
                  y = "PC2",
                  x = "PC1",
                  color = "score",
                  xlab = "PC1 (2.8% variation)",
                  ylab = "PC2 (1.1% variation)",
                  main = "PC1 and PC2 % Variation")

# by genomic selection index
ac.score <- read.csv('./score_outputs/score_genomicselection_accession.csv_updatedregions')
ac_geo_pca <- merge(geo_pca, ac.score, by = c("pi", "region", "country"))

ggpubr::ggscatter(data = r_geo_pca,
                  y = "PC2",
                  x = "PC1",
                  color = "score",
                  xlab = "PC1 (2.8% variation)",
                  ylab = "PC2 (1.1% variation)",
                  main = "PC1 and PC2 % Variation")
