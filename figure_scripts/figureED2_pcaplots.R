library(ggplot2)
library(ggpubr)

# Read PCA output and geographical/botanical race data
pca_scores3 <- read.csv("./cluster_analysis/genomic_pca_scores.csv")
envdat <- read.csv('./input_data/envdat_master_races.csv')

geodat <- envdat[,c(1:4,7)]

geo_pca <- merge(envdat, pca_scores3, by = "gen_id")

# by landrace
g_a <- ggpubr::ggscatter(data = geo_pca,
                  y = "PC2",
                  x = "PC1",
                  color = "race",
                  xlab = "PC1 (2.8% variation)",
                  ylab = "PC2 (1.1% variation)",
                  main = "PC1 and PC2 % Variation",
                  legend = 'right')

# by region
g_b <- ggpubr::ggscatter(data = geo_pca,
                  y = "PC2",
                  x = "PC1",
                  color = "region",
                  xlab = "PC1 (2.8% variation)",
                  ylab = "PC2 (1.1% variation)",
                  main = "PC1 and PC2 % Variation", legend = 'right')

# by future climate index
r.score <- read.csv('./score_outputs/score_futureclimates_sorghum_accession.csv')
r_geo_pca <- merge(geo_pca, r.score, by = c("pi", "region", "country"))

g_c <- ggpubr::ggscatter(data = r_geo_pca,
                  y = "PC2",
                  x = "PC1",
                  color = "score",
                  xlab = "PC1 (2.8% variation)",
                  ylab = "PC2 (1.1% variation)",
                  main = "PC1 and PC2 % Variation") +
  scale_color_viridis()

# by genomic selection index
ac.score <- read.csv('./score_outputs/score_genomicselection_accession.csv')
ac_geo_pca <- merge(geo_pca, ac.score, by = c("pi", "region", "country"))

g_d <- ggpubr::ggscatter(data = ac_geo_pca,
                  y = "PC2",
                  x = "PC1",
                  color = "score",
                  xlab = "PC1 (2.8% variation)",
                  ylab = "PC2 (1.1% variation)",
                  main = "PC1 and PC2 % Variation") +
  scale_color_viridis()

ggarrange(g_a, g_b, g_c, g_d)
