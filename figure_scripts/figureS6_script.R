
library(ggplot2)
library(dplyr)
library(ggpmisc)
setwd("G:/My Drive/Sorghum/R/")

r.score <- read.csv('./score_outputs/score_futureclimates_sorghum_accession.csv')
ac.score <- read.csv('./score_outputs/score_genomicselection_accession.csv')

head(ac.score)
head(r.score)

colnames(ac.score)[6] <- "ac.score"
colnames(r.score)[6] <- "r.score"

comb.score <- merge(ac.score, r.score, by = c("pi", "country", "region"))

comb.score$region[which(comb.score$region == "South Asia")] <- "Asia"

comb.score$region[which(comb.score$region == "East Asia")] <- "Asia"

a<-with(comb.score, lm(ac.score~r.score))

summary(a)

ggplot(comb.score, aes(x = ac.score, y = r.score)) +
  geom_point(aes(color = region)) +
  geom_smooth(method = "lm") + 
  stat_poly_eq(use_label(c("adj.R2","p")), label.y = "top", label.x = "right") +
  labs(x = "Adaptive Capacity Score (Genomic Prediction)", y = "Resilience Score (Future Climate)")
