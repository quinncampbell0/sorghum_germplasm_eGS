
library(ggplot2)
library(dplyr)
library(ggpmisc)
#setwd("G:/My Drive/Sorghum/R/")

r.score <- read.csv('./score_outputs/score_futureclimates_sorghum_accession_updatedregions.csv')
ac.score <- read.csv('./score_outputs/score_genomicselection_accession_updatedregions.csv')

head(ac.score)
head(r.score)

colnames(ac.score)[6] <- "ac.score"
colnames(r.score)[6] <- "r.score"

comb.score <- merge(ac.score, r.score, by = c("pi", "country", "region"))

#comb.score$region[which(comb.score$region == "South Asia")] <- "Asia"

#comb.score$region[which(comb.score$region == "East Asia")] <- "Asia"

a<-with(comb.score, lm(ac.score~r.score))

summary(a)
hist(comb.score$ac.score)
hist(comb.score$r.score)

# calculate r (pearson's)
cor(comb.score$ac.score, comb.score$r.score, use = 'complete.obs', method = 'pearson')

# Create median lines
res.mean <- mean(comb.score$r.score, na.rm = T)
ac.mean <- mean(comb.score$ac.score, na.rm = T)

## Full plot ##
p_A <- ggplot(comb.score, aes(x = ac.score, y = r.score)) +
  geom_jitter(aes(color = region), width = 0.015, height = 0.01) +
  geom_smooth(method = "lm") + 
  stat_poly_eq(use_label(c("adj.R2","p")), label.y = "top", label.x = "right") +
  labs(x = "", y = "") +
  geom_hline(yintercept = res.mean, color = 'black', linetype = 'dashed', linewidth = 1) +
  geom_vline(xintercept = ac.mean, color = 'black', linetype = 'dashed', linewidth = 1) +
  scale_color_brewer(palette = 'Dark2') +
    annotate("rect", xmin = -Inf, xmax = ac.mean, ymin = -Inf, ymax = res.mean, fill = "red", alpha = 0.1) +
    annotate("rect", xmin = ac.mean, xmax = Inf, ymin = -Inf, ymax = res.mean, fill = "yellow", alpha = 0.1) +
    annotate("rect", xmin = -Inf, xmax = ac.mean, ymin = res.mean, ymax = Inf, fill = "lightblue", alpha = 0.1) +
    annotate("rect", xmin = ac.mean, xmax = Inf, ymin = res.mean, ymax = Inf, fill = "green", alpha = 0.1) + 
  theme_minimal()

## Plot with countries pulled out

library(ggplot2)
library(ggpubr)

# Add a custom grouping variable to highlight Burundi and Ethiopia
comb.score$highlight <- ifelse(comb.score$country %in% c("Burundi", "Ethiopia"), comb.score$country, "Other")

# Plot
p_B <- ggplot(comb.score, aes(x = ac.score, y = r.score)) +
  # Use the custom grouping variable for color
  geom_jitter(aes(color = highlight), width = 0.015, height = 0.01, size = 2, alpha = ifelse(comb.score$highlight == "Other", 0.2, 1)) +
  #geom_smooth(method = "lm") + 
  #stat_poly_eq(use_label(c("adj.R2", "p")), label.y = "top", label.x = "right") +
  labs(
    x = "Adaptive Capacity Score (Genomic Prediction)",
    y = "Resilience Score (Future Climate)",
    color = "Country"
  ) +
  # Add mean lines
  geom_hline(yintercept = res.mean, color = 'black', linetype = 'dashed', linewidth = 1) +
  geom_vline(xintercept = ac.mean, color = 'black', linetype = 'dashed', linewidth = 1) +
  # Add colored quadrants
  annotate("rect", xmin = -Inf, xmax = ac.mean, ymin = -Inf, ymax = res.mean, fill = "red", alpha = 0.1) +
  annotate("rect", xmin = ac.mean, xmax = Inf, ymin = -Inf, ymax = res.mean, fill = "yellow", alpha = 0.1) +
  annotate("rect", xmin = -Inf, xmax = ac.mean, ymin = res.mean, ymax = Inf, fill = "lightblue", alpha = 0.1) +
  annotate("rect", xmin = ac.mean, xmax = Inf, ymin = res.mean, ymax = Inf, fill = "green", alpha = 0.1) +
  # Customize colors
  scale_color_manual(
    values = c("Burundi" = "#FF6F59", "Ethiopia" = "#03B5AA","Other" = "grey")
  ) +
  theme_minimal()

# Combine into one plot
ggarrange(p_A, p_B, ncol = 1)

## All countries for visualization

ggplot(comb.score, aes(x = ac.score, y = r.score)) +
  geom_point(aes(color = country)) +
  geom_smooth(method = "lm") + 
  stat_poly_eq(use_label(c("adj.R2","p")), label.y = "top", label.x = "right") +
  labs(x = "Adaptive Capacity Score (Genomic Prediction)", y = "Resilience Score (Future Climate)") +
  geom_hline(yintercept = res.mean, color = 'black', linetype = 'dashed', linewidth = 1) +
  geom_vline(xintercept = ac.mean, color = 'black', linetype = 'dashed', linewidth = 1) +
  annotate("rect", xmin = -Inf, xmax = ac.mean, ymin = -Inf, ymax = res.mean, fill = "red", alpha = 0.1) +
  annotate("rect", xmin = ac.mean, xmax = Inf, ymin = -Inf, ymax = res.mean, fill = "yellow", alpha = 0.1) +
  annotate("rect", xmin = -Inf, xmax = ac.mean, ymin = res.mean, ymax = Inf, fill = "lightblue", alpha = 0.1) +
  annotate("rect", xmin = ac.mean, xmax = Inf, ymin = res.mean, ymax = Inf, fill = "green", alpha = 0.1)
