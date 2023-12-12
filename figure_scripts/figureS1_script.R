########################################################################################
### Supplemental Figure 1

# Libraries
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(tibble)
library(tidyr)
library(rstatix)
#Load in data
setwd('G:/My Drive/Sorghum/R/')
envdat<- read.csv("./input_data/envdat_master.csv")
head(envdat)

envdat2<-envdat[,-c(1:7)]
#row.names(envdat2)<-envdat$country
head(envdat2)
envdat2_stan<-scale(envdat2)
head(envdat2_stan)

df_new <- cbind(envdat$pi, envdat$country,envdat$region, envdat2_stan)
head(df_new)
colnames(df_new)[c(1:3)] <- c("pi", "country", "region")
str(df_new)

write.csv(df_new,"stan_sorghum_env.csv", row.names = F)


data_wide<-read.csv("stan_sorghum_env.csv")
str(data_wide)

data_long <- gather(data_wide, variable, Value,topsoil_pH:prec_Q4 , factor_key=TRUE)

head(data_long)
tail(data_long)

get_box_stats <- function(y, upper_limit = max(data_long$Value) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Count =", length(y), "\n",
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2)), "\n"
  )
  )
}

# Box Plot
ggplot(data_long, aes(variable, Value, fill= variable)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  #geom_jitter(shape=16, position=position_jitter(0.05))+
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9) +
  geom_text(data=. %>% 
              group_by(variable) %>%  
              filter(Value %in% boxplot.stats(Value, coef=1.4)$out),
            aes(label=country, y=Value), nudge_x=0.1, colour="red", size=3, hjust=0.2) +
  labs(
    title = "Environmental Variation",
    x = "Range wide Variables",
    y = "Standarized variation in Environment") +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme_bw() +
  theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))