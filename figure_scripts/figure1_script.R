
# Recreate Figure 1 - Complex Heat Map

library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

Tmp.Seas_Effects <- read.csv("./input_data/Chromosomal_effects_topsoilpH.csv")
Prc.Seas_Effects <-read.csv("./input_data/Chromosomal_effects_Prc.Seas.csv")
topsoil_pH_Effects <- read.csv("./input_data/Chromosomal_effects_Tmp.Seas.csv")

envdat <- read.csv('./input_data/envdat_master.csv')
row.names(envdat) <- envdat$gen_id
envdat1 <- envdat[Tmp.Seas_Effects$gen_id,]

country <- subset(envdat1, select = c('gen_id', 'country', 'region'))

country$region[which(country$region == "East Asia" | country$region == "South Asia")] <- "Asia"

country.a <- country[with(country, order(region, country)), ]

column_order <- row.names(country.a)
regions <- country$region
countries <- country$country

c_ord <- unique(country.a$country)

out_vec <- vector()
for(i in c_ord){
  vec <- which(countries == i)
  pos <- vec[(length(vec) + 1L) %/% 2L]
  out_vec[i] <- pos
}

Tmp.Seas <- t(as.matrix(Tmp.Seas_Effects[,c(1:10)]))
Prc.Seas <- t(as.matrix(Prc.Seas_Effects[,c(1:10)]))
topsoil_pH <- t(as.matrix(topsoil_pH_Effects[,c(1:10)]))


library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 51
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)
names(mycolors) <- unique(countries)

ha = HeatmapAnnotation(Region = regions,
                       col = list(Region = c("Asia" = "pink", "East Africa" = "royalblue", "Middle East" = "forestgreen", 'South Africa' = 'red', 
                                             'West Africa' = 'purple')), show_legend = FALSE, show_annotation_name = FALSE)
ha1 = columnAnnotation(Country = countries,
                       col = list(Country = mycolors), show_legend = FALSE, 
                       foo = anno_mark(at = out_vec, labels = c_ord, side = 'bottom', labels_rot = 270), show_annotation_name = FALSE)

ht_list  = 
  Heatmap(Tmp.Seas, column_split = regions, column_gap = unit(0, "mm"), show_heatmap_legend = FALSE,
          row_order = c(1,2,3,4,5,6,7,8,9,10), top_annotation = ha, column_order = column_order, row_names_rot = 270,
          heatmap_legend_param = list(title = "Temperature Seasonality", legend_direction = 'horizontal'),
          row_names_side = "left", show_column_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE, show_parent_dend_line = FALSE) %v%
  
  Heatmap(Prc.Seas,
          row_order = c(1,2,3,4,5,6,7,8,9,10), column_order = column_order,  row_names_rot = 270, show_heatmap_legend = FALSE,
          heatmap_legend_param = list(title = "Precipitation Seasonality", legend_direction = 'horizontal'),
          row_names_side = "left", show_column_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE, show_parent_dend_line = FALSE) %v%
  
  Heatmap(topsoil_pH, 
          row_order = c(1,2,3,4,5,6,7,8,9,10), column_order = column_order, row_names_rot = 270, show_heatmap_legend = FALSE,
          heatmap_legend_param = list(title = "Topsoil pH", legend_direction = 'horizontal'), bottom_annotation = ha1, 
          row_names_side = "left", show_column_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE, show_parent_dend_line = FALSE)

pdf('chromosome_map.pdf', width = 12, height = 10)
ComplexHeatmap::draw(ht_list)
dev.off()

count.table <- envdat %>% group_by(country) %>% tally()
c_vec <- unique(country.a$country)
c_vec <- as.data.frame(c_vec)

row.names(count.table) <- count.table$country
count1 <- count.table[c_vec,]

count1$tf <- count1$n >= 20

for(i in 1:51) {
  if(count1$tf)
    decorate_annotation("foo", slice = i, {
      #grid.rect(y = 0, height = unit(2, "mm"), gp = gpar(fill = i, col = NA), just = "center")
      grid.text(paste(unique(country.a$country)[i], collapse = "\n"), x = unit(0, "mm"), just = "center", rot = 360-45)
    })
}

col_fun1 = colorRamp2(c(-1000, 0, 1000), c("blue", "white", "red"))
lgd_Tmp.Seas = Legend(col_fun = col_fun, title = "Temperature Seasonality", direction = "horizontal", legend_width=unit(1.6, 'in'))
ComplexHeatmap::draw(lgd_Tmp.Seas)

col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
lgd_Prc.Seas = Legend(col_fun = col_fun, title = "Precipitation Seasonality", direction = "horizontal", legend_width=unit(1.6, 'in'))
ComplexHeatmap::draw(lgd_Prc.Seas)

col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
lgd_topsoil_pH = Legend(c(-3, -1.5, 0, 1.5, 3), col_fun = col_fun, title = "Topsoil pH", direction = "horizontal", legend_width=unit(1.6, 'in'))
ComplexHeatmap::draw(lgd_topsoil_pH)