
setwd('G:/My Drive/Sorghum/R')

dat <- read.csv("./chromosome_map/marker_effects_chromosome_map.csv", header = T)
het <- read.csv("./input_data/envdat_master.csv")
num <- read.table("./sequence_data/sorghum_GD_numeric_rrblupformat.txt", head = T)


 # Remove the columns of the marker effects dataframe that are not well quantified
# dat <- dat[,c(1:9,11:21,34:35)]

# Prepare marker file for the heatmap
# Multiply the marker effects by the marker file
# row.names(num) <- num$taxa
gd2 <- num[num$taxa %in% het$gen_id,]
# g.in <- gd2[,-1]


#colnames(g.in)[1] <- 'gen_id'
row.names(gd2) <- gd2$taxa
g.in <- gd2[,-1]
#g.map <- as.matrix(g.in)
mat <- t(g.in)
mat_df <- as.data.frame(mat)
mat_df$Markers <- row.names(mat_df)

#colnames(mat) <- mat[1,] 
#mat <- mat[-1,]
dat <- dat[dat$rs. %in% mat_df$Markers,]

het <- subset(het, select = c('gen_id', 'topsoil_pH', 'Prc.Seas', 'Tmp.Seas'))

colnames(dat)[1:4] <- c('Markers','topsoil_pH', 'Prc.Seas', 'Tmp.Seas')

chrom_map <- dat[,c(1,5)]

mat_df1 <- merge(mat_df, chrom_map, by = 'Markers')

rownames(mat_df1) <- mat_df1$Markers
mat <- mat_df1[,-1]

#mat$chrom <- dat$chrom

#saveRDS(mat, "./chromosome_map/transposed_sorghum_genotype_matrix_for_chrommaps.rds")

comps <- colnames(dat)[2:4]

filelist <- list()
for(i in  comps){
  print(paste0(i, "chunk 1"))
  mark <- mat[,colnames(mat)!='chrom']
  mark <- as.data.frame(mark)
  mark <- lapply(mark[,1:ncol(mark)], as.numeric)
  mark <- as.data.frame(lapply(mark, `*`, dat[,i]))
  mark$chrom <- mat[,colnames(mat) == 'chrom']
  #rownames(mark) <- dat$Markers
  mark$Markers <- rownames(mark)
  
  print(paste0(i, "chunk 2"))
  #mark$chrom <- gsub("chrom", "", mark$Marker)
  #mark$chrom <- sapply(strsplit(mark$chrom, "_"), "[[", 1)
  mark <- mark[,colnames(mark)!="Markers"]

  mark$chrom <- as.factor(mark$chrom)
  print('here')
  mark2 <- aggregate(x = mark[,colnames(mark)!="chrom"],       
                     by = list(mark$chrom),
                     FUN = sum,
                     na.rm = TRUE)
  
  print(paste0(i, "chunk 3"))
  colnames(mark2)[1] <- "chrom"
  mark2$chrom <- as.numeric(mark2$chrom)
  mark2 <- mark2[order(mark2$chrom),]
  
  print(paste0(i, "chunk 4"))
  mark2 <- t(mark2)
  mark2 <- mark2[-1,]
  mark2 <- as.data.frame(mark2)
  mark2$gen_id <- rownames(mark2)
  mark2$var <- i
  #mark2 <- merge(mark2, het, by = "gen_id")
  
  print(paste0(i, "chunk 5"))
  filename <- paste(i, "Effects", sep = "_")
  new <- filename
  filelist[[length(filelist) + 1]] <- new
  assign(filename, mark2)
}

#write.csv(topsoil_pH_Effects, "Chromosomal_effects_topsoilpH.csv")
#write.csv(Prc.Seas_Effects, "Chromosomal_effects_Prc.Seas.csv")
#write.csv(Tmp.Seas_Effects, "Chromosomal_effects_Tmp.Seas.csv")

 ############################# Complex Heat Map #################################

library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

envdat <- read.csv('envdat_master.csv')
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
                       
                       #foo = anno_empty(border = FALSE,
                        #                  height = max_text_width(countries)/3))

#anno = anno_link(align_to = countries, which = "column", panel_fun = function(index){anno_text(countries[index], rot = 30, gp = gpar(fontsize = 16))}, 
#                 size = unit(2, "cm"), gap = unit(1, "cm"), width = unit(4, "cm"))

#ha2 = columnAnnotation(country = anno_text(month.name, location = 1, rot = 30, 
#                                         just = "right", gp = gpar(fontsize = 1:12+4)))

#ha3 = columnAnnotation(foo = anno_empty(border = FALSE,
#                                    width = max_text_width(countries) + unit(4, "mm")))

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



