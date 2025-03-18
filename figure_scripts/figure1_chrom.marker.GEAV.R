
############################# Complex Heat Map - Panel A #################################
library(tidyverse)
library(dplyr)
library(RColorBrewer)
#library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(SNPRelate)
library("FactoMineR")
library("factoextra")

#setwd('G:/My Drive/Sorghum/R')
pet_min_Effects <- read.csv("./input_data/Chromosomal_effects_PET_min.csv", row.names = 1)
Tmp.Seas_Effects <- read.csv("./input_data/Chromosomal_effects_Tmp.Seas.csv", row.names = 1)
Prc.Seas_Effects <-read.csv("./input_data/Chromosomal_effects_Prc.Seas.csv", row.names = 1)
topsoil_pH_Effects <- read.csv("./input_data/Chromosomal_effects_topsoilpH.csv", row.names = 1)

envdat <- read.csv('./input_data/envdat_master_updatedregions.csv')
row.names(envdat) <- envdat$gen_id
envdat1 <- envdat[pet_min_Effects$gen_id,]
envdat1$column_order <- seq(1:1937)

country <- subset(envdat1, select = c('gen_id', 'country', 'region', 'column_order'))

country.a <- country[with(country, order(country)), ]

# Re-order the levels
country.a$region <- factor(country.a$region, levels=c('West Africa', 'Central Africa', 'South Africa', 'East Africa', 'West Asia', 'South Asia', 'East Asia'))
# Re-order the data.frame
country.a <- country.a[order(country.a$region),]

column_order <- country.a$column_order

# Reorder dataframe of chromosomal effects
pet_min_Effects <- pet_min_Effects[row.names(country.a),]
Tmp.Seas_Effects <- Tmp.Seas_Effects[row.names(country.a), ]
Prc.Seas_Effects <- Prc.Seas_Effects[row.names(country.a), ]
topsoil_pH_Effects <- topsoil_pH_Effects[row.names(country.a), ]

#regions <- country$region
#countries <- country$country

c_ord <- unique(country.a$country)

out_vec <- vector()
for(i in c_ord){
  vec <- which(country.a$country == i)
  pos <- vec[(length(vec) + 1L) %/% 2L]
  out_vec[i] <- pos
}

### Custom colors 

# Define the number of colors you want
nb.cols <- 51
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)
names(mycolors) <- sample(unique(country.a$country))

# Custom region colors
rg.cols <- 7
regcols <- colorRampPalette(brewer.pal(7, "Paired"))(rg.cols)
names(regcols) <- sample(unique(country.a$region))

pet_min_m <- t(as.matrix(pet_min_Effects[,c(1:10)]))
Tmp.Seas <- t(as.matrix(Tmp.Seas_Effects[,c(1:10)]))
Prc.Seas <- t(as.matrix(Prc.Seas_Effects[,c(1:10)]))
topsoil_pH <- t(as.matrix(topsoil_pH_Effects[,c(1:10)]))
# Create overall score (summed across chromosomes)
m_list <- list(pet_min_m,Tmp.Seas, Prc.Seas, topsoil_pH)

tot_list <- lapply(m_list, function(x){
  x <- x %>% as.data.frame(.) %>% 
    bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total")) %>% as.matrix(.)
})

pet_min <- tot_list[[1]]
Prc.Seas <- tot_list[[2]]
Tmp.Seas <- tot_list[[3]]
topsoil_pH <- tot_list[[4]]

# Heat Map Annotations

ha = HeatmapAnnotation(Region = country.a$region, # region annotation
                       col = list(Region = regcols), 
                       show_legend = FALSE, show_annotation_name = FALSE)

ha1 = columnAnnotation(Country = country.a$country, # country annotation
                       col = list(Country = mycolors), show_legend = FALSE, 
                       labels = anno_mark(at = out_vec, labels = c_ord, side = 'bottom', labels_rot = 270, labels_gp = gpar(fontsize = 8)), show_annotation_name = FALSE)

####### Combined heatmap - PET plus bioclim vars

ht_list3 = 
  Heatmap(pet_min, column_split = country.a$region, column_gap = unit(0, "mm"), show_heatmap_legend = FALSE,
          row_order = c(1,2,3,4,5,6,7,8,9,10,11), top_annotation = ha, column_order = row.names(country.a$gen_id), row_names_rot = 0,
          heatmap_legend_param = list(title = "Aridity Index", legend_direction = 'horizontal'),
          row_names_side = "left", row_labels = c(1:10, "Total"), show_column_names = FALSE, show_row_dend = FALSE, 
          show_column_dend = FALSE, show_parent_dend_line = FALSE,  cluster_columns = FALSE, cluster_rows = FALSE) %v%
  
  Heatmap(Tmp.Seas,
          row_order = c(1,2,3,4,5,6,7,8,9,10,11), column_order = row.names(country.a$gen_id),  row_names_rot = 0, show_heatmap_legend = FALSE,
          heatmap_legend_param = list(title = "Heat Index", legend_direction = 'horizontal'),
          row_names_side = "left", row_labels = c(1:10, "Total"), show_column_names = FALSE, show_row_dend = FALSE, 
          show_column_dend = FALSE, show_parent_dend_line = FALSE, cluster_columns = FALSE, cluster_rows = FALSE) %v%
  
  Heatmap(Prc.Seas, 
          row_order = c(1,2,3,4,5,6,7,8,9,10,11), column_order = row.names(country.a$gen_id), row_names_rot = 0, show_heatmap_legend = FALSE,
          heatmap_legend_param = list(title = "Cold Index", legend_direction = 'horizontal'),  
          row_names_side = "left", row_labels = c(1:10, "Total"), show_column_names = FALSE, show_row_dend = FALSE, 
          show_column_dend = FALSE, show_parent_dend_line = FALSE, cluster_columns = FALSE, cluster_rows = FALSE)%v%
  
  Heatmap(topsoil_pH, 
          row_order = c(1,2,3,4,5,6,7,8,9,10,11), column_order = row.names(country.a$gen_id), row_names_rot = 0, show_heatmap_legend = FALSE,
          heatmap_legend_param = list(title = "Cold Index", legend_direction = 'horizontal'), bottom_annotation = ha1, 
          row_names_side = "left", row_labels = c(1:10, "Total"), show_column_names = FALSE, show_row_dend = FALSE, 
          show_column_dend = FALSE, show_parent_dend_line = FALSE, cluster_columns = FALSE, cluster_rows = FALSE)


pdf('./output_figures/chromosome_map_combined.pdf', width = 12, height = 10)
ComplexHeatmap::draw(ht_list3)
dev.off()

####### Plot Legends

col_fun = colorRamp2(c(-20, 0, 20), c("blue", "white", "red"))
lgd_PET.min = Legend(col_fun = col_fun, title = "PET, Driest Month", direction = "horizontal", legend_width=unit(1.6, 'in'))
ComplexHeatmap::draw(lgd_PET.min)

col_fun1 = colorRamp2(c(-1000, 0, 1000), c("blue", "white", "red"))
lgd_Tmp.Seas = Legend(col_fun = col_fun1, title = "Temperature Seasonality", direction = "horizontal", legend_width=unit(1.6, 'in'))
ComplexHeatmap::draw(lgd_Tmp.Seas)

col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
lgd_Prc.Seas = Legend(col_fun = col_fun, title = "Precipitation Seasonality", direction = "horizontal", legend_width=unit(1.6, 'in'))
ComplexHeatmap::draw(lgd_Prc.Seas)

col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
lgd_topsoil_pH = Legend(c(-3, -1.5, 0, 1.5, 3), col_fun = col_fun, title = "Topsoil pH", direction = "horizontal", legend_width=unit(1.6, 'in'))
ComplexHeatmap::draw(lgd_topsoil_pH)



############################################################################
## Circos plot - Panel B

### Prepping data for Circos plot ###
dat <- read.csv("./chromosome_map/marker_effects_chromosome_map.csv", header = T)
mat <- readRDS("./chromosome_map/transposed_sorghum_genotype_matrix_for_chrommaps.rds") # run chromosome_Complex_Heatmap.R to create this file

# Create Matrix (one chromosome)
comps <- colnames(dat)[2:4] # names of variable columns of interest
filelist <- list()
chrom_num <- 3 # set which chromosome to pull out
for(i in  comps){
  mat_chrom <- mat[which(mat$chrom == chrom_num),]
  dat_chrom <- dat[which(dat$chrom == chrom_num),]
  mark <- mat_chrom[,colnames(mat)!='chrom']
  mark <- as.data.frame(mark)
  mark <- lapply(mark[,1:ncol(mark)], as.numeric)
  mark <- as.data.frame(lapply(mark, `*`, dat_chrom[,i]))
  mark$chrom <- mat_chrom[,colnames(mat) == 'chrom']
  #rownames(mark) <- dat$Markers
  mark$position <- dat_chrom$pos
  mark$rs. <- dat_chrom$rs.
  mark$var <- i
  
  filename <- paste(i, "Marker_Effects_chrom3", sep = "_")
  new <- filename
  filelist[[length(filelist) + 1]] <- new
  assign(filename, mark)
}

# Create Matrix (full genome)
comps <- colnames(dat)[2:4] # names of variable columns of interest
filelist <- list()
for(i in  comps){
  mark <- mat[,colnames(mat)!='chrom']
  mark <- as.data.frame(mark)
  mark <- lapply(mark[,1:ncol(mark)], as.numeric)
  mark <- as.data.frame(lapply(mark, `*`, dat[,i]))
  mark$chrom <- mat[,colnames(mat) == 'chrom']
  #rownames(mark) <- dat$Markers
  mark$position <- dat$pos
  mark$rs. <- dat$rs.
  mark$var <- i
  
  filename <- paste(i, "Marker_Effects_full", sep = "_")
  new <- filename
  filelist[[length(filelist) + 1]] <- new
  assign(filename, mark)
}

# write.csv(prc_seas_Marker_Effects_chrom3, "Prc.Seas_Markereffects_Chrom3.csv")

prc.seas.chrom3.plot <- prc_seas_Marker_Effects_chrom3 %>% rename(IS22604 = IS22604.C1AB3ACXX.7.250137714, PI646283 = PI646283.D1G0BACXX.7.250129785,
                                                                  PI619920 = PI619920.D1FY1ACXX.8.250123619, IS22799 = IS22799.MRG.4.250047511,
                                                                  IS32569 = IS32569.70MUKAAXX.1.250006400, IS11473 = IS11473.D0D0HACXX.5.250046338,
                                                                  IS19859 = IS19859.MRG.4.250047487, PI586201 = PI586201.D1FY1ACXX.8.250123481,
                                                                  IS8916 = IS8916.D0D0HACXX.5.250046328, PI619829 = PI619829.D1FY1ACXX.8.250123604,
                                                                  PI586201 = PI586201.D1FY1ACXX.8.250123481, IS22701 = IS22701.C1AB3ACXX.7.250137718,
                                                                  IS22720 = IS22720.MRG.4.250047510) %>%
  select(position, rs., IS22720, IS22701, IS22799, IS8916, PI619829, PI646283, IS19859, PI586201, IS11473)

#SNPRelate
vcf.fn <- "./sequence_data/genotypes_sorghum_landraces_all.vcf" #download from FigShare link in Github

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "G:/My Drive/Sorghum/R/sequence_data/sorghum_landraces_all_snps.gds", method="copy.num.of.ref")
gds_sum <- snpgdsSummary("G:/My Drive/Sorghum/R/sequence_data/sorghum_landraces_all_snps.gds")
genofile <- snpgdsOpen("G:/My Drive/Sorghum/R/sequence_data/sorghum_landraces_all_snps.gds")
set.seed(1000)
read.gdsn(index.gdsn(genofile, "snp.id"))
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F) #ERROR: SNP on non-Autosomes***
snp.ref <- data.frame(seq_id = read.gdsn(index.gdsn(genofile, "snp.id")), 
                      rs. = read.gdsn(index.gdsn(genofile, "snp.rs.id")))

snp_chrom3 <- snpset$chr3
snps_to_keep <- snp.ref[snp.ref$seq_id %in% snp_chrom3,]

mark.dat <- prc.seas.chrom3.plot[(prc.seas.chrom3.plot$rs. %in% snps_to_keep$rs.),]
mark.dat1 <- mark.dat[order(mark.dat$position),]

col_fun1 = colorRamp2(c(-0.0055, -0.0015, 0, 0.0015, 0.0055), c("#163670", "#5d94cb", "white", "#e06e85", "#701130"))

# Write axis labels
labs<- paste0(seq(from = 0, to = 74, by = 10), "Mb")
ma <- seq(from = 1, to = 1115, by = 159)

png(filename = "./output_figures/circos_chromosome3.png", height = 1500, width = 1541, unit = "px")
circos.par(start.degree = 72, gap.after = 45)
circos.heatmap(mark.dat1[,3:11], col = col_fun1, track.height = 0.75, bg.border = "black", bg.lwd = 1, bg.lty = 1)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    cn = colnames(mark.dat)[3:11]
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                1:n - 0.5, cn, 
                cex = 3, adj = c(0, 0.5), niceFacing = T) # , facing = "bending.inside"
  }
}, bg.border = NA)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) { # labelling for megabases (outer ring)
  circos.genomicAxis(h = 'top', labels.cex = 3, labels = labs, major.at = ma)
})
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    circos.rect(CELL_META$cell.xlim[1] - convert_x(0, "mm"), 0,
                CELL_META$cell.xlim[1] - convert_x(16, "mm"), 3,
                col = "lightyellow", border = NA)
    circos.text(CELL_META$cell.xlim[1] - convert_x(9, "mm"), 1.5,
                "Somalia", cex = 3, facing = "clockwise")
    
    circos.rect(CELL_META$cell.xlim[1] - convert_x(0, "mm"), 3,
                CELL_META$cell.xlim[1] - convert_x(16, "mm"), 6,
                col = "lightpink", border = NA)
    circos.text(CELL_META$cell.xlim[1] - convert_x(9, "mm"), 4.5,
                "Uganda", cex = 3, facing = "clockwise")
    
    circos.rect(CELL_META$cell.xlim[1] - convert_x(0, "mm"), 6,
                CELL_META$cell.xlim[1] - convert_x(16, "mm"), 9,
                col = "lightgreen", border = NA)
    circos.text(CELL_META$cell.xlim[1] - convert_x(9, "mm"), 7.5,
                "Top Scores", cex = 3, facing = "clockwise")
  }
}, bg.border = NA)
circos.clear()

dev.off()


lgd_Prc.Seas = Legend(col_fun = col_fun1, title = "Marker Effect- \nPrecipitation Seasonality", direction = "horizontal", legend_width=unit(1.6, 'in'))
ComplexHeatmap::draw(lgd_Prc.Seas)
