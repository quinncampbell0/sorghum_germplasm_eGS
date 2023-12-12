
setwd("G:/My Drive/Sorghum/R/")

library(circlize)
library(ggplot2)
library(SNPRelate)
library("FactoMineR")
library("factoextra")
library(dplyr)
library(ComplexHeatmap)

### Prepping data for Circos plot ###
dat <- read.csv("./chromosome_map/marker_effects_chromosome_map.csv", header = T)
mat <- readRDS("./chromosome_map/transposed_sorghum_genotype_matrix_for_chrommaps.rds")

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

#write.csv(prc_seas_Marker_Effects_chrom3, "Prc.Seas_Markereffects_Chrom3.csv")

prc.seas.chrom3.plot <- prc_seas_Marker_Effects_chrom3 %>% rename(IS22604 = IS22604.C1AB3ACXX.7.250137714, PI646283 = PI646283.D1G0BACXX.7.250129785,
                                                             PI619920 = PI619920.D1FY1ACXX.8.250123619, IS22799 = IS22799.MRG.4.250047511,
                                                             IS32569 = IS32569.70MUKAAXX.1.250006400, IS11473 = IS11473.D0D0HACXX.5.250046338,
                                                             IS19859 = IS19859.MRG.4.250047487, PI586201 = PI586201.D1FY1ACXX.8.250123481,
                                                             IS8916 = IS8916.D0D0HACXX.5.250046328, PI619829 = PI619829.D1FY1ACXX.8.250123604,
                                                             PI586201 = PI586201.D1FY1ACXX.8.250123481, IS22701 = IS22701.C1AB3ACXX.7.250137718,
                                                             IS22720 = IS22720.MRG.4.250047510) %>%
                                                        select(position, rs., IS22720, IS22701, IS22799, IS8916, PI619829, PI646283, IS19859, PI586201, IS11473)

#Phylos n1378
#SNPRelate
vcf.fn <- "./sequence_data/genotypes_sorghum_landraces_all.vcf"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "./sequence_data/sorghum_landraces_all_snps.gds", method="copy.num.of.ref")
gds_sum <- snpgdsSummary("./sequence_data/sorghum_landraces_all_snps.gds")
genofile <- snpgdsOpen("./sequence_data/sorghum_landraces_all_snps.gds")
set.seed(1000)
read.gdsn(index.gdsn(genofile, "snp.id"))
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F) #ERROR: SNP on non-Autosomes***
snp.ref <- data.frame(seq_id = read.gdsn(index.gdsn(genofile, "snp.id")), 
                      rs. = read.gdsn(index.gdsn(genofile, "snp.rs.id")))

snp_chrom3 <- snpset$chr3
snps_to_keep <- snp.ref[snp.ref$seq_id %in% snp_chrom3,]

mark.dat <- prc.seas.chrom3.plot[(prc.seas.chrom3.plot$rs. %in% snps_to_keep$rs.),]
mark.dat1 <- mark.dat[order(mark.dat$position),]

col_fun1 = colorRamp2(c(-0.0055, -0.0015, 0, 0.0015, 0.0055), c("#701130", "#e06e85", "white",
"#5d94cb", "#163670"))

# Write axis labels
labs<- paste0(seq(from = 0, to = 74, by = 10), "Mb")
ma <- seq(from = 1, to = 1115, by = 159)


circos.par(start.degree = 72, gap.after = 40)
circos.heatmap(mark.dat1[,3:11], col = col_fun1, track.height = 0.4, bg.border = "black", bg.lwd = 1, bg.lty = 1)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    cn = colnames(mark.dat)[3:11]
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                1:n - 0.5, cn, 
                cex = 0.75, adj = c(0, 0.5), facing = "bending.inside")
  }
}, bg.border = NA)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.genomicAxis(h = 'top', labels.cex = 0.75, labels = labs, major.at = ma)
})
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    circos.rect(CELL_META$cell.xlim[1] - convert_x(1, "mm"), 0,
                CELL_META$cell.xlim[1] - convert_x(5, "mm"), 3,
                col = "orange", border = NA)
    circos.text(CELL_META$cell.xlim[1] - convert_x(3, "mm"), 1.5,
                "Somalia", cex = 0.75, facing = "clockwise")
    
    circos.rect(CELL_META$cell.xlim[1] - convert_x(1, "mm"), 3,
                CELL_META$cell.xlim[1] - convert_x(5, "mm"), 6,
                col = "pink", border = NA)
    circos.text(CELL_META$cell.xlim[1] - convert_x(3, "mm"), 4.5,
                "Uganda", cex = 0.75, facing = "clockwise")
    
    circos.rect(CELL_META$cell.xlim[1] - convert_x(1, "mm"), 6,
                CELL_META$cell.xlim[1] - convert_x(5, "mm"), 9,
                col = "lightgreen", border = NA)
    circos.text(CELL_META$cell.xlim[1] - convert_x(3, "mm"), 7.5,
                "Top Ranked", cex = 0.75, facing = "clockwise")
  }
}, bg.border = NA)
circos.clear()



lgd_Prc.Seas = Legend(col_fun = col_fun1, title = "Marker Effects", direction = "horizontal", legend_width=unit(1.6, 'in'))
ComplexHeatmap::draw(lgd_Prc.Seas)
