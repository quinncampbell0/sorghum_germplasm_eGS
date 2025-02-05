library(LDheatmap)
require(vcfR)
require(snpStats)

##### LD Heatmap- Panel A
#setwd('G:/My Drive/Sorghum/R/')
#snp <- read.vcfR("./sequence_data/genotypes_sorghum_landraces_all.vcf")
snp <- read.vcfR('C:/Users/Quinn/Desktop/plink/mydata_vcf.vcf') # VCF filtered for missingness and MAF in plink
gt <- snp@gt

snpMat <- t(gt)

#define a function to convert the value into 0,1,2
convertToNumeric <- function(x){
  gdat <- matrix(NA,nrow = nrow(x), ncol = ncol(x))
  for (m in 1:nrow(x)){
    for (n in 1:ncol(x)){
      a <-as.numeric(unlist(strsplit(x[m,n], "|"))[1]) 
      
      b <- as.numeric(unlist(strsplit(x[m,n], "|"))[3])
      gdat[m,n] <- a+b
    }
  }
  rownames(gdat) <- rownames(x)
  colnames(gdat) <- colnames(x)
  return(gdat)
}

####### Extract pairwise distance column
# Extract the SNP position and IDs
snp_positions <- as.data.frame(snp@fix[, c("ID", "CHROM", "POS")]) # Extract ID, Chromosome, and Position
colnames(snp_positions) <- c("ID", "CHROM", "POS")

# Ensure POS is numeric
snp_positions$POS <- as.numeric(snp_positions$POS)

# Calculate distances between consecutive SNPs
snp_positions <- snp_positions[order(snp_positions$CHROM, snp_positions$POS), ] # Sort by CHROM and POS
snp_positions$Distance <- c(NA, diff(snp_positions$POS)) # Compute pairwise distances

#convert to snpMatrix
gdat <- convertToNumeric(snpMat)

#load the snp_id_dist.csv, which contains the SNPs id and distance
# info <- read.csv("snp_id_dist.csv")
snpNames <-snp_positions$ID

#write.csv(snp_positions, 'snp_positions.csv', row.names = F)

snp_positions <- read.csv('snp_positions.csv')



acc <- read.csv('./input_data/envdat_master_revision.csv')

gdat1 <- gdat[acc$gen_id,]

colnames(gdat1) <- snpNames

gdat2 <- as(gdat1,"SnpMatrix")

#write.table(gdat2, 'SnpMatrix_LDHeatmap.txt')

#gdat2 <- read.table('SnpMatrix_LDHeatmap.txt', header = T)

gdat2_mat <- as.matrix(gdat2)

gdat3 <- as(gdat2_mat,"SnpMatrix")

#saveRDS(gdat3, 'SnpMatrix_LDHeatmap.RData')

#gdat3 <- readRDS('SnpMatrix_LDHeatmap.RData')

ld_acc <- ld(gdat3, depth = 10, stats="R.squared", symmetric = T)

ld_acc1 <- as.matrix(ld_acc)

LDheatmap(ld_acc1, snp_positions$POS, add.map = FALSE)

gdat_mini <- gdat3[,1:5000]
ld_mini <- ld(gdat_mini, depth = 10, stats = "R.squared", symmetrical = T)
snp_mini <- snp_positions[colnames(gdat_mini),]
ld_mini <- as.matrix(ld_mini)
#LDheatmap(gdat3,snp_positions$Distance,add.map=FALSE)
LDheatmap(ld_mini,snp_mini$POS,add.map=FALSE)



# Editing PLINK files

# Load required libraries
library(dplyr)

# Load the PLINK .map file
plink_map <- read.table("C:/Users/Quinn/Desktop/plink/myplink.map", header = FALSE, stringsAsFactors = FALSE) # map file output by PLINK
colnames(plink_map) <- c("CHR", "SNP", "Genetic_Distance", "POS")

# Load the genetic map file (assumes columns: CHR, POS, cM)
genetic_map <- read.table("./sequence_data/sorghum_GM_map.txt", header = TRUE, stringsAsFactors = FALSE) # read in from figshare linked in Github

# Ensure the columns match your data (adjust if necessary)
colnames(genetic_map) <- c("POS", "CHR", "cM")

# Interpolate genetic distances for SNPs in the .map file
plink_map <- plink_map %>%
  mutate(
    Genetic_Distance = approx(
      x = genetic_map$POS,        # Genetic map positions
      y = genetic_map$cM,         # Genetic map distances
      xout = POS,                 # SNP positions in .map file
      method = "linear",          # Linear interpolation
      rule = 2                    # Use nearest boundary for out-of-range values
    )$y
  )

# Save the updated .map file
write.table(
  plink_map,
  file = "updated_plink.map",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# LD heatmap from PLINK LD matrix

# Load LD matrix (adjust for your file format)
ld_matrix <- as.matrix(read.table("C:/Users/Quinn/Desktop/plink/ld_matrix_600.ld", header = FALSE))

#ld_1 <- colnames(ld_matrix)
#ld_2 <- sub('.', '', ld_1)

#ld_matrix1 <- rbind(ld_2, ld_matrix)

#colnames(ld_matrix1) <- NULL
#row.names(ld_matrix1) <- NULL

#ld_matrix2 <- as.numeric(ld_matrix1)

colfunc <- colorRampPalette(c('black','#E4E4E4'))
colfunc(20)

# Plot heatmap - full genome
tiff('LD_Heatmap_1.tiff', width = 1500, height = 1500)
LDheatmap(ld_matrix, add.map = FALSE, color = colfunc(20))  # Use physical distances or indices
dev.off()

### Plot by chromosome

for(chr in 1:10){
  mat_file <- paste0('C:/Users/Quinn/Desktop/plink/ld_chr',chr,'.ld')
  mat <-  as.matrix(read.table(mat_file))
  
  tiff_file <- paste0('LDHeatmap_chr',chr,'.tiff')
  
  tiff(tiff_file, width = 800, height = 800)
  LDheatmap(mat, add.map = FALSE, color = colfunc(20))  # Use physical distances or indices
  dev.off()
}
bp_matrix = as.matrix(read.table('C:/Users/Quinn/Desktop/plink/thinned_chr3.ld'))

LDheatmap(bp_matrix, add.map = FALSE, color = colfunc(20))  # Use physical distances or indices

#####################
# LD Decay Curves - Panel B

library(ggplot2)
library(ggpubr)

# Read in the LD data for chromosome 1
ld_chr1 <- read.table("C:/Users/Quinn/Desktop/plink/ld_thin_chr1.ld", header = TRUE)

ld_chr1$Distance <- abs(ld_chr1$BP_A - ld_chr1$BP_B)

# Plot for chromosome 1
ggplot(ld_chr1, aes(x = Distance, y = R2)) +
  geom_point(alpha = 0.5, color = 'lightgrey') +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  #scale_x_log10() +
  labs(title = "LD Decay Curve for Chromosome 1",
       x = "Distance (bp)", y = "r²") +
  theme_minimal()

# Loop through all chromosomes

file_names <- paste0("C:/Users/Quinn/Desktop/plink/ld_thin_chr", seq(from = 1, to = 10, by = 1), ".ld") # LD matrix for each of 10 chromosomes, output by PLINK
obj_names <- paste0("ld_chr", seq(1,10,1))

ld_list <- list()
plot_list <- list()

for(i in 1:10){
  ld_list[[i]] <- read.table(file_names[i], header = T)
  
  ld_list[[i]]$Distance <- abs(ld_list[[i]]$BP_A - ld_list[[i]]$BP_B)/1000000
  
  ld_data <- ld_list[[i]] %>%
    mutate(
      distance_bin = cut(Distance, breaks = seq(0, max(Distance), by = 0.001), include.lowest = TRUE),
      midpoint = as.numeric(sub("\\((.+),(.+)\\]", "\\1", distance_bin)) +
        as.numeric(sub("\\((.+),(.+)\\]", "\\2", distance_bin)) / 2
    )
  
  # Summarize data by bin
  ld_summary <- ld_data %>%
    group_by(distance_bin) %>%
    summarize(mean_r2 = mean(R2, na.rm = TRUE), midpoint = unique(midpoint))
  
  ld_summary$midpoint[1] <- 0.0005
  ld_summary <- ld_summary[-nrow(ld_summary),]
  
  nls_fit <- nls(mean_r2 ~ 1 / (1 + midpoint / d), data = ld_summary, start = list(d = 1000))
  curve_data <- data.frame(midpoint = seq(0, max(ld_summary$midpoint), length.out = 100))
  curve_data$fit_r2 <- predict(nls_fit, newdata = curve_data)
  
  plot_list[[i]] <-
    ggplot(ld_list[[i]], aes(x = Distance, y = R2)) +
    geom_point(alpha = 0.5, color = 'lightgrey') +
    geom_line(data = curve_data, aes(x = midpoint, y = fit_r2), color = 'blue') +
    labs(title = paste0("Chromosome ",i),
         x = "Distance (Mb)", y = "r²") +
    theme_minimal()+
    xlim(0, round(max(ld_list[[i]]$Distance)))
}

png(file = './revision_scripts/ld_decay_chrom.png', width = 2000, height = 1000, units='px')
ggarrange(plotlist=plot_list, ncol = 5, nrow = 2)
dev.off()

#### Extra code

ld_data <- ld_list[[8]] %>%
  mutate(
    distance_bin = cut(Distance, breaks = seq(0, max(Distance), by = 0.001), include.lowest = TRUE),
    midpoint = as.numeric(sub("\\((.+),(.+)\\]", "\\1", distance_bin)) +
      as.numeric(sub("\\((.+),(.+)\\]", "\\2", distance_bin)) / 2
  )

# Summarize data by bin
ld_summary <- ld_data %>%
  group_by(distance_bin) %>%
  summarize(mean_r2 = mean(R2, na.rm = TRUE), midpoint = unique(midpoint))

ld_summary$midpoint[1] <- 0.0005
ld_summary <- ld_summary[-nrow(ld_summary),]

nls_fit <- nls(mean_r2 ~ 1 / (1 + midpoint / d), data = ld_summary, start = list(d = 1000))
curve_data <- data.frame(midpoint = seq(0, max(ld_summary$midpoint), length.out = 100))
curve_data$fit_r2 <- predict(nls_fit, newdata = curve_data)

ggplot() +
  geom_point(data = ld_list[[8]], aes(x = Distance, y = R2), color = 'lightgrey') +
  geom_line(data = curve_data, aes(x = midpoint, y = fit_r2), color = 'blue') +
  theme_minimal() +
  labs(x = "Distance (bp)", y = "R2", title = "LD Decay Curve") +
  xlim(0, round(max(ld_list[[8]]$Distance)))

# Read in the LD data for chromosome 2
ld_chr2 <- read.table("C:/Users/Quinn/Desktop/plink/ld_mapthin_chr2.ld", header = TRUE)

ld_chr2$Distance <- abs(ld_chr2$BP_A-ld_chr2$BP_B)/100000

# Plot for chromosome 1
ggplot(ld_chr2, aes(x = Distance, y = R2)) +
  geom_point(alpha = 0.5, color = 'lightgrey') +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  #scale_x_log10() +
  labs(title = "LD Decay Curve for Chromosome 1",
       x = "Distance (Mb)", y = "r²") +
  theme_minimal()


####################################
# Exponential Decay curve

# Parameters for the exponential decay
A <- 1      # Initial value
b <- 0.2     # Decay constant
t <- seq(0, 20, by = 0.1)  # Time points

# Exponential decay function
y <- A * exp(-b * t)

# Create a data frame
data <- data.frame(t = t, y = y)

ggplot(data, aes(x = t, y = y)) +
  geom_point(color = "blue", size = 2) +  # Add points for data
  geom_smooth(method = "nls",
              formula = y ~ A * exp(-b * t),  # Exponential decay model
              method.args = list(start = c(A = 10, b = 0.2)),  # Initial values for A and b
              color = "darkred", size = 1.5, se = FALSE) +  # Fit line (dark red)
  labs(title = "Exponential Decay Fitting using geom_smooth()",
       subtitle = "Using nls (Nonlinear Least Squares) to fit an exponential model",
       x = "Time",
       y = "Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "italic"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

