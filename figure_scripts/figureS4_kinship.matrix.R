library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
kmat <- read.csv('./input_data/GAPIT.Genotype.Kin_Zhang.csv', header = F) # find in figshare repository linked in github

kmat[1:10,1:10]

rownames(kmat) <- kmat[,1]

colnames(kmat) <- rownames(kmat)

envdat <- read.csv('./input_data/envdat_master_revision.csv')

kmat1 <- kmat[envdat$gen_id,envdat$gen_id]

kmat1 <- as.matrix(kmat1)

### Sort by country, region

# Specify the desired order for the categories
reg_order <- c("West Africa", "Central Africa", "South Africa", 'East Africa', 'West Asia', 'South Asia', 'East Asia')

# Order the rows of the data frame based on the desired category order in column 1
envdat2 <- envdat %>% arrange(country)
ord_envdat <- envdat2[order(factor(envdat2$region, levels = reg_order)), ]


kmat2 <- kmat1[ord_envdat$gen_id,ord_envdat$gen_id]


### Create labels coordinates
c_ord <- unique(ord_envdat$country)

out_vec <- vector()
for(i in c_ord){
  vec <- which(ord_envdat$country == i)
  pos <- vec[(length(vec) + 1L) %/% 2L]
  out_vec[i] <- pos
}

### Custom colors 

# Define the number of colors you want
nb.cols <- 51
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)
names(mycolors) <- sample(unique(ord_envdat$country))

# Custom region colors
rg.cols <- 7
regcols <- colorRampPalette(brewer.pal(7, "Paired"))(rg.cols)
names(regcols) <- sample(unique(ord_envdat$region))

#### Create annotations
ca = columnAnnotation(Region = ord_envdat$region,
                      col = list(Region = regcols), 
                      #labels_gp = gpar(angle = 45),
                      show_legend = FALSE, show_annotation_name = FALSE)

ca1 = columnAnnotation(Country = ord_envdat$country,
                       col = list(Country = mycolors), show_legend = FALSE, 
                       foo = anno_mark(at = out_vec, labels = c_ord, side = 'bottom', labels_rot = 90, 
                                       labels_gp = gpar(fontsize =5)), 
                       show_annotation_name = FALSE)

ra = rowAnnotation(Region = ord_envdat$region,
                      col = list(Region = regcols), 
                      #labels_gp = gpar(angle = 45),
                      show_legend = FALSE, show_annotation_name = FALSE)

ra1 = rowAnnotation(Country = ord_envdat$country, 
                    col = list(Country = mycolors),
                    foo = anno_mark(at = out_vec, labels = c_ord, side = 'right', labels_rot = 0, 
                                       labels_gp = gpar(fontsize = 5)),
                    
                    show_legend = FALSE, show_annotation_name = FALSE)



png("./output_figures/Kinship_heatmap.png", width = 2300, height = 2000, res = 300)

Heatmap(
  kmat2,
  name = "Relatedness",                  # Legend title
  cluster_rows = FALSE,                  # Cluster rows
  cluster_columns = FALSE,               # Cluster columns
  show_row_names = FALSE,               # Hide row names
  show_column_names = FALSE,            # Hide column names
  top_annotation = ca,
  left_annotation = ra,
  bottom_annotation = ca1,
  right_annotation = ra1,
  heatmap_legend_param = list(direction = 'horizontal'),
  col = colorRampPalette(c("#f0f5f9","#0000FF"))(100) # Custom colors
)

dev.off()

