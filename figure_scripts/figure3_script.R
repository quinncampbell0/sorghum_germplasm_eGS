#### Recreate Figure 3 - Germplasm Exchange Map ####
library(dplyr)
library(maps)
library(geosphere)
library(plotfunctions)

own_accs <- read.csv('./input_data/network_map_own_accessions.csv')
full_dat <- read.csv('./input_data/network_map_point_coordinates.csv')

all_hold <- full_dat  %>% select(c(holding_country, n, hold_lat, hold_long)) %>% group_by(holding_country) %>% mutate(total_n = sum(n)) %>% mutate(n = NULL)

all_hold <- all_hold[!duplicated(all_hold$holding_country),]

map('world', fill = TRUE, col = "lightgray", mar = rep(0, 4))
## points for germplasm repositories
# need data frame with lat/long for repositories
# change size based on number of stored accessions
points(x=full_dat$hold_long, y=full_dat$hold_lat, col="slateblue", cex=3, pch=20)

# Data wrangling for pie chart
own_accs1 <- own_accs[,c(1,3)]

pie_dat <- merge(all_hold, own_accs1, by = 'holding_country')
pie_dat$all_n <- pie_dat$total_n + pie_dat$n
pie_dat$n <- as.numeric(pie_dat$n)
pie_dat$all_n <- as.numeric(pie_dat$all_n)

scaled_n <- log(pie_dat$all_n, base = 2)
#scaled_n1 <- log(pie_dat$all_n, base = 10)
#pie_dat$n[which(pie_dat$n)]

# A function to plot connections
plot_my_connection=function( dep_lon, dep_lat, arr_lon, arr_lat, ...){
  inter <- gcIntermediate(c(dep_lon, dep_lat), c(arr_lon, arr_lat), n=50, addStartEnd=TRUE, breakAtDateLine=F)             
  inter=data.frame(inter)
  diff_of_lon=abs(dep_lon) + abs(arr_lon)
  if(diff_of_lon > 180){
    lines(subset(inter, lon>=0), ...)
    lines(subset(inter, lon<0), ...)
  }else{
    lines(inter, ...)
  }
}


# Compute the connection between points
# as well as change line width based on number of sourced accessions

# Generate all pairs of coordinates
# data should be countries as rows, lat/long as columns
# think about adding 
# all_pairs <- cbind(t(combn(data$long, 2)), t(combn(data$lat, 2))) %>% as.data.frame()
# colnames(all_pairs) <- c("long1","long2","lat1","lat2")

# background map
par(mar=c(0,0,0,0))

all_full_dat <- full_dat
full_dat <- all_full_dat[which(all_full_dat$n > 5),]

# Creating a color ramp
#col_fun <- colorRampPalette(c('darkgray', 'yellow','FFB833', 'orange', "FF7733",'red'))
col_fun <- colorRampPalette(c('grey', '#DAF7A6', 'yellow','#FFC300', 'orange', "#FF5733",'#C70039'))
col_fun_points <- colorRampPalette(c('lightblue', 'blue', 'darkblue'))
col_test <- col_fun(7238)
col_points <- col_fun_points(max(all_hold$total_n))

all_hold$color <- col_points[all_hold$total_n]
#colramp <- colfun(7238)
#colramp <- hcl.colors(7238, palette = 'YlOrRd', rev = T)

width <- log(full_dat$n)/3

pdf(file = 'network_map.pdf', width = 12, height = 8)
map('world',col="black", fill=FALSE, bg="white", lwd=0.05,mar=rep(0,4), ylim=c(-90,90) )
# add every connections:
for(i in 1:nrow(full_dat)){
  plot_my_connection(full_dat$orig_long[i], full_dat$orig_lat[i], full_dat$hold_long[i], full_dat$hold_lat[i], 
                     col=col_test[full_dat$n[i]], lwd=width[i]) # adjust line color to scale by number of accessions
  #print(i)
}
#points(x=all_hold$hold_long, y=all_hold$hold_lat, col=all_hold$color, cex=log(all_hold$total_n)/2, pch=20)
draw.pie(z=cbind(pie_dat$n, pie_dat$total_n), x=pie_dat$hold_long, y=pie_dat$hold_lat, radius= scaled_n/2, col = c("darkblue", "red"), scale = F)
dev.off()

plot.new()
gradientLegend(c(0, 7500), color = col_test, side = 1, length = 0.25, n.seg = 2)
gradientLegend(c(0, 40000), color = col_points, side = 1, length = 0.25, n.seg = 4)
plot.new()
legend(x = "center", legend = c("40,000", "10,000", "1,000", "100"), pch=20, pt.cex = scale_legend/2, col = 'darkblue', 
       y.intersp = 2, x.intersp = 1, title = "# of Accessions Held by Genebank")
plot.new()
legend(x = 'center', legend = c("Sourced outside country", 'Sourced within country'), fill = c('red', 'darkblue'))
draw.pie()
# add points and names of cities
points(x=data$long, y=data$lat, col="slateblue", cex=2, pch=20)
text(rownames(data), x=data$long, y=data$lat,  col="slateblue", cex=1, pos=4)

# Add pie charts if needed to specific institutions
library(mapplots)

draw.pie(z=cbind(pie_dat$n, pie_dat$all_n), x=pie_dat$hold_long, y=pie_dat$hold_lat, radius=12, labels = NA, col = c("darkblue", "red"))

for(i in 1:nrow(pie_dat)){
  # add.pie(z=c(pie_dat$all_n[i], pie_dat$n[i]), x=pie_dat[i,][3], y=pie_dat[i,][2], radius=10, labels = NULL, col = c("darkblue", "red"))
  
}

#add.pie(z=rpois(4,10), x=0.5, y=-0.5, radius=0.3)