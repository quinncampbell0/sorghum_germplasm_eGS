#setwd('G:/My Drive/Sorghum/R/')

ld_list <- seq(0.1, 0.95, by = 0.05)
prune_set <- format(ld_list, nsmall = 2)

result_list <- lapply(prune_set, function(x){
  file <- readRDS(paste0('./xval_ldprune/xval_rrblup_LD_',x,'.RData'))
  file$xval.result$LD_prune <- x
  return(file$xval.result)
})

#full_xval <- read.csv('./revision_scripts/xval_rrblup_kfold10_allmarkers.csv')

#full_xval_tr <- full_xval %>%
#  filter(trait %in% c("topsoil_pH", "Tmp.Seas", "Prc.Seas")) %>% 
#  mutate(LD_prune = 0)

ecophys_full_xval <- readRDS('./xval_ldprune/xval_rrblup_kfold_10.RData')

ecophys_fx <- ecophys_full_xval$xval.result

ecophys_fx_tr <- ecophys_fx %>%
  filter(trait %in% c('pet_min', 'ann_heat_index', 'coldindex.mean')) %>% 
  mutate(LD_prune = 1)

xval_ld <- do.call(rbind, result_list)

#xval_ld <- xval_ld %>%
#  filter(trait %in% c('pet_min'))

xval_ld <- rbind(xval_ld, ecophys_fx_tr)

xval_ld$r.mean <- as.numeric(xval_ld$r.mean)
xval_ld$r.sd <- as.numeric(xval_ld$r.sd)
xval_ld$LD_prune <- as.numeric(xval_ld$LD_prune)

ggplot(xval_ld, aes(x = LD_prune, y = r.mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = r.mean-r.sd, ymax = r.mean+r.sd))+
  facet_wrap(~trait, scales = 'free_y')

ggplot(xval_ld, aes(x = LD_prune, y = r.mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = r.mean-r.sd, ymax = r.mean+r.sd), fill = 'blue', alpha = 0.2)+
  facet_wrap(~trait, scales = 'free_y') +
  theme_minimal()

# Only PET min

xval_ld_pet <- xval_ld %>%
  filter(trait %in% c('pet_min'))

# Add in SNP number

prune_set <- seq(0.1, 0.95, by = 0.05)
prune <- format(prune_set, nsmall = 2)

snpnums <- lapply(prune, function(x){
  
  filename <- paste0('./xval_ldprune/prune_r2_',x,'.prune.in')
  
  snps <- read.table(filename)
  
  return(nrow(snps))
})

snpnums <- unlist(snpnums)

snpnums1 <- c(snpnums, 404627)

xval_ld_pet$SNPs <- snpnums1

xval_ld_pet$LD_prune <- format(xval_ld_pet$LD_prune , nsmall = 2)

ggplot(xval_ld_pet, aes(x = LD_prune, y = r.mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = r.mean-r.sd, ymax = r.mean+r.sd), fill = 'blue', alpha = 0.2) +
  theme_minimal() + 
  labs(
    #title = "Line Plot with Standard Deviation Shaded",
    y = "r",
    x = "LD Threshold"
  )

ggplot(xval_ld_pet, aes(x = SNPs, y = r.mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), fill = 'blue', alpha = 0.2) +
  geom_point(color = 'black') + 
  geom_text(aes(label = LD_prune), vjust = -0.5, size = 3) +  # Add labels for each point
  theme_minimal() +
  scale_x_continuous(labels = scales::comma) +  # Use full numbers instead of scientific notation
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Angle x-axis labels at 45 degrees
  ) +
  labs(
    y = "r",
    x = "# of SNPs"
  )


ggplot(xval_ld_pet, aes(x = SNPs, y = r.mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = r.mean - r.sd, ymax = r.mean + r.sd), fill = 'blue', alpha = 0.2) +
  geom_point(size = 1) +
  geom_text_repel(aes(label = LD_prune), size = 3, box.padding = 0.1) +  # Automatically adjusts label positions
  theme_minimal() +
  scale_x_continuous(labels = scales::comma) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(10, 50, 10, 10)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    y = "r",
    x = "SNPs"
  )
