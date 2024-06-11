# sorghum_germplasm_eGS
Working scripts for manuscript "Climate resilience conserved in global germplasm repositories"

Genotype data in the format used in the study can be found at:
https://figshare.com/s/da4ad3fe9f749f2135b6

Original sequence data downloaded from:
https://datadryad.org/stash/dataset/doi:10.5061/dryad.jc3ht

The directories contained are described below:

chromosome_map contains scripts for creating the chromosome heatmap and circos plot contained in Figure 4 of the manuscript, as well as the data file containing the marker effects for the sorghum accessions in the study

cluster_analysis contains scripts for performing the cluster analyses described in the manuscript using genetic, climatic, spatial, and taxonomic data

cross_validation contains scripts and results for performing cross-validation on the five genomic selection models used in the manuscript

figure_scripts contains the R scripts for each of the 3 main text figures and 7 supplemental figures displayed in the manuscript

geospatial contains rasters for bioclimatic variables, temperature and precipitation variables, borders and cropland extent data stored as TIF files

input_data contains environmental and soil data, genebank data, and chromosomal effects data used as input data in the scripts contained in this github

score_calculation_scripts contains scripts to calculated the Climate Resilience and Genomic Adaptive Capacity scores described in the manuscript, as well as the script used for running environmental genomic selection using the RR-BLUP model
