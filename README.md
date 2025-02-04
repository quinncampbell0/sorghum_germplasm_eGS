# sorghum_germplasm_eGS
Working scripts for manuscript "Climate resilience conserved in global germplasm repositories: Picking the most promising parents for agile plant breeding"

Genotype data in the format used in the study can be found at:
https://figshare.com/s/da4ad3fe9f749f2135b6

Original sequence data downloaded from:
https://datadryad.org/stash/dataset/doi:10.5061/dryad.jc3ht

### **1. System Requirements**
All scripts were developed and run in R version 4.3.1.
R version 4.3.1 is compatible with Windows 7 and later versions. For macOS users, R 4.3.1 supports both Intel and ARM architectures. Ensure that you download the appropriate installer for your system.
R 4.3.1 requires a 64-bit operating system. Maximum disk space required is around 15 GB to run all scripts in this github. For optimal performance, it's recommended to have at least 256 MB of RAM available when running R.

### **2. Installation**
Previous versions of R can be installed at: https://cran.r-project.org/bin/windows/base/old/
All file paths in the provided scripts are relative to the base folder where github folders are downloaded. Set the working directory to the location where the sorghum_germplasm_eGS folders are downloaded.

## **3. Demo**
Download Github Repository into folder, i.e. C:/sorghum_EGS
To run environmental genomic selection and calculate the Genomic Adaptive Capacity Score in R:

`setwd("C:/sorghum_EGS/")`

`source("./score_calculation_scripts/GS_rrblup_full.R")`

`source("./score_calculation_scripts/score_genomicselection_rrblup.R")`

Runtime will depend on system used. EGS and cross-validation scripts may require up to 1 week of processing time on a standard laptop or desktop. 
*Please note much of the data formatting in these scripts is hard-coded, so replacing the environmental data or genotype files with other datasets will require some re-coding*


### The directories contained are described below:

  chromosome_map contains scripts for creating the chromosome heatmap and circos plot contained in Figure 4 of the manuscript, as well as the data file containing the marker effects for the sorghum accessions in the study

  cluster_analysis contains scripts for performing the cluster analyses described in the manuscript using genetic, climatic, spatial, and taxonomic data

  cross_validation contains scripts and results for performing cross-validation on the five genomic selection models used in the manuscript

  figure_scripts contains the R scripts for each of the 3 main text figures and 7 supplemental figures displayed in the manuscript

  geospatial contains rasters for bioclimatic variables, temperature and precipitation variables, borders and cropland extent data stored as TIF files

  input_data contains environmental and soil data, genebank data, and chromosomal effects data used as input data in the scripts contained in this github

  score_calculation_scripts contains scripts to calculated the Climate Resilience and Genomic Adaptive Capacity scores described in the manuscript, as well as the script used for running environmental genomic selection using the RR-BLUP model


