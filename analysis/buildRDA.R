############################################
### Creating the .rda file from raw data ###
############################################

### Load required libraries
library(vcfR)

### Load data
allLoci <- read.vcfR("G:/My Drive/Side projects/Lamprey/Lamprey_Vcfs_InitialFiltering_SRS_101018/lamprey_freebayes.targets.filtered.indep_ind.vcf")

### Save as .rda file
save(allLoci, file = "./data/WFUrbanAdaptation.rda")

