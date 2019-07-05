############################################
### Creating the .rda file from raw data ###
############################################

### Load required libraries
library(VariantAnnotation)

### Load data
allLoci <- readVcf("./extData/lamprey_freebayes.targets.filtered.indep_ind.vcf")
ageNames <- read.delim("./extData/samplenames_age_specific.csv", sep = ",", header = FALSE)

### Save as .rda file
save(allLoci, ageNames, file = "./data/SeaLampreyRapture.rda")

