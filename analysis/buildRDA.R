############################################
### Creating the .rda file from raw data ###
############################################

### Load required libraries
library(VariantAnnotation)
library(vcfR)

### Load data
allLoci.vcf <- readVcf("./extData/lamprey_freebayes.targets.filtered.indep_ind.vcf")
allLoci.vcfR <- read.vcfR("./extData/lamprey_freebayes.targets.filtered.indep_ind.vcf")
ageNames <- read.delim("./extData/samplenames_age_specific.csv", sep = ",", header = FALSE)

### Save as .rda file
save(allLoci.vcf, allLoci.vcfR, ageNames, file = "./data/SeaLampreyRapture.rda")

