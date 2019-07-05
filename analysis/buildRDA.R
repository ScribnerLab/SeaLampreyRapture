############################################
### Creating the .rda file from raw data ###
############################################

### Load required libraries
library(VariantAnnotation)
library(vcfR)

### Load data
allLoci.vcf <- readVcf("./extData/lamprey_freebayes.targets.filtered.indep_ind.vcf")
allLoci.vcfR <- read.vcfR("./extData/lamprey_freebayes.targets.filtered.indep_ind.vcf")
agePops <- read.delim("./extData/indNames.age.csv", sep = ",", header = FALSE)
indPops <- read.delim("./extData/indNames.csv", sep = ",", header = FALSE)

### Save as .rda file
save(allLoci.vcf, allLoci.vcfR, agePops, indPops, file = "./data/SeaLampreyRapture.rda")

