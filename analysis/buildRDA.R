############################################
### Creating the .rda file from raw data ###
############################################

### Load required libraries
library(VariantAnnotation)
library(vcfR)

### Load data
allLoci.vcf <- readVcf("./extData/lamprey_freebayes.targets.filtered.indep_ind.vcf")

allLoci.vcfR <- read.vcfR("./extData/lamprey_freebayes.targets.filtered.indep_ind.vcf")

agePops <-  read.delim("./extData/indNames.age.csv",
                       sep = ",",
                       header = FALSE)

indPops <- read.delim("./extData/indNames.csv", sep = ",", header = FALSE)

appendixLoci <- read.delim("./extData/Loci_Appendix_112918.txt",
                           sep = "\t",
                           header = TRUE)

onTarget_readCount <- read.delim("./extData/R_OT.txt", sep = "\t" , header = FALSE)
colnames(onTarget_readCount) <- c("name", "T", "ReadCount")

targets.chrpos <- read.table("./extData/Targets.chrpos",
                             sep = "\t" ,
                             header = TRUE)


### Save as .rda file
save(
  allLoci.vcf,
  allLoci.vcfR,
  agePops,
  indPops,
  appendixLoci,
  onTarget_readCount,
  targets.chrpos,
  file = "./data/SeaLampreyRapture.rda"
)

