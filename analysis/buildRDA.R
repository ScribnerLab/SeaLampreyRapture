############################################
### Creating the .rda file from raw data ###
############################################

### Load required libraries
library(VariantAnnotation)
library(vcfR)

### Load data
## Observed genetic data
allLoci.vcf <- readVcf("./extData/lamprey_freebayes.targetloci.filtered.subsampled.vcf")

allLoci.vcfR <- read.vcfR("./extData/lamprey_freebayes.targetloci.filtered.subsampled.vcf")

indPops <-  read.delim("./extData/indNames.csv",
                       sep = ",",
                       header = TRUE)

appendixLoci <- read.delim("./extData/Loci_Appendix_112918.txt",
                           sep = "\t",
                           header = TRUE)

onTarget_readCount <- read.delim("./extData/R_OT.txt", sep = "\t" , header = FALSE)
colnames(onTarget_readCount) <- c("name", "T", "ReadCount")

targets.chrpos <- read.table("./extData/Targets.chrpos",
                             sep = "\t" ,
                             header = TRUE)

targetDensity <- read.table("./extData/Lamprey_TargetDensity.txt",
                            sep = "\t",
                            header = TRUE)

## Pedigree analysis data
best.config.dck.age.1 <- read.table("./extData/best.config.dck.age.1.txt",
                                    sep = "\t",
                                    stringsAsFactors = F,
                                    header = TRUE)

best.config.dck.age.2 <- read.table("./extData/best.config.dck.age.2.txt",
                                    sep = "\t",
                                    stringsAsFactors = F,
                                    header = TRUE)

best.config.dck.age.3 <- read.table("./extData/best.config.dck.age.3.txt",
                                    sep = "\t",
                                    stringsAsFactors = F,
                                    header = TRUE)

best.config.scr.age.2 <- read.table("./extData/best.config.scr.age.2.txt",
                                    sep = "\t",
                                    stringsAsFactors = F,
                                    header = TRUE)

best.config.scr.age.3 <- read.table("./extData/best.config.scr.age.3.txt",
                                    sep = "\t",
                                    stringsAsFactors = F,
                                    header = TRUE)

database_v2 <- read.table("./extData/database_v2.txt",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = TRUE)

alf.dist <- read.table("./extData/alf.dist.txt",
                       sep = "\t",
                       stringsAsFactors = F,
                       header = TRUE)

### Save as .rda file
save(
  allLoci.vcf,
  allLoci.vcfR,
  indPops,
  appendixLoci,
  onTarget_readCount,
  targets.chrpos,
  targetDensity,
  best.config.dck.age.1,
  best.config.dck.age.2,
  best.config.dck.age.3,
  best.config.scr.age.2,
  best.config.scr.age.3,
  database_v2,
  alf.dist,
  file = "./data/SeaLampreyRapture.rda"
)

