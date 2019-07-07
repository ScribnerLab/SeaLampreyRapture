################################################################
######  Summary statistics for sea lamprey RAPTURE loci   ######
################################################################

##### Load data and libraries #####
library(adegenet)
library(vcfR)
library(hierfstat)
library(tidyverse)
library(SeaLampreyRapture)
data(SeaLampreyRapture)

##### Convert data into required formats #####
dat.genind <- vcfR2genind(allLoci.vcfR, sep = "[|/]") ## Create genind object

### Assign pop names to genind file
dat.genind@pop <- indPops$V1

### Create hierfstat object
dat.hierfstat.pops <- genind2hierfstat(dat.genind)

###### Analyze genetic diversity data #####
### Assemble genetic diversity tibble, including subpopulation assignment
geneticDiversity <- basic.stats(dat.hierfstat.pops)
MAF.DCJ <- as.data.frame(minorAllele(dat.genind[pop="DCJ"]))
MAF.SC <- as.data.frame(minorAllele(dat.genind[pop="SC"]))
MAF.BR <- as.data.frame(minorAllele(dat.genind[pop="BR"]))
MAF.CARP <- as.data.frame(minorAllele(dat.genind[pop="CARP"]))
MAF.SM <- as.data.frame(minorAllele(dat.genind[pop="SM"]))
MAF.all <- as.data.frame(minorAllele(dat.genind))

### Overall locus-specific measures
geneticDiversity.loci <- cbind(
  MAF.all,
  geneticDiversity$perloc$Ho,
  geneticDiversity$perloc$Hs,
  geneticDiversity$perloc$Fis
)

names(geneticDiversity.loci) <- c("MAF", "Ho", "He", "Fis")

#write.csv(geneticDiversity.loci, "geneticDiversity.loci.csv")

## By population locus-specific measures
Locus <- names(dat.genind@all.names) ## Names of loci
Hobs <- as.data.frame(geneticDiversity$Ho) ## Obs heterozygosity per locus
Hobs$Locus <- Locus
Hexp <- as.data.frame(geneticDiversity$Hs) ## Exp heterozygosity per locus
Hexp$Locus <- Locus
Fis <- as.data.frame(geneticDiversity$Fis) ## Fis per locus
Fis$Locus <- Locus

names(MAF.DCJ) <- "DCJ"
names(MAF.SC) <- "SC"
names(MAF.BR) <- "BR"
names(MAF.CARP) <- "CARP"
names(MAF.SM) <- "SM"

MAF <- cbind(MAF.DCJ,
             MAF.SC$SC,
             MAF.BR$BR,
             MAF.CARP$CARP,
             MAF.SM$SM)
names(MAF) <- c("DCJ", "SC", "BR", "CARP", "SM")
MAF$Locus <- Locus

## Assemble measure and population specific dataframes into a dataframe
popLocusGD.df <-  cbind(Hobs, Hexp)

## Assemble measure and population specific dataframes into a tibble
Hobs.tib <- Hobs %>%
  gather("Subpop", "Hobs", 1:5)

Hexp.tib <- Hexp %>%
  gather("Subpop", "Hexp", 1:5)

Fis.tib <- Fis %>%
  gather("Subpop", "Fis", 1:5)

MAF.tib <- MAF %>%
  gather("Subpop", "MAF", 1:5)

MAF.tib %>%
  mutate(Subpop = as.factor(Subpop)) %>%
  group_by(Subpop) %>%
  summarize(Prop. = 1 - (sum(MAF <= 0.05)/n()), N = n() - sum(MAF <= 0.05), total.N = sum(MAF > 0))

## Assumble overall tibble
gd.tib <- Hobs.tib %>%
  left_join(Hexp.tib) %>%
  left_join(Fis.tib) %>%
  left_join(MAF.tib)

#### Pairwise Fsts for all unrelated lamprey ####
pairwise.fsts <- pairwise.fst(dat.genind)

