################################################################
######  Summary statistics for sea lamprey RAPTURE loci   ######
################################################################

##### Load data and libraries #####
library(adegenet)
library(vcfR)
library(hierfstat)
library(tidyverse)
library(gridExtra)
library(ggthemes)
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
geneticDiversity <- basic.stats(dat.hierfstat.pops) ## Calculate diversity measures per locus
MAF.DCJ <- as.data.frame(minorAllele(dat.genind[pop="DCJ"])) ## Calculate minor allele frequencies per locus for each subpopulation
MAF.SC <- as.data.frame(minorAllele(dat.genind[pop="SC"]))
MAF.BR <- as.data.frame(minorAllele(dat.genind[pop="BR"]))
MAF.CARP <- as.data.frame(minorAllele(dat.genind[pop="CARP"]))
MAF.SM <- as.data.frame(minorAllele(dat.genind[pop="SM"]))
MAF.all <- as.data.frame(minorAllele(dat.genind))

## Locus-specific measures
geneticDiversity.loci <- cbind(
  MAF.all,
  geneticDiversity$perloc$Ho,
  geneticDiversity$perloc$Hs,
  geneticDiversity$perloc$Fis
)

names(geneticDiversity.loci) <- c("MAF", "Ho", "He", "Fis")

#write.csv(geneticDiversity.loci, "geneticDiversity.loci.csv")

## Measure and population specific dataframes
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
             MAF.SM$SM) ## MAF per locus for each subpop
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

PAs <- MAF.tib %>%
  mutate(Subpop = as.factor(Subpop)) %>%
  filter(MAF > 0) %>%
  group_by(Locus) %>%
  summarize(Ntemp = n_distinct(Subpop)) %>%
  filter(Ntemp == 1)

temp2 <- MAF.tib %>%
  filter(Locus %in% PAs$Locus) %>%
  filter(MAF > 0) %>%
  group_by(Subpop) %>%
  summarize(PAs = n())

downsampled_pairwise <- pairwise.fst(genind.downsampled)

## Assumble overall tibble
gd.tib <- Hobs.tib %>%
  left_join(Hexp.tib) %>%
  left_join(Fis.tib) %>%
  left_join(MAF.tib)

# ### Locus-specific summary data
# Hobs.loci <- gd.tib %>%
#   group_by(Locus) %>%
#   summarize(mean.Hobs = mean(Hobs))
#
#
# temp <- gd.tib %>%
#   group_by(Locus) %>%
#   summarize(mean.MAF = mean(MAF))






### Plot data
## Observed heterozygosity
Hobs.figure <- gd.tib %>%   ## Histogram of Hobs
  ggplot(aes(x = Hobs)) +
  geom_histogram(breaks = seq(0,1.0, by = 0.02)) +
  xlab("Observed heterozygosity") +
  ylab("Count") +
  theme_bw(base_size = 22) +
  NULL

Hobs.facetted <- Hobs.figure + facet_wrap( ~ Subpop, ncol = 1, scales = "fixed", strip.position ="right")

## Expected heterozygosity
Hexp.figure <- gd.tib %>%   ## Histogram of Hobs
  ggplot(aes(x = Hexp)) +
  geom_histogram(breaks = seq(0,0.5, by = 0.05)) +
  xlab("Expected heterozygosity") +
  ylab("Count") +
  theme_bw(base_size = 22) +
  NULL

Hexp.facetted <- Hexp.figure + facet_wrap( ~ Subpop, ncol = 1, scales = "fixed", strip.position ="right")

## Fis

Fis.figure <- gd.tib %>%   ## Histogram of Hobs
  ggplot(aes(x = Fis)) +
  geom_histogram(breaks = seq(-1.0,1.0, by = 0.05)) +
  xlab("Fis") +
  ylab("Count") +
  theme_bw(base_size = 22) +
  NULL

Fis.facetted <- Fis.figure + facet_wrap( ~ Subpop, ncol = 1, scales = "fixed", strip.position ="right")



## MAF

MAF.figure2 <- MAF.tib %>%   ## Histogram of Hobs
  ggplot(aes(x = MAF)) +
  geom_histogram(breaks = seq(0,0.5, by = 0.01)) +
  xlab("MAF") +
  ylab("Count") +
  theme_linedraw() +
  theme(strip.text.x = element_text(size = 14, angle = 0)) +
  theme(legend.position="none",
        text = element_text(size=14),
        axis.text.x = element_text(angle=0, hjust=1))

MAF.facetted2 <- MAF.figure2 + facet_wrap( ~ Subpop, ncol = 1, scales = "fixed", strip.position ="right")

p <- ggplot(MAF.tib, aes(x=MAF, fill = Subpop)) +
  geom_histogram(binwidth = 0.02, color = "black") +
  xlim(-0.01, 0.5) +
  facet_wrap(~Subpop,nrow = 5)


p + ggtitle("") +
  xlab("MAF") + ylab("Locus Count") +
  theme_linedraw() +
  theme(strip.text.x = element_text(size = 14, angle = 0)) +
  theme(legend.position="none",
        text = element_text(size=14),
        axis.text.x = element_text(angle=0, hjust=1))


#### Pairwise Fsts for all unrelated lamprey ####
pairwise.fsts <- pairwise.fst(dat.genind) ## Fsts using subsetted genind object

