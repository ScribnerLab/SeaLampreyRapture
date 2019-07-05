#######################################################################
####  Outlier analysis of sea lamprey RAPTURE data using OutFLANK  ####
#######################################################################

### Load packages and data
library(OutFLANK)
library(tidyverse)
library(ggrepel)
library(VariantAnnotation)
library(snpStats)
library(SeaLampreyRapture)

data(SeaLampreyRapture)

### Convert files from VCF to 0,1,2,9 format
dat <- genotypeToSnpMatrix(allLoci, uncertain=FALSE)

dat1 <- as.data.frame(dat$genotypes@.Data)

write.csv(dat1, "./tmp/data1temp.csv")
dat1 <- read.csv("./tmp/data1temp.csv") ## Writing then re-reading the csv converts from raw --> int

## Convert from VCF allele representation to the count of reference alleles (9 = no data)
dat1[dat1 == 0] <- 9
dat1[dat1 == 1] <- 0
dat1[dat1 == 2] <- 1
dat1[dat1 == 3] <- 2

## Read in age-specific population names
dat1$pop <- ageNames$V1 ## Assign populations to each individual, stored within genind object

## Subset to only retain largest year class from DCJ and SC
dat2 <- filter(dat1, pop %in% c("BR", "CARP", "DCJ-2", "SC-2", "SM"))

##### Build input files ######
genotype <- dat2[, 2:(ncol(dat2)-1)]
ind <- paste("pop", dat2$pop) # vector with the name of population
locinames <- as.character(colnames(genotype)) # vector of name of loci (actually just index)
FstDataFrame1 <- MakeDiploidFSTMat(genotype, locinames, ind)

##### Perform test #####
out_trim <- OutFLANK(FstDataFrame1, NumberOfSamples=111, qthreshold = 0.05, Hmin = 0.1)

P1 <- pOutlierFinderChiSqNoCorr(FstDataFrame1, Fstbar = out_trim$FSTNoCorrbar,
                                dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)

## Organize dataframe
P1$Index <- as.integer(rownames(P1))
P1$label <- substr(P1$LocusName, 1, 17) # The names of outlier loci
P1$transformed.p <- -log(P1$pvalues)
P1$scaffold <- as.integer(substr(P1$LocusName, 7, 10)) # The names of outlier loci
P1$scaffold[P1$scaffold > 44] <- 9999
P1$scaffold <- factor(P1$scaffold)

## Drop unneeded columns
drops <- c("LocusName","He", "T1", "T2", "T1NoCorr", "T2NoCorr", "meanAlleleFreq", "pvalues",
           "pvaluesRightTail", "qvalues", "OutlierFlag")
plotting.df <- P1[ , !(names(P1) %in% drops)]

axis.df <- plotting.df %>% group_by(scaffold) %>%
  summarize(center = mean(Index)) %>%
  mutate(label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", " ", "12", " ", "14", " ",
                   "16", " ", "18", " ", "20", " ", "22", " ", "24", " ", "26", " ", "28", " ",
                   "30", " ", " ", " ", " ", "35", " ", " ", " ", " ", "40", " ", " ",
                   " ", " ", ">44"))

### Create Manhattan plot
ggplot(P1, aes(Index, transformed.p, label = label)) +
  geom_point(aes(color = P1$scaffold)) +
  ylab(expression(-log[10](p))) +
  scale_colour_manual(values = rep(c("gray80", "gray60"), 48), aesthetics = "color") +
  geom_text_repel(
    data = subset(P1, OutlierFlag == "TRUE")
  ) +
  geom_point(data = P1[P1$OutlierFlag == "TRUE", ], color = "red", shape = 18) +
  scale_x_continuous(label = axis.df$label, breaks= axis.df$center, name = "Scaffold") +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  theme_classic(base_size = 18) +
  theme(legend.position="none") +
  NULL

#### Create plot of full Fst distribution
ggplot(P1, aes(x = FSTNoCorr)) +
  geom_histogram(binwidth = 0.005) +
  geom_vline(xintercept = 0.2674457, size = 1, lty = 5) +
  ylab("Frequency") +
  xlab(expression(F["ST"])) +
  theme_classic(base_size = 24) +
  NULL


### "OF" flags homogeneous data as outliers, therefore, we do not depend on the T or F flag
### Calculated P values are then transformed using qvalue
qval <- qvalue(P1$pvaluesRightTail)$qvalues
alpha <- 0.05
outliers <- which(qval<alpha)
length(outliers) # Number of outlier loci

## Find outlier loci names, rather than numbers
outliers.num <- as.numeric(outliers)
outlierLoci <- locinames[outliers.num]
outlierLoci <- substr(outlierLoci,1,nchar(outlierLoci)-4) # The names of outlier loci
outlierLoci.df <- as.data.frame(outlierLoci)

#write.csv(outlierLoci, "outlierLoci.csv")

