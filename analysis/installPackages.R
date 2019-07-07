#########################################################
#####         Install packages needed for           #####
#####          WFUrbanAdaptation analyses           #####
#########################################################

list.of.packages <- c("tidyverse",
                      "vcfR",
                      "ggrepel",
                      "devtools",
                      "egg",
                      "adegenet",
                      "hierfstat",
                      "gridExtra",
                      "ggthemes")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

#### Install qvalue from BioConductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
biocLite("VariantAnnotation")
biocLite("snpStats")
biocLite("quantsmooth")
install_github("whitlock/outFLANK", force = TRUE)
