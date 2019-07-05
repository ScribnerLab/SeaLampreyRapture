################################################################
######      Sea lamprey DAPC based on RAPTURE panel       ######
################################################################

##### Load data and libraries #####
library(adegenet)
library(vcfR)
library(tidyverse)
library(SeaLampreyRapture)

data(SeaLampreyRapture)

##### Convert vcfR data to genind #####
dat.genind <- vcfR2genind(allLoci.vcfR, sep = "[|/]") ## Create genind object

### Assign population and age class names to genind file
dat.genind@pop <- agePops$V1

###### Perform discriminant analysis of principal components ######
dapc.all <- dapc(dat.genind, var.contrib = TRUE, scale = FALSE, n.pca = 100, n.da = nPop(dat.genind) - 1)

### ID loci with 10% of highest loadings for LD1 & LD2
loadings <- as.data.frame(dapc.all$var.contr)[1:2]
#loadings$LD1.std <- loadings$LD1/0.00147533
#loadings$LD2.std <- loadings$LD2/0.0009887752

## Ranking loci by eigenvalue
perc.rank <- function(x) trunc(rank(x))/length(x)
loadings.LD1 <- within(loadings, LD1.rank <- perc.rank(LD1))
loadings.LD1 <- loadings.LD1[order(loadings.LD1$LD1, decreasing = TRUE),]
loadings.LD1.subset <- subset(loadings.LD1, LD1.rank > 0.9)
tail(loadings.LD1.subset)

loadings.LD2 <- within(loadings, LD2.rank <- perc.rank(LD2))
loadings.LD2 <- loadings.LD2[order(loadings.LD2$LD2, decreasing = TRUE),]
loadings.LD2.subset <- subset(loadings.LD2, LD2.rank > 0.9)
#write.csv(loadings, "dapcLoadings.loci.csv")


### Plot DAPC results
## Core plot
dapc.df <- cbind(as.data.frame(dapc.all$ind.coord[,1:4]), dat.genind@pop)
names(dapc.df) <- c("A1", "A2", "A3", "A4", "pop")
centroids <- aggregate(cbind(A1, A2)~pop, dapc.df, mean)

dapc.plot <- ggplot(dapc.df, aes(x = A1, y = A2, color = pop, shape = pop)) +
  geom_point(alpha = 0.25, size = 3) +
  geom_point(data = centroids, size = 5) +
  ylab("Discriminant function 2") +
  xlab("Discriminant function 1") +
  stat_ellipse() +
  scale_shape_manual(name = "Population & Age",
                     labels = c("BR", "CARP", "DCJ, Age 1", "DCJ, Age 2", "DCJ, Age 3", "SC, Age 2", "SC, Age 3", "SM"),
                     values = c(8, 19, 17, 17, 17, 18, 18, 15)) +
  scale_color_manual(name = "Population & Age",
                     labels = c("BR", "CARP", "DCJ, Age 1", "DCJ, Age 2", "DCJ, Age 3", "SC, Age 2", "SC, Age 3", "SM"),
                     values = c("#000000", "#98a93c", "#33ffcc", "#00b359", "#336600", "#ff4242", "#d00000", "#CC79A7")) +
  ylim(-13, 8.0) +
  theme_classic(base_size = 22) +
  theme(legend.position="none") +
  guides(size = FALSE,
         shape = guide_legend(override.aes = list(size = 5))) +
  NULL


## Plot variation explained by retained principal components
# Set up data frame
dapcPCAvar <- dapc.all$pca.eig
dapcPCAvar <- 100 * cumsum(dapcPCAvar) / sum(dapcPCAvar)
dapcPCAvar <- as.data.frame(dapcPCAvar)
dapcPCAvar$comp <- seq(1, nrow(dapcPCAvar))
dapcPCAvar$color <- cut(
  dapcPCAvar$comp,
  breaks = c(-Inf, 100, Inf),
  labels = c("retained", "not retained")
)

# Perform plotting
pca.plot <- ggplot(dapcPCAvar, aes(x = comp, y = dapcPCAvar, fill = color)) +
  scale_fill_manual(values = c("gray28", "gray90")) +
  geom_bar(stat = "identity") +
  ylab("Variation explained (%)") +
  xlab("Principal component") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  NULL

#### Create final figure
final.dapc <- dapc.plot + annotation_custom(ggplotGrob(pca.plot),
                                            xmin = 3.5, ymin = -13.5,
                                            xmax = 13, ymax = -4)

#ggsave("lampreyDAPC.png", plot = final.dapc, width = 10, height = 7, units = "in")
