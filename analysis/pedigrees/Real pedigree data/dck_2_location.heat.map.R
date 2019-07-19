#originally created on: Nov. 15th, 2018
#by: Nick Sard

#ABOUT: This script was written to evaluate sibship relationships within a stream network
#There is a set of known locations were lamprey were sampled that can be evaluated across
#all pairwise comparisons. The question that we are trying to answer is how many siblings
#were sampled in the same location relative to other sampling locations

#loading libraries
library(tidyverse)
library(RColorBrewer)
library(SeaLampreyRapture)
data("SeaLampreyRapture")

#picking pretty colors
my.colors <- brewer.pal(7,"OrRd")

### Function ###
#first defining a function I will use later
unique.pairs <- function(ids){
  x <- expand.grid(ids,ids)
  x <- x[ifelse(x$Var1==x$Var2,T,F)==F,]
  x$both <- NULL
  row.names(x) <- NULL
  return(x)
} # end of unique.pairs

#############
### Age 1 ###
#############

#getting and location type, lat, long, df info for merging later
df2 <- database_v2 %>%
  mutate(loc = as.factor(loc), type = as.factor(type)) %>%
  filter(loc == "DCK") %>%
  group_by(type) %>%
  dplyr::slice(1) %>%
  select(type,lat,long)
head(df2)

#selecting cols of interest for merge later
df1 <- database_v2 %>%
  select(id,lat,long,type)
head(df1)

##making sure names are correct in pedigree
#colnames(best.config.dck.age.1) <- c("OffspringID","MotherID","FatherID")
#head(df)

#fixing offspring ID thing - removing a space
best.config.dck.age.1$OffspringID <- gsub(pattern = " ",replacement = "",x = best.config.dck.age.1$OffspringID)
best.config.dck.age.1$ClusterIndex <- NULL
head(best.config.dck.age.1)

#now merging df with df1 to connect IDs with sampling location
df <-
  merge(x = best.config.dck.age.1,
        y = df1,
        by.x = "OffspringID",
        by.y = "id")
head(df)
#note that the "new" should be the same number of rows as the "original"

#getting offspring names and making a df of unique pairs
offs <- df[,1]
tmp <- unique.pairs(offs)
head(tmp)

#making a df that I can use logic to determine sibship
#basically merging the pedigree that has location info attached twice to the
#unique pairs DF

#first one on the first offs ID
tmp1 <- merge(x = df,y = tmp, by.x = "OffspringID",by.y = "Var1",all.y = T)
tmp1$mps.y <- paste(tmp1$FatherID,tmp1$MotherID,sep="_")
names(tmp1)[1] <- "OffspringID.y"
head(tmp1)

#now on the second offs ID
tmp2 <- suppressMessages(merge(x = df,y = tmp1, by.x = "OffspringID",by.y = "Var2",all.y = T))
tmp2$mps.x <- paste(tmp2$FatherID.x,tmp2$MotherID.x,sep="_")
names(tmp2)[1] <- "OffspringID.x"
head(tmp2)

#making three columns to fill in with 1 for if they are full sibs or half sibs
tmp2$fs <- 0
tmp2$mat_hs <- 0
tmp2$pat_hs <- 0

#using loci to identify FS and HS
tmp2$fs[tmp2$mps.x == tmp2$mps.y] <- 1
tmp2$mat_hs[tmp2$MotherID.x == tmp2$MotherID.y] <- 1
tmp2$pat_hs[tmp2$FatherID.x ==  tmp2$FatherID.y] <- 1
#getting read of full-sib "picks ups" with the above logic
tmp2$mat_hs[tmp2$mps.x == tmp2$mps.y] <- 0
tmp2$pat_hs[tmp2$mps.x == tmp2$mps.y] <- 0
head(tmp2)

#getting rid of all un_related dyads
tmp2$all <-  tmp2$fs + tmp2$mat_hs + tmp2$pat_hs
table(tmp2$all)

#filtering on that related (all) column
tmp3 <- filter(tmp2, all == 1)
head(tmp3)

#doing a bit of remaning
names(tmp3)[1] <- "off1"
names(tmp3)[7] <- "off2"
head(tmp3)

#the next part is a bit involved...
#the issue at hand is that the unique.pairs functions provides a
#row for both "Off1 - Off2" and the "Off2 - Off1" dyad. They have the same relationship
#(see example below), but I only want one of those comparisons
tmp3 %>% filter(off1 == "PMDCKJ1_10" & off2 == "PMDCKJ1_14")
tmp3 %>% filter(off2 == "PMDCKJ1_10" & off1 == "PMDCKJ1_14")

#making a column to identify duplicate dyad information
tmp3$both <- paste(pmin(tmp3$off1,tmp3$off2),pmax(tmp3$off1,tmp3$off2))
table(duplicated(tmp3$both))

#filtering out just one
tmp3 <- tmp3[duplicated(tmp3$both)==F,]
head(tmp3)

#now I want to get counts of related dyads within each *location* pairise comparison
tmp4 <- tmp3 %>%
  group_by(type.x, type.y) %>%
  dplyr::count() %>%
  as.data.frame()

head(tmp4)

#the above does not contain locations with 0 observations, so use expand grid
#to get all comparisons - making names character strings and making a unique ID with paste again
df3 <- expand.grid(df2$type,df2$type)
df3$Var1 <- as.character(df3$Var1)
df3$Var2 <- as.character(df3$Var2)
df3$both <- paste0(pmin(df3$Var1,df3$Var2),pmax(df3$Var1,df3$Var2))
head(df3)

#getting the same format for names in tmp4 for a merge
tmp4$type.x <- as.character(tmp4$type.x)
tmp4$type.y <- as.character(tmp4$type.y)
tmp4$both <- paste0(pmin(tmp4$type.x,tmp4$type.y),pmax(tmp4$type.x,tmp4$type.y))
head(tmp4)

#now combining
df3 <- merge(df3,tmp4,by="both",all.x = T)
df3 <- df3 %>% select(Var1,Var2,n)

#filling NAs with zero
df3$n[is.na(df3$n)] <- 0
head(df3)

#making into wideform
df4 <- df3 %>%
  spread(Var2, n)

head(df4)

#the problem is that some values are on the top of the diagonal and others are on
#the bottom. I need to get them on the same side using some base R functions
row.names(df4) <- df4$Var1
df4$Var1 <- NULL
df4 <- replace(df4, !lower.tri(df4,diag = T), NA)

#adding locations names back in
df5 <- as.data.frame(cbind(rownames(df4),as.data.frame(df4)))
names(df5)[1] <- "locA"
rownames(df5) <- NULL
df5

#Changing it to long form and getting rid of NAs
df5 <- gather(df5, key = "locB",value = "value",-locA)
df5 <- df5[!is.na(df5$value),]
df5$value2 <- df5$value
df5$value2[df5$value == 0]<-NA
df5$age <- "Age - 1"
age1 <- df5
head(age1)

#############
### Age 2 ###
#############

#getting and location type, lat, long, df info for merging later
df2 <- database_v2 %>%
  mutate(loc = as.factor(loc), type = as.factor(type)) %>%
  filter(loc == "DCK") %>%
  group_by(type) %>%
  dplyr::slice(1) %>%
  select(type,lat,long)
head(df2)

#selecting cols of interest for merge later
df1 <- database_v2 %>%
  select(id,lat,long,type)
head(df1)

##making sure names are correct in pedigree
#colnames(best.config.dck.age.2) <- c("OffspringID","MotherID","FatherID")
#head(df)

#fixing offspring ID thing - removing a space
best.config.dck.age.2$OffspringID <- gsub(pattern = " ",replacement = "",x = best.config.dck.age.2$OffspringID)
best.config.dck.age.2$ClusterIndex <- NULL
head(best.config.dck.age.2)

#now merging df with df1 to connect IDs with sampling location
df <-
  merge(x = best.config.dck.age.2,
        y = df1,
        by.x = "OffspringID",
        by.y = "id")
head(df)
#note that the "new" should be the same number of rows as the "original"

#getting offspring names and making a df of unique pairs
offs <- df[,1]
tmp <- unique.pairs(offs)
head(tmp)

#making a df that I can use logic to determine sibship
#basically merging the pedigree that has location info attached twice to the
#unique pairs DF

#first one on the first offs ID
tmp1 <- merge(x = df,y = tmp, by.x = "OffspringID",by.y = "Var1",all.y = T)
tmp1$mps.y <- paste(tmp1$FatherID,tmp1$MotherID,sep="_")
names(tmp1)[1] <- "OffspringID.y"
head(tmp1)

#now on the second offs ID
tmp2 <- suppressMessages(merge(x = df,y = tmp1, by.x = "OffspringID",by.y = "Var2",all.y = T))
tmp2$mps.x <- paste(tmp2$FatherID.x,tmp2$MotherID.x,sep="_")
names(tmp2)[1] <- "OffspringID.x"
head(tmp2)

#making three columns to fill in with 1 for if they are full sibs or half sibs
tmp2$fs <- 0
tmp2$mat_hs <- 0
tmp2$pat_hs <- 0

#using loci to identify FS and HS
tmp2$fs[tmp2$mps.x == tmp2$mps.y] <- 1
tmp2$mat_hs[tmp2$MotherID.x == tmp2$MotherID.y] <- 1
tmp2$pat_hs[tmp2$FatherID.x ==  tmp2$FatherID.y] <- 1
#getting read of full-sib "picks ups" with the above logic
tmp2$mat_hs[tmp2$mps.x == tmp2$mps.y] <- 0
tmp2$pat_hs[tmp2$mps.x == tmp2$mps.y] <- 0
head(tmp2)

#getting rid of all un_related dyads
tmp2$all <-  tmp2$fs + tmp2$mat_hs + tmp2$pat_hs
table(tmp2$all)

#filtering on that related (all) column
tmp3 <- filter(tmp2, all == 1)
head(tmp3)

#doing a bit of remaning
names(tmp3)[1] <- "off1"
names(tmp3)[7] <- "off2"
head(tmp3)

#the next part is a bit involved...
#the issue at hand is that the unique.pairs functions provides a
#row for both "Off1 - Off2" and the "Off2 - Off1" dyad. They have the same relationship
#(see example below), but I only want one of those comparisons
tmp3 %>% filter(off1 == "PMDCKJ1_10" & off2 == "PMDCKJ1_14")
tmp3 %>% filter(off2 == "PMDCKJ1_10" & off1 == "PMDCKJ1_14")

#making a column to identify duplicate dyad information
tmp3$both <- paste(pmin(tmp3$off1,tmp3$off2),pmax(tmp3$off1,tmp3$off2))
table(duplicated(tmp3$both))

#filtering out just one
tmp3 <- tmp3[duplicated(tmp3$both)==F,]
head(tmp3)

#now I want to get counts of related dyads within each *location* pairise comparison
tmp4 <- tmp3 %>%
  group_by(type.x, type.y) %>%
  dplyr::count() %>%
  as.data.frame()

head(tmp4)

#the above does not contain locations with 0 observations, so use expand grid
#to get all comparisons - making names character strings and making a unique ID with paste again
df3 <- expand.grid(df2$type,df2$type)
df3$Var1 <- as.character(df3$Var1)
df3$Var2 <- as.character(df3$Var2)
df3$both <- paste0(pmin(df3$Var1,df3$Var2),pmax(df3$Var1,df3$Var2))
head(df3)

#getting the same format for names in tmp4 for a merge
tmp4$type.x <- as.character(tmp4$type.x)
tmp4$type.y <- as.character(tmp4$type.y)
tmp4$both <- paste0(pmin(tmp4$type.x,tmp4$type.y),pmax(tmp4$type.x,tmp4$type.y))
head(tmp4)

#now combining
df3 <- merge(df3,tmp4,by="both",all.x = T)
df3 <- df3 %>% select(Var1,Var2,n)

#filling NAs with zero
df3$n[is.na(df3$n)] <- 0
head(df3)

#making into wideform
df4 <- df3 %>%
  spread(Var2, n)

head(df4)

#the problem is that some values are on the top of the diagonal and others are on
#the bottom. I need to get them on the same side using some base R functions
row.names(df4) <- df4$Var1
df4$Var1 <- NULL
df4 <- replace(df4, !lower.tri(df4,diag = T), NA)

#adding locations names back in
df5 <- as.data.frame(cbind(rownames(df4),as.data.frame(df4)))
names(df5)[1] <- "locA"
rownames(df5) <- NULL
df5

#Changing it to long form and getting rid of NAs
df5 <- gather(df5, key = "locB",value = "value",-locA)
df5 <- df5[!is.na(df5$value),]
df5$value2 <- df5$value
df5$value2[df5$value == 0]<-NA
df5$age <- "Age - 2"
age2 <- df5
head(age2)

#############
### Age 3 ###
#############

#getting and location type, lat, long, df info for merging later
df2 <- database_v2 %>%
  mutate(loc = as.factor(loc), type = as.factor(type)) %>%
  filter(loc == "DCK") %>%
  group_by(type) %>%
  dplyr::slice(1) %>%
  select(type,lat,long)
head(df2)

#selecting cols of interest for merge later
df1 <- database_v2 %>%
  select(id,lat,long,type)
head(df1)

##making sure names are correct in pedigree
#colnames(best.config.dck.age.3) <- c("OffspringID","MotherID","FatherID")
#head(df)

#fixing offspring ID thing - removing a space
best.config.dck.age.3$OffspringID <- gsub(pattern = " ",replacement = "",x = best.config.dck.age.3$OffspringID)
best.config.dck.age.3$ClusterIndex <- NULL
head(best.config.dck.age.3)

#now merging df with df1 to connect IDs with sampling location
df <-
  merge(x = best.config.dck.age.3,
        y = df1,
        by.x = "OffspringID",
        by.y = "id")
head(df)
#note that the "new" should be the same number of rows as the "original"

#getting offspring names and making a df of unique pairs
offs <- df[,1]
tmp <- unique.pairs(offs)
head(tmp)

#making a df that I can use logic to determine sibship
#basically merging the pedigree that has location info attached twice to the
#unique pairs DF

#first one on the first offs ID
tmp1 <- merge(x = df,y = tmp, by.x = "OffspringID",by.y = "Var1",all.y = T)
tmp1$mps.y <- paste(tmp1$FatherID,tmp1$MotherID,sep="_")
names(tmp1)[1] <- "OffspringID.y"
head(tmp1)

#now on the second offs ID
tmp2 <- suppressMessages(merge(x = df,y = tmp1, by.x = "OffspringID",by.y = "Var2",all.y = T))
tmp2$mps.x <- paste(tmp2$FatherID.x,tmp2$MotherID.x,sep="_")
names(tmp2)[1] <- "OffspringID.x"
head(tmp2)

#making three columns to fill in with 1 for if they are full sibs or half sibs
tmp2$fs <- 0
tmp2$mat_hs <- 0
tmp2$pat_hs <- 0

#using loci to identify FS and HS
tmp2$fs[tmp2$mps.x == tmp2$mps.y] <- 1
tmp2$mat_hs[tmp2$MotherID.x == tmp2$MotherID.y] <- 1
tmp2$pat_hs[tmp2$FatherID.x ==  tmp2$FatherID.y] <- 1
#getting read of full-sib "picks ups" with the above logic
tmp2$mat_hs[tmp2$mps.x == tmp2$mps.y] <- 0
tmp2$pat_hs[tmp2$mps.x == tmp2$mps.y] <- 0
head(tmp2)

#getting rid of all un_related dyads
tmp2$all <-  tmp2$fs + tmp2$mat_hs + tmp2$pat_hs
table(tmp2$all)

#filtering on that related (all) column
tmp3 <- filter(tmp2, all == 1)
head(tmp3)

#doing a bit of remaning
names(tmp3)[1] <- "off1"
names(tmp3)[7] <- "off2"
head(tmp3)

#the next part is a bit involved...
#the issue at hand is that the unique.pairs functions provides a
#row for both "Off1 - Off2" and the "Off2 - Off1" dyad. They have the same relationship
#(see example below), but I only want one of those comparisons
tmp3 %>% filter(off1 == "PMDCKJ1_10" & off2 == "PMDCKJ1_14")
tmp3 %>% filter(off2 == "PMDCKJ1_10" & off1 == "PMDCKJ1_14")

#making a column to identify duplicate dyad information
tmp3$both <- paste(pmin(tmp3$off1,tmp3$off2),pmax(tmp3$off1,tmp3$off2))
table(duplicated(tmp3$both))

#filtering out just one
tmp3 <- tmp3[duplicated(tmp3$both)==F,]
head(tmp3)

#now I want to get counts of related dyads within each *location* pairise comparison
tmp4 <- tmp3 %>%
  group_by(type.x, type.y) %>%
  dplyr::count() %>%
  as.data.frame()

head(tmp4)

#the above does not contain locations with 0 observations, so use expand grid
#to get all comparisons - making names character strings and making a unique ID with paste again
df3 <- expand.grid(df2$type,df2$type)
df3$Var1 <- as.character(df3$Var1)
df3$Var2 <- as.character(df3$Var2)
df3$both <- paste0(pmin(df3$Var1,df3$Var2),pmax(df3$Var1,df3$Var2))
head(df3)

#getting the same format for names in tmp4 for a merge
tmp4$type.x <- as.character(tmp4$type.x)
tmp4$type.y <- as.character(tmp4$type.y)
tmp4$both <- paste0(pmin(tmp4$type.x,tmp4$type.y),pmax(tmp4$type.x,tmp4$type.y))
head(tmp4)

#now combining
df3 <- merge(df3,tmp4,by="both",all.x = T)
df3 <- df3 %>% select(Var1,Var2,n)

#filling NAs with zero
df3$n[is.na(df3$n)] <- 0
head(df3)

#making into wideform
df4 <- df3 %>%
  spread(Var2, n)

head(df4)

#the problem is that some values are on the top of the diagonal and others are on
#the bottom. I need to get them on the same side using some base R functions
row.names(df4) <- df4$Var1
df4$Var1 <- NULL
df4 <- replace(df4, !lower.tri(df4,diag = T), NA)

#adding locations names back in
df5 <- as.data.frame(cbind(rownames(df4),as.data.frame(df4)))
names(df5)[1] <- "locA"
rownames(df5) <- NULL
df5

#Changing it to long form and getting rid of NAs
df5 <- gather(df5, key = "locB",value = "value",-locA)
df5 <- df5[!is.na(df5$value),]
df5$value2 <- df5$value
df5$value2[df5$value == 0]<-NA
df5$age <- "Age - 3"
age3 <- df5
head(age3)

#combining
all.ages <- rbind(age1,age2,age3)
head(all.ages)

#now the figure
tiff(filename = "./analysis/pedigrees/Real pedigree data/Output/all.ages.dyads.heatmap.tiff",width = 14,height = 8,units = "in",res = 300)
ggplot(all.ages, aes(x = locA, y=locB, fill=value2,label=value2))+
  facet_wrap(~age)+
  geom_tile(color="black")+
  geom_text(size=5)+
  scale_fill_gradientn(colours = my.colors)+
  scale_y_discrete(position = "right")+
  theme_bw()+
  labs(x = "",y="", fill="Related dyads")+
  theme(legend.title= element_text(size=18),
        legend.text = element_text(size =12),
        strip.text = element_text(size=14),
        axis.text = element_text(size=18,colour = "black"),
        axis.text.x = element_text(angle = 90,vjust = 0.5),
        legend.position = c(0.10, 0.75))
dev.off()

#fin!
