#originally created on: Nov. 15th, 2018
#by: Nick Sard

#ABOUT: This script was written to test whether clusters of siblings within a locations
#were expected by random change within a stream network

#loading libraries
library(tidyverse)
library(SeaLampreyRapture)
data("SeaLampreyRapture")

### Function ###
#first defining a function I will use later
unique.pairs <- function(ids){
  x <- expand.grid(ids,ids)
  x <- x[ifelse(x$Var1==x$Var2,T,F)==F,]
  x$both <- NULL
  row.names(x) <- NULL
  return(x)
} # end of unique.pairs

#making buckets to be filled
real.data <- NULL
sim.out <- NULL

#############
### age 1 ###
#############

#reading in data
df <- best.config.dck.age.1
df$ClusterIndex <- NULL
head(df)

#reading in loc data
df1 <- database_v2
head(df1)

#making some lat, long, df info for merging later
df2 <- df1 %>% filter(loc == "DCK") %>% group_by(type) %>% slice(1) %>% select(type,lat,long)
head(df2)

#selecting cols of interest for merge with pedigree
df1 <- df1 %>% select(id,lat,long,type)
head(df1)

#making sure names are correct
colnames(df) <- c("OffspringID","MotherID","FatherID")
head(df)

#fixing a name thing
df$OffspringID <- gsub(pattern = " ",replacement = "",x = df$OffspringID)
head(df)

#now merging df with df1
df <-
  merge(x = df,
        y = df1,
        by.x = "OffspringID",
        by.y = "id")
head(df)

#getting offspring names and making a df of unique pairs
offs <- df[,1]
tmp <- unique.pairs(offs)
head(tmp)

#making a df that I can use logic to determine sibship
tmp1 <- merge(x = df,y = tmp, by.x = "OffspringID",by.y = "Var1",all.y = T)
tmp1$mps.y <- paste(tmp1$FatherID,tmp1$MotherID,sep="_")
names(tmp1)[1] <- "OffspringID.y"
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
tmp2$mat_hs[tmp2$mps.x == tmp2$mps.y] <- 0
tmp2$pat_hs[tmp2$mps.x == tmp2$mps.y] <- 0
head(tmp2)

#getting rid of all un_related dyads
tmp2$all <- tmp2$fs + tmp2$mat_hs + tmp2$pat_hs
head(tmp2)

#getting cols I want for the randomization prep
tmp3 <- tmp2 %>% select(OffspringID.x, OffspringID.y, type.x, type.y,all)
tmp3$within <- ifelse(tmp3$type.x == tmp3$type.y, "Y","N")
head(tmp3)

#counts of related dyads, related within a single location, and total number of dyads overall
dyads.related <- (tmp3 %>% filter(all == 1) %>% nrow())
dyads.related.within <- (tmp3 %>% filter(within == "Y" & all == 1) %>% nrow())
dyads.total <- nrow(tmp3)

#now the proportion of within location dyads relative to all related dyads
real.prop <- dyads.related.within/dyads.total
real.prop

#now for a randomization, which asks the question given the sampling distribution of
#larvae, would we expext that the number of related individuals within sampling location
#relative to all related individuals be expected by chance

#the for-loop genartes a null distribution to test that question
tmp4 <- tmp3
out <- NULL
i <- 1
i <- NULL
for(i in 1:1000){
  print(i)

  #resetting "all" to zero
  tmp4$all <- 0

  #shuffling the rows randomly
  tmp4$rand <- runif(n = nrow(tmp4))
  tmp4 <- tmp4 %>% arrange(rand)

  #making 1:dyads.related rows "related" - note that these are shuffled
  tmp4$all[1:dyads.related] <- 1

  #calculating the within location related expectation
  within.rand <- (tmp4 %>% filter(within == "Y" & all == 1) %>% nrow())
  total <- nrow(tmp4)
  prop <- within.rand/total
  sim <- i

  #saving this information for later
  out1 <- data.frame(sim,within.rand,dyads.related,total,prop,stringsAsFactors = F)
  out1
  out <- rbind(out,out1)
}
out$age <- "Age - 1"
head(out)

#saving for later
sim.out <- rbind(sim.out,out)
head(sim.out)

#doing the stats
pvalue <- length(out$prop[out$prop>real.prop])/nrow(out)
real.data1 <- tibble(dyads.total,dyads.related.within,dyads.related,real.prop, age= "Age - 1",pvalue)
real.data <- rbind(real.data,real.data1)
real.data

#############
### age 2 ###
#############

#reading in data
df <- best.config.dck.age.2
df$ClusterIndex <- NULL
head(df)

#reading in loc data
df1 <- database_v2
head(df1)

#making some lat, long, df info for merging later
df2 <- df1 %>% filter(loc == "DCK") %>% group_by(type) %>% slice(1) %>% select(type,lat,long)
head(df2)

#selecting cols of interest for merge with pedigree
df1 <- df1 %>% select(id,lat,long,type)
head(df1)

#making sure names are correct
colnames(df) <- c("OffspringID","MotherID","FatherID")
head(df)

#fixing a name thing
df$OffspringID <- gsub(pattern = " ",replacement = "",x = df$OffspringID)
head(df)

#now merging df with df1
df <-
  merge(x = df,
        y = df1,
        by.x = "OffspringID",
        by.y = "id")
head(df)

#getting offspring names and making a df of unique pairs
offs <- df[,1]
tmp <- unique.pairs(offs)
head(tmp)

#making a df that I can use logic to determine sibship
tmp1 <- merge(x = df,y = tmp, by.x = "OffspringID",by.y = "Var1",all.y = T)
tmp1$mps.y <- paste(tmp1$FatherID,tmp1$MotherID,sep="_")
names(tmp1)[1] <- "OffspringID.y"
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
tmp2$mat_hs[tmp2$mps.x == tmp2$mps.y] <- 0
tmp2$pat_hs[tmp2$mps.x == tmp2$mps.y] <- 0
head(tmp2)

#getting rid of all un_related dyads
tmp2$all <- tmp2$fs + tmp2$mat_hs + tmp2$pat_hs
head(tmp2)

#getting cols I want for the randomization prep
tmp3 <- tmp2 %>% select(OffspringID.x, OffspringID.y, type.x, type.y,all)
tmp3$within <- ifelse(tmp3$type.x == tmp3$type.y, "Y","N")
head(tmp3)

#counts of related dyads, related within a single location, and total number of dyads overall
dyads.related <- (tmp3 %>% filter(all == 1) %>% nrow())
dyads.related.within <- (tmp3 %>% filter(within == "Y" & all == 1) %>% nrow())
dyads.total <- nrow(tmp3)

#now the proportion of within location dyads relative to all related dyads
real.prop <- dyads.related.within/dyads.total
real.prop

#now for a randomization, which asks the question given the sampling distribution of
#larvae, would we expext that the number of related individuals within sampling location
#relative to all related individuals be expected by chance

#the for-loop genartes a null distribution to test that question
tmp4 <- tmp3
out <- NULL
i <- 1
i <- NULL
for(i in 1:1000){
  print(i)

  #resetting "all" to zero
  tmp4$all <- 0

  #shuffling the rows randomly
  tmp4$rand <- runif(n = nrow(tmp4))
  tmp4 <- tmp4 %>% arrange(rand)

  #making 1:dyads.related rows "related" - note that these are shuffled
  tmp4$all[1:dyads.related] <- 1

  #calculating the within location related expectation
  within.rand <- (tmp4 %>% filter(within == "Y" & all == 1) %>% nrow())
  total <- nrow(tmp4)
  prop <- within.rand/total
  sim <- i

  #saving this information for later
  out1 <- data.frame(sim,within.rand,dyads.related,total,prop,stringsAsFactors = F)
  out1
  out <- rbind(out,out1)
}
out$age <- "Age - 2"
head(out)

#saving for later
sim.out <- rbind(sim.out,out)
head(sim.out)

#doing the stats
pvalue <- length(out$prop[out$prop>real.prop])/nrow(out)
real.data1 <- tibble(dyads.total,dyads.related.within,dyads.related,real.prop, age= "Age - 2",pvalue)
real.data <- rbind(real.data,real.data1)
real.data

#############
### age 3 ###
#############

#reading in data
df <- best.config.dck.age.3
df$ClusterIndex <- NULL
head(df)

#reading in loc data
df1 <- database_v2
head(df1)

#making some lat, long, df info for merging later
df2 <- df1 %>% filter(loc == "DCK") %>% group_by(type) %>% slice(1) %>% select(type,lat,long)
head(df2)

#selecting cols of interest for merge with pedigree
df1 <- df1 %>% select(id,lat,long,type)
head(df1)

#making sure names are correct
colnames(df) <- c("OffspringID","MotherID","FatherID")
head(df)

#fixing a name thing
df$OffspringID <- gsub(pattern = " ",replacement = "",x = df$OffspringID)
head(df)

#now merging df with df1
df <-
  merge(x = df,
        y = df1,
        by.x = "OffspringID",
        by.y = "id")
head(df)

#getting offspring names and making a df of unique pairs
offs <- df[,1]
tmp <- unique.pairs(offs)
head(tmp)

#making a df that I can use logic to determine sibship
tmp1 <- merge(x = df,y = tmp, by.x = "OffspringID",by.y = "Var1",all.y = T)
tmp1$mps.y <- paste(tmp1$FatherID,tmp1$MotherID,sep="_")
names(tmp1)[1] <- "OffspringID.y"
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
tmp2$mat_hs[tmp2$mps.x == tmp2$mps.y] <- 0
tmp2$pat_hs[tmp2$mps.x == tmp2$mps.y] <- 0
head(tmp2)

#getting rid of all un_related dyads
tmp2$all <- tmp2$fs + tmp2$mat_hs + tmp2$pat_hs
head(tmp2)

#getting cols I want for the randomization prep
tmp3 <- tmp2 %>% select(OffspringID.x, OffspringID.y, type.x, type.y,all)
tmp3$within <- ifelse(tmp3$type.x == tmp3$type.y, "Y","N")
head(tmp3)

#counts of related dyads, related within a single location, and total number of dyads overall
dyads.related <- (tmp3 %>% filter(all == 1) %>% nrow())
dyads.related.within <- (tmp3 %>% filter(within == "Y" & all == 1) %>% nrow())
dyads.total <- nrow(tmp3)

#now the proportion of within location dyads relative to all related dyads
real.prop <- dyads.related.within/dyads.total
real.prop

#now for a randomization, which asks the question given the sampling distribution of
#larvae, would we expext that the number of related individuals within sampling location
#relative to all related individuals be expected by chance

#the for-loop genartes a null distribution to test that question
tmp4 <- tmp3
out <- NULL
i <- 1
i <- NULL
for(i in 1:1000){
  print(i)

  #resetting "all" to zero
  tmp4$all <- 0

  #shuffling the rows randomly
  tmp4$rand <- runif(n = nrow(tmp4))
  tmp4 <- tmp4 %>% arrange(rand)

  #making 1:dyads.related rows "related" - note that these are shuffled
  tmp4$all[1:dyads.related] <- 1

  #calculating the within location related expectation
  within.rand <- (tmp4 %>% filter(within == "Y" & all == 1) %>% nrow())
  total <- nrow(tmp4)
  prop <- within.rand/total
  sim <- i

  #saving this information for later
  out1 <- data.frame(sim,within.rand,dyads.related,total,prop,stringsAsFactors = F)
  out1
  out <- rbind(out,out1)
}
out$age <- "Age - 3"
head(out)

#saving for later
sim.out <- rbind(sim.out,out)
head(sim.out)

#doing the stats
pvalue <- length(out$prop[out$prop>real.prop])/nrow(out)
real.data1 <- tibble(dyads.total,dyads.related.within,dyads.related,real.prop, age= "Age - 3",pvalue)
real.data <- rbind(real.data,real.data1)
real.data

#writing to file
write.table(x = real.data,file = "./analysis/pedigrees/Read pedigree data/Output/randomization.results.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

#fin!
