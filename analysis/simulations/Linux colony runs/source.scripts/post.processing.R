#!/usr/bin/env Rscript

#originally created on: August 29, 2017
#by: Nick Sard

#ABOUT: This script was written to compare a known pedigree to the inferred from COLONY

#loading in libraries
library(tidyverse)

#loading my own functions
unique.pairs <- function(ids){
  x <- expand.grid(ids,ids)
  x <- x[ifelse(x$Var1==x$Var2,T,F)==F,] #removing dyads comparing one individual to itself
  x$both <- NULL
  row.names(x) <- NULL
  return(x)
} # end of unique.pairs 

#reading in the pedigree
df <- read.table(file = "best.config.txt",header=T,stringsAsFactors = F)
colnames(df) <- c("off","par1","par2","clus")
#head(df)

#identifying unique assigned mate pairs
df$mp1 <- paste(df$par1,df$par2,sep = "_")
#head(df)

#separating the off id into the true parent ids and off kid id
df$kid2 <- df$off
df$kid2 <- gsub("F","_F",df$kid2)
df$kid2 <- gsub("C","_C",df$kid2)
df <- separate(df, kid2 ,into = c("male","female","kid"), sep="_")
#head(df)

#identifying true mate pairs
df$mp2 <- paste(df$male,df$female,sep="_")
#head(df)

#writing to file for hand annotation
#write.table(x = df,file = "Output/pedigree.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

#I am calculating a suite things to look at 
#differnces beween estimated and known values to look a bias
#full sibling familes estimates and known
est.mp <- length(unique(df$mp1))
real.mp <- length(unique(df$mp2))

#number of parents estimated and known
est.par <- length(unique(df$par1))+length(unique(df$par2))
real.par <- length(unique(df$male))+length(unique(df$female))

#for inferred parents
#comparing rs information into a df to calculation reproductive success (rs)
my.par1 <- unique(df$par1)
my.par2 <- unique(df$par2)
rs.dist.est <- data.frame(par = c(my.par1,my.par2),
                          type = c(rep("par1",times=length(my.par1)), rep("par2",times=length(my.par2))),
                          rs = 0, stringsAsFactors = F)

#using a for loop to count RS
i <- NULL
for(i in 1:nrow(rs.dist.est)){
  if(rs.dist.est$type[i] == "par1"){
    rs.dist.est$rs[i] <- nrow(df[df$par1 == rs.dist.est$par[i],])
  } else {
    rs.dist.est$rs[i] <- nrow(df[df$par2 == rs.dist.est$par[i],])
  }
}
est.mean.rs <- round(mean(rs.dist.est$rs),2)
est.var.rs <- round(var(rs.dist.est$rs),2)

#for known parents...doing the same thing 
my.par1 <- unique(df$male)
my.par2 <- unique(df$female)
rs.dist.real <- data.frame(par = c(my.par1,my.par2),
                          type = c(rep("par1",times=length(my.par1)), rep("par2",times=length(my.par2))),
                          rs = 0, stringsAsFactors = F)
i <- NULL
for(i in 1:nrow(rs.dist.real)){
  if(rs.dist.real$type[i] == "par1"){
    rs.dist.real$rs[i] <- nrow(df[df$male == rs.dist.real$par[i],])
  } else {
    rs.dist.real$rs[i] <- nrow(df[df$female == rs.dist.real$par[i],])
  }
}
real.mean.rs <- round(mean(rs.dist.real$rs),2)
real.var.rs <- round(var(rs.dist.real$rs),2)

#######################################
### Inferred data set dyad analysis ###
#######################################

#head(df)
est <- df[,1:3]
colnames(est) <- c("OffspringID","FatherID","MotherID")
#head(est)

#getting offspring names and making a est of unique pairs
offs <- est[,1]
tmp <- unique.pairs(offs)
#tmp

#making a est that I can use logic to determine sibship
tmp1 <- merge(x = est,y = tmp, by.x = "OffspringID",by.y = "Var1",all.y = T)
tmp1$mps.y <- paste(tmp1$FatherID,tmp1$MotherID,sep="_")
names(tmp1)[1] <- "OffspringID.y"
tmp2 <- suppressMessages(merge(x = est,y = tmp1, by.x = "OffspringID",by.y = "Var2",all.y = T))
tmp2$mps.x <- paste(tmp2$FatherID.x,tmp2$MotherID.x,sep="_")
names(tmp2)[1] <- "OffspringID.x"

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
#head(tmp2)
	
#saving for later
est.ped <- tmp2

#these are the number of dyads
est.FS <- sum(tmp2$fs)
est.mHS <- sum(tmp2$mat_hs)
est.pHS <- sum(tmp2$pat_hs)

###################################
### Real data set dyad analysis ###
###################################
#head(df)
real <- df[,c(1,6,7)]
colnames(real) <- c("OffspringID","FatherID","MotherID")
#head(real)

#getting offspring names and making a real of unique pairs
offs <- real[,1]
tmp <- unique.pairs(offs)
#tmp

#making a real that I can use logic to determine sibship
tmp1 <- merge(x = real,y = tmp, by.x = "OffspringID",by.y = "Var1",all.y = T)
tmp1$mps.y <- paste(tmp1$FatherID,tmp1$MotherID,sep="_")
names(tmp1)[1] <- "OffspringID.y"
tmp2 <- suppressMessages(merge(x = real,y = tmp1, by.x = "OffspringID",by.y = "Var2",all.y = T))
tmp2$mps.x <- paste(tmp2$FatherID.x,tmp2$MotherID.x,sep="_")
names(tmp2)[1] <- "OffspringID.x"

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
#head(tmp2)

#saving for later 
real.ped <- tmp2

#these are the number of dyads
real.FS <- sum(tmp2$fs)
real.mHS <- sum(tmp2$mat_hs)
real.pHS <- sum(tmp2$pat_hs)

#summarizing the data
real <- c(real.mp,real.FS,real.mHS,real.pHS,real.par,real.mean.rs,real.mean.rs)
est <- c(est.mp,est.FS,est.mHS,est.pHS,est.par,est.mean.rs,est.mean.rs)
paramater <- c("mp","nFS","nMHS","nPHS","npar","mean.rs","var.rs")
df1 <- data.frame(paramater,est,real,stringsAsFactors = F)
#df1

#combining the two
#head(est.ped)
#head(real.ped)

#getting the types for estimate ped
est.ped <- est.ped[,c(1,9:11)]
est.ped$type <- "UR"
est.ped$type[est.ped$fs == 1] <- "FS"
est.ped$type[est.ped$mat_hs == 1] <- "HS"
est.ped$type[est.ped$pat_hs == 1] <- "HS"
#table(est.ped$type)
#head(est.ped)

#getting the types for realimate ped
real.ped <- real.ped[,c(1,9:11)]
real.ped$type <- "UR"
real.ped$type[real.ped$fs == 1] <- "FS"
real.ped$type[real.ped$mat_hs == 1] <- "HS"
real.ped$type[real.ped$pat_hs == 1] <- "HS"
#table(real.ped$type)
#head(real.ped)

#table(real.ped$type,est.ped$type)
est.ped <- est.ped[,c(1,5)]
real.ped <- real.ped[,c(1,5)]
colnames(est.ped) <- c("off.est","est.type")
colnames(real.ped) <- c("off.real","real.type")
peds <- cbind(est.ped,real.ped)
#head(peds)

#table(peds$off.est == peds$off.real)
if(unique(peds$off.est == peds$off.real) != T){stop("Names don't match up in pedigree comparison")}

x <- peds %>%
  group_by(est.type,real.type) %>%
  summarize(ns = n())  %>% ungroup() %>% as.data.frame()
#x

#making a matrix to understand if and when errors in dyadic relationships occured
out <- data.frame(real.type = c("FS","HS","UR"), FS = 0, HS = 0 , UR = 0, stringsAsFactors = F)
#out

#filling that table in with the appropriate information
i <- 1
i <- NULL
for(i in 1:nrow(out)){
  
  if(nrow(x[x$real.type == out$real.type[i],])>0){
    my.est.type <-  x[x$real.type == out$real.type[i],1]
    my.est.type <- as.character(my.est.type)
    if(length(my.est.type[my.est.type == "FS"]) == 1){
      out$FS[i] <- x$ns[x$real.type == out$real.type[i] & x$est.type == "FS"]
    }
    if(length(my.est.type[my.est.type == "HS"]) == 1){
      out$HS[i] <- x$ns[x$real.type == out$real.type[i] & x$est.type == "HS"]
    }
    if(length(my.est.type[my.est.type == "UR"])){
      out$UR[i] <- x$ns[x$real.type == out$real.type[i] & x$est.type == "UR"]      
    }
  }
}
out$total = out$FS + out$HS + out$UR
out$pFS  = out$FS/out$total
out$pHS = out$HS/out$total
out$pUR = out$UR/out$total
#out

#this is called a confusion matrix 
conf.mat <- out

#writing to out
write.table(x = df1,file = "output.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
write.table(x = conf.mat,file = "conf.mat.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

#fin!
