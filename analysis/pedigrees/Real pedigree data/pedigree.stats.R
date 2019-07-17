#originally created on: Nov. 14th, 2018
#by: Nick Sard

#ABOUT: This script was written to calculation the Number of successful parents,
#half-sib dyads, and full-sib dyads within a given pedigree, along with reproductive success and
#variation in reproductive success, and coancestry

#loading libraries
library(tidyverse)
library(SeaLampreyRapture)
data("SeaLampreyRapture")

#loading my own functions
source("./R/coancestry.R")

#setting working directory and reading in data
setwd("C:/Users/sard/Google Drive/R/Data analysis/2018/Sea Lamprey/Github/Real pedigree data/")

#reading in data
my.files <- list.files(path = "./extData/",pattern = "best.config")
my.files

#for loop to do the work in all files
i <- 1
i <- NULL
out <- NULL
for(i in 1:length(my.files)){
  print(i)

  #reading in the best configuration file
  df <- read.table(file = paste0("./extData/",my.files[i]),header = T,sep = "\t",stringsAsFactors = F)
  head(df)

  #some other stats

  #first manipulating the file ID to get Age and location information
  my.id <- my.files[i]
  my.id <- gsub(pattern = "best.config.",replacement = "",x = my.id)
  my.id <- gsub(pattern = "age.",replacement = "",x = my.id)
  my.id <- gsub(pattern = ".txt",replacement = "",x = my.id)
  my.id <- gsub(pattern = "\\.",replacement = "_",x = my.id)
  my.id

  #using my in-house coancestry script to calculate coancestery
  tmp <- coancestry(df = df[,-4],info = T)
  head(tmp)

  #getting the number of unique parents observed in the pedigree
  #as a measure of successful parents
  n.par <- length(unique(gsub(pattern = " ",replacement = "",x = c(df$FatherID,df$MotherID))))

  #getting the number of genotyped offspring
  n.gt <- nrow(df)

  #placing that parent and genotyped offspring information into a DF
  tmp <- rbind(tmp,data.frame(type = "npar",value = n.par))
  tmp <- rbind(tmp,data.frame(type = "ngt",value = n.gt))

  #keeping the information I want
  tmp <- tmp[tmp$type %in% c("n_hs","n_fs","ngt","coancestry","npar","n.gt"),]

  #identifying which location and year these data are associated with
  tmp$id <- my.id

  #now going to calculate reproductive success (RS) and variance in RS
  #short and way to do this get a table of parent IDs, do this twice (one for each sex)
  #take the second column whish is vector of RS values
  my.rs <- rbind(as.data.frame(table(df$FatherID)),as.data.frame(table(df$MotherID)))[,2]
  my.mean <- mean(my.rs)
  my.var <- var(my.rs)

  #placing that RS and var(RS) information into a DF
  tmp <- rbind(tmp,data.frame(type = c("mean","var"),value = c(my.mean,my.var),id = my.id))

  #saving to out
  out <- rbind(out,tmp)
}
out$value <- round(out$value,4)

#spreading the data to make it wideform
out <- out %>% spread(key = id,value = value)
out

#writing to file
write.table(x = out,file = "Output/pedigree.summary.k.vk.stats.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

#fin!
