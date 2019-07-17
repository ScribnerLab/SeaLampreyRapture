#originally created on: Nov. 30th, 2018
#by: Nick Sard

#ABOUT: This script was written to test if the proportion of within sampling location related dyads differed among
#age classes

#loading libraries
library(tidyverse)

#setting working directory and reading in data
setwd("C:/Users/sard/Google Drive/R/Data analysis/2018/Sea Lamprey/Real data pedigrees/")

df <-
  matrix(c(370, 1406,2670,21170),
         nrow = 2,
         dimnames = list(related = c("Related", "total"),
                         within = c("age1", "age2")))
df
fisher.test(df)

df <-
  matrix(c(370, 1406,52,870),
         nrow = 2,
         dimnames = list(related = c("Related", "total"),
                         within = c("age1", "age2")))
df
fisher.test(df)

df <-
  matrix(c(2670,21170,52,870),
         nrow = 2,
         dimnames = list(related = c("Related", "total"),
                         within = c("age2", "age3")))
df
fisher.test(df)
