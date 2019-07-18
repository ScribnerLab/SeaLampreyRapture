#!/usr/bin/env Rscript  

#written by: Nick Sard
#written on: September 5 2017

#ABOUT: This script was written to gather all output files from colony simulations and
#collate them all into a single file - identified by the simulation ID

#identifying files with output in the names
my.files <- list.files(pattern="output")

#using a for loop to go through all the files and collate them into a single summary file
#i <- 1
i <- NULL
out <- NULL
for(i in 1:length(my.files)){
	print(i)
	df <- read.table(my.files[i],head=T,stringsAsFactors=F)
	df$sim <- gsub("output\\.","",my.files[i])
	out <- rbind(out,df)
} #end of i for loop

#writing out to file
write.table(out, "output.summary.txt",sep="\t",quote=F,col.names=T,row.names=F)

#fin!
