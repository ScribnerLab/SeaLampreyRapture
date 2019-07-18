#originally created on: August 21, 2017
#by: Nick Sard

#ABOUT: This script was written to summarize information in confusion matrices calculated for
#a set of simulated genotypic datasets run through COLONY

#loading in libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(SeaLampreyRapture)
data("SeaLampreyRapture")

#setting working directory and reading in data
myfiles <- list.files(path = "./analysis/simulations/Post Linux analyis/Input",pattern = "conf")
myfiles

#reading in files
i <- 5
i <- NULL
df <- NULL
for(i in 1:length(myfiles)){
  print(myfiles[i])

  #reading  in file and getting cols I wan
  df1 <- read.table(file = paste0("./analysis/simulations/Post Linux analyis/Input/",myfiles[i]),header = T,sep = "\t",stringsAsFactors = F)
  df1 <- df1[,c("real.type","FS","HS","UR","sim")]
  head(df1)
  names(df1)[5] <- "rep"

  #getting some naming done
  my.loci <- gsub(pattern = "_conf*.+",replacement = "",x = myfiles[i])
  my.loci <- gsub(pattern = "*.+loci",replacement = "",x = my.loci)
  my.loci

  my.par <- gsub(pattern = "_conf*.+",replacement = "",x = myfiles[i])
  my.par <- gsub(pattern = "npar",replacement = "",x = my.par)
  my.par <- gsub(pattern = "_nloci*.+",replacement = "",x = my.par)
  my.par

  my.sim <- gsub(pattern = "_conf*.+",replacement = "",x = myfiles[i])
  my.sim <- gsub(pattern = "npar",replacement = "Parents = ",x = my.sim)
  my.sim <- gsub(pattern = "_nloci",replacement = " & Loci = ",x = my.sim)
  my.sim

  df1$par <- my.par
  df1$loci <- my.loci
  df1$rep <- gsub(pattern = "conf.mat.colony2_sim_",replacement = "",x = df1$rep)
  df1$rep <- gsub(pattern = "conf.mat.",replacement = "",x = df1$rep)
  df1$rep <- gsub(pattern = ".Dat",replacement = "",x = df1$rep)
  df1$sim <- my.sim
  head(df1)

  #quick stats
  df1 <- df1 %>%
    mutate(total = FS + HS + UR) %>%
    mutate(pFS = round(FS/total,3), pHS = round(HS/total,3), pUR = round(UR/total,3))
  head(df1)

  #combining for larger df
  df <- rbind(df,df1)
}
head(df)
table(df$sim)

#writing combined table to file
write.table(x = df, file = "./analysis/simulations/Post Linux analyis/Output/simulation.output.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

#now for a graphic looking at the false negative assignment rate
tmp <- df %>%
  select(sim,rep,par,loci,real.type:UR,pFS:pUR) %>%
  gather("est.type","val",FS:UR) %>%
  filter(est.type == real.type) %>%
  gather("prop.type","val2",pFS:pUR) %>%
  mutate(prop2 = gsub(pattern = "^p",replacement = "",x = prop.type)) %>%
  filter(est.type == prop2)
head(tmp)

#a bit more naming
head(tmp)
tmp$real.type2[tmp$real.type == "FS"] <- "Full-sibs"
tmp$real.type2[tmp$real.type == "HS"] <- "Half-sibs"
tmp$real.type2[tmp$real.type == "UR"] <- "Unrelated"
tmp$loci <- paste(tmp$loci,"loci")
tmp$par <- paste(tmp$par,"parents")
head(tmp)
tmp$prop2

#calculating full negative rate
tmp$fneg = 1-tmp$val2
head(tmp)

tiff(filename = "./analysis/simulations/Post Linux analyis/Output/false.negative.relationship.tiff",width = 8,height = 8,units = "in",res = 300)
ggplot(tmp, aes(x=real.type2,y=fneg))+
  facet_grid(par~loci)+
  geom_boxplot()+
  theme_bw()+
  labs(x="Dyad relationship",y="False negative assignment rate\n")+
  theme(axis.text.x = element_text(angle=90,vjust=0.5),
        axis.text = element_text(size=18,color = "black"),
        strip.text = element_text(size=16,color = "black"),
        axis.title = element_text(size=22))
dev.off()

head(tmp)

#back to df to flip things  on their head, and looking at accuracy
head(df)
tmp1 <- df %>%
  select(real.type:sim) %>%
  gather("est.type","val",FS:UR) %>%
  spread("real.type","val") %>%
  mutate(total = FS + HS + UR) %>%
  mutate(pFS = round(FS/total,3), pHS = round(HS/total,3), pUR = round(UR/total,3)) %>%
  select(sim,rep,par,loci,est.type:UR,pFS:pUR) %>%
  gather("real.type","val",FS:UR) %>%
  filter(est.type == real.type) %>%
  gather("prop.type","val2",pFS:pUR) %>%
  mutate(prop2 = gsub(pattern = "^p",replacement = "",x = prop.type)) %>%
  filter(est.type == prop2)

#a bit more naming
head(tmp1)
tmp1$real.type2[tmp1$real.type == "FS"] <- "Full-sibs"
tmp1$real.type2[tmp1$real.type == "HS"] <- "Half-sibs"
tmp1$real.type2[tmp1$real.type == "UR"] <- "Unrelated"
tmp1$loci <- paste(tmp1$loci,"loci")
tmp1$par <- paste(tmp1$par,"parents")

#calculating false positive rate
tmp1$fpos <- 1- tmp1$val2
head(tmp1)

tiff(filename = "./analysis/simulations/Post Linux analyis/Output/accuracy.relationship.tiff",width = 8,height = 8,units = "in",res = 300)
ggplot(tmp1, aes(x=real.type2,y=1-fpos))+
  facet_grid(par~loci)+
  geom_boxplot()+
  theme_bw()+
  labs(x="Dyad relationship",y="Accuracy rate\n")+
  theme(axis.text.x = element_text(angle=90,vjust=0.5),
        axis.text = element_text(size=18,color = "black"),
        strip.text = element_text(size=16,color = "black"),
        axis.title = element_text(size=22))
dev.off()

#looking are distribution of known dyads type (HS, FS, UR) sizes among simulations types
head(df)

#some more naming
df$real.type2[df$real.type == "FS"] <- "Full-sibs"
df$real.type2[df$real.type == "HS"] <- "Half-sibs"
df$real.type2[df$real.type == "UR"] <- "Unrelated"
head(df)

tiff(filename = "./analysis/simulations/Post Linux analyis/Output/family.sizes.tiff",width = 8,height = 8,units = "in",res = 300)
ggplot(df, aes(x=par,y=total,color=loci))+
  facet_wrap(~real.type2)+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  labs(x="Number of parents in original breeding matrix",y="Number of known dyads (log10 scale)",color="Number of\nloci")+
  theme(axis.text.x = element_text(angle=90,vjust=0.5),
        axis.text = element_text(size=18,color = "black"),
        strip.text = element_text(size=16,color = "black"),
        axis.title = element_text(size=22))
dev.off()

#fin!
