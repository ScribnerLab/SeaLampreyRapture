#originally created on: August 21, 2017
#by: Nick Sard

#ABOUT: This script was written to create a input file (Input3.Par) for the simulation executable within COLONY.

# This script is written to allow the user to simulate a broad range of simulation parameters. In the context of the
# sea lamprey paper, we focused on a few discrete set of parameters for simplicity. We held the number of parents constant
# to 10, 100, or 1000. In addition, we held the number of loci  genotyped constant to 100, 200, or 500. Finally, we held
# the number of offspring constant to 100 throughout all 9 combinations of the above parameters.
# We also allowed the sex ratio to vary uniformally from 1 to 2. The number of offspring produced by a female was also
# allowed to vary from 25000 to 10000.
# See the "Setting variables" section below for more information.

# The script below assumes a new set of output directories are created each time the script is run for each combination of
# parents and loci simulated. It provides 100 simulated datasets per combination of parents and loci. There are four files
# created for each simulation. The actual allele frequencies used to simulate offspring (allele.dist), the original breeding
# matrix (brd.mat), the simulation input file (Input3.Par), and the actual colony input file (colony.dat). Each file is saved
# with the simulation number (1-100).

#Assumptions:
#1) Source scripts with personalized functions are loaded. The breeding matrix functions (brd.mat.function.R)
#   enable the creation of "full" breeding matrix and subsampling of the breeding matrix to some number of genotyped
#   individuals. There are functions to move from a breeding matrix to a pedigree (something similar to a Best Configuration
#   file from COLONY). Finally, there is a function to calculate a suite of summary statistics from a breeding matrix
#2) The script is really about assembling a Input3.Par file, the standard input file for the simulation model of COLONY. The
#   script than provides the assembled InputPar3.Par file to the COLONY simulation executable. Currently, it assumes that
#   the the simulation executable is in the default directory for COLONY "C:/ZSL/Colony/simu2.exe"; however, the user can
#   change the directory by changing [my.colony.path] in the colony.sim.create.R source file.

#######################################
### Creating generic file structure ###
#######################################

#creating the set of directories necessary for simulations to save files to
#again, must be changed each time a new set of parameters below is changed
dir.create("./analysis/pedigrees/Simulations/Create simulations/Output/")
dir.create("./analysis/pedigrees/Simulations/Create simulations/Output/sims")
dir.create("./analysis/pedigrees/Simulations/Create simulations/Output/sims/Npar10_Nloci100")
dir.create("./analysis/pedigrees/Simulations/Create simulations/Output/sims/Npar10_Nloci100/breeding_matrices")
dir.create("./analysis/pedigrees/Simulations/Create simulations/Output/sims/Npar10_Nloci100/alf_dist")
dir.create("./analysis/pedigrees/Simulations/Create simulations/Output/sims/Npar10_Nloci100/colony_dat")
dir.create("./analysis/pedigrees/Simulations/Create simulations/Output/sims/Npar10_Nloci100/input_par")

#making some basic folder names
#to be used later when saving files to specific output folders
sim.folder <- "./analysis/pedigrees/Simulations/Create simulations/Output/sims/Npar10_Nloci100/"
mat.folder <- "./analysis/pedigrees/Simulations/Create simulations/Output/sims/Npar10_Nloci100/breeding_matrices/"
alf.folder <- "./analysis/pedigrees/Simulations/Create simulations/Output/sims/Npar10_Nloci100/alf_dist/"
dat.folder <- "./analysis/pedigrees/Simulations/Create simulations/Output/sims/Npar10_Nloci100/colony_dat/"
par.folder <- "./analysis/pedigrees/Simulations/Create simulations/Output/sims/Npar10_Nloci100/input_par/"

#loading in libraries
library(tidyverse)
library(SeaLampreyRapture)
data("SeaLampreyRapture")

#to stop scientific notation in COLONY input files
options(scipen=999)

#source scripts
source("./R/brd.mat.functions.R")
source("./R/colony.sim.create.R")

#I want to simulate a set of offspring genotypes based on allele frequencies observed in my data, so
#I am going to read in a file containing the allele frequencies and randomly draw some number (100, 200, or 500)
#of allele frequencies from it. Mine is from the discovery RAD VCF from STACKs. The minor allele fruency (q) is
#used to calculate p from the 3460 rad_tags

#view the allele frequency distrubution file
head(alf.dist)

#########################
### setting variables ###
#########################

# To simulate data using the COLONY simulation module the user needs several pieces of information; an essential component is
# the breeding matrix. Literally, a matrix containing some number of columns and rows representing mothers and fathers,
# respectively. The cells within the matrix represent the number of offspring produced by all possiblem mother/father
# combinations. Most of what occurs below is about making that breeding matrix and tracking associated information

#Minimum and maximum expected number of parents spawning in study area
min.parents <- 10
max.parents <- 10
#for the sea lamprey paper we evaluated offspring produced by 10, 100, 1000 parents to evaluate the power to infer
#sibship across a broad range of parameters.

#Maximum expected number of offspring collected using specific gear
maxoff <- 100
minoff <- 100
# for the sea lamprey paper, we collected small number of offspring ~100 for some age classes. Here, I set
#all simulations to 100 genotyped offspring

#Minimum and maximum expected sex ratio - Based on Duong et al. 2013 MEC
min.sex.ratio <- 1
max.sex.ratio <- 2

#Minimum and maximum expected number of fertilized eggs per female
min.fert <- 25000
max.fert <- 100000
#each female has some number of eggs that are fertilized. Based on XXX, we chose to set the minimum and maximum values to
#the written above.

#setting minimun and maximum bounds for mean number of mate pairs
min.mates <- 3
max.mates <- 3
#center of the poission distribution below is 3 mates per female, based on XXX.

#defining the number of loci that the user wants to be genotyped
nloci <- 100

#creating an object to store information about parameters specific for each simulation
sim.info <- NULL

#now for the simulations

#setting the desired number of simulated datasets
nsims <- 100
k <- 1
k <- NULL
for(k in 1:nsims){

  ############################
  ### setting output names ###
  ############################

  mating.str.name <- paste0(mat.folder,"mating_matrix_sim_",k,".txt")
  sim.input.name <- paste0(par.folder,"input3_sim_",k,".Par")
  colony.input.name <- paste0(dat.folder,"colony2_sim_",k,'.Dat')
  alf.dist.name <- paste0(alf.folder,"alf.dist_sim_",k,'.txt')

  ####################################################
  ### Clearing tmp varaibles before each simuation ###
  ####################################################
  n.mom <- NULL
  n.dad <- NULL
  n.par <- NULL
  n.off <- NULL
  sex.ratio <- NULL
  sex.ratio2 <- NULL
  picks <- NULL
  mat.str <- NULL
  ms1 <- NULL
  mat.sub <- NULL
  mp.lost <- NULL
  ped2 <-  NULL
  mat2 <- NULL
  ms2 <- NULL

  ######################################################
  ### Create a breeding matrix for entire study area ###
  ######################################################

  print(k)

  #defining the number of moms, dads, and rs distribution (poisson in this case) (20 and 5 mean, sd)

  #picking the number of parents from a uniform distrbution
  n.par <- round(runif(n = 1,min = min.parents, max = max.parents))
  n.par

  # assuming the the sex ratio of successful males to females is NOT 1:1
  sex.ratio <- round(runif(n = 1,min = min.sex.ratio,max = max.sex.ratio),1)
  sex.ratio

  #needed a way to pick of combination of males and females that would exactly (or as close as possible) to
  #the sex ratio that was choosen, here is my solution

  #enumerate all possible combinations of males and females, calculate sex ratio (sr) between them,
  #get the absoluate difference between their sr and the real sr, and pick the one with smallest difference
  picks <- expand.grid(1:n.par,1:n.par) %>%
    mutate(sr =  Var1/Var2) %>%
    mutate(npar = Var1 + Var2) %>%
    mutate(diffs = abs(sr - sex.ratio)) %>%
    mutate(diffnp = abs(npar - n.par)) %>%
    filter(diffnp == min(diffnp)) %>%
    filter(diffs == min(diffs))
  picks

  #sometimes its possible to pick more than one 'pick' - e.g. when the sex ratio is 2,
  #so I put in an if statement to randomly select one of the rows in that case
  if(nrow(picks)== 1){
    n.mom <- picks$Var2
    n.dad <- picks$Var1
    sex.ratio2 <- round(picks$sr,2)
    diff <- picks$diff
  } else {
    warning(paste("There are",nrow(picks),"sex ratios that work. Picking one randomly."))
    my.row <- as.numeric(sample(x = row.names(picks),size = 1))
    picks <- picks[my.row,]
    n.mom <- picks$Var2
    n.dad <- picks$Var1
    sex.ratio2 <- round(picks$sr,2)
    diff <- picks$diff
  }

  #using breding matrix function to create the breeding matrix
  mat.str <- brd.mat(moms = n.mom,dads = n.dad,min.mates = min.mates,max.mates = max.mates,min.fert = min.fert,max.fert = max.fert)
  head(mat.str)

  #getting summary stats
  ms1 <- mat.stats(mat = mat.str)
  ms1$type <- "Before"
  ms1$lost <- 0
  ms1$sim <- k
  ms1

  #######################################
  ### Subsamping full breeding matrix ###
  #######################################

  #randomly selecting the number of offspring that will be sampled in a given sampling effort
  n.off = round(runif(n = 1,min = minoff, max = maxoff))
  n.off

  #creating full pedigree with that breeding matrix, and sub-sampling each mate pair randomly
  mat.sub <- mat.sub.sample(mat = mat.str,noff = n.off)
  head(mat.sub)
  mat.sub

  #recording how many are lost
  mplost <- nrow(mat.sub[mat.sub$off1 == 0,])
  mplost

  #getting rid of the zeros
  mat.sub <- mat.sub[mat.sub$off1 != 0,]

  #converting to small "complete" genetic pedigree
  ped2 <- mat2ped(ped = mat.sub)
  head(ped2)
  ped2

  #converting back to a breeding matrix to calculate a range of stats
  mat2 <- ped2mat(ped = ped2)
  head(mat2)

  #getting numbers of mates and rs per sex
  ms2 <- mat.stats(mat = mat2)
  ms2$type <- "After"
  ms2$lost <- mplost
  ms2$sim <- k
  ms2

  #saving information associated with each breeding matrix for graphics
  sim.info1 <- rbind(ms1,ms2)
  sim.info <- rbind(sim.info,sim.info1)
  sim.info

  ##################################
  ### Selecting loci of interest ###
  ##################################

  #for this first set of simulations, going for 500 loci - one per chromsome
  head(alf.dist)

  #now selecting 500 loci at random
  alf.dist1 <- alf.dist %>% mutate(picks = runif(n = n())) %>% arrange(picks) %>% head(n=nloci) %>% select(ID:q)
  head(alf.dist1)

  ##########################
  ### Making Colony file ###
  ##########################

  #making the file
  colony.sim.create(mat2, update.alfs = 0, spp.type = 2, selfing.rate = 0, inbreeding = 0,
                  alf.dist = alf.dist1, ploidy = 2, fem.gamy = 0, mal.gamy = 0, clone = 0, sib.scale = 0,
                  sib.prior = 0, known.alfs = 0, n.reps = 1, run.reps = 1, run.length = 3,
                  monitor = 0, windows.version = 0, likelihood.type = 1, likelihood.precision = 3,
                  prob.mom = 0.01, prob.dad = 0.01, n.mating.str = 1, nmarkers = nrow(alf.dist1), gt.err = 0.02,
                  mut.err = 0.0001, pr.miss.gt = 0, marker.type = 0, n.alleles = 2, allele.dist = 3,
                  map.length = -1, inbrd.coef = 0, pat.sib.size = 1, mat.sib.size = 1,n.mom = NA, n.dad = NA,
                  delete.sup.files = T)

  #writing mating matrix to file and modifying some names
  write.table(x = mat2,file = mating.str.name,append = F,quote = F,sep = " ",row.names = F,col.names = F)
  write.table(x = alf.dist1,file =alf.dist.name,append = F,quote = F,sep = " ",row.names = F,col.names = F)
  file.rename(from = "input3.Par",to = sim.input.name)
  file.rename(from = "COLONY2_1.DAT",to = colony.input.name)
} #end of creating simulated datasets and colony input files


#writing each to file for evulation via graphics
write.table(x = sim.info,file = paste0(sim.folder,"simulation.information.txt"),append = T,quote = F,sep = "\t",row.names = F,col.names = T)

#fin!
