###########################
### colony.sim.create() ###
###########################

#Description

#Input Parameters:

#defining path to COLONY on local machine
my.colony.path <- "C:/ZSL/Colony/simu2.exe"

#many

colony.sim.create <- function(mat.str, update.alfs = 0, spp.type = 2, selfing.rate = 0, inbreeding = 0,
                              alf.dist, ploidy = 2, fem.gamy = 0, mal.gamy = 0, clone = 0, sib.scale = 0,
                              sib.prior = 0, known.alfs = 0, n.reps = 1, run.reps = 1, run.length = 3,
                              monitor = 0, windows.version = 0, likelihood.type = 1, likelihood.precision = 3,
                              prob.mom = 0, prob.dad = 0, n.mating.str = 1, nmarkers = 10, gt.err = 0.02,
                              mut.err = 0.0001, pr.miss.gt = 0, marker.type = 0, n.alleles = 2, allele.dist = 3,
                              map.length = -1, inbrd.coef = 0, pat.sib.size = 1, mat.sib.size = 1,n.mom = NA, n.dad = NA,
                              delete.sup.files = T){
  
  #############################################################
  ### saving a bunch of objects for lines in the input file ###
  #############################################################
  
  #getting current working directory and fixing slashes for running on linux
  my.project <- paste("'sim'","! Output file name",sep = "\t")
  run.reps <- paste(run.reps,"! Number of replicates",sep = "\t")
  likelihood.type <- paste(likelihood.type,"! 0/1/2 Pair-Likelihood(PLS)/Full likelihood(FL)/FL-PLS",sep = "\t")
  likelihood.precision <- paste(likelihood.precision,"! 0/1/2/3=low/medium/high/very high precision for likelihood",sep = "\t")
  spp.type <- paste(spp.type,selfing.rate,sep = " ")
  spp.type <- paste(spp.type, "! 2/1=Dioecious/Monoecious, Selfing rate for monoecious",sep = "\t")
  n.mating.str <- paste(n.mating.str, "! Seed for random number generator",sep = "\t")
  n.mom <- ncol(mat.str)
  n.dad <- nrow(mat.str)
  npar <- paste(paste(n.dad,n.mom),"! # of males & females in mating structure",sep = "\t")
  mat.str1 <- matrix(data = 0, nrow = n.dad, ncol = n.mom)
  
  #at this time I am stating that fathers and mothers genotypes are not to be included in the simulation
  n.dad.gts <- rep(0,times = n.dad)
  n.mom.gts <- rep(0,times = n.mom)
  ncand <- paste(paste(0,0),"! # of unrelated males & females",sep = "\t")
  prcand <- paste(paste(prob.dad,prob.mom),"! Prob. of males & females are true parents",sep = "\t")
  nloci <- paste(nmarkers, "! Number of loci",sep = "\t") 
  pr.miss.gt <- paste(pr.miss.gt, "! Prob. of missing genotypes",sep = "\t") 
  gt.err <- rep(gt.err,times=nmarkers)
  gt.err <- paste(gt.err,collapse = " ")
  gt.err <- paste(gt.err, "! Allelic dropout rates",sep = "\t") 
  mut.err <- rep(mut.err,times=nmarkers)
  mut.err <- paste(mut.err,collapse = " ")
  mut.err <- paste(mut.err, "! Non-allelic error rates",sep = "\t") 
  marker.type <- rep(marker.type,times=nmarkers)
  marker.type <- paste(marker.type,collapse = " ")
  marker.type <- paste(marker.type, "! Codominant/Dominant (0/1) markers",sep = "\t") 
  n.alleles <- rep(n.alleles,times=nmarkers)
  n.alleles <- paste(n.alleles,collapse = " ")
  n.alleles <- paste(n.alleles, "! # of alleles per locus",sep = "\t") 
  ploidy <- paste(ploidy, "! 1/n=HaploDiploid/n-ploid species",sep="\t")
  gamy <- paste(paste(mal.gamy,fem.gamy),"! 0/1=Monogamy/Polygamy for males & females",sep = "\t")
  random.seed <- round(runif(n = 1,min = 1,max = 9999))
  random.seed <- paste(random.seed, "! Seed for random number generator",sep = "\t")
  sib.prior <- paste(sib.prior,pat.sib.size,mat.sib.size)
  sib.prior <- paste(sib.prior,"! 0/1/2/3=No/Weak/Medium/Strong, R & R (Pat/Mat sibship size)",sep = "\t")
  clone <- paste(clone,"! 0/1=Clone inference =No/Yes",sep="\t")
  sib.scale <- paste(sib.scale,"! 0/1=Scale full sibship=No/Yes",sep = "\t")
  known.alfs <- paste(known.alfs,"! 0/1=Unknown/Known population allele frequency",sep = "\t")
  update.alfs <- paste(update.alfs,"! Not update/update allele frequency",sep = "\t")
  n.reps <- paste(n.reps,"! Number of runs",sep = "\t")
  run.length <- paste(run.length,"! 0/1/2/3/4 VeryShort/Short/Medium/Long/VeryLong Length of run",sep = "\t")
  map.length <- paste(map.length,"! Map length in Morgans",sep = "\t")
  inbrd.coef <- paste(inbrd.coef,"! Inbreeding coefficient of parents", sep = "\t")
  inbreeding <- paste(inbreeding,"! 0/1=No inbreeding/inbreeding", sep = "\t")
  monitor <- paste(monitor,"! 0/1=Monitor method by Iterate#/Time in second",sep = "\t")
  windows.version <- paste(windows.version,"! Windows version",sep = "\t")
  monitor.interval <- paste(100000,"! Monitor interval in Iterate# / in seconds",sep = "\t")
  
  #making sure a input3.Par file doesnt exist here
  my.files <- list.files(pattern = "input3.Par")
  if(length(my.files)>0){stop("input3.Par already exists in this directory. Delete or move it")}
  
  if(allele.dist != 3 ){
    allele.dist <- paste(allele.dist, "! Type of allelic distribution(s)",sep = "\t") 
    #making the actual file
    cat(my.project, run.reps, likelihood.type, likelihood.precision, spp.type, n.mating.str,npar,
        file = "input3.Par",sep = "\n",append = T)
    write.table(mat.str,file = "input3.Par",sep = " ",append = T,row.names = F,col.names = F)
    cat("",file = "input3.Par",sep = "\n",append = T)
    write.table(mat.str1,file = "input3.Par",sep = " ",append = T,row.names = F,col.names = F)
    cat("",file = "input3.Par",sep = "\n",append = T)
    write.table(mat.str1,file = "input3.Par",sep = " ",append = T,row.names = F,col.names = F)
    cat("",file = "input3.Par",sep = "\n",append = T)
    write.table(mat.str1,file = "input3.Par",sep = " ",append = T,row.names = F,col.names = F)
    cat("",file = "input3.Par",sep = "\n",append = T)
    write.table(matrix(data = n.dad.gts,nrow = 1),file = "input3.Par",sep = " ",append = T,row.names = F,col.names = F)
    write.table(matrix(data = n.mom.gts,nrow = 1),file = "input3.Par",sep = " ",append = T,row.names = F,col.names = F)
    cat("",file = "input3.Par",sep = "\n",append = T)
    cat(ncand,prcand, nloci, pr.miss.gt, gt.err, mut.err, marker.type, n.alleles,allele.dist, ploidy,
        gamy, random.seed, sib.prior, clone, sib.scale, known.alfs, update.alfs, n.reps,
        run.length, map.length, inbrd.coef,inbreeding, windows.version, monitor.interval,
        file = "input3.Par",sep = "\n",append = T)
  } else {
    if(is.data.frame(alf.dist)!=T){stop("Allele frequency distrbution is not a data frame.")}
    allele.dist <- paste(allele.dist, "! Type of allelic distrubution(s)",sep = "\t")
    #making the actual file
    cat(my.project, run.reps, likelihood.type, likelihood.precision, spp.type, n.mating.str,npar,
        file = "input3.Par",sep = "\n",append = T)
    write.table(mat.str,file = "input3.Par",sep = " ",append = T,row.names = F,col.names = F)
    cat("",file = "input3.Par",sep = "\n",append = T)
    write.table(mat.str1,file = "input3.Par",sep = " ",append = T,row.names = F,col.names = F)
    cat("",file = "input3.Par",sep = "\n",append = T)
    write.table(mat.str1,file = "input3.Par",sep = " ",append = T,row.names = F,col.names = F)
    cat("",file = "input3.Par",sep = "\n",append = T)
    write.table(mat.str1,file = "input3.Par",sep = " ",append = T,row.names = F,col.names = F)
    cat("",file = "input3.Par",sep = "\n",append = T)
    write.table(matrix(data = n.dad.gts,nrow = 1),file = "input3.Par",sep = " ",append = T,row.names = F,col.names = F)
    write.table(matrix(data = n.mom.gts,nrow = 1),file = "input3.Par",sep = " ",append = T,row.names = F,col.names = F)
    cat("",file = "input3.Par",sep = "\n",append = T)
    cat(ncand,prcand, nloci, pr.miss.gt, gt.err, mut.err, marker.type, n.alleles,allele.dist,
        file = "input3.Par",sep = "\n",append = T)
    write.table(alf.dist[,-1],file = "input3.Par",sep = " ",append = T,row.names = F,col.names = F)
    cat(ploidy, gamy, random.seed, sib.prior, clone, sib.scale, known.alfs, update.alfs, n.reps,
        run.length, map.length, inbrd.coef,inbreeding, windows.version, monitor.interval,
        file = "input3.Par",sep = "\n",append = T)
  }
  
  
  #making sure a COLONY2.dat file doesnt exist here
  my.files <- list.files(pattern = "COLONY2_1.DAT")
  if(length(my.files)>0){stop("input3.Par already exists in this directory. Delete or move it")}
  
  #running the simulation module of COLONY
  suppressWarnings(system(my.colony.path,show.output.on.console = F))
  
  if(delete.sup.files == T){
    #removing extra files
    my.files <- list.files(pattern = "_1.txt")
    suppressWarnings(for(i in 1:length(my.files)){file.remove(my.files)})
    my.files <- list.files(pattern = ".BestConfig")
    for(i in 1:length(my.files)){file.remove(my.files)}
  }

  
  #print
  #print("Input3.Par has been created")
  
}

