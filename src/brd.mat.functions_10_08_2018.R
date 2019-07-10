
brd.mat <- function(moms = 10,dads = 10,min.mates = 1, min.fert = 2500, max.fert = 6500,max.mates = 3){
  
  #making a blank matrix to be filled in
  mat <- matrix(data = 0,nrow = dads,ncol = moms)
  mat
  
  #picking a middle mates
  middle.mates <- sample(min.mates:max.mates,size=1)
  
  #defining mate.probs for females
  possible.mates <- 1:dads
  df <- data.frame(possible.mates)
  df$prob <- (middle.mates^df$possible.mates)*(exp(-middle.mates))/factorial(df$possible.mates)
  df
  
  #defining values for mate parameters
  min.mates <- 1
  max.mates <- dads
  mat.probs <- df$prob
  
  #for loop for mom and dad to ID which mates
  
  ##################
  ### moms first ###
  ##################
  i <- 1
  i <- NULL
  black.listed.rows <- NA
  for(i in 1:ncol(mat)){
    
    #identifying the number of mates for a given mom
    mates <- sample(min.mates:max.mates,size = 1,prob = mat.probs)
    
    #determining which males she can still mate with 
    #that if male as the maximum number of mates allowed 
    #then he can no longer be mated with
    
    #first defining all possible males to mate
    my.rows <- 1:nrow(mat)
    
    #now determining if the blast.listed.rows has information
    #to account for, if so those males are removed
    suppressWarnings( if(!is.na(black.listed.rows)){my.rows <- my.rows[!(my.rows %in% black.listed.rows)]})
    
    #if a female can still mate with at least 1 male
    if(length(my.rows)>0){
      
      #then mates are selected and filled with a 1
      my.rows <- sample(x = my.rows,size = mates)
      mat[my.rows,i] <- 1
      mat
      #checking to see if more males should be added to the black list
      
      #first checking to see which at the max mates
      blr <- which(rowSums(mat)>=max.mates)
    
      if(length(blr) > 0){
        black.listed.rows <- blr
      }
    }
  }
  mat
  
  
  #resetting to original values for mate parameters
  min.mates <- NULL
  max.mates <- NULL
  mat.probs <- NULL
  
  #defining mate.probs for males
  possible.mates <- 1:moms
  df <- data.frame(possible.mates)
  df$prob <- (middle.mates^df$possible.mates)*(exp(-middle.mates))/factorial(df$possible.mates)
  df
  
  #defining original values for stuff
  min.mates <- 1
  max.mates <- moms
  mat.probs <- df$prob
  
  #################
  ### dads next ###
  #################
  i <- 1
  i <- NULL
  for(i in 1:nrow(mat)){
    current.mates <- sum(mat[i,])
    
    if(current.mates == 0){
      mates <- sample(min.mates:max.mates,size = 1,prob = mat.probs)
      my.cols <- sample(x = 1:ncol(mat),size = mates)
      mat[i,my.cols] <- 1
      mat
    } 
  }  
  mat
  #now need to fill in those successful mate pairs with 
  #some number of viable offpsring
  i <- 1
  i <- NULL
  for(i in 1:ncol(mat)){
    my.mates <- NULL
    my.mates <- which(mat[,i]>0)
    my.mates
    my.off <- sample(x = min.fert:max.fert,size = 1)
    mate.props <- (rev(1:length(my.mates))/(length(my.mates)+1)) /  sum(rev(1:length(my.mates))/(length(my.mates)+1))
    my.mates <- sample(x = my.mates,size = length(my.mates),replace = F)
    mat[my.mates,i] <- round(mate.props*my.off)
    mat
  }
  mat
  return(mat)
}

round(rev(1:2)/(2+1)/sum(rev(1:2)/(2+1)) * 5524)
mat.sub.sample <- function(mat,noff = 250){
  

 #making some generic dad and mom names and converting mat to DF
 dads <- paste0("dads",1:nrow(mat)) 
 moms <- paste0("moms",1:ncol(mat))
 mat <- as.data.frame(mat)
 colnames(mat) <- moms
 mat <- data.frame(cbind(dads,mat),stringsAsFactors=F)
 mat$dads <- as.character(mat$dads)
 
 #making the breeding matrix into mate pair RS data frame (long form)
 ped1 <- gather(data = mat,key = "moms",value = "off",-dads) %>% filter(off != 0)
 
 #identifying each mate pair (mp) my number
 ped1$mp <- 1:nrow(ped1)
 
 #using a for loop to make a vecter the size of the sum(ped1$off)
 #with each mp number replicated the number of times they have fertilized eggs
 i <- 1
 ped <- NULL
 for(i in 1:nrow(ped1)){
   ped2 <- rep(x=ped1$mp[i],times = ped1$off[i])
   ped <- c(ped,ped2)
 }
 
 #randomly sub-sampling ped the size of noff
 my.picks <- sample(x = ped,size = noff,replace = F)
 #getting counts of each mp left
 my.picks <- as.data.frame(table(my.picks),stringsAsFactors = F)
 names(my.picks) <- c("mp","off1")
 
 #merging with ped1 and making sure to fill in all mp lost with zero
 ped1 <- merge(x = ped1,y = my.picks,by = "mp",all.x =T)
 ped1$off1[is.na(ped1$off1)] <- 0
 ped1 <- ped1 %>% arrange(mp)
 return(ped1)
}

mat.stats <- function(mat){
  
  #calculating the stats
  if(ncol(mat)>2){
    female.rs <- round(mean(colSums(mat[,-1])),2)
  } else{
    female.rs <- sum(mat[,-1])
  }
 
  if(nrow(mat)>1){
    male.rs <- round(mean(rowSums(mat[,-1])),2)    
  } else {
    male.rs <- sum(mat[1,])
  }

  mp.rs <- round(mean(as.matrix(mat[,-1])),3)
  mat.str2 <- mat[,-1]
  mat.str2[mat.str2 > 0] <- 1
  mp.count <- sum(mat.str2)
  females.mates <- round(mean(colSums(mat.str2)),2)
  male.mates <- round(mean(rowSums(mat.str2)),2)
  max.female.mates <- max(colSums(mat.str2))
  max.male.mates <- max(rowSums(mat.str2))
  n.mom <- ncol(mat)
  n.dad <- nrow(mat)
  n.par <- n.mom+n.dad
  sex.ratio <- round(n.dad/n.mom,2)
  npar <- n.mom + n.dad
  noff <- sum(mat)
  
  #making a DF to return back
  mat.info <- data.frame(n.par,n.mom,n.dad,sex.ratio, noff,
                          female.rs,females.mates,
                          male.rs,male.mates,mp.count,mp.rs)
  return(mat.info)
}

mat2ped <- function(ped){
  
  #makin sure there are not zeros
  ped <- ped[ped$off1 != 0,]
  
  #making a generic pedigree with the remaining off
  i <- 1
  i <- NULL
  ped.out <- NULL
  for(i in 1:nrow(ped)){
    off1 <- paste(paste0("o",1:ped$off1[i]),ped$moms[i],ped$dads[i],sep="_")  
    off1 <- gsub(pattern = "moms",replacement = "m",x = off1)
    off1 <- gsub(pattern = "dads",replacement = "d",x = off1)
    ped2 <- data.frame(off = off1, mom = ped$moms[i], dad = ped$dads[i])
    ped.out <- rbind(ped.out,ped2)
  }
  
  return(ped.out)
}


ped2mat <- function(ped){
  
  #identifying matepairs
  ped$mp <- paste(ped$dad,ped$mom,sep = "_")
  
  #getting counts for each successful pair
  ped1 <- as.data.frame(table(ped$mp))
  
  #making the matrix
  mat <- ped1 %>% separate(col = Var1,into = c("dads","moms"),sep = "_") %>% spread(key = moms,value = Freq,fill = 0)

  #want to recreate original breeding mat
  
  #re-ordering cols
  cols <- colnames(mat)[-1]
  cols <- data.frame(cols = cols, order = as.numeric(gsub(pattern = "moms",replacement = "",x = cols)))
  cols <- cols %>% arrange(order) %>% select(cols) %>% pull() %>% as.character()
  cols <- c("dads",cols)
  mat <- mat[,cols]
  
  #re-ordering rows
  rows <- data.frame(rows = row.names(mat), order = as.numeric(gsub(pattern = "dads",replacement = "",x = mat$dads)))
  rows <- rows %>% arrange(order) %>% select(rows) %>% pull() %>% as.character() %>%as.numeric()
  mat <- mat[rows,]
  
  #making it a matrix
  row.names(mat) <- mat$dads
  mat$dads <- NULL
  return(mat)
}