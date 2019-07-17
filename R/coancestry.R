#################
### coancestry
#################

#Description
#This fuction takes a pedigree and calculates coancestry based on the number of half- and full-sibs.

#Input Parameters:
#df - a dataframe with three columns in this order offspring IDs, mother IDs,and father IDs 
#info - Default - F, A TRUE(T)/False(F) variable that determines if just the coancestry value is returned (F) or
#       a table with counts of half- and full-sibs are returned as well
#graphic - Default - F, A TRUE(T)/False(F) variable that determines if a graphic of coancestry values among all dyads is
#          displayed. Note: Requires that ggplot2 is installed as a package


## parentage functions
coancestry <- function(df,info=F,graphic=F){
  
  #making sure names are correct
  colnames(df) <- c("OffspringID","MotherID","FatherID")
  
  #first defining a function I will use later
  unique.pairs <- function(ids){
    x <- expand.grid(ids,ids)
    x <- x[ifelse(x$Var1==x$Var2,T,F)==F,]
    x$both <- NULL
    row.names(x) <- NULL
    return(x)
  } # end of unique.pairs
  
  #getting offspring names and making a df of unique pairs
  offs <- df[,1]
  tmp <- unique.pairs(offs)
  tmp
  
  #making a df that I can use logic to determine sibship
  tmp1 <- merge(x = df,y = tmp, by.x = "OffspringID",by.y = "Var1",all.y = T)
  tmp1$mps.y <- paste(tmp1$FatherID,tmp1$MotherID,sep="_")
  names(tmp1)[1] <- "OffspringID.y"
  tmp2 <- suppressMessages(merge(x = df,y = tmp1, by.x = "OffspringID",by.y = "Var2",all.y = T))
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
  head(tmp2)
  
  if(graphic == T){
    require(ggplot2)
    tmp2$sibship <- "UR"
    tmp2$sibship[tmp2$fs == 1]<- "FS"
    tmp2$sibship[tmp2$mat_hs == 1]<- "mHS"
    tmp2$sibship[tmp2$pat_hs == 1]<- "pHS"
    head(tmp2)
    ggplot(tmp2, aes(x = OffspringID.x,y=OffspringID.y,fill=sibship))+
      geom_bin2d(color="black")+
      theme_bw()+
      labs(x="",y="",fill="Sibship")
  }
  
  #calculating coancestry
  full_sibs_n <- sum(tmp2$fs)
  half_sibs_n <- sum(tmp2$mat_hs,tmp2$pat_hs)
  n <- nrow(df)
  nt <- (n^2)-n
  nt2 <- nrow(tmp2)
  output <- ((0.25*full_sibs_n)+(0.125*half_sibs_n))/nt
  output2 <- ((0.25*full_sibs_n)+(0.125*half_sibs_n))/nt2
  
  #making summary
  output1 <- data.frame(type = c("nt","nt2","n_hs","n_fs","coancestry","coancestry2"),
                       value = c(nt,nt2,half_sibs_n,full_sibs_n,output,output2))
  
  #now returning apprportiate information                     
  if(info == F){
    #returning coancestry value
    return(output)
  } 
  
  if(info == T){
    #returning coancestry value
    return(output1)
  }

  } # end of coancestry function