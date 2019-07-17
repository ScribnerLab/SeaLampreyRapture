#originally created on: Feb 6th, 2019
#by: Nick Sard

#ABOUT: This script was written to visualize mating behavior among sea lamprey mate pairs
#detected in best configuration pedigrees in Duffin's Creek

#loading libraries
library(tidyverse)
library(SeaLampreyRapture)
data("SeaLampreyRapture")

#adding age informatioon
best.config.dck.age.1$age <- "Age - 1"
best.config.dck.age.2$age <- "Age - 2"
best.config.dck.age.3$age <- "Age - 3"

#combining all together
df <- rbind(best.config.dck.age.1,
            best.config.dck.age.2,
            best.config.dck.age.3)
head(df)

#getting a mate pair column by pasting together the IDs in both parent columns
#age info needs to be in there to count MPs per age class
df$mp <- paste(df$FatherID,df$MotherID,df$age,sep = "@")
df$mp <- gsub(pattern = "par[0-9]_",replacement = "",x = df$mp)
head(df)

#grouping by mp and getting counts - i.e. mate pair reproductive success per age
df1 <- df %>%
  group_by(mp) %>%
  count() %>%
  separate(col = mp,into = c("dad","mom","age"),sep = "@")
head(df1)

#now releveling mom and dad ids - for the figure

#moms
class(df1$mom)
df1$mom1 <- factor(as.numeric(df1$mom),levels = 1:max(as.numeric(df1$mom)))

#dads
class(df1$dad)
df1$dad1 <- factor(as.numeric(df1$dad),levels = 1:max(as.numeric(df1$dad)))

ggplot(df1,aes(x=dad1,y=mom1,fill=n,label=n))+
  facet_wrap(~age,scales = "free",ncol=3) +
  geom_bin2d(color="black")+
  geom_text(size=3)+
  theme_bw()+
  scale_fill_gradient(low = "white", high = "red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=90,vjust=.5,size=10,color = "black"),
        axis.text.y = element_text(size=10,color = "black"),
        strip.text = element_text(size=18),
        axis.title = element_text(size=22))+
  labs(x="Unknown sex 1",y="Unknown sex 2",fill="Number of\noffspring")

#note that sex is not know, so I labelled the axes as such.
#fin!
