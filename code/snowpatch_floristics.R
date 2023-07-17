###
### Code for manuscript
###
### Title: Early-melting snowpatch plant communities are transitioning into novel ecosystems
### Authors: John Morgan and Zac Walker 
### Year: 2023
### Script: snowpatch_floristics
### About: This script is for nmds ordinations
###

rm(list = ls()) # clear environment

### Load packages ------------------------------------

library(dplyr)
library(vegan)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(cowplot)
library(ecodist)

### Load and clean data ------------------------------------

#floristics
scov <- read.csv("data/snowp_floristics.csv",
                      header=T)
names(scov)[1] <- "site"
names(scov)[2] <- "date"
names(scov)[3] <- "year"
names(scov)[4] <- "quad"

#env data
sadd <- read.csv("data/snowp_additional.csv",
                 header=T)

#join datasets
scov <- as_tibble(left_join(sadd, scov))

#check
head(scov)

#remove singletons
scov <- scov %>% select(!c("Gentianella.muelleriana", "Geranium.spp.","Poa.helmsii")) 

#make a presence absence version of the floristic data
spa <- ifelse(scov[17:87]>0,1,0)
spa2 <- scov[1:16]
spa <- cbind(spa2, spa)

#standardise floristic abundance data
scover <- decostand(scov[17:87], method = "hel")

#run ordinations - one for pres ab and one for abundance. k = 3 becaus k = 2 stress is too high.
ord.pa <- metaMDS(spa[17:87], distance = "jaccard", k=3, 
                  trymax=9999, autotransform = TRUE, expand=FALSE, 
                  trace=FALSE, plot=FALSE, binary=T)

ord.cov <- metaMDS(scover, distance = "bray", k=3,
                   trymax=9999, autotransform = TRUE, expand=FALSE,
                   trace=FALSE, plot=FALSE)


## check stress of ordination
ord.pa #stress = 0.18
ord.cov #stress = 0.16

## produce shepard plots to check fit of ordination
stressplot(ord.pa) #non-metric fit R2 = 0.965; linear fit R2 = 0.741
title(main = "Shepard plot of presence-absence NMDS")
stressplot(ord.cov) #non-metric fit R2 = 0.974; linear fit R2 = 0.809
title(main = "Shepard plot of cover NMDS")



#PRESENCE ABSENCE
#extract NMDS scores (x and y coordinates) and add columns to data frame 
data.scores.pa = as.data.frame(scores(ord.pa$points))
data.scores.pa$site = spa$site
data.scores.pa$year = spa$year
data.scores.pa$shrub = spa$shrubcov
data.scores.pa$drygrass = spa$drygrass
data.scores.pa$year = as.factor(data.scores.pa$year) # make year factor
head(data.scores.pa)

#COVER

#extract NMDS scores (x and y coordinates) and add columns to data frame 
data.scores.cov = as.data.frame(scores(ord.cov$points))
data.scores.cov$site = scov$site
data.scores.cov$year = scov$year
data.scores.cov$shrub = scov$shrubcov
data.scores.cov$drygrass = scov$drygrass
data.scores.cov$year = as.factor(data.scores.cov$year) # make year factor
head(data.scores.cov)

################
#vector fitting
vec.sp.pa<-envfit(ord.pa$points, spa[4:5], perm=9999, choices=c(1,2,3))
vec.sp.df.pa<-as.data.frame(vec.sp.pa$vectors$arrows*sqrt(vec.sp.pa$vectors$r))
vec.sp.df.pa$species<-rownames(vec.sp.df.pa)

vec.sp.cov<-envfit(ord.cov$points, scov[4:5], perm=9999,choices=c(1,2,3))
vec.sp.df.cov<-as.data.frame(vec.sp.cov$vectors$arrows*sqrt(vec.sp.cov$vectors$r))
vec.sp.df.cov$species<-rownames(vec.sp.df.cov)

###############
#plotting

a = ggplot(data.scores.pa, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 4, alpha = 0.6, aes(colour = year))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "white", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())+
  geom_segment(data=vec.sp.df.pa,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black") + 
  geom_text(data=vec.sp.df.pa,aes(x=MDS1,y=MDS2,label=species),size=5,
            nudge_x = -0.005, nudge_y = -0.05)+
  ggtitle("Presence-absence data")+
  labs(x = "NMDS1", y = "NMDS2")+
  xlim(-0.85,0.85) 

a

b = ggplot(data.scores.pa, aes(x = MDS1, y = MDS3)) + 
  geom_point(size = 4, alpha = 0.6, aes(colour = year))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  geom_segment(data=vec.sp.df.pa,aes(x=0,xend=MDS1,y=0,yend=MDS3),
               arrow = arrow(length = unit(0.5, "cm")),colour="black") + 
  geom_text(data=vec.sp.df.pa,aes(x=MDS1,y=MDS3,label=species),size=5,
            nudge_x = -0.005, nudge_y = 0.07)+
  ggtitle("Presence-absence data")+
  labs(x = "NMDS1", y = "NMDS3")+
  xlim(-0.85,0.85)

b

c = ggplot(data.scores.cov, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 4, alpha = 0.6, aes(colour = year))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())+
  geom_segment(data=vec.sp.df.cov,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black")+
  geom_text(data=vec.sp.df.cov,aes(x=MDS1,y=MDS2,label=species),size=5,
            nudge_x = -0.005, nudge_y = -0.05)+
  ggtitle("Cover data")+
  labs(x = "NMDS1", y = "NMDS2")+
  xlim(-0.85,0.85)
  
c

d = ggplot(data.scores.cov, aes(x = MDS1, y = MDS3)) + 
  geom_point(size = 4, alpha = 0.6, aes(colour = year))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  geom_segment(data=vec.sp.df.cov,aes(x=0,xend=MDS1,y=0,yend=MDS3),
               arrow = arrow(length = unit(0.5, "cm")),colour="black")+
  geom_text(data=vec.sp.df.cov,aes(x=MDS1,y=MDS3,label=species),size=5,
            nudge_x = -0.005, nudge_y = -0.05)+
  ggtitle("Cover data")+
  labs(x = "NMDS1", y = "NMDS3")+
  xlim(-0.85,0.85)

d

#combine all four into one figure
ggarrange(a,b,c,d, nrow = 2, ncol = 2, labels = "AUTO",
          common.legend = TRUE, legend = "right")

#####################################
  
  
  
