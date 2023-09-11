###
### Code for manuscript
###
### Title: Early-melting snowpatch plant communities are transitioning into novel ecosystems
### Authors: John Morgan and Zac Walker 
### Year: 2023
### Script: snowpatch_shrubs
### About: This script is for producing size class distributions
###

rm(list=ls()) #clear environment

### Load packages ------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggpubr)

### SIZE CLASS DISTRIBUTIONS  ------------------------------------
###
###

### Load and clean data ------------------------------------


shrb <- read.csv("data/snowpatch_shrubdata.csv") #cleaned data including: 7 sites and 6 genera
shrb <- shrb[-1] #remove first column

#create new object for each taxa of interest
Acrothamnus <- shrb %>% 
  filter(sp2 == "Acrothamnus")
Asterolasia <- shrb %>% 
  filter(sp2 == "Asterolasia")
Melicytus <- shrb %>% 
  filter(sp2 == "Melicytus")
OleariaF <- shrb %>% 
  filter(sp == "Olearia frostii")
OleariaB <- shrb %>% 
  filter(sp == "Olearia brevipedunculata")
OleariaP <- shrb %>% 
  filter(sp == "Olearia phlogopappa")
Phebalium <- shrb %>% 
  filter(sp2 == "Phebalium")
Pimelea <- shrb %>% 
  filter(sp2 == "Pimelea")


#plot by diameter for each species
a <- ggplot(Acrothamnus, aes(x = wdt_1)) +
  facet_wrap(vars(year), ncol = 2, scales = "fixed")+
  geom_histogram(binwidth=5)+
  theme_cowplot()+
  xlab("")+
  ylab("# of individuals")+
  labs(title = "Acrothamnus montanus")+
  theme(plot.title = element_text(face = "italic"))

b <- ggplot(Asterolasia, aes(x = wdt_1)) +
  facet_wrap(vars(year), ncol = 2, scales = "fixed")+
  geom_histogram(binwidth=5)+
  theme_cowplot()+
  xlab("")+
  ylab("")+
  labs(title = "Asterolasia trymalioides")+
  theme(plot.title = element_text(face = "italic"))


ff <- ggplot(OleariaF, aes(x = wdt_1)) +
  facet_wrap(vars(year), ncol = 2, scales = "fixed")+
  geom_histogram(binwidth=5)+
  theme_cowplot()+
  xlab("")+
  ylab("# of individuals")+
  labs(title = "Olearia frostii")+
  theme(plot.title = element_text(face = "italic"))

fb <- ggplot(OleariaF, aes(x = wdt_1)) +
  facet_wrap(vars(year), ncol = 2, scales = "fixed")+
  geom_histogram(binwidth=5)+
  theme_cowplot()+
  xlab("")+
  ylab("# of individuals")+
  labs(title = "Olearia brevipedunculata")+
  theme(plot.title = element_text(face = "italic"))


fp <- ggplot(OleariaP, aes(x = wdt_1)) +
  facet_wrap(vars(year), ncol = 2, scales = "fixed")+
  geom_histogram(binwidth=5)+
  theme_cowplot()+
  xlab("")+
  ylab("")+
  labs(title = "Olearia phlogopappa")+
  theme(plot.title = element_text(face = "italic"))

g <- ggplot(Phebalium, aes(x = wdt_1)) +
  facet_wrap(vars(year), ncol = 2, scales = "fixed")+
  geom_histogram(binwidth=5)+
  theme_cowplot()+
  xlab("Shrub diameter (cm)")+
  ylab("# of individuals")+
  labs(title = "Phebalium squamulosum")+
  theme(plot.title = element_text(face = "italic"))

h <- ggplot(Pimelea, aes(x = wdt_1)) +
  facet_wrap(vars(year), ncol = 2, scales = "fixed")+
  geom_histogram(binwidth=5)+
  theme_cowplot()+
  xlab("Shrub diameter (cm)")+
  ylab("")+
  labs(title = "Pimelea axiflora")+
  theme(plot.title = element_text(face = "italic"))

i <- ggplot(Melicytus, aes(x = wdt_1)) +
  facet_wrap(vars(year), ncol = 2, scales = "fixed")+
  geom_histogram(binwidth=5)+
  theme_cowplot()+
  xlab("Shrub diameter (cm)")+
  ylab("# of individuals")+
  labs(title = "Melicytus angustifolius")+
  theme(plot.title = element_text(face = "italic"))

#out all 8 plots into one figure
shrub1 <- ggarrange(a,b,i,fb,ff,fp,g,h, ncol = 3, nrow = 3)

ggsave("shrub1.png", plot = shrub1, width = 170, height = 220, units = "mm", bg = "white")




