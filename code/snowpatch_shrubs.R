###
### Code for manuscript
###
### Title: Early-melting snowpatch plant communities are transitioning into novel ecosystems
### Authors: John Morgan and Zac Walker 
### Year: 2023
### Script: snowpatch_shrubs
### About: This script is for producing size class distributions and condit curves
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


### CONDIT CURVES  ------------------------------------
###
###

shrubs <- read.csv("data/snowpatch_condit.csv") #max shrub size 99
shrubs <- shrubs[-1]

xy90 <- 
  shrubs %>% group_by(year, sp) %>% 
  summarise(n = n()) %>% 
  filter(year == "2019") %>% 
  mutate(percent = (n/6352)*100) %>% 
  filter(percent > 1)

specie90 <- xy90$sp #get list of species
df90 <- shrubs %>% filter(sp %in% specie90) %>% 
  filter(year == "1990")

df <- shrubs %>% filter(sp %in% specie90) %>% 
  filter(year == "2019")

#filter shrubs to just sp of interest
Acrothamnus90 <- df90 %>% filter(sp2 == "Acrothamnus")
Asterolasia90 <- df90 %>% filter(sp2 == "Asterolasia")
Melicytus90 <-  df90 %>% filter(sp2 == "Melicytus")
OleariaB90  <- df90 %>% filter(sp == "Olearia brevipedunculata")
OleariaF90  <- df90 %>% filter(sp == "Olearia frostii")
OleariaP90 <- df90 %>% filter(sp == "Olearia phlogopappa")
Phebalium90 <- df90 %>% filter(sp2 == "Phebalium")
Pimelea90 <- df90 %>% filter(sp2 == "Pimelea")

Acrothamnus19 <- df %>% filter(sp2 == "Acrothamnus")
Asterolasia19 <- df %>% filter(sp2 == "Asterolasia")
Melicytus19 <-  df %>% filter(sp2 == "Melicytus")
OleariaB19  <- df %>% filter(sp == "Olearia brevipedunculata")
OleariaF19  <- df %>% filter(sp == "Olearia frostii")
OleariaP19 <- df %>% filter(sp == "Olearia phlogopappa")
Phebalium19 <- df %>% filter(sp2 == "Phebalium")
Pimelea19 <- df %>% filter(sp2 == "Pimelea")

breaks <- c(1,5,10,15,20,30,40,50,100)# set up cut-off values 
tags <- c("[1-5)", "[5-10)","[10-15)",
          "[15-20)", "[20-30)","[30-40)","[40-50)", "[50-100)")# specify interval/bin labels

# bucketing values into bins for each species
#2019
bin_ac <- cut(Acrothamnus19$wdt_1, 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE, 
              labels=tags)

bin_as <- cut(Asterolasia19$wdt_1, 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE, 
              labels=tags)

bin_m <- cut(Melicytus19$wdt_1, 
             breaks=breaks, 
             include.lowest=TRUE, 
             right=FALSE, 
             labels=tags)

bin_ob <- cut(OleariaB19$wdt_1, 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE, 
              labels=tags)

bin_of <- cut(OleariaF19$wdt_1, 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE, 
              labels=tags)

bin_op <- cut(OleariaP19$wdt_1, 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE, 
              labels=tags)

bin_ph <- cut(Phebalium19$wdt_1, 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE, 
              labels=tags)


bin_pim <- cut(Pimelea19$wdt_1, 
               breaks=breaks, 
               include.lowest=TRUE, 
               right=FALSE, 
               labels=tags)

###
#1990
# bucketing values into bins for each species
bin_ac90 <- cut(Acrothamnus90$wdt_1, 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE, 
              labels=tags)

bin_as90 <- cut(Asterolasia90$wdt_1, 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE, 
              labels=tags)

bin_m90 <- cut(Melicytus90$wdt_1, 
             breaks=breaks, 
             include.lowest=TRUE, 
             right=FALSE, 
             labels=tags)

bin_ob90 <- cut(OleariaB90$wdt_1, 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE, 
              labels=tags)

bin_of90 <- cut(OleariaF90$wdt_1, 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE, 
              labels=tags)

bin_op90 <- cut(OleariaP90$wdt_1, 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE, 
              labels=tags)

bin_ph90 <- cut(Phebalium90$wdt_1, 
              breaks=breaks, 
              include.lowest=TRUE, 
              right=FALSE, 
              labels=tags)

bin_pim90 <- cut(Pimelea90$wdt_1, 
               breaks=breaks, 
               include.lowest=TRUE, 
               right=FALSE, 
               labels=tags)

bin_combine <- data.frame(Acrothamnus_montanus.90 = summary(bin_ac90),
                          Asterolasia_trymalioides.90 = summary(bin_as90),
                          Melicytus_dentatus.90 = summary(bin_m90),
                          Olearia_brevipedunculata.90 = summary(bin_ob90),
                          Olearia_frostii.90 = summary(bin_of90),
                          Olearia_phlogopappa.90 = summary(bin_op90),
                          Phebalium_squamulosum.90 = summary(bin_ph90),
                          Pimelea_axiflora.90 = summary(bin_pim90),
                          Acrothamnus_montanus.19 = summary(bin_ac),
                          Asterolasia_trymalioides.19 = summary(bin_as),
                          Melicytus_dentatus.19 = summary(bin_m),
                          Olearia_brevipedunculata.19 = summary(bin_ob),
                          Olearia_frostii.19 = summary(bin_of),
                          Olearia_phlogopappa.19 = summary(bin_op),
                          Phebalium_squamulosum.19 = summary(bin_ph),
                          Pimelea_axiflora.19 = summary(bin_pim),
                          wdth = c(4,5,5,5,10,10,10,50))

bins <- bin_combine %>% 
  mutate(across(everything(), ~ . / wdth)+1) %>% 
  select(!wdth) %>% 
  mutate(across(everything(), log10))

bins$di = log10(c(3,7.5,12.5,17.5,25,35,45,75))

#make it a long dataframe for easier plotting
bins_long <- bins %>% 
  pivot_longer(1:16,names_to = "species", values_to = "number") %>% 
  separate_wider_delim(cols = "species", delim = ".", names = c("species","year"))

#### get lm slope coef's for plotting
df <- bins_long #make new object for function

#create function
regression=function(df){ 
  reg_fun<-lm(formula=df$number~df$di)
  slope<-round(coef(reg_fun)[2],2)
}
regressions_data<-plyr::ddply(bins_long,c("species", "year"),regression) #apply regression across species to get slope for each
colnames(regressions_data)<-c ("species","year","slope")#fix colnames
arrange(regressions_data, year, slope)

bins_long$year[bins_long$year == "19"] <- "2019"
bins_long$year[bins_long$year == "90"] <- "1990"

#plot regression lines (sorting facet_wrap from most negative slope to most positive)
condit <- ggplot(bins_long, aes(x = di, y = number, group = year))+
  facet_wrap(~factor(species, levels=c('Acrothamnus_montanus',
                                       'Asterolasia_trymalioides',
                                       'Melicytus_dentatus',
                                       'Olearia_brevipedunculata',
                                       'Olearia_frostii',
                                       'Olearia_phlogopappa',
                                       'Phebalium_squamulosum',
                                       'Pimelea_axiflora'),
                     labels = c('Acrothamnus montanus',
                                'Asterolasia trymalioides',
                                'Melicytus angustifolius',
                                'Olearia brevipedunculata',
                                'Olearia frostii',
                                'Olearia phlogopappa',
                                'Phebalium squamulosum',
                                'Pimelea axiflora')))+
  ylim(-.5,3.5)+
  geom_point(aes(colour = year))+
  geom_smooth(method = lm, se =F, aes(colour = year))+ 
  #geom_label(data=regressions_data,inherit.aes=FALSE,
  #           aes(x = 2, y = 0,label=paste("slope:",round(slope,2))))+
  geom_text(data=filter(regressions_data, year == "90"),inherit.aes=FALSE, size = 5, colour = "#F8766D",
             aes(x = 1.5, y = 3.3, label=paste("slope:",slope)))+
  geom_text(data=filter(regressions_data, year == "19"),inherit.aes=FALSE, size = 5, colour = "#00BFC4",
            aes(x = 1.5, y = 2.6, label=paste("slope:",slope)))+
  xlab("Adjusted size class - log10(di)")+
  ylab("Adjusted number of individuals - log10(ni + 1)")+
  theme_cowplot(16) +
  theme(strip.text = element_text(face = "italic"))

condit

ggsave("shrub2.png", plot = condit, width = 240, height = 180, units = "mm", bg = "white")






