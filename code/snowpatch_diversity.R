###
### Code for manuscript
###
### Title: Early-melting snowpatch plant communities are transitioning into novel ecosystems
### Authors: John Morgan and Zac Walker 
### Year: 2023
### Script: snowpatch_diversity
### About: This script is for the diversity analyses
###
###

rm(list = ls()) # clear environment

### Load packages ------------------------------------

library(vegan)
library(tidyr)
library(dplyr)
library(betapart)
library(ggplot2)
library(cowplot)

### Load data ------------------------------------
zzz <- read.csv("data/snowpatch_df.csv")
zzz <- as_tibble(zzz[-1]) #remove index column and make tibble

### Hill diversity indicies ------------------------------------

#calculate hill numbers
hill <- renyi(zzz[-1], hill = TRUE)

#new dataframe with site, species richness, hill-shannon and hill simpsons
ee <- data.frame(
  site = zzz$site,
  hill0 = hill[1],
  hill1 = hill[4],
  hill2 = hill[5])

names(ee)[2] <- "hill0"
names(ee)[3] <- "hill1"
names(ee)[4] <- "hill2"

#seperate site into site name and type (grouping)
ee <- separate_wider_delim(ee, site, delim = "_", names = c("site","type"))

#make type a factor
ee$type <- factor(as.factor(ee$type), levels = c("old","min","max"),
                  labels = c("1982", "2022 (low)", "2022 (high)"))

#check significance with anova and TukeyHSD
ahill0 <- aov(hill0 ~ type, data = ee) %>% TukeyHSD()
ahill1 <- aov(hill1 ~ type, data = ee) %>% TukeyHSD()
ahill2 <- aov(hill2 ~ type, data = ee) %>% TukeyHSD()

#create labels for significance based on TukeyHSD - for plotting in ggplot
ann_text <- data.frame(hill = c("hill0","hill0","hill0",
                                "hill1","hill1","hill1",
                                "hill2","hill2","hill2"),
                       type = as.factor(c("2022 (high)", "2022 (low)", "1982",
                                          "2022 (high)", "2022 (low)", "1982",
                                          "2022 (high)", "2022 (low)", "1982")),
                       value = c(32,32,32,
                                 32,32,32,
                                 32,32,32),
                       lab = c("b","b","a",
                               "b","b","a",
                               "a","a","a"))

#make longer for plotting
eee <- ee %>% gather(3:5, key = "hill", value = "value")

#labels for facet wrap
div_lab <- c(
  `hill0` = "Alpha diversity",
  `hill1` = "Hill-Shannon",
  `hill2` = "Hill-Simpson")

#plot
gg_div <- ggplot(eee, aes(x = type, y = value))+
  facet_wrap(vars(hill), labeller = as_labeller(div_lab))+
  geom_jitter(size = 3, width = 0.1, alpha = 0.5, aes(colour = type))+
  geom_boxplot(fill = "#00000000")+
  theme_cowplot(14)+
  theme(legend.position = "none")+
  ylab("Effective species number")+
  xlab("Survey year")+
  geom_text(data = ann_text, aes(x = type, label = lab))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

gg_div

ggsave("div.png", plot = gg_div, width = 183, units = "mm", bg = "white")


#### Beta diversity ------------------------------------

beta_z<-ifelse(zzz[2:73]>0,1,0)#convert zzz to presence absence
distpair<-beta.pair(beta_z)
bd<-betadisper(distpair[[3]],ee$type) #note: uses ee dataframe from above

anova(bd) # Sig diff between groups
TukeyHSD(bd) # check pairwise
bd$group.distances #check mean group distances

plot(bd,
     hull = FALSE,
     cex = 1.2,
     pch = c(19,19,19),
     ellipse = TRUE,
     main = "",
     sub = "",
     #col = c("red", "green4", "blue"))
     col = c("#DF536B","#61D04F","#2297E6"))



#make dataframe for plotting
beta_distc <- data.frame(group = bd$group,
                         distances = bd$distances)

#make dataframe for annotating Tukey differences
ann_text2 <- data.frame(type = as.factor(c("2022 (high)", "2022 (low)", "1982")),
                        lab = c("b","ab","a"),
                        distances = as.numeric(c("0.5","0.5","0.5")))

#plot
gg_beta <- ggplot(beta_distc, aes(x=group, y = distances))+
  geom_jitter(aes(colour = group), alpha = 0.5,
              size = 3, width = 0.02, height = 0)+
  geom_boxplot(width = 0.3, fill = "#00000000")+
  theme_cowplot()+
  theme(legend.position = "none")+
  xlab("Survey year")+
  ylab("Distance from centroid")+
  geom_text(data = ann_text2, aes(x = type, label = lab))

gg_beta

ggsave("beta.png", plot = gg_beta, width = 92, units = "mm", bg = "white")


#### Spatial Autocorrelation ------------------------------------
library(ade4)
#read in file of coordinates
coords <- read.csv("data/coords.csv")

#filter community dataset (from earlier in this file) to be only 2022 dataset
zzzy <- zzz %>% filter(!grepl("old", site))

#create distance matrix for community data and for coordinate data
comm <- dist(zzzy[2:73])
diff <- dist(coords[2:3])

#run mantel test
man <- mantel.rtest(comm, diff, nrepet = 9999)

#check file
man 

#No autocorrelation: Simulated p-value: 0.2504
