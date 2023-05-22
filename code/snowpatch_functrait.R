###
### Code for manuscript
###
### Title: Early-melting snowpatch plant communities are transitioning into novel ecosystems
### Authors: John Morgan and Zac Walker 
### Year: 2023
### Script: snowpatch_functrait
### About: This script is for functional trait analysis and community weighted means
###

rm(list = ls()) # clear environment

### Load packages ------------------------------------

library(geometry)
library(ape)
library(dplyr)
library(tidyr)

#
#
#
###### FUNCTIONAL TRAIT BY GROUP (e.g. 1982, 2022 low shrub, 2022 high shrub) ######### 
#
#
#


### Load and clean data ------------------------------------
abundances = read.csv("data/snowp_floristics.csv",
                      header=T) #ordinal cover

#fix names
names(abundances)[1] <- "site"
names(abundances)[2] <- "date"
names(abundances)[3] <- "year"
names(abundances)[4] <- "quad"

#load additional data
shrubcover <- read.csv("data/snowp_additional.csv",
                       header=T)
#fix name
names(shrubcover)[1] <- "site"

#join
abund <- left_join(abundances, shrubcover) %>% 
  relocate(site, year, quad, shrubcov, drygrass)

#dataframe of just site [1] and each species [2:74]
abund2 <- abund[,-c(2:17)]

### Split data into three dataframes ------------------------------------
### ..... 1982, 2022 low shrub cover, 2022 high shrub cover....

ab_sh_max <- abund %>% 
  group_by(site) %>% 
  slice_max(order_by = SHRUB, with_ties = FALSE) %>% 
  filter(year == "2022")

ab_sh_min <- abund %>% 
  group_by(site) %>% 
  slice_min(order_by = SHRUB, with_ties = FALSE)

abund2max <- ab_sh_max[,-c(2:16)]
abund2old <- ab_sh_min %>% 
  filter(year != "2022")
abund2old <- abund2old[,-c(2:16)]
abund2min <- ab_sh_min %>% 
  filter(year == "2022")
abund2min <- abund2min[,-c(2:16)]
abund2max$type = "max"
abund2min$type = "min"
abund2old$type = "old"
abund2_all <- rbind(abund2max, abund2min, abund2old)
abund2_all <- abund2_all %>% unite("site", site, type)
head(abund2_all)

### Load trait data  ------------------------------------

traits = read.csv("data/snowp_traits.csv", 
                  header=T, row.names=1)

#reduce to four traits: sla, ldmc, seed mass and height
traits <- traits[1:4]

#names in order to fit with func. code
ab <- abund2_all
ab <- as.data.frame(ab)
rownames(ab) <- ab[,1]
ab <- ab[-1]
obs <- ab
abundances <- ab

#initiate binary version of obs, get size of data matrix
o.bin = obs
samples = length(abundances[,1])
taxa = length(obs[1,])
trait.s = length(traits[1,])
samples #check
taxa #check
trait.s #check

# z transform each trait
for (k in 1:trait.s) {
  mean.t = mean(traits[,k])
  sd.t = sd(traits[,k])
  for (j in 1:taxa) {
    traits[j, k] = (traits[j,k]-mean.t)/(sd.t)
  }
}


# T = number of traits
T<-dim(traits)[2]

# c = number of communities
C<-dim(abundances)[1]

# check coherence of number of species in 'traits' and 'abundances'
if (dim(abundances)[2]!=dim(traits)[1]) stop(" error : different number of species in 'traits' and 'abundances' matrices ")

# check absence of NA in 'traits'
if (length(which(is.na(traits)==T))!=0) stop(" error : NA in 'traits' matrix ")

# replacement of NA in 'abundances' by '0'
abundances[which(is.na(abundances))]<- 0

# definition of vector for results, with communities'names as given in 'abundances'
nbsp<-rep(NA,C) ; names(nbsp)<-row.names(abundances)
FRic<-rep(NA,C) ; names(FRic)<-row.names(abundances)
FEve<-rep(NA,C) ; names(FEve)<-row.names(abundances)
FDiv<-rep(NA,C) ; names(FDiv)<-row.names(abundances)

# scaling and centering of each trait according to all species values
traitsCS<-scale(traits, center=TRUE, scale=TRUE)


############################################################################################################
# loop to compute on each community the three Functional Diversity Indices, plus the number of species

for (i in 1:C)
{
  # selection of species present in the community
  esppres<-which(abundances[i,]>0)
  
  #  number of species in the community
  S<-length(esppres) ; nbsp[i]<-S
  
  # check if
  if ( S<=T) stop(" number of species must be higher than number of traits ")
  
  # filter on 'traits' and 'abundances' to keep only values of species present in the community
  tr<-traitsCS[esppres,] ;  ab<-as.matrix(abundances[i,esppres] )
  
  # scaling of abundances
  abondrel<-ab/sum(ab)
  
  # functional diversity indices
  
  # FRIC
  
  # use of convhulln function
  
  # volume
  FRic[i]<-round(convhulln(tr,"FA")$vol,6)
  
  # identity of vertices
  vert0<-convhulln(tr,"Fx TO 'vert.txt'")
  vert1<-scan("vert.txt",quiet=T)
  vert2<-vert1+1
  vertices<-vert2[-1]
  
  
  # FEve
  
  # computation of inter-species euclidian distances
  distT<-dist(tr, method="euclidian")
  
  # computation of Minimum Spanning Tree and conversion of the 'mst' matrix into 'dist' class
  linkmst<-mst(distT) ; mstvect<-as.dist(linkmst)
  
  # computation of the pairwise cumulative relative abundances and conversion into 'dist' class
  abond2<-matrix(0,nrow=S,ncol=S)
  for (q in 1:S)
    for (r in 1:S)
      abond2[q,r]<-abondrel[q]+abondrel[r]
  abond2vect<-as.dist(abond2)
  
  # computation of EW for the (S-1) segments pour relier les S points
  EW<-rep(0,S-1)
  
  flag<-1
  for (m in 1:((S-1)*S/2))
  {if (mstvect[m]!=0) {EW[flag]<-distT[m]/(abond2vect[m]) ; flag<-flag+1}}
  
  # computation of the PEW and comparison with 1/S-1, finally computation of FEve
  minPEW<-rep(0,S-1)  ;  OdSmO<-1/(S-1)
  for (l in 1:(S-1))
    minPEW[l]<-min( (EW[l]/sum(EW)) , OdSmO )
  FEve[i]<-round( ( (sum(minPEW))- OdSmO) / (1-OdSmO ) ,6)
  
  
  # FDiv
  
  # traits values of vertices
  trvertices<-tr[vertices,]
  
  # coordinates of the center of gravity of the vertices (Gv)
  baryv<-apply(trvertices,2,mean)
  
  # euclidian dstances to Gv (dB)
  distbaryv<-rep(0,S)
  for (j in 1:S)
    distbaryv[j]<-( sum((tr[j,]-baryv)^2) )^0.5
  
  # mean of dB values
  meandB<-mean(distbaryv)
  
  # deviations to mean of db
  devdB<-distbaryv-meandB
  
  # relative abundances-weighted mean deviation
  abdev<-abondrel*devdB
  
  # relative abundances-weighted mean of absolute deviations
  ababsdev<-abondrel*abs(devdB)
  
  # computation of FDiv
  FDiv[i]<-round( (sum(abdev)+meandB) / (sum(ababsdev)+meandB) ,6)
  
} 

# end of i


res<-list(nbsp=nbsp, FRic=FRic,FEve=FEve, FDiv=FDiv)
invisible(res)

res2 <- as.data.frame(res)
res2$site_type <- rownames(res2)

res3 <- res2 %>% 
  separate(site_type, c("site", "type"), sep = "_") %>% 
  relocate(site, type)

ab_sh_max2 <- ab_sh_max %>% 
  select(site, quad, shrubcov, drygrass)
ab_sh_max2$type = "max"

ab_sh_min2 <- ab_sh_min %>% 
  filter(year == "2022") %>% 
  select(site, quad, shrubcov, drygrass)
ab_sh_min2$type = "min"

ab_sh_old2 <- ab_sh_min %>% 
  filter(year != "2022") %>% 
  select(site, quad, shrubcov, drygrass)
ab_sh_old2$type = "old"

abshh <- rbind(ab_sh_old2, ab_sh_max2, ab_sh_min2)

res4 <- left_join(abshh, res3) %>% 
  gather(key = "trait", value = "value", 7:9)

res4$type <- factor(as.factor(res4$type), levels = c("old","min","max"),
             labels = c("1982", "2022 (low)", "2022 (high)"))


#check if any func. are significantly different between groups - none are
summary(aov(value ~ type, data = filter(res4, trait == "FEve")))
summary(aov(value ~ type, data = filter(res4, trait == "FRic")))
aov(value ~ type, data = filter(res4, trait == "FRic")) %>% TukeyHSD()
summary(aov(value ~ type, data = filter(res4, trait == "FDiv")))



ann_text4 <- data.frame(type = as.factor(c("2022 (high)", "2022 (low)", "1982",
                                           "2022 (high)", "2022 (low)", "1982",
                                           "2022 (high)", "2022 (low)", "1982")),
                        value = c(0.95,0.95,0.95,
                                  0.95,0.95,0.95,
                                  21,21,21),
                        trait = c("FDiv", "FDiv", "FDiv",
                                  "FEve","FEve", "FEve",
                                  "FRic","FRic","FRic"),
                        lab = c("a","a","a",
                                "a","a","a",
                                "b","ab","a"))


#create boxplots of functional trait metrics by group type --------------
ft_1 <- ggplot(res4, aes(x = type, y = value)) +
  facet_wrap(vars(trait), scales = "free_y",
             strip.position = "left")+
  geom_jitter(aes(colour = type), size = 3, alpha = 0.5, width = 0.1)+
  geom_boxplot(fill = "#00000000")+
  theme_cowplot()+
  theme(legend.position="none",
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_text(size = 16))+
  xlab("Survey year")+
  ylab(NULL)+
  ggtitle("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_text(data = ann_text4, aes(x = type, y = value, label = lab))

ft_1

ggsave("ft1.png", plot = ft_1, width = 183, height = 80, units = "mm", bg = "white")


############
############
############
############

#
#
#
###### FUNCTIONAL TRAIT BY QUADRAT for plotting ######### 
#
#
#


abundances = read.csv("data/snowp_floristics.csv",
                      header=T)

names(abundances)[1] <- "site"
names(abundances)[2] <- "date"
names(abundances)[3] <- "year"
names(abundances)[4] <- "quad"

shrubcover <- read.csv("data/snowp_additional.csv",
                       header=T)
names(shrubcover)[1] <- "site"

abund <- left_join(abundances, shrubcover) %>% 
  relocate(site, year, quad, shrubcov, drygrass) %>% 
  unite("site_quad", site, quad)

abund2 <- abund[,-c(2:15)]

traits = read.csv("data/snowp_traits.csv", 
                  header=T, row.names=1)

traits <- traits[1:4]
head(traits)

#attach(traits)
ab <- abund2
ab <- as.data.frame(ab)
rownames(ab) <- ab[,1]
ab <- ab[-1]
obs <- ab

abundances <- ab

#initiate binary version of obs, get size of data matrix
o.bin = obs
samples = length(abundances[,1])
taxa = length(obs[1,])
trait.s = length(traits[1,])
samples
taxa
trait.s

#if unsure, test for same number of taxa in traits:
taxa2 = length(traits[,1])
taxa2
taxa

# z transform each trait

for (k in 1:trait.s) {
  mean.t = mean(traits[,k])
  sd.t = sd(traits[,k])
  for (j in 1:taxa) {
    traits[j, k] = (traits[j,k]-mean.t)/(sd.t)
  }
}


# T = number of traits
T<-dim(traits)[2]

# c = number of communities
C<-dim(abundances)[1]

# check coherence of number of species in 'traits' and 'abundances'
if (dim(abundances)[2]!=dim(traits)[1]) stop(" error : different number of species in 'traits' and 'abundances' matrices ")

# check absence of NA in 'traits'
if (length(which(is.na(traits)==T))!=0) stop(" error : NA in 'traits' matrix ")

# replacement of NA in 'abundances' by '0'
abundances[which(is.na(abundances))]<- 0

# definition of vector for results, with communities'names as given in 'abundances'
nbsp<-rep(NA,C) ; names(nbsp)<-row.names(abundances)
FRic<-rep(NA,C) ; names(FRic)<-row.names(abundances)
FEve<-rep(NA,C) ; names(FEve)<-row.names(abundances)
FDiv<-rep(NA,C) ; names(FDiv)<-row.names(abundances)

# scaling and centering of each trait according to all species values
traitsCS<-scale(traits, center=TRUE, scale=TRUE)


############################################################################################################
# loop to compute on each community the three Functional Diversity Indices, plus the number of species

for (i in 1:C)
{
  # selection of species present in the community
  esppres<-which(abundances[i,]>0)
  
  #  number of species in the community
  S<-length(esppres) ; nbsp[i]<-S
  
  # check if
  if ( S<=T) stop(" number of species must be higher than number of traits ")
  
  # filter on 'traits' and 'abundances' to keep only values of species present in the community
  tr<-traitsCS[esppres,] ;  ab<-as.matrix(abundances[i,esppres] )
  
  # scaling of abundances
  abondrel<-ab/sum(ab)
  
  # functional diversity indices
  
  # FRIC
  
  # use of convhulln function
  
  # volume
  FRic[i]<-round(convhulln(tr,"FA")$vol,6)
  
  # identity of vertices
  vert0<-convhulln(tr,"Fx TO 'vert.txt'")
  vert1<-scan("vert.txt",quiet=T)
  vert2<-vert1+1
  vertices<-vert2[-1]
  
  
  # FEve
  
  # computation of inter-species euclidian distances
  distT<-dist(tr, method="euclidian")
  
  # computation of Minimum Spanning Tree and conversion of the 'mst' matrix into 'dist' class
  linkmst<-mst(distT) ; mstvect<-as.dist(linkmst)
  
  # computation of the pairwise cumulative relative abundances and conversion into 'dist' class
  abond2<-matrix(0,nrow=S,ncol=S)
  for (q in 1:S)
    for (r in 1:S)
      abond2[q,r]<-abondrel[q]+abondrel[r]
  abond2vect<-as.dist(abond2)
  
  # computation of EW for the (S-1) segments pour relier les S points
  EW<-rep(0,S-1)
  
  flag<-1
  for (m in 1:((S-1)*S/2))
  {if (mstvect[m]!=0) {EW[flag]<-distT[m]/(abond2vect[m]) ; flag<-flag+1}}
  
  # computation of the PEW and comparison with 1/S-1, finally computation of FEve
  minPEW<-rep(0,S-1)  ;  OdSmO<-1/(S-1)
  for (l in 1:(S-1))
    minPEW[l]<-min( (EW[l]/sum(EW)) , OdSmO )
  FEve[i]<-round( ( (sum(minPEW))- OdSmO) / (1-OdSmO ) ,6)
  
  
  # FDiv
  
  # traits values of vertices
  trvertices<-tr[vertices,]
  
  # coordinates of the center of gravity of the vertices (Gv)
  baryv<-apply(trvertices,2,mean)
  
  # euclidian dstances to Gv (dB)
  distbaryv<-rep(0,S)
  for (j in 1:S)
    distbaryv[j]<-( sum((tr[j,]-baryv)^2) )^0.5
  
  # mean of dB values
  meandB<-mean(distbaryv)
  
  # deviations to mean of db
  devdB<-distbaryv-meandB
  
  # relative abundances-weighted mean deviation
  abdev<-abondrel*devdB
  
  # relative abundances-weighted mean of absolute deviations
  ababsdev<-abondrel*abs(devdB)
  
  # computation of FDiv
  FDiv[i]<-round( (sum(abdev)+meandB) / (sum(ababsdev)+meandB) ,6)
  
} 

# end of i


res<-list(nbsp=nbsp, FRic=FRic,FEve=FEve, FDiv=FDiv)
invisible(res)

res2 <- as.data.frame(res)
res2$site_quad <- rownames(res2)

res3 <- res2 %>% 
  separate(site_quad, c("site", "quad"), sep = "_") %>% 
  relocate(site, quad)

res4 <- abund %>% 
  separate(site_quad, c("site", "quad"), sep = "_") %>% 
  relocate(site, quad) %>% 
  select(site, quad, year, shrubcov, drygrass) %>% 
  left_join(res3) %>% 
  gather(key = "trait", value = "value", 7:9)

res4$year <- as.factor(res4$year)

#create plot of functional trait metrics by shrub cover in quadrat --------------
ft_2 <- ggplot(res4, aes(x = shrubcov, y = value)) +
  facet_wrap(vars(trait), scales = "free_y",
             strip.position = "left")+

  geom_jitter(aes(colour = year), size = 3, alpha = 0.5, width = 0.1)+
  geom_smooth(method = lm)+
  theme_cowplot()+
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y = element_text(size = 16))+
  xlab("Shrub cover (%)")+
  ylab(NULL)

ft_2

ggsave("ft2.png", plot = ft_2, width = 183, height = 80, units = "mm", bg = "white")


###############################
###############################
###############################
###############################

#
#
#
# Community Weighted Mean Traits -----------------------------
#
#
#

library(BAT)

##
### Load and clean data ------------------------------------
abundances = read.csv("data/snowp_floristics.csv",
                      header=T) #ordinal cover

#fix names
names(abundances)[1] <- "site"
names(abundances)[2] <- "date"
names(abundances)[3] <- "year"
names(abundances)[4] <- "quad"

#load additional data
shrubcover <- read.csv("data/snowp_additional.csv",
                       header=T)
#fix name
names(shrubcover)[1] <- "site"

#join
abund <- left_join(abundances, shrubcover) %>% 
  relocate(site, year, quad, shrubcov, drygrass)

#dataframe of just site [1] and each species [2:74]
abund2 <- abund[,-c(2:17)]

### Split data into three dataframes ------------------------------------
### ..... 1982, 2022 low shrub cover, 2022 high shrub cover....

ab_sh_max <- abund %>% 
  group_by(site) %>% 
  slice_max(order_by = SHRUB, with_ties = FALSE) %>% 
  filter(year == "2022")

ab_sh_min <- abund %>% 
  group_by(site) %>% 
  slice_min(order_by = SHRUB, with_ties = FALSE)

abund2max <- ab_sh_max[,-c(2:16)]
abund2old <- ab_sh_min %>% 
  filter(year != "2022")
abund2old <- abund2old[,-c(2:16)]
abund2min <- ab_sh_min %>% 
  filter(year == "2022")
abund2min <- abund2min[,-c(2:16)]
abund2max$type = "max"
abund2min$type = "min"
abund2old$type = "old"
abund2_all <- rbind(abund2max, abund2min, abund2old)
abund2_all <- abund2_all %>% unite("site", site, type)
head(abund2_all)

######


ab <- abund2_all[-1]
head(ab)

traits = read.csv("data/snowp_traits.csv", header=T, row.names=1)

traits <- scale(traits[1:4])
head(traits)


cwms <- as.data.frame(cwm(ab, traits))

cwms$site_type <- abund2_all$site


cwms2 <- cwms %>% 
  separate(site_type, c("site", "type"), sep = "_") %>% 
  relocate(site, type) %>% 
  gather(key = "trait", value = "value", 3:6)


cwms2$type <- factor(as.factor(cwms2$type), levels = c("old","min","max"),
                    labels = c("1982", "2022 (low)", "2022 (high)"))

cwms2$trait <- factor(as.factor(cwms2$trait), levels = c("HT","LDMC","SLA","SEED"),
                     labels = c("Height", "LDMC", "SLA","Seed mass"))


summary(aov(value ~ type, data = filter(cwms2, trait == "Height")))
aov(value ~ type, data = filter(cwms2, trait == "Height")) %>% TukeyHSD()
summary(aov(value ~ type, data = filter(cwms2, trait == "LDMC")))
summary(aov(value ~ type, data = filter(cwms2, trait == "Seed mass")))
summary(aov(value ~ type, data = filter(cwms2, trait == "SLA")))

ann_text3 <- data.frame(type = as.factor(c("2022 (high)", "2022 (low)", "1982",
                        "2022 (high)", "2022 (low)", "1982",
                        "2022 (high)", "2022 (low)", "1982",
                        "2022 (high)", "2022 (low)", "1982")),
                         value = c(1.1,1.1,1.1,
                                   1.1,1.1,1.1,
                                   1.1,1.1,1.1,
                                   1.1,1.1,1.1),
                        trait = c("Height", "Height", "Height",
                                  "LDMC","LDMC", "LDMC",
                                  "Seed mass","Seed mass","Seed mass",
                                  "SLA","SLA","SLA"),
                         lab = c("b","b","a",
                                 "a","a","a",
                                 "a","a","a",
                                 "a","a","a"))



ft_3 <- ggplot(cwms2, aes(x = type, y = value)) +
  facet_wrap(vars(trait))+
  geom_boxplot()+
  geom_jitter(aes(colour = type), size = 3, alpha = 0.5, width = 0.1)+
  theme_cowplot()+
  theme(legend.position="none",
        strip.placement = "inside",
        strip.text.x = element_text(size = 16))+
  xlab("Survey year")+
  ylab("Standardised trait value")+
  geom_text(data = ann_text3, aes(x = type, y = value, label = lab))

ft_3

ggsave("ft3.png", plot = ft_3, width = 183, height = 160, units = "mm", bg = "white")

##############################
##############################
##############################
##############################



