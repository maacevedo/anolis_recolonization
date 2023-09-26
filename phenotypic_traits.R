##################################################
#Analysis of Phenotypic traits in chrono        #
#SVL, Weight, Body Condition                    #
#Limbs (f1,f2,b1, and b2)                       #
#################################################

#############
#Libraries  #
#############

library(ggplot2)
library(sjPlot)
library(sjmisc)
library(viridis)
library(cowplot)
library(extrafont)
library(lme4)
library(RLRsim)
library(pbkrtest)
library(effects)
library(tidyr)
library(tidyverse)
library(plyr)
library(dplyr)
library(knitr)
library(MuMIn)
library(gtExtras)
library(broom)
library(dplyr)
library(ggeffects)
library(gt)
library(htmlTable)
library(ggpubr)
library(rstatix)
library(ggprism)
library(magick)
library(grid)

#######################################################################
#Upload data                                                          #
#######################################################################

data.g <- read.csv("Data/data_phenotypic_traits.csv",
                 header=TRUE)
dim(data.g)
head(data.g)


#set levels of forest_age, making old forest the baseline
data.g$forest_age <- factor(data.g$forest_age, levels=c("Old","Mid","Young"))

#set levels of site making El Verde (place with the most samples) the baseline
data.g$site <- factor(data.g$site, levels=c("El Verde","Carite1","Carite2"))

var_select=c("site")
count_freq=count(data.g,var_select)

var_select2=c("site","sex")
count_freq2=count(data.g,var_select2)
count_freq2 %>% htmlTable
N=sum(count_freq[,2])
N

#######################################################################
#SVL                                                                  #
#######################################################################

#Fixed effects model
mod.svl.fm <- lm(log(svl)~forest_age*site+sex,data=data.g)
summary(mod.svl.fm)

#random effects model
mod.svl.re <- lmer(log(svl)~forest_age*site+sex+(1|season),data=data.g,
                   REML = FALSE) #singular fit

#For supporting information
tab_model(mod.svl.fm,show.se=TRUE,show.stat=TRUE)

#predictions for plot
ee_svl <- data.frame(Effect(c("forest_age","site","sex"),mod.svl.fm))


#######################################################################
#SVL Fig                                                              #
#######################################################################

ee_svl$forest_age=factor(ee_svl$forest_age,levels=c("Young","Mid","Old"))
ee_svl$sex=factor(ee_svl$sex,levels=c("Males","Females","Juveniles"))

psvl <- ggplot(ee_svl,
               aes(forest_age,exp(fit)))+
  geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper), color = forest_age),
                position = position_dodge(0.0), width = 0.2)+
  facet_grid(sex~site,scales="free")+
  geom_point(aes(color = forest_age,shape=forest_age), position = position_dodge(0.0),size=2)+
  scale_color_manual(values = c("gold2","darkolivegreen3","darkolivegreen4"))+
  ylab("Snout-to-vent length (mm)")+xlab("Forest Age")+theme_bw()+theme(legend.position = "",axis.text=element_text(size=12),
                                                                        axis.title=element_text(size=12))

sign.differences <- tibble::tribble(
  ~group1, ~group2, ~p.signif, ~y.position, ~site, ~sex,
  "Young", "Old",   "*",       58,          "El Verde", "Males",
  "Young", "Old",   "*",       46,          "El Verde", "Females",
  "Young", "Old",   "*",       38,          "El Verde", "Juveniles",
  "Young", "Old",   "*",       68,          "Carite1",  "Males",
  "Young", "Old",   "*",       54,          "Carite1",  "Females",
  "Young", "Old",   "*",       44,          "Carite1",  "Juveniles"
)

sign.differences$sex=factor(sign.differences$sex,levels=c("Males","Females","Juveniles"))
sign.differences$site=factor(sign.differences$site,levels=c("El Verde","Carite1","Carite2"))

psvl <- psvl+add_pvalue(sign.differences)

svl <- image_read("svl.png")
svl <- image_trim(svl)
psvl
grid.raster(svl,x=0.23,y=0.87, width=0.2)


#######################################################################
#Weight                                                               #
#######################################################################

#Fixed effects model
mod.w.fm <- lm(log(weight)~forest_age*site+sex,data=data.g)
summary(mod.w.fm)

mod.w.re <- lmer(log(weight)~forest_age+site+sex+(1|season),data=data.g) 
#singular fit

#for supporting information
tab_model(mod.w.fm,show.se=TRUE,show.stat=TRUE)

#######################################################################
#Body Condition                                                       #
#######################################################################

#filter males
data.gm <- data.g %>% select(uid,site,forest_age,svl,weight,sex,season) %>% 
  filter(sex == "Males")

data.gm=data.gm[complete.cases(data.gm[,4:5]),] 

data.gm$residuals <- residuals(lm(log10(weight)~log10(svl),data=data.gm))

#Fixed effects model
mod.bc.m <- lm(residuals~forest_age*site,data=data.gm)
summary(mod.bc.m)


mod.bc.m.re <- lmer(residuals~forest_age*site+(1|season),data=data.gm, REML=FALSE)
tab_model(mod.bc.m.re,show.se=TRUE,show.stat=TRUE)

#LRT to assess if RE are needed
bc_eLRT <- exactLRT(mod.bc.m.re,mod.bc.m)

mod.bc.m.re <- lmer(residuals~forest_age*site+(1|season),data=data.gm) 
mod.bc.m.re2 <- lmer(residuals~site+(1|season),data=data.gm) 

#LRT testing for differences by forest age
kr.bc <-KRmodcomp(mod.bc.m.re,mod.bc.m.re2) 

#######################################################################
#Limbs Comparisons                                                    #
#######################################################################

##################
#f1 (radius/ulna)# 
##################

mod.f1 <- lm(log(f1)~log(svl)+forest_age*site+sex,data=data.g) #fixed effects model

#fixed effect of sample with random effect of within sample sampling date
mod.f1.re <- lmer(log(f1)~log(svl)+forest_age*site+sex+(1|season),data=data.g, REML=FALSE) #fixed effects model

#LRT to assess if RE are needed
f1.eLRT <- exactLRT(mod.f1.re,mod.f1)
f1.eLRT

mod.f1.re <- lmer(log(f1)~log(svl)+forest_age*site+sex+(1|season),data=data.g) 
mod.f1.re.2 <- lmer(log(f1)~log(svl)+site+sex+(1|season),data=data.g) 

#Kenward-Roger method
kr.f1 <-KRmodcomp(mod.f1.re,mod.f1.re.2) 
summary(kr.f1)

#Make inference
summary(mod.f1.re)
tab_model(mod.f1.re,show.se=TRUE,show.stat=TRUE)
mod.f1.re_r2 <- r.squaredGLMM(mod.f1.re) #r2 r2c is the conditional (fixed and RE)

summary(mod.f1.re)

ee_f1 <- data.frame(Effect(c("forest_age","site","sex"),mod.f1.re))

#########################
#Figure f1 (radius/ulna)#
#########################

ee_f1$forest_age=factor(ee_f1$forest_age,levels=c("Young","Mid","Old"))
ee_f1$sex=factor(ee_f1$sex,levels=c("Males","Females","Juveniles"))

pf1 <- ggplot(ee_f1,
              aes(forest_age,exp(fit)))+
  geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper), color = forest_age),
                position = position_dodge(0.0), width = 0.2)+
  facet_grid(sex~site,scale="free")+
  geom_point(aes(color = forest_age,shape=forest_age), position = position_dodge(0.0),size=5)+
  scale_color_manual(values = c("gold2","darkolivegreen3","darkolivegreen4"))+
  ylab("Radius/ulna (mm)")+xlab("Forest Age")+theme_bw()+theme(legend.position = "",axis.text=element_text(size=12),
                                                               axis.title=element_text(size=12))
pf1.sign.differences <- tibble::tribble(
  ~group1, ~group2, ~p.signif, ~y.position, ~site, ~sex,
  "Young", "Old",   "*",       9,          "El Verde", "Males",
  "Young", "Old",   "*",       8.7,          "El Verde", "Females",
  "Young", "Old",   "*",       9,          "El Verde", "Juveniles",
  "Young", "Old",   "*",       10.5,          "Carite2",  "Males",
  "Young", "Old",   "*",       10,          "Carite2",  "Females",
  "Young", "Old",   "*",       10.5,          "Carite2",  "Juveniles",
  "Old", "Mid",   "**",       9,          "Carite2",  "Males",
  "Old", "Mid",   "**",       8.5,          "Carite2",  "Females",
  "Old", "Mid",   "**",       9,          "Carite2",  "Juveniles"
)

pf1.sign.differences$sex=factor(pf1.sign.differences$sex,levels=c("Males","Females","Juveniles"))
pf1.sign.differences$site=factor(pf1.sign.differences$site,levels=c("El Verde","Carite1","Carite2"))

pf1 <- pf1+add_pvalue(pf1.sign.differences)

f1 <- image_read("f1.png")
f1 <- image_trim(f1)
pf1
grid.raster(f1,x=0.13,y=0.91, width=0.1)


###############
#f2           #
###############

mod.f2 <- lm(log(f2)~log(svl)+forest_age*site+sex,data=data.g) #fixed effects model

#fixed effect of sample with random effect of within sample sampling date
mod.f2.re <- lmer(log(f2)~log(svl)+forest_age*site+sex+(1|season),data=data.g, REML=FALSE) #fixed effects model

#LRT to assess if RE are needed
f2.eLRT=exactLRT(mod.f2.re,mod.f2)

mod.f2.re <- lmer(log(f2)~log(svl)+forest_age*site+sex+(1|season),data=data.g) 
mod.f2.re.2 <- lmer(log(f2)~log(svl)+site+sex+(1|season),data=data.g) 

#Kenward-Roger method
kr.f2 <-KRmodcomp(mod.f2.re,mod.f2.re.2) 
summary(kr.f2)

#Make inference
summary(mod.f2.re)
tab_model(mod.f2.re,show.se=TRUE,show.stat=TRUE)
r2.f2=round(r.squaredGLMM(mod.f2.re),2)

ee_f2 <- data.frame(Effect(c("forest_age","site","sex"),mod.f2.re))

#########################
#Figure f2              #
#########################

ee_f2$forest_age=factor(ee_f2$forest_age,levels=c("Young","Mid","Old"))
ee_f2$sex=factor(ee_f2$sex,levels=c("Males","Females","Juveniles"))

pf2 <- ggplot(ee_f2,
              aes(forest_age,exp(fit)))+
  geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper), color = forest_age),
                position = position_dodge(0.0), width = 0.2)+
  facet_grid(sex~site,scale="free")+
  coord_cartesian(ylim = c(7, 13))+
  geom_point(aes(color = forest_age,shape=forest_age), position = position_dodge(0.0),size=5)+
  scale_color_manual(values = c("gold2","darkolivegreen3","darkolivegreen4"))+
  ylab("Humerus (mm)")+xlab("Forest Age")+theme_bw()+
  theme(legend.position = "",axis.text=element_text(size=12),
        axis.title=element_text(size=12))

pf2.sign.differences <- tibble::tribble(
  ~group1, ~group2, ~p.signif, ~y.position, ~site, ~sex,
  "Young", "Old",   "*",       11,          "El Verde", "Males",
  "Young", "Old",   "*",       10.5,          "El Verde", "Females",
  "Young", "Old",   "*",       11,          "El Verde", "Juveniles",
  "Mid", "Old",   "**",       10.5,          "El Verde", "Males",
  "Mid", "Old",   "**",       10,          "El Verde", "Females",
  "Mid", "Old",   "**",       10.5,          "El Verde", "Juveniles",
  "Young", "Old",   "**",       11.3,          "Carite2",  "Males",
  "Young", "Old",   "**",       10.7,          "Carite2",  "Females",
  "Young", "Old",   "**",       11,          "Carite2",  "Juveniles",
  "Old", "Mid",   "**",       10.5,          "Carite2",  "Males",
  "Old", "Mid",   "**",       10,          "Carite2",  "Females",
  "Old", "Mid",   "**",       10.3,          "Carite2",  "Juveniles"
)

pf2.sign.differences$sex=factor(pf2.sign.differences$sex,levels=c("Males","Females","Juveniles"))
pf2.sign.differences$site=factor(pf2.sign.differences$site,levels=c("El Verde","Carite1","Carite2"))

pf2 <- pf2+add_pvalue(pf2.sign.differences)

f2 <- image_read("f2.png")
f2 <- image_trim(f2)
pf2
grid.raster(f2,x=0.13,y=0.91, width=0.1)

###############
#b1           #
###############

mod.b1 <- lm(log(b1)~log(svl)+forest_age*site+sex,data=data.g) #fixed effects model
summary(mod.b1)

#fixed effect of sample with random effect of within sample sampling date
mod.b1.re <- lmer(log(b1)~log(svl)+forest_age*site+sex+(1|season),data=data.g, REML=FALSE) #fixed effects model

#LRT to assess if RE are needed
b1.eLRT <- exactLRT(mod.b1.re,mod.b1)

mod.b1.re <- lmer(log(b1)~log(svl)+forest_age*site+sex+(1|season),data=data.g) 
mod.b1.re.2 <- lmer(log(b1)~log(svl)+site+sex+(1|season),data=data.g) 

#Kenward-Roger method
kr.b1 <-KRmodcomp(mod.b1.re,mod.b1.re.2) 
summary(kr.b1)

#Make inference
summary(mod.b1.re)
tab_model(mod.b1.re,show.se=TRUE,show.stat=TRUE)

r2.b1=round(r.squaredGLMM(mod.b1.re),2)

###############
#b2           #
###############

mod.b2 <- lm(log(b2)~log(svl)+forest_age*site+sex,data=data.g) #fixed effects model

#fixed effect of sample with random effect of within sample sampling date
mod.b2.re <- lmer(log(b2)~log(svl)+forest_age*site+sex+(1|season),data=data.g, REML=FALSE) #fixed effects model

#LRT to assess if RE are needed
b2.eLRT <- exactLRT(mod.b2.re,mod.b2)

mod.b2.re <- lmer(log(b2)~log(svl)+forest_age*site+sex+(1|season),data=data.g) 
mod.b2.re.2 <- lmer(log(b2)~log(svl)+site+sex+(1|season),data=data.g) 

#Kenward-Roger method
kr.b2 <-KRmodcomp(mod.b2.re,mod.b2.re.2) 
summary(kr.b2)

#Make inference
summary(mod.b2.re)
tab_model(mod.b2.re,show.se=TRUE,show.stat=TRUE)

r2.b2=round(r.squaredGLMM(mod.b2.re),2)


##########################################
#Ad-hoc analysis on broken tails
#as proxy for predation
##########################################

data.g <- data.g %>% mutate(tail_broken = recode(tail,
                                                   'tc'=0,
                                                   'tb'=1,
                                                   'tr'=1))

mod_tail <- glm(tail_broken~site*forest_age+sex,family="binomial",data=data.g)
summary(mod_tail)
tab_model(mod_tail,show.se=TRUE,show.stat=TRUE)
