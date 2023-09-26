##################################################
#Analysis of malaria parasitism in chrono       #
#################################################

library(ggplot2)
library(sjPlot)
library(tidyverse)
library(ggeffects)
library(broom)
library(ggprism)
library(RLRsim)
library(pbkrtest)
library(lme4)
library(pscl)

###################
#Upload data
###################
infection <- read.csv("Data/data_parasitism.csv",header=TRUE)
head(infection)

###################
#Format data
###################

infection$forest_age <- factor(infection$forest_age,levels=c("Old","Mid","Young"))
infection$site <- factor(infection$site,levels=c("El Verde","Carite1","Carite2"))

####################
#Model
####################

mod.inf <- glm(infection~forest_age+site+sex,family="binomial", data=infection)
summary(mod.inf)
tab_model(mod.inf)

pR2(mod.inf)

#singular fit
mod.inf.re <- glmer(infection~forest_age+site+sex+(1|season),data=infection,family=binomial)

#Calculate odds ratios with their SE
model.df <- tidy(mod.inf)  # Convert model to dataframe for easy manipulation
model.df <- model.df %>% 
  mutate(or = exp(estimate),  # Odds ratio/gradient
         var.diag = diag(vcov(mod.inf)),  # Variance of each coefficient
         or.se = sqrt(or^2 * var.diag)) 

#Prepare data for plotting
eff_inf <- ggpredict(mod.inf,terms=c("forest_age","site"))
eff_inf$x <- factor(eff_inf$x,levels=c("Young","Mid","Old"))

sign.differences <- tibble::tribble(
  ~group1, ~group2, ~p.signif, ~y.position, ~group,
  "Young", "Old",   "*",       0.75,           "El Verde",
  "Young", "Old",   "*",       0.85,           "Carite1",
  "Young", "Old",   "*",       0.85,           "Carite2"
  )

sign.differences$group <- factor(sign.differences$group,levels=c("El Verde","Carite1","Carite2"))

pinf <- ggplot(eff_inf,
               aes(x=x,y=predicted))+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, color = x),
                position = position_dodge(0.0), width = 0.2)+
  facet_grid(~group)+
  geom_point(aes(color = x,shape=x), position = position_dodge(0.0),size=5)+
  scale_color_manual(values = c("gold2","darkolivegreen3","darkolivegreen4"))+
  ylab("Probability of infection")+xlab("Forest Age")+theme_bw()+
  theme(legend.position = "",axis.text=element_text(size=12),axis.title=element_text(size=12))

pinf <- pinf+add_pvalue(sign.differences)

ggsave(pinf, filename = "Fig_inf_pred.pdf", device = cairo_pdf, 
       width = 6, height = 3, units = "in") 
