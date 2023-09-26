###########################
#Upload data
############################

library(Rmisc)

#transect data
transect.ev.19 <- read.csv("Data/anolis_transect_data_ev_19.csv",header=TRUE)
transect.ev.21 <- read.csv("Data/anolis_transect_data_ev_21.csv",header=TRUE)

#site data
anolisSiteData_ev_21 <- read.csv("Data/anolissitedata_ev_21.csv",header = TRUE)
anolisSiteData_ev_19 <- read.csv("Data/anolissitedata_ev_19.csv",header = TRUE)

#########################################
#Source distance function
##########################################

source("distance_function.R")

##############################
#Bootstrap to estimate
#lambda at ev
##############################

ev_21 <- list()
ev_19 <- list()

#create data sets -1 row systematically 
for (i in 1:dim(transect.ev.19)[1]){

  ev_19[[i]] <- transect.ev.19[-i,]
  
}

for (i in 1:dim(transect.ev.21)[1]){
  
  ev_21[[i]] <- transect.ev.21[-i,]
  
}

#create data frame with all combinations of indices of 19 and 21
combinations <- expand.grid(1:dim(transect.ev.19)[1],1:dim(transect.ev.21)[1])

#initiate vectors to store results
lambdas_mid <- numeric()
lambdas_old <- numeric()
lambdas_young <- numeric()


for (j in 1:dim(combinations)[1]){
  
  temp19=N_est(ev_19[[combinations[j,]$Var1]],anolisSiteData_ev_19)
  temp21=N_est(ev_21[[combinations[j,]$Var2]],anolisSiteData_ev_21)
  
  lambdas_mid[j]=temp21$Estimate[1]/temp19$Estimate[1]
  lambdas_old[j]=temp21$Estimate[2]/temp19$Estimate[2]
  lambdas_young[j]=temp21$Estimate[3]/temp19$Estimate[3]

  print(paste0("Current iteration: ", j))
}

save.image(file="bootstrap.RData")

#Make dataframe for figure
lambdas <- (rbind(CI(lambdas_young),CI(lambdas_mid),CI(lambdas_old)))
lambdas$forest_age <- c("young","mid","old")  
lambdas$forest_age <- factor(lambdas$forest_age,levels=c("young","mid","old"))

#######
#Figure

density_plot <- ggplot(lambdas,aes(x=forest_age, y=mean, colour=forest_age,fill=forest_age))+
  geom_point(aes(color=forest_age,shape=forest_age,size=1))+
  geom_linerange(aes(ymin=lower,ymax=upper,color=forest_age))+
  scale_color_manual(values = c("gold2","darkolivegreen3","darkolivegreen4"))+
  theme_bw()+ylab(expression(hat(lambda)))+
  xlab("Forest Age")+theme(
    plot.tag.position=c(0.5,-0.1),plot.margin=margin(t = 10, r = 10, b = 40, l = 10))+
  theme(legend.position = "",axis.text=element_text(size=12),
        axis.title=element_text(size=12))

ggsave(density_plot, filename = "Fig_lambdas.pdf", device = cairo_pdf, 
       width = 5, height = 5, units = "in") 
