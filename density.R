##################################################
#Analysis of density and population dynamics    #
#Density a lambda                               #
#################################################

#######################
#Distance sampling    #
#######################

#transect data
transect.ev.19 <- read.csv("Data/anolis_transect_data_ev_19.csv",
                           header=TRUE)
transect.ev.21 <- read.csv("Data/anolis_transect_data_ev_21.csv",
                           header=TRUE)
transect.carite.21 <- read.csv("Data/anolis_transect_data_carite_21.csv",
                               header=TRUE)

#site data
anolisSiteData_carite_21 <- read.csv("Data/anolissitedata_carite_21.csv",
                                     header = TRUE)
anolisSiteData_ev_21 <- read.csv("Data/anolissitedata_ev_21.csv",
                                 header = TRUE)
anolisSiteData_ev_19 <- read.csv("Data/anolissitedata_ev_19.csv",
                                 header = TRUE)

#########################################
#Estimate N
##########################################

source("distance_function.R")

N_est_carite_21=N_est(transect.carite.21,anolisSiteData_carite_21)
N_est_ev_21=N_est(transect.ev.21,anolisSiteData_ev_21)
N_est_ev_19=N_est(transect.ev.19,anolisSiteData_ev_19)

N_ests <- rbind(N_est_ev_21,N_est_carite_21)
N_ests$sites <- c(rep("El Verde",3),rep("Carite2",3))
N_ests$sites=factor(N_ests$sites,levels=c("El Verde", "Carite2"))
N_ests$Forest_age=fct_recode(N_ests$Forest_age, Young="young",Mid="mid",Old="old")

###################################
#Density Figure                   #
###################################

density_plot <- ggplot(N_ests,aes(x=Forest_age, y=Estimate, colour=Forest_age,fill=Forest_age))+
  geom_errorbar(aes(ymin=lcl,ymax=ucl,color=Forest_age),linewidth=1)+
  geom_point(aes(color=Forest_age,shape=Forest_age,size=5))+
  facet_grid(~sites)+
  scale_color_manual(values = c("gold2","darkolivegreen3","darkolivegreen4"))+
  theme_bw()+ylab(expression(hat(N)))+
xlab("Forest Age")+theme(
    plot.tag.position=c(0.5,-0.1),plot.margin=margin(t = 10, r = 10, b = 40, l = 10))+
  theme(legend.position = "",axis.text=element_text(size=12),
        axis.title=element_text(size=12))
density_plot

ggsave(density_plot, filename = "Fig_abundance.pdf", device = cairo_pdf, 
       width = 8, height = 4, units = "in") 
