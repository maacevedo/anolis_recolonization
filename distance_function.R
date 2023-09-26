library(REdaS)
library(Distance)
library(ggplot2)
library(tidyverse)


N_est=function(transect,anolisSiteData){

relative_bearing <- numeric()
for (i in 1:dim(transect)[1]){
  if(transect$angle[i] < transect$main_heading[i]){ 
    relative_bearing[i] <- transect$angle[i]+360-transect$main_heading[i]}else
    { 
      relative_bearing[i] <- transect$angle[i]-transect$main_heading[i]}
}  

########################################
#Calculate relative bearing adjusted
#i.e. calculate angle to horizontal line
#########################################

relative_bearing_adjusted <- numeric()  
for (i in 1:dim(transect)[1]){  
  if(relative_bearing[i] >= 0 & relative_bearing[i] < 91){
    relative_bearing_adjusted[i] <- relative_bearing[i]
  }else if(relative_bearing[i] >= 91 & relative_bearing[i] < 181){
    relative_bearing_adjusted[i] <- 180 - relative_bearing[i]
  }else if(relative_bearing[i] >= 181 & relative_bearing[i] < 270){
    relative_bearing_adjusted[i] <- relative_bearing[i]-180
  }else
    relative_bearing_adjusted[i] <- 360 - relative_bearing[i]
}

#convert to radians to use sin function
relative_bearing_adjusted_rad <- deg2rad(relative_bearing_adjusted)

#Calculate the oposite length (perpendicular distance to the transect)
dist <- sin(relative_bearing_adjusted_rad)*transect$distance_observer

dist.df <- data.frame(Region.Label=transect$forest_age,Area=50,Effort=1,distance=dist,Sample.Label=transect$site)
estimate <- ds(dist.df,truncation=6)
gof.e <- gof_ds(estimate)

#save N
N<- estimate$dht$individuals$N[1:3,]
N$Forest_age <- c("old","mid","young")
N$Forest_age <- factor(N$Forest_age, levels = c("young","mid","old"))

return(N)
}