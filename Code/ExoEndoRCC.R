#############################################################
#Author: Wesley Maddox
#perform merged exo-endo rich-club analysis
#both exocytosis and endocytosis are included in this dataset
#############################################################
source("C:/Users/Wesley/Documents/14-15/Spring15/EECS459/Project/Code/GeneratePPI.R")
source("C:/Users/Wesley/Documents/14-15/Spring15/EECS459/Project/Code/RCCFns.R")

#ppi is called exoEndo.graph
phi <- c()
for (x in 1:max(degree(exoEndo.graph))){
  phi <- c(phi,RCC(x,exoEndo.graph))
}

rho <- c()
for (x in 1:max(degree(exoEndo.graph))){
  rho <- c(rho, normalize.RCC(exoEndo.graph,10,x))
}

plot(phi)
plot(rho)