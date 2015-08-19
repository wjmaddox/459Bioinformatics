######################################
#FullVall.R
#Author: Wesley Maddox
#Date: 5/3/15
#Rich-Club Validation on Full Network
######################################

library(igraph)
library(foreach)
library(doParallel)
library(parallel)

print("loading Full Network")
load("~/EECS459/fullNetwork.RData")
print("making clusters")
cl <- makeCluster(16)
registerDoParallel(cl)
print("loading Rich-Club Functions")
source("~/EECS459/RichClubFns.R")

rcc.vals <- rcc.validate(exoEndo.graph)
rcc.val <- as.numeric(colSums(rcc.vals)/length(V(exoEndo.graph)))
save(rcc.val,file="~/EECS459/FullValRichClub.RData")
