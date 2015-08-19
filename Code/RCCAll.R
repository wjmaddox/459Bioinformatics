#################################################
#requires RCCFns.R & GeneratePPI.R
#Name: RunRCCAll.R
#################################################
#also foreach & doparallel
library(foreach)

cl <- makeCluster(3,type="SOCK")
registerDoSNOW(cl)

n <- V(exoEndo.graph)$name
foreach(i=1:length(n),.combine=c,.export=c('searchInts','edgeSwitching2','RCC','normalize.RCC'),.packages='igraph') %dopar% (normalize.RCC(exoEndo.graph,10,i))
