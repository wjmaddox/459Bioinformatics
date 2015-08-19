############################################################################
#GeneratePPI.R
#Author: Wesley Maddox
#Date: 5/3/15
#Purpose: read in both Biogrid and GeneOntology datasets and convert them to
#a merged PPI-network with both exocytosis and endocytosis genes
############################################################################
library(igraph)

#BIOGRID-fixed removes the first 35 lines from the entire biogrid file
biogrid <- read.delim("BIOGRID-fixed.txt",header=TRUE,fill=TRUE,stringsAsFactors=FALSE)

exoEndo <- read.delim("GeneOntologyExoEndo.txt",header=FALSE,fill=TRUE,stringsAsFactors=FALSE)
exoEndo.tmp <- matrix(unlist(strsplit(exoEndo$V1,":")),byrow=TRUE,ncol=2)

exoEndo$TypeOfInt <- exoEndo.tmp[,1]
exoEndo$Name <- exoEndo.tmp[,2]

#only want PPIs, so remove reactome 
exoEndo2 <- exoEndo[-which(exoEndo$TypeOfInt=="Reactome"),]

proteins <- unique(exoEndo2$V2)
proteins <- matrix(unlist(strsplit(proteins,"_")),ncol=1,byrow=TRUE)
proteins <- proteins[-which(proteins[,1]=="HUMAN"),]


exoEndo.whichOnly <- intersect(which(biogrid$OFFICIAL_SYMBOL_A%in%proteins),which(biogrid$OFFICIAL_SYMBOL_B%in%proteins))
exoEndo.edgelist <- as.matrix(biogrid[exoEndo.whichOnly,c(3,4)])
#want to remove self-interactions & double interactions
exoEndo.edgelist <- unique(exoEndo.edgelist[-which(exoEndo.edgelist[,1]==exoEndo.edgelist[,2]),])
exoEndo.graph <- graph.edgelist(exoEndo.edgelist,directed=FALSE)

save(exoEndo.graph,file="fullNetwork.RData")
