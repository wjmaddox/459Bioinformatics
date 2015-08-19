##########################################################################
#Name: GeneratePartialPPI.R
#Author: Wesley Maddox
#Purpose: Generate PPI list for exocytosis and endocytosis from BIOGRID dbf
#NOT a Script
##########################################################################

includeType <- function(row){
	if (row%in%both){
		return("BOTH")
	}
	else if(row%in%endo){
		return("Endocytosis")
	}
	else{return("Exocytosis")}
}
g1 <- grep("endocytosis",exoEndo$V6)
g2 <- grep("exocytosis",exoEndo$V6)
g3 <- intersect(g1,g2)

both <- exoEndo[g3,]
proteinsB <- unique(both$V2)
proteinsB <- matrix(unlist(strsplit(proteinsB,"_")),ncol=1,byrow=T)
proteinsB <- proteinsB[-which(proteinsB[,1]=="HUMAN"),]

Exo <- exoEndo2[-which(exoEndo2$Type=="Endocytosis"),]

proteinsX <- unique(Exo$V2)
proteinsX <- matrix(unlist(strsplit(proteinsX,"_")),ncol=1,byrow=TRUE)
proteinsX <- proteinsX[-which(proteinsX[,1]=="HUMAN"),]
both.whichOnly <- intersect(which(biogrid$OFFICIAL_SYMBOL_A%in%proteinsB),which(biogrid$OFFICIAL_SYMBOL_B%in%proteinsB))
both.edgelist <- as.matrix(biogrid[val.whichOnly,c(3,4)])
both.edgelist <- unique(both.edgelist[-which(both.edgelist[,1]==both.edgelist[,2]),])
both.graph <- graph.edgelist(both.edgelist,directed=F)
#exoEndo.which <- union(which(biogrid$OFFICIAL_SYMBOL_A%in%proteins),which(biogrid$OFFICIAL_SYMBOL_B%in%proteins))
#exoEndo.edgelist <- as.matrix(biogrid[exoEndo.which,c(3,4)])

#exoEndo.graph <- graph.edgelist(exoEndo.edgelist,directed=FALSE)

exo.whichOnly <- intersect(which(biogrid$OFFICIAL_SYMBOL_A%in%proteinsX),which(biogrid$OFFICIAL_SYMBOL_B%in%proteinsX))
exo.edgelist <- as.matrix(biogrid[exo.whichOnly,c(3,4)])

exo.edgelist <- unique(exo.edgelist[-which(exo.edgelist[,1]==exo.edgelist[,2]),])
exo.graph <- graph.edgelist(Endo.edgelist,directed=FALSE)

Endo <- exoEndo2[-which(exoEndo2$Type=="Exocytosis"),]

proteinsE <- unique(Endo$V2)
proteinsE <- matrix(unlist(strsplit(proteinsE,"_")),ncol=1,byrow=TRUE)
proteinsE <- proteinsE[-which(proteinsE[,1]=="HUMAN"),]

endo.whichOnly <- intersect(which(biogrid$OFFICIAL_SYMBOL_A%in%proteinsE),which(biogrid$OFFICIAL_SYMBOL_B%in%proteinsE))
endo.edgelist <- as.matrix(biogrid[endo.whichOnly,c(3,4)])

endo.edgelist <- unique(endo.edgelist[-which(endo.edgelist[,1]==endo.edgelist[,2]),])
endo.graph <- graph.edgelist(endo.edgelist,directed=FALSE)