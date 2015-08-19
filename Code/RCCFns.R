################################################################
#Name: RCCFns.R
#Date: 3/23/15
#Author: Wesley Maddox
#Purpose: Rich-Club Coefficient calculation functions
################################################################
library(igraph)
searchInts <- function(g.edge,n1,n2){
  #returns 1 if there is an interaction, 0 else
  c1 <- which(g.edge[,1]==n1)
  c2 <- which(g.edge[,2]==n1)
  
  n1.edge <- cbind(c(g.edge[c1,1],g.edge[c2,2]),c(g.edge[c1,2],g.edge[c2,1]))
  m<-which(n1.edge==n2)
  if(length(m)>0){return(1)}
  else{return(0)}
}

edgeSwitching2<-function(g.edge,nr){#g is edgelist
  #output is new edgelist
  cols <- 1:2
  col <- sample(cols,1)
  ocol <- cols[! cols %in% col]
  
  nints <- dim(g.edge)[1]
  
  for (x in 1:(nr*nints)){
    int1 <- sample(1:nints,1)
    int2 <- sample(1:nints,1)
    
    if (length(unique(c(g.edge[int1,col],g.edge[int2,col],g.edge[int1,ocol],g.edge[int2,ocol])))==4){#ensures different nodes
      #include search for prior interactions
      #ex: A-B and C-D
      #print(x)
      s1 <- searchInts(g.edge,g.edge[int1,col],g.edge[int2,ocol]) #A and D
      s2 <- searchInts(g.edge,g.edge[int2,col],g.edge[int1,ocol]) #C and B
      
      if ((s1==0)&&(s2==0)){
        #print(x)
        edge1 <- c(g.edge[int1,col],g.edge[int2,ocol])
        edge2 <- c(g.edge[int2,col],g.edge[int1,ocol])
        g.edge <- g.edge[-c(int1,int2),]
        g.edge <- rbind(g.edge,rbind(edge1,edge2))      
      }
    }
  }
  return(g.edge)
}

RCC <- function(k, graph){
  #k is degree
  #graph is graph object
  graph.edge <- get.edgelist(graph)
  rc <- names(which(degree(graph)>k))
  
  nGTk <- length(rc) #nodes with degree greater than k
  
  edgesGTk <- graph.edge[intersect(which(graph.edge[,1]%in%rc),
                                   which(graph.edge[,2]%in%rc)),]
  eGTk <- dim(edgesGTk)[1] #edges among nodes with degree greater than k
  
  phi.k <- 2*eGTk/(nGTk*(nGTk-1))
  return(phi.k)
}


normalize.RCC <- function(graph,nr,k){
  graph.edge <- get.edgelist(graph)
  
  RCC2 <- c()
  for(q in 1:50){
    tmp.edge <- edgeSwitching2(graph.edge,nr)
    tmp.graph <- graph.edgelist(tmp.edge,directed=FALSE)
    RCC2 <- c(RCC2, RCC(k,tmp.graph))
  }
  
  rho.k <- RCC(k,graph)/mean(RCC2)
  return(rho.k)
}