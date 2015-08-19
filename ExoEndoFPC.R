#source("C:/Users/Wesley/Documents/14-15/Spring15/EECS459/Project/Code/GeneratePPI.R")
#ppi is called ExoEndo.graph
library(igraph)
load("fullNetwork.RData")
semilc <- function(assg){
  qu <- rep(0,vcount(assg))
  for (n in 1:vcount(assg)){
    nn <- neighbors(assg,n) #set of nearest neighbors of n
    w2 <- 0
    w3 <- 0
    for (w in nn){
      w2 <- neighbors(assg,n)  #nearest neighbors of w
      for (w.n in w2){
        w3 <- neighbors(assg,w.n) #next nearest neighbors of w2
      }
    }
    qu[n] <- sum(length(w2),length(w3))
  }
  clv <- rep(0,vcount(assg))
  for (n in 1:vcount(assg)){
    clv[n] <- sum(qu[neighbors(assg,n)])
  }
  names(clv) <- names(degree(assg))
  return(clv)
}

nmc <- function(assg){
  nmc.mat <- matrix(nrow=vcount(assg),ncol=16) #16 is all motifs
  for (n in 1:vcount(assg)){
    n.neighbors <- neighbors(assg,n)
    #n.Nneighbors <- c()
    #for (b in n.neighbors){ }
    assg.sub <- subgraph(assg,c(n,n.neighbors))
    motifs.2 <- ecount(assg.sub)
    motifs.3 <- graph.motifs(assg.sub,size=3)
    motifs.4 <- graph.motifs(assg.sub,size=4)
    nmc.mat[n,] <- c(motifs.2,motifs.3,motifs.4)
  }
  nmc.mat[is.na(nmc.mat)]<-0
  #unweighted matrix for motifs
  s <- cov(nmc.mat)
  
  s.eigen <- eigen(s)
  l <- which(s.eigen$values==max(s.eigen$values))
  a <- s.eigen$vectors[l,]
  
  i.score <- c()
  for(n in 1:nrow(nmc.mat)){
    i.score[n] <- sum(a*nmc.mat[n,])
  }
  names(i.score) <- names(degree(assg))
  return(i.score)
}


fpc <- function(assg){
  c1 <- degree(assg)
  c2 <- betweenness(assg)
  c3 <- closeness(assg)
  c4 <- transitivity(assg,type="localundirected")
  c5 <- graph.coreness(assg)
  c6 <- evcent(assg)$vector
  #c7 <- page.rank(assg)$vector #using pageRank instead of semi-local centrality
  c7 <- semilc(assg)
  c8 <- nmc(assg)
  
  c <- cbind(c1,c2,c3,c4,c5,c6,c7,c8)
  c[is.na(c)]<-0
  cv <- cov(c)
  w <- eigen(cv)$values
  
  fpc <- c()
  for (x in 1:dim(c)[1]){
    fpc <- c(fpc,sum(w*c[x,]))
  }
  names(fpc) <- names(degree(assg))
  return(fpc)
  #names(which(fpc==max(fpc)))
}

fpc.score <- fpc(exoEndo.graph)
save(fpc.score,file="FullFPC.RData")
#plot(fpc.score)
print("Saved Score to FullFPC")
#t <- tail(names(fpc.score[order(fpc.score)]),n=25)
#assg.edge <- get.edgelist(assg)

#t.g <- graph.edgelist(exoEndo.edgelist[intersect(which(exoEndo.edgelist[,1]%in%t),which(exoEndo.edgelist[,2]%in%t)),],
                      directed=F)
#tkplot(t.g)


fpc.validate <- function(assg){
  n <- V(assg)$name
  f.rem <- function(assg, i){
    return(delete.vertices(assg,n[i]))
  }
  #dopar here
  fpc.valsT <- foreach(i=1:length(n),.combine=rbind,.export=c('fpc','nmc','semilc'),.packages='igraph') %do% (fpc(f.rem(assg,i)))
  fpc.vals <- c(0,fpc.valsT[1,])
  for(j in 2:nrow(fpc.valsT)){
    fpc.tmp <- as.numeric(fpc.valsT[j,])
    fpc.vals <- rbind(fpc.vals,c(fpc.tmp[1:(j-1)],0,fpc.tmp[j:length(fpc.tmp)]))
  }
  return(fpc.vals)
}

fpc.val <- as.numeric(colSums(fpc.vals)/22)
save(fpc.val,file="FullValFPC.RData")
print("Saved Validated Score to FullValFPC")