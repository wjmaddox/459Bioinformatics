library(igraph)

semilc <- function(assg){
  qu <- rep(0,vcount(assg))
  for (n in 1:vcount(assg)){
    nn <- neighbors(assg,n) #set of nearest neighbors of n
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