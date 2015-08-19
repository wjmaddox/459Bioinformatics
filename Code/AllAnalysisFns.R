###################################################3
#AllAnalysisFns.R
#Author: Wesley Maddox
#5/3/15
#contains all functions to be sourced
#####################################################

#####################################################
#fpc scoring
#####################################################
#semi-local centrality measure
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

#nearest neighbors centrality measure
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

#fpc score
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
  
  c <- cbind(c1,c2,c3,c4,c5,c6,c7)
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

##############################
#rich-club analysis
##############################
#deltaS function
deltaS <- function(a,b,m,n){
  return(1/(1+exp(-(m-n*b)*tan(pi*a/2))))
}

#core-Score function
runCS <- function(a,b,assg){
  require(igraph)
  print(c(a,b))
  assg.mat <- get.adjacency(assg)
  c.StarU <- deltaS(a,b,1:vcount(assg),vcount(assg))
  c.Star <- c.StarU/sum(c.StarU)
  c.i <- sample(c.Star,vcount(assg))
  #nest function in order to only have to pass c.i
  sa.optim <- function(c.i){
    c <- c.i%*%t(c.i)
    r <- sum(as.matrix(assg.mat%*%c))
    return(-r)
  }
  out <- optim(fn=sa.optim,par=c.i,method="SANN")
  cs <- out$par*(-out$value)
  return(cs)
}

#richClub scoring wrapper
runRichClub <- function(graph2){
  graph2.mat <- get.adjacency(graph2)
  
  c.StarU <- deltaS(aT,bT,1:vcount(graph2),vcount(graph2))
  c.Star <- c.StarU/sum(c.StarU)
  c.i <- sample(c.Star,vcount(graph2))
  c.j <- sample(c.Star,vcount(graph2))
  
  #cs.ab <- matrix(nrow=10,ncol=10)
  cs.ab <- data.frame(matrix(nrow=100,ncol=vcount(graph2)))
  #cs.ab <- c()
  z<-1
  l <- matrix(nrow=100,ncol=2)
  s <- seq(.1,1,.1)
  for (a in 1:10){
    for (b in 1:10){
      #cs.ab[z,] <- runCS(s[a],s[b],graph2)
      #print(paste(s[a],s[b]))
      l[z,] <- c(s[a],s[b])
      z<-z+1
    }
  }
  #rewrite of this command
  cs.ab <- foreach(z=1:nrow(l),.combine=rbind,.packages()) %dopar% (runCS(l[z,1],l[z,2],graph2))
  
  cs.score <- colSums(cs.ab)
  cs.norm <- cs.score/max(cs.score)
  names(cs.norm) <- names(degree(graph2))
  #turn into function
  #also match the scores into a contour plot for alpha and beta
  
  maxScore <- c()
  for (x in 1:100){
    maxScore <- c(maxScore,which(cs.ab[x,]==max(cs.ab[x,])))
  }
  
  maxScore.mat <- matrix(maxScore,ncol=10,byrow=TRUE)
  rownames(maxScore.mat) <- s
  colnames(maxScore.mat) <- s
  return(list(maxScore.mat,cs.norm))
}