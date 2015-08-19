#########################################################
#RichClubFns.R
#########################################################
deltaS <- function(a,b,m,n){
  return(1/(1+exp(-(m-n*b)*tan(pi*a/2))))
}

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

runRichClub <- function(graph2){
  graph2.mat <- get.adjacency(graph2)
  
  #c.StarU <- deltaS(aT,bT,1:vcount(graph2),vcount(graph2))
  #c.Star <- c.StarU/sum(c.StarU)
  #c.i <- sample(c.Star,vcount(graph2))
  #c.j <- sample(c.Star,vcount(graph2))
  
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
  cs.ab <- foreach(z=1:nrow(l),.combine=rbind,.packages='igraph',.export=c('deltaS','runCS','runRichClub')) %dopar% (runCS(l[z,1],l[z,2],graph2))
  
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

rcc.validate <- function(assg){
  n <- V(assg)$name
  f.rem <- function(assg, i){
    return(delete.vertices(assg,n[i]))
  }
  #dopar here
  fpc.valsT <- foreach(i=1:length(10),.combine=rbind,.export=c('deltaS','runCS','runRichClub'),.packages='igraph') %do% (runRichClub(f.rem(assg,i))[[2]])
  fpc.vals <- c(0,fpc.valsT[1,])
  for(j in 2:nrow(fpc.valsT)){
    fpc.tmp <- as.numeric(fpc.valsT[j,])
    fpc.vals <- rbind(fpc.vals,c(fpc.tmp[1:(j-1)],0,fpc.tmp[j:length(fpc.tmp)]))
  }
  return(fpc.vals)
}
