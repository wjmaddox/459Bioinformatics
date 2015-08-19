#source("C:/Users/Wesley/Documents/14-15/Spring15/EECS459/Project/Code/GeneratePPI.R")
library(igraph,foreach)
library(doParallel)
cl <- makeCluster(3,type="SOCK")
registerDoSNOW(cl)
#need to save these functions in separate files and then source them
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

#finalList <- runRichClub(exoEndo.graph)
#save(finalList,file="~/EECS459/FullRichClub.RData")

rccW <- function(assg){
  listT <- runRichClub(assg)
  return(as.numeric(listT[[2]]))
}
rcc.validate <- function(assg){
  n <- V(assg)$name
  f.rem <- function(assg, i){
    return(delete.vertices(assg,n[i]))
  }
  #dopar here
  fpc.valsT <- foreach(i=1:length(n),.export=c('deltaS','runCS','runRichClub','rccW'),.packages=c('igraph','doParallel')) %dopar% (rccW(f.rem(assg,i)))
  #fpc.vals <- c(0,fpc.valsT[1,])
  
  return(fpc.valsT)
}

rcc.vals <- c(0,rcc.valT2[1,])
for(j in 2:(nrow(rcc.valT2)-1)){
  rcc.tmp <- as.numeric(rcc.valT2[j,])
  #print(length(rcc.tmp))
  rv <- c(rcc.tmp[1:(j-1)],0,rcc.tmp[j:length(rcc.tmp)])
  print(length(rv))
  rcc.vals <- rbind(rcc.vals, rv)
}
rcc.vals <- rbind(rcc.vals,c(as.numeric(rcc.valT2[j,]),0))

rcc.val <- as.numeric(colSums(rcc.vals)/length(V(exoEndo.graph)))
save(rcc.val,file="~/EECS459/FullValRichClub.RData")
#ms.mat <- finalList[[1]]
#cs.norm <- finalList[[2]]
#library(corrplot)
#col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F","cyan", "#007FFF", "blue", "#00007F"))
#corrplot(maxScore.mat,is.corr=FALSE,method="color",col = col4(34),cl.length=18,title="Max for Different Combinations of Alpha&Beta", mar=c(1,1,1,1))

#make these ggplot
#measures <- data.frame(t(rbind(degree(exoEndo.graph),cs,fpc)))

#g1 <- ggplot(measures,aes(x=))
#corrplot(maxScore.mat,is.corr=FALSE,method="color",col = col4(10))
#add in title and change legend
#finalList <- list(cs.norm=cs.norm,ms=maxScore.mat)
  #return(finalList)
################################################
#full network validation functions
################################################
rcc.validate1 <- function(assg){
  n <- V(assg)$name
  f.rem <- function(assg, i){
    return(delete.vertices(assg,n[i]))
  }
  #dopar here
  fpc.valsT <- foreach(i=1:56,.export=c('deltaS','runCS','runRichClub','rccW'),.packages=c('igraph','doParallel')) %dopar% (rccW(f.rem(assg,i)))
  #fpc.vals <- c(0,fpc.valsT[1,])
  
  return(fpc.valsT)
}
rcc.val1 <- rcc.validate1(exoEndo.graph)
save(rcc.val1,file="EECS459/FullValRC1.RData")

rcc.validate2 <- function(assg){
  n <- V(assg)$name
  f.rem <- function(assg, i){
    return(delete.vertices(assg,n[i]))
  }
  #dopar here
  fpc.valsT <- foreach(i=57:112,.export=c('deltaS','runCS','runRichClub','rccW'),.packages=c('igraph','doParallel')) %dopar% (rccW(f.rem(assg,i)))
  #fpc.vals <- c(0,fpc.valsT[1,])
  
  return(fpc.valsT)
}
rcc.val2 <- rcc.validate2(exoEndo.graph)
save(rcc.val2, file="EECS459/FullValRC2.RData")

rcc.validate3 <- function(assg){
  n <- V(assg)$name
  f.rem <- function(assg, i){
    return(delete.vertices(assg,n[i]))
  }
  #dopar here
  fpc.valsT <- foreach(i=113:170,.export=c('deltaS','runCS','runRichClub','rccW'),.packages=c('igraph','doParallel')) %dopar% (rccW(f.rem(assg,i)))
  #fpc.vals <- c(0,fpc.valsT[1,])
  
  return(fpc.valsT)
}
rcc.val3 <- rcc.validate3(exoEndo.graph)
save(rcc.val3, file="EECS459/FullValRC3.RData")

rcc.validate4 <- function(assg){
  n <- V(assg)$name
  f.rem <- function(assg, i){
    return(delete.vertices(assg,n[i]))
  }
  #dopar here
  fpc.valsT <- foreach(i=171:length(n),.export=c('deltaS','runCS','runRichClub','rccW'),.packages=c('igraph','doParallel')) %dopar% (rccW(f.rem(assg,i)))
  #fpc.vals <- c(0,fpc.valsT[1,])
  
  return(fpc.valsT)
}
rcc.val4 <- rcc.validate4(exoEndo.graph)
save(rcc.val4, file="EECS459/FullValRC4.RData")
