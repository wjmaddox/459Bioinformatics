########################################################
#Name: RCCExample.R
#Author: Wesley Maddox
#Date: 3/27/15
#Purpose: Test example of rich-club coefficient metrics
########################################################

library(igraph)
source("C://Users//Wesley//Documents//14-15//Spring15//EECS459//Project//Code//RCCFns.R")
assg <- graph.formula(a-f,a-b-c-d-e-f-g,g-a-h-i-j,g-b-k-l-m,
                      g-c-n-o-p,g-d-q-r-s,g-e-t-v,g-f-w-x,b-e-s)
#degree ranges from 1to6
phi <- c()
for (x in 1:5){
  phi <- c(phi,RCC(x,assg))
}

rho <- c()
for (x in 1:5){
  rho <- c(rho, normalize.RCC(assg,10,x))
}

plot(phi)
plot(rho)

#######################################################
#FPC
#######################################################
c1 <- degree(assg)
c2 <- betweenness(assg)
c3 <- closeness(assg)
c4 <- transitivity(assg,type="localundirected")
c5 <- graph.coreness(assg)
c6 <- evcent(assg)$vector
c7 <- page.rank(assg)$vector #using pageRank instead of semi-local centrality
#c8 <- 

c <- cbind(c1,c2,c3,c4,c5,c6,c7)
c[is.na(c)]<-0
cv <- cov(c)
w <- eigen(cv)$values

fpc <- c()
for (x in 1:dim(c)[1]){
fpc <- c(fpc,sum(w*c[x,]))
}
names(fpc) <- names(degree(assg))
names(which(fpc==max(fpc)))
plot(fpc)

t <- tail(names(fpc[order(fpc)]))
assg.edge <- get.edgelist(assg)

t.g <- graph.edgelist(assg.edge[intersect(which(assg.edge[,1]%in%t),which(assg.edge[,2]%in%t)),],
                      directed=F)
tkplot(t.g)
