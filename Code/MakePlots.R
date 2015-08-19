##############################################################
#MakePlots.R
#Purpose: script that makes plots for results of all analysis
#############################################################
library(ggplot2)
setwd("Z:/EECS459")

load("FullRCC.RData")
dgF <- names(rccEE)
rc <- data.frame(cbind(rccEE,dgF))
rc$rccEE <- as.numeric(as.character(rc$rccEE))
rc$dgF <- as.numeric(as.character(rc$dgF))
ggplot(rc,aes(x=dgF,y=rccEE))+geom_point()+xlim(0,30)+ylim(0,3.2)+labs(x="Degree",y="Normalized Rich-Club Coefficient",title="Exocytosis and EndoCytosis Rich-Club Coefficients")

load("ExoRCC.RData")
dgF2 <- names(rccEXO)
rc2 <- data.frame(cbind(rccEXO,dgF2))
rc2$rccEXO <- as.numeric(as.character(rc2$rccEXO))
rc2$dgF <- as.numeric(as.character(rc2$dgF))
ggplot(rc2,aes(x=dgF,y=rccEXO))+geom_point()+xlim(0,10)+ylim(0,3.2)+labs(x="Degree",y="Normalized Rich-Club Coefficient",title="Exocytosis Rich-Club Coefficients")

load("EndoRCC.RData")
dgF3 <- names(rccENDO)
rc3 <- data.frame(cbind(rccENDO,dgF3))
rc3$rccENDO <- as.numeric(as.character(rc3$rccENDO))
rc3$dgF <- as.numeric(as.character(rc3$dgF))
ggplot(rc3,aes(x=dgF,y=rccENDO))+geom_point()+xlim(0,25)+ylim(0,3)+labs(x="Degree",y="Normalized Rich-Club Coefficient",title="Endocytosis Rich-Club Coefficients")

rc$Type <- "Full"
rc2$Type <- "Exocytosis"
rc3$Type <- "Endocytosis"
colnames(rc3) <- c("Normalized Rich-Club Coefficient","Degree","Network")
load("FullFPC.RData")
load("fullNetwork.RData")
names(fpc.score) <- V(exoEndo.graph)$name
fpc.score2 <- fpc.score[order(fpc.score)]
write.csv(as.data.frame(fpc.score2),"~/14-15//Spring15//EECS459//Project//FullFPC.csv",quote=FALSE,row.names=TRUE)

load("FullValFPC.RData")
fpcdf <- data.frame(cbind(fpc.score,fpc.val))
ggplot(fpcdf,aes(x=fpc.score,y=fpc.val))+geom_point()+labs(x="FPC Score",y="FPC Validation")

load("ExoFPC.RData")
load("EXONetwork.RData")
fpc.EXO2 <- fpc.EXO
names(fpc.EXO2) <- V(exo.graph)$name
fpc.EXO2 <- fpc.EXO2[order(fpc.EXO2)]
write.csv(as.data.frame(fpc.EXO2),"~/14-15/Spring15/EECS459/Project/ExoFPC.csv",quote=FALSE,row.names=TRUE)

load("ExoValFPC.RData")
fpcdf2 <- data.frame(cbind(fpc.EXO,fpc.valX))
ggplot(fpcdf2,aes(x=fpc.EXO,fpc.valX))+geom_point()+labs(x="FPC Score",y="FPC Validation")

load("EndoFPC.RData")
load("ENDONetwork.RData")
fpc.ENDO2 <- fpc.ENDO
names(fpc.ENDO2) <- V(endo.graph)$name
fpc.ENDO2 <- fpc.ENDO2[order(fpc.ENDO2)]
write.csv(as.data.frame(fpc.ENDO2),"~/14-15/Spring15/EECS459/Project/EndoFPC.csv",quote=FALSE,row.names=TRUE)

load("EndoValFPC.RData")
fpcdf3 <- data.frame(cbind(fpc.ENDO,fpc.valENDO))
ggplot(fpcdf3,aes(x=fpc.ENDO,y=fpc.valENDO))+geom_point()+labs(x="FPC Score",y="FPC Validation")

load("FullRichClub.RData")
library(corrplot)
corrplot(finalList[[1]],is.corr=FALSE)
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F","cyan", "#007FFF", "blue", "#00007F"))
corrplot(finalList[[1]],is.corr=FALSE,method="color",col = col4(210),cl.length=210,mar=c(1,1,1,1),cl.pos="n")
dev.off()
rc.full <- finalList[[2]][order(finalList[[2]])]
write.csv(as.data.frame(rc.full),"~/14-15/Spring15/EECS459/Project/FullRichClub.csv",quote=F,row.names=T)

load("ExoRichClub.RData")
rcE <- finalListX[[2]][order(finalListX[[2]])]
rcE <- rcE/max(rcE)
write.csv(as.data.frame(rcE),"~/14-15/Spring15/EECS459/Project/ExoRichClub.csv",quote=F,row.names=T)

load("EndoRichClub2.RData")
rcE2 <- finalListE[[2]][order(finalListE[[2]])]

write.csv(as.data.frame(rcE2),"~/14-15/Spring15/EECS459/Project/EndoRichClub.csv",quote=F,row.names=T)

load("EndoValRC.RData")
rcdf <- data.frame(cbind(finalListE[[2]],rcc.valENDO))
ggplot(rcdf,aes(x=V1,y=rcc.valENDO))+geom_point()+labs(x="Core/Periphery Score",y="Core/Periphery Validation")

load("ExoValRichClub.RData")
#rcc.valX <- rcc.valX/max(rcc.valX)
finalListX[[2]]<- finalListX[[2]]/max(finalListX[[2]])
rcdf2 <- data.frame(cbind(finalListX[[2]],rcc.valX))
ggplot(rcdf2,aes(x=V1,y=rcc.valX))+geom_point()+labs(x="Core/Periphery Score",y="Core/Periphery Validation")

corrs <- matrix(nrow=0,ncol=2)
colnames(corrs) <- c("Correlation","Comparison")
corrs <- rbind(corrs,c(cor(finalList[[2]],fpc.score),"FPC C/P Full"))
corrs <- rbind(corrs,c(cor(finalListE[[2]],fpc.ENDO),"FPC C/P Endocytosis"))
corrs <- rbind(corrs,c(cor(finalListX[[2]],fpc.EXO),"FPC C/P Exocytosis"))
corrs <- rbind(corrs,c(cor(fpc.score,fpc.val),"FPC FPC.Val Full"))
corrs <- rbind(corrs,c(cor(fpc.ENDO,fpc.valENDO),"FPC FPC.Val Endocytosis"))
corrs <- rbind(corrs,c(cor(fpc.EXO,fpc.valX),"FPC FPC.Val Exocytosis"))
corrs <- rbind(corrs,c(cor(finalListX[[2]],rcc.valX), "C/P C/P.Val Exocytosis"))
corrs <- rbind(corrs,c(cor(finalListE[[2]],rcc.valENDO),"C/P C/P.Val Endocytosis"))

write.csv(corrs, "~/14-15/Spring15/EECS459/Project/CorrelationResults.csv",quote=FALSE,row.names=FALSE)
ggplot(data.frame(cbind(finalList[[2]],fpc.score)),aes(x=V1,y=fpc.score))+geom_point()+labs(x="Core/Periphery Score",y="FPC Score",title="Full Network")
ggplot(data.frame(cbind(finalListE[[2]],fpc.ENDO)),aes(x=V1,y=fpc.ENDO))+geom_point()+labs(x="Core/Periphery Score",y="FPC Score",title="Endocytosis")
ggplot(data.frame(cbind(finalListX[[2]]/max(finalListX[[2]]),fpc.EXO)),aes(x=V1,y=fpc.EXO))+geom_point()+labs(x="Core/Periphery Score",y="FPC Score",title="Exocytosis")


V(exoEndo.graph)[c("COPS5","LYN","SRC","GRB2","EGFR")]$color <- "red"
V(exo.graph)$color <- "lightblue"
V(exo.graph)[c("RPH3AL","STXBP1","STX1A","SYTL4","RAB27B")]$color <- "red"
V(endo.graph)$color <- "lightblue"
V(endo.graph)[c("TNK2","CDC42","SRC","GRB2","EGFR")]$color <- "red"

fullVal <- c(rcc.val1,rcc.val2)
rcc.valT2 <- do.call(rbind,fullVal)
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
