#!/usr/bin/Rscript

ROC=read.csv("./ROC",sep="\t",header=F)
colnames(ROC)=c("x","y")

jpeg(file = "plot.jpg")

#plot(ROC,xlab="1-Sp",ylab="Sn",cex=1,pch=16,col="red",xlim=c(0:1),ylim=c(0:1))
#把cex在设置成0就隐藏了所有的点
plot(ROC[,1:2],xlab="1-Sp",ylab="Sn",cex=0.4,pch=16,col="yellow",xlim=c(0:1),ylim=c(0:1))

My.loess=loess(x ~ y, data=ROC)
Fit=fitted(My.loess)
lines(Fit,ROC$y,col="red")
arrows(0.1,0,0.1,1,angle=0,col="green");
arrows(0,0,1,1,angle=0,col="blue")

dev.off()