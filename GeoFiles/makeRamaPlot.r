
setwd("~/GitHub/RSA-normalization-values/GeoFiles/")
code<-"Asp"
fileOver<-paste(code, "_SA_Over",sep='')
fileUnder<-paste(code, "_SA_Under",sep='')

Over<-read.delim(fileOver)
Under<-read.delim(fileUnder)

fileOut<-paste(code, "_RamaPlotRSA.pdf", sep='')

pdf(fileOut, width=10.5, height=10)
##png(fileOut, width=1000, height=1000, units= "px")
##BOTTOM, LEFT, TOP RIGHT
par(mar=c(6,7,2,6))
##LABEL, AXIS, TICKS
par(mgp=c(0,1.5,0))
rcol<-rgb(150/255, 150/255, 150/255)
rcol2<-rgb(225/255, 225/255, 225/255)
black<-rgb(0,0,0)
plot(Under$Phi*180/pi, Under$Psi*180/pi, pch=16, col=rcol, xlab="", ylab="", cex.axis=2.75, las=1) 
points(Over$Phi*180/pi, Over$Psi*180/pi, pch=21, col=rcol2, bg="black", cex=1.5)
mtext(expression(psi), side=2, line=4.5, outer= FALSE, las=1, cex=2.75)
mtext(expression(phi), side=1, line=4.5, outer=FALSE, cex=2.75)

legend("topright", c("RSA < 1", "RSA > 1"), bg=rgb(1,1,1), cex=2.75, col=c(rcol,black), pch=c(16,16), box.lwd = 1)
#box.lty=1,box.col="white", bg="white"
dev.off()

