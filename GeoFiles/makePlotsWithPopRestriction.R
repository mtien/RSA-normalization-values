setwd("~/GitHub/RSA-normalization-values/GeoFiles/")
code<- "Asp"
fileC<-paste(code, "_max_theoretical_bins_Again",sep='')
fileE<-paste(code, "_max_emperical_bins_pop_restriction",sep='')
##fileV<-paste("EmpericalVCalculated_diff_pop_nonZeroed_with_pop_restriction_", code,sep='')

##tC<- paste(nom, " Calculated Weighted", sep='')
##fC<- paste(tC, "POP.png", sep='')

##tE<- paste(nom, " Emperical Weighted", sep='')
##fE<- paste(tE, "POP.png", sep='')

##tV<- paste(nom, " Calculated and Emperical Difference against Emperical Population", sep='')
##fV<- paste(tV, "POP.png", sep='')

##png(fC)
##source("makeHeatplotCalc.r")
##dev.off()

##png(fE)
##source("makeHeatplotEmp.r")
##dev.off()

##png(fV)
##par(mar=c(7,7,5,5))
##par(mgp=c(3,1,0))
##source("makeEmpCalVpop.r")
##dev.off()

dataC<-read.delim(fileC)
mnum<-max(dataC$maxSA)
zC<-matrix(dataC$maxSA, 72, byrow=T)

dataE<-read.delim(fileE)
zE<-matrix(dataE$maxSA, 72, byrow=T)

fA<- paste(code, "_Rama_HSV.svg", sep='')

##rspec<-c(0)
##for(r in 1:mnum)
##{
##rspec<-c(rspec, rgb(r/mnum,0, (mnum-r)/mnum))
##}

rspec<-c(rgb(1,1,1), hsv(seq(.9, 1.0, length=56),1,1), hsv(seq(0, .15, length=40),1,1) )

##png(fA, width=1600, height=600, units="px")
svg(fA, width=16, height=6)
##BOTTOM, LEFT, TOP RIGHT
par(mar=c(8,9,8,5))
##LABEL, AXIS, TICKS
par(mgp=c(0,1.5,0))
##par(adj=c(10,10,10,10))
par(mfrow=c(1,2))
par(xpd=TRUE, mar=par()$mar+c(0,0,0,10))
##image(seq(-180,180,5), seq(-180,180,5), zE, zlim=c(0,mnum), col=rev(gray(0:200.0/200.0)), xlab="", ylab="", cex.lab=4.0, cex.axis=2.0, main="A", cex.main=4.0, las=1)
image(seq(-180,180,5), seq(-180,180,5), zE, zlim=c(0,mnum), col=rspec, xlab="", ylab="", cex.lab=2.0, cex.axis=1.5, cex.main=2.0, las=1)
mtext(expression(psi), side=2, line=3.75, outer= FALSE, las=1, cex=1.5)
mtext(expression(phi), side=1, line=3.25, outer=FALSE, cex=1.5)
mtext("A", side=3, line=2, outer= FALSE, adj=0, cex=1.5, font=2)
##gray(log(1:10)/max(log(1:10)))
##rspec2<-c(0)
##for(r in 1:10)
##{
##rspec2<-c(rspec2, heat.colors(r))
##rspec2<-c(rspec2, rgb( log(r)/max(log(1:10)),0,log(11-r)/max(log(1:10)) ) )
##}




par(mar=c(8,9,8,5))
par(mgp=c(0,1.5,0))
par(xpd=TRUE, mar=par()$mar+c(0,0,0,10))
##image(seq(-180,180,1), seq(-180,180,1), zC, zlim=c(0,mnum), col=rev(gray(0:200.0/200.0)), xlab="", ylab=", cex.lab=4.0, cex.axis=2.0, main="B", cex.main=4.0, las=1)
image(seq(-180,180,5), seq(-180,180,5), zC, zlim=c(0,mnum), col=rspec, xlab="", ylab="", cex.lab=2.0, cex.axis=1.5, cex.main=2.0, las=1)
mtext(expression(psi), side=2, line=3.75, outer= FALSE, las=1, cex=1.5)
mtext(expression(phi), side=1, line=3.0, outer=FALSE, cex=1.5)
mtext("B", side=3, line=2, outer= FALSE, adj=0, cex=1.5, font=2)
##colorlegend(posx=c(.91, .96), col=rev, zlim=c(0,mnum), zval=c(0, mnum))

#345
rspec2<-c()
for(r in 0:8)
  {
  rspec2<-c(rspec2, rspec[(12*r) +1])
  }

legend(345,
       200, 
       y.intersp=.5,
       legend=c("", "", "", "", "", "", "", "", ""), 
       fill=rev(rspec2), 
       border="black",
       cex=3,
       bty="n")

text(378,-172,"0", cex=1.5, xpd = TRUE)
##text(260,0,"30", cex=1.5, xpd = TRUE)
##text(260,20,"30", cex=1.5, xpd = TRUE)
text(370,-15,as.character(mnum/2), cex=1.5, xpd = TRUE)
##text(260,50,"90", cex=1.5, xpd = TRUE)
text(360,136,as.character(mnum), cex=1.5, xpd = TRUE)
text(398,185,"SA Scale", cex=1.5, xpd = TRUE)


dev.off()
