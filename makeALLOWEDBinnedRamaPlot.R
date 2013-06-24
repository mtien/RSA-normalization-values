code<- "Ala"
fileA<-paste(code, "_max_bins_ALLOWED",sep='')

setwd("~/GitHub/RSA-normalization-values/AllowedBins/")

dataA<-read.delim(fileA)
mnum<-max(dataA$max_theo_SA)
zC<-matrix(dataA$max_theo_SA, 72, byrow=T)
zE<-matrix(dataA$max_obs_SA, 72, byrow=T)

setwd("~/GitHub/RSA-normalization-values/Figures/")

fA<- paste(code, "_Rama_HSV_ALLOWED.svg", sep='')

rspec<-c(rgb(1,1,1), hsv(seq(.9, 1.0, length=56),1,1), hsv(seq(0, .15, length=40),1,1) )

svg(fA, width=16, height=6)
par(mar=c(8,9,8,5))
par(mgp=c(0,1.5,0))
par(mfrow=c(1,2))
par(xpd=TRUE, mar=par()$mar+c(0,0,0,10))
image(seq(-180,180,5), seq(-180,180,5), zE, zlim=c(0,mnum), col=rspec, xlab="", ylab="", cex.lab=2.0, cex.axis=1.5, cex.main=2.0, las=1)
mtext(expression(psi), side=2, line=3.75, outer= FALSE, las=1, cex=1.5)
mtext(expression(phi), side=1, line=3.25, outer=FALSE, cex=1.5)
mtext("A", side=3, line=2, outer= FALSE, adj=0, cex=1.5, font=2)

par(mar=c(8,9,8,5))
par(mgp=c(0,1.5,0))
par(xpd=TRUE, mar=par()$mar+c(0,0,0,10))
image(seq(-180,180,5), seq(-180,180,5), zC, zlim=c(0,mnum), col=rspec, xlab="", ylab="", cex.lab=2.0, cex.axis=1.5, cex.main=2.0, las=1)
mtext(expression(psi), side=2, line=3.75, outer= FALSE, las=1, cex=1.5)
mtext(expression(phi), side=1, line=3.0, outer=FALSE, cex=1.5)
mtext("B", side=3, line=2, outer= FALSE, adj=0, cex=1.5, font=2)


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
text(370,-15,as.character(mnum/2), cex=1.5, xpd = TRUE)
text(360,136,as.character(mnum), cex=1.5, xpd = TRUE)
text(398,185,"SA Scale", cex=1.5, xpd = TRUE)


dev.off()
