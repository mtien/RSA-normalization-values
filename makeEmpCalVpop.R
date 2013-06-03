code<-"Asp"
tV<- paste(code, "_DifferenceVPopulation", sep='')
fV<- paste( tV, ".pdf", sep='')

setwd("~/GitHub/RSA-normalization-values/EmpiricalVTheoretical/")
fileV<-paste("EmpericalVCalculated_diff_pop_nonZeroed_", code,sep='')
d<-read.delim(fileV)

setwd("~/GitHub/RSA-normalization-values/Figures/")
pdf(fV, width=11.5, height=10)
##png(fV)
par(mar=c(7,9,5,5))
par(mgp=c(3,2,0))
plot(x=d$pop, y=d$SA_diff, log="x", xlab= " ", ylab=" ", las=1, pch=16, cex.axis=3.0, cex=2, bty="l")
##points(x=d$pop, y=d$SA_diff, pch=16, col="gray"
axis(1, c(1,2,5,10,50,200), tck=0, cex.axis=3.0) 
mtext("Difference Between Maximum SA", side=2, line=6.5, cex=3.0)
mtext("Number of RSA Values per bin", side=1, line=4.5, cex=3.0)
dev.off()
