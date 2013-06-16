data.path <- "../Wilke/EmpVCalc/"

code <- "ALA"
tV<- paste(code, "_DifferenceVPopulation", sep='')

fV<- paste(tV, ".pdf", sep='')

pdf(fV, width=10, height=8)
par(mar=c(7,7,5,5))
par(mgp=c(3,1,0))

fileV<-paste(data.path, "EmpericalVCalculated_diff_pop_nonZeroed_with_pop_restriction_", code, sep='')
d<-read.delim(fileV)
q<-plot(x=d$pop, y=d$SA_diff, log="x", xlab= " ", ylab=" ", las=1, pch=16, frame=FALSE, cex.axis=1.8, cex=1.2)
##points(x=d$pop, y=d$SA_diff, pch=16, col="gray")
mtext(expression(paste(Delta, "(Maximum SA per Bin)" )), side=2, line=4.5, outer= FALSE, cex=1.8)
mtext("Number of Residues per Bin", side=1, line=4.5, outer=FALSE, cex=1.8)
dev.off()
