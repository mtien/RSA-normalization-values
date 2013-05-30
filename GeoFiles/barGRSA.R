
setwd("~/GitHub/RSA-normalization-values/GeoFiles/")
AA<-c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Pro", "Phe", "Ser", "Thr", "Trp", "Tyr", "Val")
aa<-c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "P", "F", "S", "T", "W", "Y", "V")
maxRSA1<-c()
maxRSA2<-c()

for(a in 1:length(AA))
{

theo1<-paste(AA[a], "_geo", sep='')
theo2<-paste(AA[a], "_Rose_RSA", sep='')
d1<-read.delim(theo1)
d2<-read.delim(theo2)
temp1<-d1$RSA
temp2<-d2$RSA

count1<-0.0
count2<-0.0

for( i in 1:length(temp1))
{
if(temp1[i]>1.0)
{
count1=count1+1.0
}
}

for( i in 1:length(temp2))
{
if(temp2[i]>1.0)
{
count2=count2+1.0
}
}

maxRSA1<-c(maxRSA1, (count1/length(temp1)))
maxRSA2<-c(maxRSA2, (count2/length(temp2)))

}

color1<-colors()[23]
color2<-colors()[24]

colors<-c(color1, color2)
maxRSA<-c(maxRSA1, maxRSA2)*100

RSAmat<- matrix(maxRSA, ncol=20, byrow=FALSE)
num<-2.75

pdf("BarGraphRSA.pdf", width=12, height=8)
##png("BarGraphRoseRSA.png", width=1280, height=1024, units= "px" )
##BOTTOM, LEFT, TOP RIGHT
par(mar=c(6,6,2,2))
##LABEL, AXIS, TICKSpar(mar=c(4,0,0,0))
par(mgp=c(4,1,0))

barplot(RSAmat, xlab= "Amino Acids", ylab= "Percent of Residues with RSA > 1", beside=TRUE, names.arg=aa ,col=colors, cex.names=1.77, cex.lab=num, ylim=c(0,3.0), cex.axis=2.0, las=1)
legend("topright", c("Miller (1987)", "Rose (1985)"), fill=colors, cex=num)
dev.off()

##color1<-colors()[23]
##color2<-colors()[24]
##colors<-c(color1, color2)
##png("BarGraphRSA.png", width=1600, height=1000, units= "px")
##BOTTOM, LEFT, TOP RIGHT
##par(mar=c(12,12,4,2))
##par(mar=c(4,0,0,0))
##LABEL, AXIS, TICKS
##par(mgp=c(7,2,0))

##barplot(RSAmat, xlab= "Amino Acids", ylab= "Percent of Residues with RSA > 1", beside=TRUE, names.arg=aa ,col=colors, cex.names=3.0, cex.lab=4.0, ylim=c(0,3.0), cex.axis=3.0, las=1)
##legend("topright", c("Miller (1987)", "Rose (1985)"), fill=colors, cex=4.0)
##
##dev.off()
