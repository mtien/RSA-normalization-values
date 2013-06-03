
AA<-c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Pro", "Phe", "Ser", "Thr", "Trp", "Tyr", "Val")
aa<-c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "P", "F", "S", "T", "W", "Y", "V")
maxRSA1<-c()
maxRSA2<-c()

for(a in 1:length(AA))
{
theo1<-paste(AA[a], "_geo", sep='')
theo2<-paste(AA[a], "_Rose_RSA", sep='')

setwd("~/GitHub/RSA-normalization-values/GeoFiles/")
d1<-read.delim(theo1)

setwd("~/GitHub/RSA-normalization-values/RoseRSA/")
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

pdf("~/GitHub/RSA-normalization-values/Figures/BarGraphRSA.pdf", width=12, height=8)
par(mar=c(6,6,2,2))
##LABEL, AXIS, TICKSpar(mar=c(4,0,0,0))
par(mgp=c(4,1,0))

barplot(RSAmat, xlab= "Amino Acids", ylab= "Percent of Residues with RSA > 1", beside=TRUE, names.arg=aa ,col=colors, cex.names=1.77, cex.lab=num, ylim=c(0,3.0), cex.axis=2.0, las=1)
legend("topright", c("Miller (1987)", "Rose (1985)"), fill=colors, cex=num)
dev.off()
