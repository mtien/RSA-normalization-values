code<-"Ala"

setwd("~/GitHub/RSA-normalization-values/GeoFiles/ALLOWED/")
f_all<-paste(code, "_ALLOWED_geo", sep='')
d<-read.delim(f_all)
h<-hist(d$SA/129.0, plot=FALSE)

setwd("~/GitHub/RSA-normalization-values/Figures")
pdf("Alanine_RSA_distribution.pdf", width=11.5, height=9)
par(mar=c(7,9.5,5,5))
par(mgp=c(3,2,0))

plot(h, main="", ylim=c(0,20000),xlab= " ", ylab=" ", cex.axis=3.0, xaxt='n', xlim= c(0,1.0), bty="l")
axis(1, c(0,.2, .4, .6, .8, 1.0), tck=-.01, cex.axis=3.0)
mtext("Frequency", side=2, line=6.5, cex=3.0)
mtext("RSA (normalized by max theoretical)", side=1, line=4.5, cex=3.0)
dev.off()
