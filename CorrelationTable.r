setwd("~/GitHub/RSA-normalization-values/")
wolfenden<- read.delim("Wolfden.txt")
wolf<- wolfenden$Hydro

Rose<-read.delim("rose.txt")
rose<-Rose$Hydro

Kyte<-read.delim("Kite_Doolittle.txt")
kyte<-Kyte$Hydro

Scales<- read.csv("Hydrophobicity_Scales_Updated.txt", row.names=1)

cor.wolf.rose<-cor.test(wolf, rose)
cor.wolf.empMean<-cor.test(wolf, Scales$meanEmp)
cor.wolf.theoMean<-cor.test(wolf, Scales$meanTheo)
cor.wolf.empMedian<-cor.test(wolf, Scales$medianEmp)
cor.wolf.theoMedian<-cor.test(wolf, Scales$medianTheo)
cor.wolf.empSqrt<-cor.test(wolf, Scales$meanEmpSqrt)
cor.wolf.theoSqrt<-cor.test(wolf, Scales$meanTheoSqrt)
cor.wolf.empBC<-cor.test(wolf, meanEmpBC)
cor.wolf.theoBC<-cor.test(wolf, meanTheoBC)
cor.wolf.fraction100<-cor.test(wolf, Scales$ScaleFrac)
cor.wolf.fraction95<-cor.test(wolf, Scales$ScaleFrac95)

wolfCorr<- c(cor.wolf.rose$estimate,
cor.wolf.empMean$estimate,
cor.wolf.theoMean$estimate,
cor.wolf.empMedian$estimate,
cor.wolf.theoMedian$estimate,
cor.wolf.empSqrt$estimate,
cor.wolf.theoSqrt$estimate,
cor.wolf.empBC$estimate,
cor.wolf.theoBC$estimate,
cor.wolf.fraction100$estimate,
cor.wolf.fraction95$estimate)

wolfPval<- c(cor.wolf.rose$p.value,
cor.wolf.empMean$p.value,
cor.wolf.theoMean$p.value,
cor.wolf.empMedian$p.value,
cor.wolf.theoMedian$p.value,
cor.wolf.empSqrt$p.value,
cor.wolf.theoSqrt$p.value,
cor.wolf.empBC$p.value,
cor.wolf.theoBC$p.value,
cor.wolf.fraction100$p.value,
cor.wolf.fraction95$p.value)

cor.kyte.rose<-cor.test(kyte, rose)
cor.kyte.empMean<-cor.test(kyte, Scales$meanEmp)
cor.kyte.theoMean<-cor.test(kyte, Scales$meanTheo)
cor.kyte.empMedian<-cor.test(kyte, Scales$medianEmp)
cor.kyte.theoMedian<-cor.test(kyte, Scales$medianTheo)
cor.kyte.empSqrt<-cor.test(kyte, Scales$meanEmpSqrt)
cor.kyte.theoSqrt<-cor.test(kyte, Scales$meanTheoSqrt)
cor.kyte.empBC<-cor.test(kyte, meanEmpBC)
cor.kyte.theoBC<-cor.test(kyte, meanTheoBC)
cor.kyte.fraction100<-cor.test(kyte, Scales$ScaleFrac)
cor.kyte.fraction95<-cor.test(kyte, Scales$ScaleFrac95)

kyteCorr<- c(cor.kyte.rose$estimate,
cor.kyte.empMean$estimate,
cor.kyte.theoMean$estimate,
cor.kyte.empMedian$estimate,
cor.kyte.theoMedian$estimate,
cor.kyte.empSqrt$estimate,
cor.kyte.theoSqrt$estimate,
cor.kyte.empBC$estimate,
cor.kyte.theoBC$estimate,
cor.kyte.fraction100$estimate,
cor.kyte.fraction95$estimate)

kytePval<- c(cor.kyte.rose$p.value,
cor.kyte.empMean$p.value,
cor.kyte.theoMean$p.value,
cor.kyte.empMedian$p.value,
cor.kyte.theoMedian$p.value,
cor.kyte.empSqrt$p.value,
cor.kyte.theoSqrt$p.value,
cor.kyte.empBC$p.value,
cor.kyte.theoBC$p.value,
cor.kyte.fraction100$p.value,
cor.kyte.fraction95$p.value)

labels<-c("Rose", "Empirical Mean", "Theoretical Mean", "Empirical Median", "Theoretical Median", "Empirical Sqrt Mean", "Theoretical Sqrt Mean", "Empirical Box-Cox Mean", "Theoretical Box-Cox Mean", "100% Buried", "95% Buried")

Correlations<- data.frame(labels=labels, wolfCorr= wolfCorr, wolfPval=wolfPval, kyteCorr=kyteCorr, kytePval=kytePval)

