setwd("~/GitHub/RSA-normalization-values/Correlation/")
wolfenden<- read.delim("Wolfden.txt")
wolf<- wolfenden$Hydro

Kyte<-read.delim("Kite_Doolittle.txt")
kyte<-Kyte$Hydro

Radzicka<- read.delim("Radzicka.txt")
radz<- Radzicka$Hydro

MacCallum<-read.delim("MacCallum.txt")
mac<-MacCallum$Hydro

Moon<-read.delim("Moon.txt")
moon<-Moon$Hydro

Wimley<-read.delim("Wimley.txt")
wim<-Wimley$Hydro

Fauchere<- read.delim("Fauchere.txt")
fau<-Fauchere$Hydro

Rose<-read.delim("rose.txt")
rose<-Rose$Hydro
rose19<-c(rose[1:14], rose[16:20])
rose17<-c(rose[1:7], rose[10:14], rose[16:20])
Scales<- read.table("Hydrophobicity_Scales_Updated.txt", row.names=1, header=TRUE)
Scales19<- read.table("Hydrophobicity_Scales_Updated19.txt", row.names=1, header=TRUE)
Scales17<- read.table("Hydrophobicity_Scales_Updated17.txt", row.names=1, header=TRUE)

cor.wolf.rose<-cor.test(wolf, rose)
cor.wolf.empMean<-cor.test(wolf, Scales$meanEmp)
cor.wolf.theoMean<-cor.test(wolf, Scales$meanTheo)
cor.wolf.empMedian<-cor.test(wolf, Scales$medianEmp)
cor.wolf.theoMedian<-cor.test(wolf, Scales$medianTheo)
cor.wolf.empSqrt<-cor.test(wolf, Scales$meanEmpSqrt)
cor.wolf.theoSqrt<-cor.test(wolf, Scales$meanTheoSqrt)
cor.wolf.empBC<-cor.test(wolf, Scales$meanEmpBC)
cor.wolf.theoBC<-cor.test(wolf, Scales$meanTheoBC)
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
cor.kyte.empBC<-cor.test(kyte, Scales$meanEmpBC)
cor.kyte.theoBC<-cor.test(kyte, Scales$meanTheoBC)
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

cor.fau.rose<-cor.test(fau, rose)
cor.fau.empMean<-cor.test(fau, Scales$meanEmp)
cor.fau.theoMean<-cor.test(fau, Scales$meanTheo)
cor.fau.empMedian<-cor.test(fau, Scales$medianEmp)
cor.fau.theoMedian<-cor.test(fau, Scales$medianTheo)
cor.fau.empSqrt<-cor.test(fau, Scales$meanEmpSqrt)
cor.fau.theoSqrt<-cor.test(fau, Scales$meanTheoSqrt)
cor.fau.empBC<-cor.test(fau, Scales$meanEmpBC)
cor.fau.theoBC<-cor.test(fau, Scales$meanTheoBC)
cor.fau.fraction100<-cor.test(fau, Scales$ScaleFrac)
cor.fau.fraction95<-cor.test(fau, Scales$ScaleFrac95)

fauCorr<- c(cor.fau.rose$estimate,
cor.fau.empMean$estimate,
cor.fau.theoMean$estimate,
cor.fau.empMedian$estimate,
cor.fau.theoMedian$estimate,
cor.fau.empSqrt$estimate,
cor.fau.theoSqrt$estimate,
cor.fau.empBC$estimate,
cor.fau.theoBC$estimate,
cor.fau.fraction100$estimate,
cor.fau.fraction95$estimate)

fauPval<- c(cor.fau.rose$p.value,
cor.fau.empMean$p.value,
cor.fau.theoMean$p.value,
cor.fau.empMedian$p.value,
cor.fau.theoMedian$p.value,
cor.fau.empSqrt$p.value,
cor.fau.theoSqrt$p.value,
cor.fau.empBC$p.value,
cor.fau.theoBC$p.value,
cor.fau.fraction100$p.value,
cor.fau.fraction95$p.value)


cor.radz.rose<-cor.test(radz, rose19)
cor.radz.empMean<-cor.test(radz, Scales19$meanEmp)
cor.radz.theoMean<-cor.test(radz, Scales19$meanTheo)
cor.radz.empMedian<-cor.test(radz, Scales19$medianEmp)
cor.radz.theoMedian<-cor.test(radz, Scales19$medianTheo)
cor.radz.empSqrt<-cor.test(radz, Scales19$meanEmpSqrt)
cor.radz.theoSqrt<-cor.test(radz, Scales19$meanTheoSqrt)
cor.radz.empBC<-cor.test(radz, Scales19$meanEmpBC)
cor.radz.theoBC<-cor.test(radz, Scales19$meanTheoBC)
cor.radz.fraction100<-cor.test(radz, Scales19$ScaleFrac)
cor.radz.fraction95<-cor.test(radz, Scales19$ScaleFrac95)

radzCorr<- c(cor.radz.rose$estimate,
cor.radz.empMean$estimate,
cor.radz.theoMean$estimate,
cor.radz.empMedian$estimate,
cor.radz.theoMedian$estimate,
cor.radz.empSqrt$estimate,
cor.radz.theoSqrt$estimate,
cor.radz.empBC$estimate,
cor.radz.theoBC$estimate,
cor.radz.fraction100$estimate,
cor.radz.fraction95$estimate)

radzPval<- c(cor.radz.rose$p.value,
cor.radz.empMean$p.value,
cor.radz.theoMean$p.value,
cor.radz.empMedian$p.value,
cor.radz.theoMedian$p.value,
cor.radz.empSqrt$p.value,
cor.radz.theoSqrt$p.value,
cor.radz.empBC$p.value,
cor.radz.theoBC$p.value,
cor.radz.fraction100$p.value,
cor.radz.fraction95$p.value)

cor.moon.rose<-cor.test(moon, rose)
cor.moon.empMean<-cor.test(moon, Scales$meanEmp)
cor.moon.theoMean<-cor.test(moon, Scales$meanTheo)
cor.moon.empMedian<-cor.test(moon, Scales$medianEmp)
cor.moon.theoMedian<-cor.test(moon, Scales$medianTheo)
cor.moon.empSqrt<-cor.test(moon, Scales$meanEmpSqrt)
cor.moon.theoSqrt<-cor.test(moon, Scales$meanTheoSqrt)
cor.moon.empBC<-cor.test(moon, Scales$meanEmpBC)
cor.moon.theoBC<-cor.test(moon, Scales$meanTheoBC)
cor.moon.fraction100<-cor.test(moon, Scales$ScaleFrac)
cor.moon.fraction95<-cor.test(moon, Scales$ScaleFrac95)

moonCorr<- c(cor.moon.rose$estimate,
cor.moon.empMean$estimate,
cor.moon.theoMean$estimate,
cor.moon.empMedian$estimate,
cor.moon.theoMedian$estimate,
cor.moon.empSqrt$estimate,
cor.moon.theoSqrt$estimate,
cor.moon.empBC$estimate,
cor.moon.theoBC$estimate,
cor.moon.fraction100$estimate,
cor.moon.fraction95$estimate)

moonPval<- c(cor.moon.rose$p.value,
cor.moon.empMean$p.value,
cor.moon.theoMean$p.value,
cor.moon.empMedian$p.value,
cor.moon.theoMedian$p.value,
cor.moon.empSqrt$p.value,
cor.moon.theoSqrt$p.value,
cor.moon.empBC$p.value,
cor.moon.theoBC$p.value,
cor.moon.fraction100$p.value,
cor.moon.fraction95$p.value)

cor.mac.rose<-cor.test(mac, rose17)
cor.mac.empMean<-cor.test(mac, Scales17$meanEmp)
cor.mac.theoMean<-cor.test(mac, Scales17$meanTheo)
cor.mac.empMedian<-cor.test(mac, Scales17$medianEmp)
cor.mac.theoMedian<-cor.test(mac, Scales17$medianTheo)
cor.mac.empSqrt<-cor.test(mac, Scales17$meanEmpSqrt)
cor.mac.theoSqrt<-cor.test(mac, Scales17$meanTheoSqrt)
cor.mac.empBC<-cor.test(mac, Scales17$meanEmpBC)
cor.mac.theoBC<-cor.test(mac, Scales17$meanTheoBC)
cor.mac.fraction100<-cor.test(mac, Scales17$ScaleFrac)
cor.mac.fraction95<-cor.test(mac, Scales17$ScaleFrac95)

macCorr<- c(cor.mac.rose$estimate,
cor.mac.empMean$estimate,
cor.mac.theoMean$estimate,
cor.mac.empMedian$estimate,
cor.mac.theoMedian$estimate,
cor.mac.empSqrt$estimate,
cor.mac.theoSqrt$estimate,
cor.mac.empBC$estimate,
cor.mac.theoBC$estimate,
cor.mac.fraction100$estimate,
cor.mac.fraction95$estimate)

macPval<- c(cor.mac.rose$p.value,
cor.mac.empMean$p.value,
cor.mac.theoMean$p.value,
cor.mac.empMedian$p.value,
cor.mac.theoMedian$p.value,
cor.mac.empSqrt$p.value,
cor.mac.theoSqrt$p.value,
cor.mac.empBC$p.value,
cor.mac.theoBC$p.value,
cor.mac.fraction100$p.value,
cor.mac.fraction95$p.value)

cor.wim.rose<-cor.test(wim, rose)
cor.wim.empMean<-cor.test(wim, Scales$meanEmp)
cor.wim.theoMean<-cor.test(wim, Scales$meanTheo)
cor.wim.empMedian<-cor.test(wim, Scales$medianEmp)
cor.wim.theoMedian<-cor.test(wim, Scales$medianTheo)
cor.wim.empSqrt<-cor.test(wim, Scales$meanEmpSqrt)
cor.wim.theoSqrt<-cor.test(wim, Scales$meanTheoSqrt)
cor.wim.empBC<-cor.test(wim, Scales$meanEmpBC)
cor.wim.theoBC<-cor.test(wim, Scales$meanTheoBC)
cor.wim.fraction100<-cor.test(wim, Scales$ScaleFrac)
cor.wim.fraction95<-cor.test(wim, Scales$ScaleFrac95)

wimCorr<- c(cor.wim.rose$estimate,
cor.wim.empMean$estimate,
cor.wim.theoMean$estimate,
cor.wim.empMedian$estimate,
cor.wim.theoMedian$estimate,
cor.wim.empSqrt$estimate,
cor.wim.theoSqrt$estimate,
cor.wim.empBC$estimate,
cor.wim.theoBC$estimate,
cor.wim.fraction100$estimate,
cor.wim.fraction95$estimate)

wimPval<- c(cor.wim.rose$p.value,
cor.wim.empMean$p.value,
cor.wim.theoMean$p.value,
cor.wim.empMedian$p.value,
cor.wim.theoMedian$p.value,
cor.wim.empSqrt$p.value,
cor.wim.theoSqrt$p.value,
cor.wim.empBC$p.value,
cor.wim.theoBC$p.value,
cor.wim.fraction100$p.value,
cor.wim.fraction95$p.value)

labels<-c("Rose", "Empirical Mean", "Theoretical Mean", "Empirical Median", "Theoretical Median", "Empirical Sqrt Mean", "Theoretical Sqrt Mean", "Empirical Box-Cox Mean", "Theoretical Box-Cox Mean", "100_Buried", "95_Buried")

Correlations<- data.frame(labels=labels, wolfCorr= wolfCorr, wolfPval=wolfPval, kyteCorr=kyteCorr, kytePval=kytePval, radzCorr= radzCorr, radzPval=radzPval, moonCorr= moonCorr, moonPval=radzPval, macCorr= macCorr, macPval=macPval, wimCorr=wimCorr, wimPval=wimPval, fauCorr=fauCorr, fauPval=fauPval)
write.table(Correlations, "Hydrophobicity_Correlations.txt",quote=FALSE, row.names=FALSE, sep='\t')


Corr_trans<-data.frame(rbind(wolfCorr= wolfCorr, wolfPval=wolfPval, kyteCorr=kyteCorr, kytePval=kytePval, radzCorr= radzCorr, radzPval=radzPval, moonCorr= moonCorr, moonPval=radzPval, macCorr= macCorr, macPval=macPval, wimCorr=wimCorr, wimPval=wimPval, fauCorr=fauCorr, fauPval=fauPval))
colnames(Corr_trans)<-labels
write.table(Corr_trans, "Hydrophobicity_Correlations_Transposed.txt",quote=FALSE, row.names=FALSE, sep='\t')

C_labels<-c("EXP","Rose", "Theoretical Mean","Empirical Mean","100% Buried", "95% Buried")
C_table<- data.frame( Rose=Corr_trans$"Rose", Theoretical_Mean= Corr_trans$"Theoretical Mean", Empirical_Mean=Corr_trans$"Empirical Mean", buried_100= Corr_trans$"100_Buried",  buried_95=Corr_trans$"95_Buried")
row.names(C_table)<-names(Correlations)[2:length(names(Correlations))]
write.table(C_table, "Hydrophobicity_Correlations_Figure.txt",quote=FALSE, row.names=TRUE, sep='\t')
