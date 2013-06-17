get19<-function(tempARR){
  temp<-c(tempARR[1:14], tempARR[16:20])
  return(temp)
}

get17<-function(tempARR){
  temp<-c(tempARR[1:7], tempARR[10:14], tempARR[16:20])
  return(temp)
}

library("TeachingDemos")
library("MASS")

AA<-c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")
aa<-c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
setwd("~/GitHub/RSA-normalization-values/Correlation/")
NormValues<- read.table("NormalizationValuesByPercentDataCoverageAndGenerous.txt", row.names=1, header=TRUE)

Theoretical<-NormValues$max_theo_97
Empirical<-NormValues$max_emp_97

meanTheo<-c()
##SDTheo<-c()
meanTheoSqrt<-c()
medianTheo<-c()
meanTheoBC<-c()
TheoRSA<-c()

meanEmp<-c()
##SDEmp<-c()
meanEmpSqrt<-c()
medianEmp<-c()
meanEmpBC<-c()
EmpRSA<-c()

ScaleFrac<-c()
ScaleFrac95<-c()

for(a in 1:length(AA))
{
setwd("~/GitHub/RSA-normalization-values/GeoFiles/")
SA_AA<- paste(aa[a], "_geo", sep='')
f<- read.delim(SA_AA)

TheoRSA<-c(TheoRSA, (f$SA+1)/(Theoretical[a]+1))

EmpRSA<-c(EmpRSA, (f$SA+1)/(Empirical[a]+1))
}

######################################################################

start_max=1
test<- boxcox(TheoRSA~1, lambda=seq(0,start_max, length=100))
lambdaT<-test$x[which.max(test$y)]

while(lambdaT > 0.95*max(test$x))
{
start_max= start_max+1
this_seq= seq(0, start_max, length=100)
test<- boxcox(TheoRSA~1, this_seq)
lambdaT<-test$x[which.max(test$y)]
}

start_max=1
test<- boxcox(EmpRSA~1, lambda=seq(0,start_max, length=100))
lambdaE<-test$x[which.max(test$y)]

while(lambdaE > 0.95*max(test$x))
{
start_max= start_max+1
this_seq= seq(0, start_max, length=100)
test<- boxcox(EmpRSA~1, this_seq)
lambdaE<-test$x[which.max(test$y)]
}

#######################################################################

for(a in 1:length(AA))
{
setwd("~/GitHub/RSA-normalization-values/GeoFiles/")
SA_AA<- paste(aa[a], "_geo", sep='')
f<- read.delim(SA_AA)

denom<- length(f$SA)
variableTemp<-f$SA
number<- length(variableTemp[variableTemp<1])
ScaleFrac<- append( ScaleFrac, number/denom)

TRSA<-(f$SA)/Theoretical[a]
number<- length(TRSA[TRSA<.05])
ScaleFrac95<- c(ScaleFrac95, number/denom)


TRSA<-(f$SA)/Theoretical[a]

meanTheo<-append(meanTheo, mean(1-TRSA) )
##SDTheo<-c(SDTheo, sd(TheoRSA))
meanTheoSqrt<-append(meanTheoSqrt, mean(1-sqrt(TRSA)))
medianTheo<- append(medianTheo, median(1-TRSA))

TRSA_transformed<-(f$SA+1)/(Theoretical[a]+1)
meanTheoBC<-c(meanTheoBC, mean( bct(TRSA_transformed, lambdaT) ) )

ERSA<-(f$SA)/Empirical[a]

meanEmp<- append(meanEmp, mean(1-ERSA) )
##SDEmp<- c(SDEmp, sd(EmpRSA))
meanEmpSqrt<-append(meanEmpSqrt, mean(1-sqrt(ERSA)))
medianEmp<-append(medianEmp, median(1-ERSA))

ERSA_transformed<-(f$SA+1)/(Empirical[a]+1)
meanEmpBC<-c(meanEmpBC, mean( bct(ERSA_transformed, lambdaE) ))
}

##out<-data.frame(names=AA, Theoretical_Mean=meanTheo, Theoretical_Standard_Deviation= SDTheo, Empirical_Mean= meanEmp, Empirical_Standard_Deviation=SDEmp)
##out<-data.frame(names=AA, Empirical_Mean= meanEmp, Empirical_Log_mean=meanEmpLog, Theoretical_Mean=meanTheo, Theoretical_LOG_Mean= meanTheoLog)
##mean_normal_SQRT<-data.frame(names=AA, Empirical_Mean= meanEmp, Empirical_Sqrt_mean=meanEmpLog, Theoretical_Mean=meanTheo, Theoretical_Sqrt_Mean= meanTheoLog)

Scales<-data.frame(AminoAcid=AA, meanEmp=meanEmp, medianEmp=medianEmp,
meanEmpSqrt=meanEmpSqrt, meanEmpBC= meanEmpBC, meanTheo=meanTheo, 
medianTheo=medianTheo, meanTheoSqrt=meanTheoSqrt, meanTheoBC= meanTheoBC,
ScaleFrac=ScaleFrac, ScaleFrac95= ScaleFrac95)

MeanScales<-data.frame(AminoAcid=AA, meanEmp=meanEmp, medianEmp=medianEmp,
meanEmpSqrt=meanEmpSqrt, meanEmpBC= meanEmpBC, meanTheo=meanTheo, 
medianTheo=medianTheo, meanTheoSqrt=meanTheoSqrt, meanTheoBC= meanTheoBC)

BuriedScales<-data.frame(AminoAcid=AA,
ScaleFrac=ScaleFrac, ScaleFrac95= ScaleFrac95)

setwd("~/GitHub/RSA-normalization-values/Correlation/")
##write.table(MeanScales, file="MeanHydrophobicityScales.txt", quote=FALSE, row.names=FALSE, sep='\t')
##write.table(BuriedScales, file="BuriedHydrophobicityScales.txt", quote=FALSE, row.names=FALSE, sep='\t')

write.table(Scales, file="Hydrophobicity_Scales_updated.txt",quote=FALSE, row.names=FALSE, sep='\t')

Scales19<-data.frame(AminoAcid=get19(AA), meanEmp=get19(meanEmp), 
                     medianEmp=get19(medianEmp),meanEmpSqrt=get19(meanEmpSqrt), 
                     meanEmpBC=get19(meanEmpBC), meanTheo=get19(meanTheo),
                     medianTheo=get19(medianTheo), meanTheoSqrt=get19(meanTheoSqrt), 
                     meanTheoBC= get19(meanTheoBC),
                   ScaleFrac=get19(ScaleFrac), ScaleFrac95= get19(ScaleFrac95))

Scales17<-data.frame(AminoAcid=get17(AA), meanEmp=get17(meanEmp), 
                     medianEmp=get17(medianEmp),meanEmpSqrt=get17(meanEmpSqrt), 
                     meanEmpBC=get17(meanEmpBC), meanTheo=get17(meanTheo),
                     medianTheo=get17(medianTheo), meanTheoSqrt=get17(meanTheoSqrt), 
                     meanTheoBC= get17(meanTheoBC),
                     ScaleFrac=get17(ScaleFrac), ScaleFrac95= get17(ScaleFrac95))
write.table(Scales19, file="Hydrophobicity_Scales_updated19.txt",quote=FALSE, row.names=FALSE, sep='\t')
write.table(Scales17, file="Hydrophobicity_Scales_updated17.txt",quote=FALSE, row.names=FALSE, sep='\t')
