
AA<-c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")
Aa<-c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")

maxTheo<-c()
maxEmp<-c()

maxTheo_pop<-c()
maxEmp_pop<-c()

for(a in 1:length(Aa))
{
setwd("~/GitHub/RSA-normalization-values/TheoreticalBins/no_pop_restriction/")
theo<-paste(Aa[a], "_max_theoretical_bins_Again", sep='')
d<-read.delim(theo)
maxTheo<-c(maxTheo,max(d$maxSA))

setwd("~/GitHub/RSA-normalization-values/EmpiricalBins/no_pop_restriction/")
emp<-paste(Aa[a],"_max_emperical_bins",sep='')
e<-read.delim(emp)
maxEmp<- c(maxEmp, max(e$maxSA))

}

for(a in 1:length(Aa))
{
  setwd("~/GitHub/RSA-normalization-values/TheoreticalBins/pop_restriction/")
  theo<-paste(Aa[a], "_max_theoretical_bins_Again_pop_restriction", sep='')
  d<-read.delim(theo)
  maxTheo_pop<-c(maxTheo_pop,max(d$maxSA))
  
  setwd("~/GitHub/RSA-normalization-values/EmpiricalBins/pop_restriction/")
  emp<-paste(Aa[a],"_max_emperical_bins_pop_restriction",sep='')
  e<-read.delim(emp)
  maxEmp_pop<- c(maxEmp_pop, max(e$maxSA))
  
}

out<-data.frame(names=Aa, Theoretical=maxTheo, Empirical= maxEmp)
out_pop<-data.frame(names=Aa, Theoretical=maxTheo_pop, Empirical= maxEmp_pop)
setwd("~/GitHub/RSA-normalization-values/")
write.csv(out, file="NormalizationValues_no_pop.txt")
write.csv(out_pop, file="NormalizationValues_pop.txt")