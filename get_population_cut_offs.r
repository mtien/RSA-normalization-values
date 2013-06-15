AA<-c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")
Aa<-c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Pro", "Phe", "Ser", "Thr", "Trp", "Tyr", "Val")

pop_80<-c()
max_theo_80<-c()
max_emp_80<-c()

pop_97<-c()
max_theo_97<-c()
max_emp_97<-c()

max_theo_100<-c()
max_emp_100<-c()

for(a in 1:length(Aa))
{
  setwd("~/GitHub/RSA-normalization-values/AllBins")
  f_all<-paste(Aa[a], "_max_bins_all", sep='')
  d<-read.delim(f_all)
  max_theo_100<-c(max_theo_100, max(d$max_theo_SA))
  max_emp_100<-c(max_emp_100, max(d$max_obs_SA))
  
  counter<-1
  dtemp<- d[unlist(d$obs_bin_pop> counter),]
  while( sum(dtemp$obs_bin_pop)/sum(d$obs_bin_pop) > .97 ){
    counter<- counter+1
    dtemp<- d[unlist(d$obs_bin_pop> counter),]
  }
  
  counter<- counter-1
  pop_97<-c(pop_97, counter)
  
  dtemp<- d[unlist(d$obs_bin_pop> counter),]
  max_theo_97<-c(max_theo_97, max(dtemp$max_theo_SA))
  max_emp_97<-c(max_emp_97, max(dtemp$max_obs_SA))
  
  while( sum(dtemp$obs_bin_pop)/sum(d$obs_bin_pop) > .80 ){
    counter<- counter+1
    dtemp<- d[unlist(d$obs_bin_pop> counter),]
  }
  
  counter<- counter-1
  pop_80<-c(pop_80, counter)
  
  dtemp<- d[unlist(d$obs_bin_pop> counter),]
  max_theo_80<-c(max_theo_80, max(dtemp$max_theo_SA))
  max_emp_80<-c(max_emp_80, max(dtemp$max_obs_SA))
}

data<-data.frame(AminoAcid=Aa, 
                   max_theo_100=max_theo_100,
                   max_emp_100= max_emp_100,
                   pop_97= pop_97,
                   max_theo_97=max_theo_97,
                   max_emp_97=max_emp_97,
                   pop_80=pop_80,
                   max_theo_80=max_theo_80,
                   max_emp_80=max_emp_80)
setwd("~/GitHub/RSA-normalization-values/")
write.table(data, file="NormalizationValuesByPercentDataCoverage.txt", quote=FALSE, row.names=FALSE, sep='\t')
