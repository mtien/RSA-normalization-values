AA<-c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")
aa<-c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")

write(c("Emp/Theo","AA", "SA", "phi", "psi"),file = "max_exps_RamaVals1.txt", ncolumns = 5, append=FALSE)

for (i in 1:length(AA)){
  emp<-read.delim(paste("GeoFiles/ALLOWED/",aa[i],"_ALLOWED_geo",sep = ""))
  aa_emp<- emp[unlist(emp$SA==max(emp$SA)),]

  write(paste("Emp", AA[i], aa_emp$SA, aa_emp$phi*180.0/pi, aa_emp$psi*180.0/pi, sep = " "), file = "max_exps_RamaVals1.txt", append = TRUE)

  theo<-read.delim(paste("AnglesIteratedThroughAgain/ALLOWED/", "AnglesIteratedThroughAgain_ALLOWED_",AA[i],sep = ""))
  AA_theo<- theo[unlist(theo$SA==max(theo$SA)),]
  
  write(paste("Theo", AA[i], AA_theo$SA, AA_theo$Phi, AA_theo$Psi, sep = " "), file = "max_exps_RamaVals1.txt", append = TRUE)

}
