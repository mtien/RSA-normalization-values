setwd("~/GitHub/RSA-normalization-values/GeoFiles/ALLOWED/")
emp<-read.delim("Ala_ALLOWED_geo")
ala_emp<- emp[unlist(emp$SA==max(emp$SA)),]

ala_emp$psi*180.0/pi
ala_emp$phi*180.0/pi

setwd("~/GitHub/RSA-normalization-values/AnglesIteratedThroughAgain/ALLOWED/")
theo<-read.delim("AnglesIteratedThroughAgain_ALLOWED_ALA")
ala_theo<- theo[unlist(theo$SA==max(theo$SA)),]

ala_theo$Psi
ala_theo$Phi