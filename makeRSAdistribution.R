code<-"Ala"

setwd("~/GitHub/RSA-normalization-values/GeoFiles/ALLOWED/")
f_all<-paste(code, "_ALLOWED_geo", sep='')
d<-read.delim(f_all)
hist(d$SA/128.0)