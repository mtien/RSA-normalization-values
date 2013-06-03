#!/usr/bin/python

import sys, os
directory_geo="/Users/Matthew/Documents/GitHub/RSA-normalization-values/GeoFiles"
directory_over="/Users/Matthew/Documents/GitHub/RSA-normalization-values/SA_Over_Under/Over"
directory_under="/Users/Matthew/Documents/GitHub/RSA-normalization-values/SA_Over_Under/Under"

AA=sys.argv[1]
Aa= AA.capitalize()

fileIn=file(os.path.join(directory_geo,Aa+"_geo"))
fileIn.readline()

fileOutUnder=file(os.path.join(directory_under,AA+"_SA_Under"), "w")
fileOutOver=file(os.path.join(directory_over,AA+"_SA_Over"), "w")

fileOutUnder.write("RSA\tPsi\tPhi\n")
fileOutOver.write("RSA\tPsi\tPhi\n")

for line in fileIn:
    data=line.split()
    if( float(data[2])<=1):
        fileOutUnder.write(data[2] + "\t" +data[6] + "\t" + data[7] + "\n")
    else:
        fileOutOver.write(data[2] + "\t" +data[6] + "\t" + data[7] + "\n")

fileIn.close()
fileOutUnder.close()
fileOutOver.close()
