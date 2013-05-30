#!/usr/bin/python

import sys

AA=sys.argv[1]
Aa= AA.capitalize()

fileIn=file(Aa+"_geo")
fileIn.readline()

fileOutUnder=file(AA+"_SA_Under", "w")
fileOutOver=file(AA+"_SA_Over", "w")

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
