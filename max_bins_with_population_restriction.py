#!usr/bin/python
import sys, os
import math

directory_geo="/Users/Matthew/Documents/GitHub/RSA-normalization-values/GeoFiles"
directory_theo="/Users/Matthew/Documents/GitHub/RSA-normalization-values/AnglesIteratedThroughAgain"
##directory_bin_theo="/Users/Matthew/Documents/GitHub/RSA-normalization-values/TheoreticalBins/no_pop_restriction"
##directory_bin_emp="/Users/Matthew/Documents/GitHub/RSA-normalization-values/EmpiricalBins/no_pop_restriction"

directory_bin_theo="/Users/Matthew/Documents/GitHub/RSA-normalization-values/TheoreticalBins/pop_restriction"
directory_bin_emp="/Users/Matthew/Documents/GitHub/RSA-normalization-values/EmpiricalBins/pop_restriction"


inp=sys.argv[1]
tempNom= inp[:3]
nom=tempNom.capitalize()
fname=os.path.join(directory_geo,nom+"_geo")
filein= open(fname)

filein.readline()
key=[]
bins={}
for i in range(-36,36):
    for j in range(-36,36):
        bins[(i*5,j*5)]=[0,0]
        key.append((i*5,j*5))

for line in filein:
    info=line.split()
    psi= (float(info[6])*180.0)/math.pi
    phi= (float(info[7])*180.0)/math.pi
    rsa= float(info[1])
    b=((int(phi)/5)*5, (int(psi)/5)*5)
    bins[b][1]+=1
    if(rsa > bins[b][0]):
        bins[b][0]=rsa

filein.close()
fout=os.path.join(directory_bin_emp, inp+"_max_emperical_bins_pop_restriction")
##fout=os.path.join(directory_bin_emp, inp+"_max_emperical_bins")
fileout=open(fout, 'w')
fileout.write("Phi\tPsi\tmaxSA\tpopulation\n")

for k in key:
    if( bins[k][1] >= 5):
        fileout.write(str(k[0])+"\t"+str(k[1])+"\t"+ str(bins[k][0]) +"\t" + str(bins[k][1])+ "\n")
    else:
        fileout.write(str(k[0])+"\t"+str(k[1])+"\t"+ str(0) +"\t" + str(bins[k][1])+ "\n")


fileout.close()
    
    
