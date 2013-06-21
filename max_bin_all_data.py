#!usr/bin/python
import sys, os
import math

directory_geo="/Users/Matthew/Documents/GitHub/RSA-normalization-values/GeoFiles"
directory_theo="/Users/Matthew/Documents/GitHub/RSA-normalization-values/AnglesIteratedThroughAgain"

directory_bin_all="/Users/Matthew/Documents/GitHub/RSA-normalization-values/AllBins"

key=[]
bins={}
for i in range(-36,36):
    for j in range(-36,36):
        ##max observed, max theoretical, population
        bins[(i*5,j*5)]=[0,0,0]
        key.append((i*5,j*5))

inp=sys.argv[1]
tempNom= inp[:3]
nom=tempNom.capitalize()
fname=os.path.join(directory_geo,nom+"_geo")
filein= open(fname)
filein.readline()

for line in filein:
    info=line.split()
    psi= (float(info[6])*180.0)/math.pi
    phi= (float(info[7])*180.0)/math.pi
    sa= float(info[1])
    b=((int(phi)/5)*5, (int(psi)/5)*5)
    bins[b][2]+=1
    if(sa > bins[b][0]):
        bins[b][0]=sa

filein.close()

cap= tempNom.upper()

fname=os.path.join(directory_theo ,"AnglesIteratedThroughAgain" + cap)
filein= open(fname)
filein.readline()

for line in filein:
    info=line.split()
    psi= float(info[2])
    phi= float(info[1])
    sa= float(info[0])
    b=( (int(phi)/5)*5, (int(psi)/5)*5 )
    if(sa > bins[b][1]):
        bins[b][1]=sa

filein.close()

fout=os.path.join(directory_bin_all , inp+"_max_bins_all")
fileout=open(fout, 'w')
fileout.write("Phi\tPsi\tmax_obs_SA\tmax_theo_SA\tobs_bin_pop\n")

for k in key:
    fileout.write(str(k[0])+"\t"+str(k[1]) +"\t"+ str(bins[k][0]) +"\t"+ str(bins[k][1]) +"\t"+ str(bins[k][2]) +"\n")


fileout.close()

