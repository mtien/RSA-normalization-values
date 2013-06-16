#!usr/bin/python
import sys, os
import math

directory_bin_all="/Users/Matthew/Documents/GitHub/RSA-normalization-values/AllBins"
directory_bin_gen="/Users/Matthew/Documents/GitHub/RSA-normalization-values/GenerousBins"


def inbounds(pair):
    phi=pair[0]
    psi=pair[1]
    if( phi>= -180 and phi <180 and psi >=-180 and psi < 180 ):
        return True
    else:
        return False

def search_generous(position, counted):
    phi= position[0]
    psi= position[1]
    
    for i in range(0, 21, 5):
        for j in range(0, 21, 5):
            if(i+j <=20):
                neg= (phi-i, psi-j)
                pos= (phi+i, psi+j)
                if( neg not in counted and inbounds(neg) ):
                    counted.append(neg)
                if( pos not in counted and inbounds(pos) ):
                    counted.append(pos)
    return counted
                

fileIn= open("NormalizationValuesByPercentDataCoverage.txt")
fileIn.readline()

for line in fileIn:
    info=line.strip().split("\t")
    AA=info[0]
    ALLOWED_cutoff= int(info[3])
    
    f_nom=os.path.join(directory_bin_all , AA+"_max_bins_all")
    fdic=open(f_nom)
    fdic.readline()

    rama=[]

    for line2 in fdic:
        data= line2.strip().split("\t")
        phi= int(data[0])
        psi= int(data[1])
        obs_bin_pop= int(data[4])
        if(obs_bin_pop>= ALLOWED_cutoff):
            rama.append((phi,psi))
    fdic.close()
    
    print AA
    seen=[]
    for key in rama:
        seen= search_generous(key, seen)

    fdic=open(f_nom)
    fdic.readline()

    f_out= os.path.join(directory_bin_gen , AA+"_max_bins_GENEROUS")
    fileout= open(f_out, 'w')
    fileout.write("Phi\tPsi\tmax_obs_SA\tmax_theo_SA\tobs_bin_pop\n")

    print "GENEROUS"
    for line2 in fdic:
        data= line2.strip().split("\t")
        phi= int(data[0])
        psi= int(data[1])
        data_point= (phi, psi)
        if(data_point in seen):
            fileout.write(line2)
    fileout.close()
    fdic.close()
    
fileIn.close()
