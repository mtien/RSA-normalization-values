#!usr/bin/python
import sys, os
import math

directory_bin_all="/Users/Matthew/Documents/GitHub/RSA-normalization-values/AllBins"
directory_bin_core="/Users/Matthew/Documents/GitHub/RSA-normalization-values/CoreBins"
directory_bin_corr="/Users/Matthew/Documents/GitHub/RSA-normalization-values/Correlation/"

fileIn= open(os.path.join(directory_bin_corr, "NormalizationValuesByPercentDataCoverageAndGenerous.txt"))
fileIn.readline()

for line in fileIn:
    info=line.strip().split("\t")
    AA=info[0]
    CORE_cutoff= int(info[6])
    print AA
    f_nom=os.path.join(directory_bin_all , AA+"_max_bins_all")
    fdic=open(f_nom)
    header=fdic.readline()
    
    f_out= os.path.join(directory_bin_core , AA+"_max_bins_CORE")
    fileOut= open(f_out, 'w')
    fileOut.write(header)
    for line2 in fdic:
        data= line2.strip().split("\t")
        obs_bin_pop= int(data[4])
        if(obs_bin_pop > CORE_cutoff):
            fileOut.write(line2)
        else:
            phi= data[0]
            psi= data[1]
            fileOut.write(phi + "\t" + psi + "\t0\t0\t" + str(obs_bin_pop) + "\n")

    fdic.close()
    fileOut.close()
    
fileIn.close()
