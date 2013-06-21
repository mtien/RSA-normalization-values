#!usr/bin/python
import sys, os
import math

directory_bin_all="/Users/Matthew/Documents/GitHub/RSA-normalization-values/AllBins"
directory_geo="/Users/Matthew/Documents/GitHub/RSA-normalization-values/GeoFiles"
directory_allowed_geo="/Users/Matthew/Documents/GitHub/RSA-normalization-values/GeoFiles/ALLOWED"
directory_bin_corr="/Users/Matthew/Documents/GitHub/RSA-normalization-values/Correlation"
directory_theo="/Users/Matthew/Documents/GitHub/RSA-normalization-values/AnglesIteratedThroughAgain"
directory_allowed_theo="/Users/Matthew/Documents/GitHub/RSA-normalization-values/AnglesIteratedThroughAgain/ALLOWED"

fileIn= open(os.path.join(directory_bin_corr, "NormalizationValuesByPercentDataCoverageAndGenerous.txt"))
fileIn.readline()

for line in fileIn:
    info=line.strip().split("\t")
    AA=info[0]
    ALLOWED_cutoff= int(info[3])
    print AA
    f_nom=os.path.join(directory_bin_all , AA+"_max_bins_all")
    fdic=open(f_nom)
    fdic.readline()

    rama=[]

    for line2 in fdic:
        data= line2.strip().split("\t")
        phi= int(data[0])
        psi= int(data[1])
        obs_bin_pop= int(data[4])
        if(obs_bin_pop> ALLOWED_cutoff):
            rama.append((phi,psi))
    fdic.close()
    
    f_geo=os.path.join(directory_geo, AA+"_geo")
    fileGeo=open(f_geo)
    header=fileGeo.readline()

    f_out= os.path.join(directory_allowed_geo , AA+"_ALLOWED_geo")
    fileOut= open(f_out, 'w')
    fileOut.write(header)

    print "make"
    for line2 in fileGeo:
        data= line2.strip().split()
        psi= (int((float(data[6])*180.0)/math.pi)/5)*5
        phi= (int((float(data[7])*180.0)/math.pi)/5)*5
        data_point= (phi, psi)
        if(data_point in rama):
            fileOut.write(line2)
    fileOut.close()
    fileGeo.close()



    cap=AA.upper()
    f_theo=os.path.join(directory_theo, "AnglesIteratedThroughAgain" + cap)
    fileTheo=open(f_theo)
    header=fileTheo.readline()
    
    f_out= os.path.join(directory_allowed_theo , "AnglesIteratedThroughAgain_ALLOWED_" + cap)
    fileOut= open(f_out, 'w')
    fileOut.write(header)

    print "theo"
    for line2 in fileTheo:
        data= line2.strip().split()
        psi= float(data[2])
        phi= float(data[1])
        b=((int(phi)/5)*5, (int(psi)/5)*5)
        if(b in rama):
            fileOut.write(line2)
    fileOut.close()
    fileTheo.close()
    
fileIn.close()
