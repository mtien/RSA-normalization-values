#!usr/bin/python
import sys
import math

inp=sys.argv[1]
tempNom= inp[:3]
nom=tempNom.capitalize()
fname=nom+"_geo"
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
##    print rsa
##    print phi
##    print psi
    b=((int(phi)/5)*5, (int(psi)/5)*5)
##    print b
    bins[b][1]+=1
    if(rsa > bins[b][0]):
        bins[b][0]=rsa

filein.close()
fout=inp+"_max_emperical_bins_pop_restriction"
fileout=open(fout, 'w')
fileout.write("Phi\tPsi\tmaxSA\tpopulation\n")

for k in key:
    if( bins[k][1] >= 1):
        fileout.write(str(k[0])+"\t"+str(k[1])+"\t"+ str(bins[k][0]) +"\t" + str(bins[k][1])+ "\n")
    else:
        fileout.write(str(k[0])+"\t"+str(k[1])+"\t"+ str(0) +"\t" + str(bins[k][1])+ "\n")


fileout.close()
    
    
