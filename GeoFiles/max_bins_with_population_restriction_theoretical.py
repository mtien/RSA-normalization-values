#!usr/bin/python
import sys
import math

inp=sys.argv[1]
fname=inp + "_max_emperical_bins_pop_restriction"
filein= open(fname)
filein.readline()
key=[]
bins={}
for line in filein:
    data=line.split()
    phi=int(data[0])
    psi=int(data[1])
    sa= int(float(data[2]))
    pop= int(data[3])
    bins[(phi,psi)]=[sa, pop, 0]
    key.append((phi,psi))

filein.close()

fname="AnglesIteratedThroughAgain" + inp
filein= open(fname)
filein.readline()
flag=0
temp=0
for line in filein:
    info=line.split()
    psi= float(info[2])
    phi= float(info[1])
    sa= float(info[0])
##    if(flag==0):
##        temp=phi
##        flag=1
##    if(phi==temp):
##        phi=-phi
    b=((int(phi)/5)*5, (int(psi)/5)*5)
    if(sa > bins[b][2]):
        bins[b][2]=sa
	##print bins[b][2]

filein.close()

fout=inp+"_max_theoretical_bins_Again"
fileout=open(fout, 'w')
fileout.write("Phi\tPsi\tmaxSA\n")

for k in key:
    if( bins[k][0]==0):
	fileout.write(str(k[0])+"\t"+str(k[1])+"\t"+ str(0)+"\n")
    else:
	fileout.write(str(k[0])+"\t"+str(k[1])+"\t"+ str(bins[k][2])+"\n")


fileout.close()
