#!usr/bin/python
import sys, os
import math
directory_geo="/Users/Matthew/Documents/GitHub/RSA-normalization-values/GeoFiles/"

inp=sys.argv[1]
tempNom= inp[:3]
nom=tempNom.capitalize()
cap= inp.upper()

fin=os.path.join(directory_geo,nom+ "_geo")
filein= open(fin)

filein.readline()
key=[]
bins={}
for i in range(-36,36):
    for j in range(-36,36):
        bins[(i*5,j*5)]=[0,0,0]
        key.append((i*5,j*5))

for line in filein:
    info=line.split()
    psi= (float(info[6])*180.0)/math.pi
    phi= (float(info[7])*180.0)/math.pi
    sa= float(info[1])
    b=((int(phi)/5)*5, (int(psi)/5)*5)
    bins[b][1]+=1
    if(sa > bins[b][0]):
        bins[b][0]=sa

filein.close()

directory_theo="/Users/Matthew/Documents/GitHub/RSA-normalization-values/AnglesIteratedThroughAgain"
fname=os.path.join(directory_theo ,"AnglesIteratedThroughAgain" + cap)
filein= open(fname)
filein.readline()
flag=0
temp=0
for line in filein:
    info=line.split()
    psi= float(info[2])
    phi= float(info[1])
    sa= float(info[0])
    b=((int(phi)/5)*5, (int(psi)/5)*5)
    if(sa > bins[b][2]):
        bins[b][2]=sa

filein.close()

directory_evt="/Users/Matthew/Documents/GitHub/RSA-normalization-values/EmpiricalVTheoretical"
fout=os.path.join(directory_evt ,"EmpericalVCalculated_diff_pop_nonZeroed_"+ nom)


fileout=open(fout, 'w')
fileout.write("SA_diff\tpop\n")
for k in key:
    ##if(emp_bins[k][1]>=5):
    fileout.write( str(bins[k][2]-bins[k][0]) + "\t" + str(bins[k][1])+"\n")
    
fileout.close()
