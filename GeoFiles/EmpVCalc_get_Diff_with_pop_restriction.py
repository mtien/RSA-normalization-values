#!usr/bin/python
import sys
import math

inp=sys.argv[1]
tempNom= inp[:3]
nom=tempNom.capitalize()

fin=nom+ "_geo"
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
    psi= (float(info[5])*180.0)/math.pi
    phi= (float(info[6])*180.0)/math.pi
    sa= float(info[1])
##    print rsa
##    print phi
##    print psi
    b=((int(phi)/5)*5, (int(psi)/5)*5)
##    print b
    bins[b][1]+=1
    if(sa > bins[b][0]):
        bins[b][0]=sa

filein.close()
##
##calc_bins={}
##
##for i in range(-36,36):
##    for j in range(-36,36):
##        calc_bins[(i*5,j*5)]=0
##fr="AnglesIteratedThrough" +inp
##fileR=open(fr)
##fileR.readline()
##flag=0
##temp=0
##for line in fileR:
##    info=line.split()
##    rsa=float(info[0])
##    phi=float(info[1])
##    psi=float(info[2])
##    if(flag==0):
##        temp=phi
##        flag=1
##    if(phi==temp):
##        phi=-phi
##    x=(int(phi)/5)*5
##    y=(int(psi)/5)*5
##    if( rsa> calc_bins[(x,y)]):
##        calc_bins[(x,y)]=rsa
##fileR.close()

##inp=sys.argv[1]
##fname=inp + "_max_emperical_bins_pop_restriction"
##filein= open(fname)
##filein.readline()
##key=[]
##bins={}
##for line in filein:
##    data=line.split()
##    phi=int(data[0])
##    psi=int(data[1])
##    sa= int(float(data[2]))
##    pop= int(data[3])
##    bins[(phi,psi)]=[sa, pop, 0]
##    key.append((phi,psi))
##
##filein.close()

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


fout="EmpericalVCalculated_diff_pop_nonZeroed_with_pop_restriction_"+ inp
fileout=open(fout, 'w')
fileout.write("SA_diff\tpop\n")
for k in key:
##    if(emp_bins[k][1]>=5):
    fileout.write( str(bins[k][2]-bins[k][0]) + "\t" + str(bins[k][1])+"\n")
    
fileout.close()
