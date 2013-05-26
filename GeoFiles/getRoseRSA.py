#!usr/bin/python
import sys
import math

AA=["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Pro", "Phe", "Ser", "Thr", "Trp", "Tyr", "Val"]
residue_max_acc = {'A': 118.1, 'R': 256.0, 'N': 165.5, 'D': 158.7, \
		   'C': 146.1, 'Q': 193.2, 'E': 186.2, 'G': 88.1,  \
		   'H': 202.5, 'I': 181.0, 'L': 193.1, 'K': 225.8, \
		   'M': 203.4, 'F': 222.8, 'P': 146.8, 'S': 129.8, \
		   'T': 152.5, 'W': 266.3, 'Y': 236.8, 'V': 164.5}
resdict = { 'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe': 'F', \
	    'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Lys': 'K', 'Leu': 'L', \
	    'Met': 'M', 'Asn': 'N', 'Pro': 'P', 'Gln': 'Q', 'Arg': 'R', \
	    'Ser': 'S', 'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y' }

for aa in AA:
    f=open(aa+"_geo")
    first= f.readline()
    lab="AA\tRSA\n"
    fout=open(aa+ "_Rose_RSA", "w")
    fout.write(lab)
    for line in f:
        stuff=line.split()
        num=float(stuff[1])/residue_max_acc[resdict[aa]]
        fout.write(stuff[0]+ "\t" + str(num) +"\n")
    fout.close()
    f.close()
