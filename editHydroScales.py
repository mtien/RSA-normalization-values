#!/usr/bin/python
import os
from math import log10, floor
directory_bin_corr="/Users/Matthew/Documents/GitHub/RSA-normalization-values/Correlation"

def round_sig(x, sig=3):
    ##print x
    if(x>0):
        return round(x, sig-int(floor(log10(x)))-1)
    else:
        return round(x, sig-int(floor(log10(abs(x))))-1)

fileIn= open(os.path.join(directory_bin_corr , "Hydrophobicity_Scales_updated_Figures.txt"))

code= {"ALA": "Alanine", "ASN": "Asparagine", "ARG": "Arginine", "ASP": "Asparate",
 "CYS": "Cysteine", "GLN": "Glutamine","GLU": "Glutamate","GLY": "Glycine",
 "HIS": "Histidine", "ILE":"Isoleucine", "LEU": "Leucine", "LYS": "Lysine",
 "MET": "Methionine", "PHE": "Phenylalanine", "PRO": "Proline", "SER": "Serine",
 "THR": "Threonine", "TRP": "Tryptophan", "TYR": "Tyrosine", "VAL": "Valine"}


fout=open(os.path.join(directory_bin_corr ,"HydrophobicityScalesTex.txt"),"w")
fout.write("\hline ")

fileIn.readline()
for line in fileIn:
    tempLine=line.strip()
    data=tempLine.split("\t")
    last=len(data)-1
    fout.write(code[data[0]]+ " & ")
    for i in range(1,len(data)-1):
        temp=round_sig(float(data[i]))
        fout.write(str(temp)+ " & ")
        
    fout.write(str(round_sig(float(data[last])))  + "\\\\ \n")


    
fileIn.close()
fout.close()    
