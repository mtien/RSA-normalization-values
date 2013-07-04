#!usr/bin/python

import sys, os, math, string, re, gzip, urllib, shutil, Bio
import gzip, time
import cStringIO
from Bio.PDB import *
from Bio.PDB import PDBParser
import Bio.PDB.Dice
from Bio.PDB.DSSP import *
from numpy import *
from Bio.PDB.Polypeptide import *
import DSSPData as dd
import subprocess
from Geometry import *
from PeptideBuilder import *
import random 


Rotamer_Variance_0=['G', 'A', 'P']
Rotamer_Variance_1=['S', 'C']
Rotamer_Variance_Branched=['V', 'T', 'I', 'L']
Rotamer_Variance_2=['H', 'V', 'T', 'Y', 'W', 'N', 'D', 'F']
Rotamer_Variance_3=['I', 'L', 'Q', 'E', 'M', 'K', 'R']

resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
	    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
	    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
	    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

def checkTripeptide(structure):
    res_list= structure.get_residues()
    XX=res_list[1]
    flag=True
    
    O= XX['O']
    for atom in XX:
        if((atom.get_name() != 'CA') or (atom.get_name() != 'C') or (atom.get_name() != 'O') or (atom.get_name() != 'N')):
            flag= (atom-O)>1.96
        if(flag==False):
            return flag

    return flag
            

def makeTripeptide(geoGly1, geoXX, geoGly2):
    
    tripeptide= initialize_res(geoGly1)
    tripeptide=add_residue(tripeptide, geoXX)
    tripeptide=add_residue(tripeptide, geoGly2)

    return tripeptide

def getSA(structure, AA):
    ft="temp" + AA+ ".pdb"
    w= PDBIO()
    w.set_structure(structure)
    w.save(ft)
    DSSPOUTFILE="tempDSSPout"+ AA
    command = "dssp "+ ft +" > " + DSSPOUTFILE
    subprocess.call( command, shell="True")
    dd_ob=dd.DSSPData()
    dd_ob.parseDSSP(DSSPOUTFILE)
    if dd_ob.isEmpty():
        print "could not get the secondary structure for", pdb_id + chain_id
        return ( [], [], [], "", "" )

    RESIDUES= dd_ob.getAAs()
    SOLACC=dd_ob.getACC()
    return float( SOLACC[1] )

def changeRotamers(geoGly1, geoXX, geoGly2):
    acc_list=[]
    AA= geoXX.residue_name
    rotamers=[60, -60, 180]
    if(AA in Rotamer_Variance_1):
        for i in rotamers:
            geoXX.inputRotamers([i])
            structure= makeTripeptide(geoGly1, geoXX, geoGly2)
            if(checkTripeptide):
                acc_list.append( getSA(structure, AA))
            else:
                acc_list.append(0)
    elif(AA in Rotamer_Variance_2):
        if(AA in Rotamer_Variance_Branched):
            for i in rotamers:
                rotamers2= [60, -60, 180]
                rotamers2.pop(i)
                for j in rotamers2:
                    geoXX.inputRotamers([i,j])
                    structure= makeTripeptide(geoGly1, geoXX, geoGly2)
                    if(checkTripeptide):
                        acc_list.append( getSA(structure, AA))
                    else:
                        acc_list.append(0)
        else:
            for i in rotamers:
                for j in rotamers:
                    geoXX.inputRotamers([i,j])
                    structure= makeTripeptide(geoGly1, geoXX, geoGly2)
                    if(checkTripeptide):
                        acc_list.append( getSA(structure, AA))
                    else:
                        acc_list.append(0)
    else:
        if(AA =="I"):
            for i in range(0,10):
                rotamers2=[60,-60,180]
                rotamer_1=random.choice(rotamers)
                rotamers2.pop(temp1)
                rotamer_2=random.choice(rotamers2)
                rotamer_3=random.choice(rotamers)
                geoXX.inputRotamers([rotamer_1, rotamer_2, rotamer_3])
                structure= makeTripeptide(geoGly1, geoXX, geoGly2)
                if(checkTripeptide):
                    acc_list.append(getSA(structure, AA))
                else:
                    acc_list.append(0)
            
        elif(AA== "L"):
            for i in range(0,10):
                rotamers2=[60,-60,180]
                rotamer_1=random.choice(rotamers)
                rotamers2.pop(temp1)
                rotamer_2=random.choice(rotamers2)
                rotamer_3=random.choice(rotamers)
                geoXX.inputRotamers([rotamer_3, rotamer_1, rotamer_2])
                structure= makeTripeptide(geoGly1, geoXX, geoGly2)
                if(checkTripeptide):
                    acc_list.append(getSA(structure, AA))
                else:
                    acc_list.append(0)
        else:
            for i in range(0,10):
                geoXX.generateRandomRotamers()
                structure= makeTripeptide(geoGly1, geoXX, geoGly2)
                if(checkTripeptide):
                    acc_list.append(getSA(structure, AA))
                else:
                    acc_list.append(0)

    return max(acc_list)

AA=sys.argv[1]
ftemp="AnglesIteratedThroughAgain"+ AA

fileR=open(ftemp, 'w')
fileR.write("SA\tPhi\tPsi\n")    

geoGly1= geometry("G")
geoXX= geometry(resdict[AA])
geoGly2= geometry("G")
for i in range(-180 , 180):
    print str(i)
    for j in range(-180 , 180):
        acc=0
        geoXX.phi=i
        geoGly2.psi=j
        if(resdict[AA] not in Rotamer_Variance_0):
            acc=changeRotamers(geoGly1, geoXX, geoGly2)
        else:
            acc=getSA( makeTripeptide(geoGly1, geoXX, geoGly2), resdict[AA])
        fileR.write(str(acc)+ "\t"+ str(geoXX.phi) + "\t" + str(geoGly2.psi) +"\n")

fileR.close()
