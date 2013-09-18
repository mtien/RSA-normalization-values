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
import Geometry
import PeptideBuilder
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

def make_pdb_file(struct, file_nom):
    outfile = PDBIO()
    outfile.set_structure(struct)
    outfile.save( file_nom)
    return file_nom 

def build_exact_model(structure, chain_index):
    model=structure[0]
    chain=model[chain_index]
    model_structure_geo=[]
    prev="0"
    N_prev="0"
    CA_prev="0"
    C_prev="0"
    O_prev="0"
    prev_res=""
    rad=180.0/math.pi
    for res in chain:
        name=res.get_resname()
        if(name !="HOH"):
            geo=Geometry.geometry(resdict[name])
            if(prev=="0"):
                N_prev= res['N']
                CA_prev= res['CA']
                C_prev= res['C']
                O_prev= res['O']
                prev="1"
            else:   
                n1=N_prev.get_vector()
                ca1=CA_prev.get_vector()
                c1=C_prev.get_vector()
                o1=O_prev.get_vector()
                
                O_curr=res['O']
                C_curr=res['C']
                N_curr=res['N']
                CA_curr=res['CA']
                                
                o=O_curr.get_vector()
                c=C_curr.get_vector()
                n=N_curr.get_vector()
                ca=CA_curr.get_vector()

                geo.CA_C_N_angle=calc_angle(ca1, c1, n)*rad
                geo.C_N_CA_angle=calc_angle(c1, n, ca)*rad

                geo.peptide= N_curr-C_prev

                psi= calc_dihedral(n1, ca1, c1, n) ##goes to current res
                omega= calc_dihedral(ca1, c1, n, ca) ##goes to current res
                phi= calc_dihedral(c1, n, ca, c) ##goes to current res

                geo.psi_im1=psi*rad
                geo.omega=omega*rad
                geo.phi=phi*rad

                geo.CA_N_length= CA_curr - N_curr
                geo.CA_C_length= CA_curr - C_curr 
                geo.C_O_length= C_curr - O_curr

                geo.N_CA_C_angle= calc_angle(n, ca, c)*rad
                geo.CA_C_O_angle= calc_angle(ca, c, o)*rad

                geo.N_CA_C_O= calc_dihedral(n, ca, c, o)*rad

                N_prev= res['N']
                CA_prev= res['CA']
                C_prev= res['C']
                O_prev= res['O']
                
                if(name=='ALA'):
                    geo.CA_CB_length=CA_curr-res['CB']
                    
                    geo.C_CA_CB_angle= calc_angle(c, ca, res['CB'].get_vector())*rad
                    
                    geo.N_C_CA_CB_diangle= calc_dihedral(n, c, ca, res['CB'].get_vector())*rad

                elif(name=='GLU'):
                    geo.CA_CB_length=CA_curr-res['CB']
                    geo.C_CA_CB_angle=calc_angle(c, ca, res['CB'].get_vector())*rad
                    geo.N_C_CA_CB_diangle=calc_dihedral(n, c, ca, res['CB'].get_vector() )*rad

                    geo.CB_CG_length=res['CB']-res['CG']
                    geo.CA_CB_CG_angle=calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector())*rad
                    geo.N_CA_CB_CG_diangle=calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector())*rad

                    geo.CG_CD_length=res['CG']-res['CD']
                    geo.CB_CG_CD_angle=calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector() )*rad
                    geo.CA_CB_CG_CD_diangle=calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector())*rad

                    geo.CD_OE1_length=res['CD']-res['OE1']
                    geo.CG_CD_OE1_angle=calc_angle(res['CG'].get_vector(), res['CD'].get_vector(), res['OE1'].get_vector() )*rad
                    geo.CB_CG_CD_OE1_diangle= calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(),res['CD'].get_vector(), res['OE1'].get_vector())*rad

                    geo.CD_OE2_length=res['CD']-res['OE2']
                    geo.CG_CD_OE2_angle=calc_angle(res['CG'].get_vector(), res['CD'].get_vector(), res['OE2'].get_vector() )*rad
                    geo.CB_CG_CD_OE2_diangle= calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(),res['CD'].get_vector(), res['OE2'].get_vector())*rad
                    
                elif(name=='HIS'):
                    geo.CA_CB_length=CA_curr-res['CB']
                    geo.CB_CG_length=(res['CB']-res['CG'])
                    geo.CG_ND1_length=(res['CG']-res['ND1'])
                    geo.CG_CD2_length=(res['CG']-res['CD2'])
                    geo.ND1_CE1_length=(res['ND1']-res['CE1'])
                    geo.CD2_NE2_length=(res['CD2']-res['NE2'])

                    geo.C_CA_CB_angle=(calc_angle(c, ca, res['CB'].get_vector()))*rad
                    geo.CA_CB_CG_angle=(calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector() ) )*rad
                    geo.CB_CG_ND1_angle= (calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['ND1'].get_vector() ) )*rad
                    geo.CB_CG_CD2_angle=(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector() ) )*rad
                    geo.CG_ND1_CE1_angle=(calc_angle(res['CG'].get_vector(), res['ND1'].get_vector(), res['CE1'].get_vector() ) )*rad
                    geo.CG_CD2_NE2_angle=(calc_angle(res['CG'].get_vector(), res['CD2'].get_vector(), res['NE2'].get_vector() ) )*rad
                    
                    geo.N_C_CA_CB_diangle=(calc_dihedral(n, c, ca, res['CB'].get_vector()))*rad
                    geo.N_CA_CB_CG_diangle=(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector()))*rad
                    geo.CA_CB_CG_ND1_diangle=(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['ND1'].get_vector()))*rad
                    geo.CA_CB_CG_CD2_diangle=(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector()))*rad
                    geo.CB_CG_ND1_CE1_diangle=(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['ND1'].get_vector(), res['CE1'].get_vector()))*rad
                    geo.CB_CG_CD2_NE2_diangle=(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector(), res['NE2'].get_vector()))*rad
            
            model_structure_geo.append(geo)
    model_structure=PeptideBuilder.make_structure_from_geos(model_structure_geo)
    return model_structure
            
def makeTripeptide(geoGly1, geoXX, geoGly2):
    
    tripeptide= initialize_res(geoGly1)
    tripeptide=add_residue(tripeptide, geoXX)
    tripeptide=add_residue(tripeptide, geoGly2)

    return tripeptide

def getSA(structure, AA):
    ft="temp" + AA +".pdb"
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

if __name__=="__main__":
    pdb_file="3ssbC_GHE_tripeptide.pdb"
    parser=PDBParser()
    emp_structure=parser.get_structure('sample', pdb_file)
    emp_model=emp_structure[0]
    emp_chain=emp_model['C']
    print getSA(emp_structure, "EMPIRICAL")

    theo_structure=build_exact_model(emp_structure, "C")
    print getSA(theo_structure, "THEORETICAL")
    

    
