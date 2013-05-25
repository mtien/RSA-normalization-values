#!/usr/bin/python
import sys, os, math, string, re, gzip, urllib, shutil, Bio
import gzip, time
import cStringIO
from Bio.PDB import PDBParser
import Bio.PDB.Dice
from Bio.PDB.DSSP import *
from numpy import *
from Bio.PDB.Polypeptide import *

import DSSPData as dd
import subprocess



resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
	    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
	    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
	    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }


residue_max_acc = {'A': 113.0, 'R': 241.0, 'N': 158.0, 'D': 151.0, \
		   'C': 140.0, 'Q': 189.0, 'E': 183.0, 'G': 85.0,  \
		   'H': 194.0, 'I': 182.0, 'L': 180.0, 'K': 211.0, \
		   'M': 204.0, 'F': 218.0, 'P': 143.0, 'S': 122.0, \
		   'T': 146.0, 'W': 259.0, 'Y': 229.0, 'V': 160.0}


def buildLocalPDBName( pdb_id ):
	"""Constructs the name of the pdb file in the local
	pdb cache from the pdb identifier. The path to the
	cache is currently hardcoded."""
	dir = "structures"
	# we create the path if it doesn't exist yet
	if not os.access( dir, os.R_OK ):
		os.makedirs( dir )
	pdb_id = string.lower( pdb_id )
	return os.path.join( dir, pdb_id + ".pdb.gz" )

def downloadStatus( count, bsize, tsize ):
	"""Reports the download status of a file"""
	if count*bsize > tsize:
		print "\t", tsize, "/", tsize, "done"
	else:
		if count/8. == count/8:
			print "\t", count*bsize, "/", tsize

def retrievePDBStructure( pdb_id ):
	pdb_id = string.lower( pdb_id )
	middle = pdb_id[1:3]
	source = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/" + middle + "/pdb" + pdb_id + ".ent.gz"
	target = buildLocalPDBName( pdb_id )
	print "   downloading structure", pdb_id
	print source
	try:
		urllib.urlretrieve( source, target, downloadStatus )
		"""Here is to determine whether the target file is empty"""
		f = open( target, "r" )
		f.seek( 0, 2 )
		if f.tell() == 0:
			print "      structure not available. trying archive of obsolete structures..."
			source = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/obsolete/pdb/" + middle + "/pdb" + pdb_id + ".ent.gz"
			print source
			try:
				urllib.urlretrieve( source, target, downloadStatus )
				f = open( target, "r")
				f.seek( 0, 2 )
				if f.tell() == 0:
					print "      doesn't work. structure not available."
					return False
				else:
					return True
			except:
				print "      doesn't work. structure not available."
				return False
		else:
			return True
	
	except:
		print "      structure not available. trying archive of all structures..."
		source = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/" + "pdb" + pdb_id + ".ent.gz"
		print source
		try:
			urllib.urlretrieve( source, target, downloadStatus )
			f = open( target, "r")
			f.seek( 0, 2 )
			if f.tell() == 0:
				print "      doesn't work. structure not available."
				return False
			else:
				return True
		except:
			print "      doesn't work. structure not available."
			return False

def parsePDBStructure( pdb_id ):
	filename = buildLocalPDBName( pdb_id )
	if not os.access( filename, os.R_OK ):
		found = retrievePDBStructure( pdb_id )
	else:
		found = True
	
	if not found:
		try:
			os.unlink(filename)
			print "   structure", pdb_id, "cannot be retrieved"
			return False
		except WindowsError:
			return False
		
	
	parser=PDBParser()
	tmpfname = "tmpfile1.pdb"
	f_in = gzip.open(filename, 'rb')
	f_out = open(tmpfname, 'wb')
	f_out.writelines(f_in)
	f_out.close()
	f_in.close()

	print pdb_id, filename, tmpfname
	
	time_count = 0
	while not os.access( tmpfname, os.R_OK ):
		time.sleep(1)
		time_count += 1
		if time_count > 10:
			return False
	try:
		structure = parser.get_structure( pdb_id, tmpfname )
	except:
		time.sleep(1)
		time_count += 1
		#os.system("rm " + filename)
		os.unlink(filename)
		print "   PDBException: 'No parent' for", pdb_id
		os.unlink( tmpfname )
		return False
	return structure

def getStruInfo( pdb_id, chain_id, DSSP ):
	"""Default setting: DSSP executable file should be in the local directory."""

	if chain_id == '-':
		chain_id = ' '

	structure = parsePDBStructure( pdb_id )
	
	if not structure:
		print "could not get structure", pdb_id + chain_id
		return  ( [], [], [], "", "" )
	print "   calculating secondary structure for", pdb_id + chain_id
	
	model = structure[0]
	chain = model[chain_id]

	Bio.PDB.Dice.extract( structure, chain_id, 0, 50000, 'tempfile.pdb' )
	command = "dssp tempfile.pdb > tempDSSPout"
	subprocess.call( command, shell="True")
	dd_ob=dd.DSSPData()
	dd_ob.parseDSSP("tempDSSPout")
        
	if dd_ob.isEmpty():
		print "could not get the secondary structure for", pdb_id + chain_id
		return ( [], [], [], "", "" )

	RESIDUES= dd_ob.getAAs()
	RESNUMS= dd_ob.getResnums()
	SS= dd_ob.getSecStruc()
	SOLACC=dd_ob.getACC()
	
	ss_list= []
	res_list=[]
	seq_list = []
	acc_list = []
	rsa_list = []
	
        distance_list=[]
        psi_list=[]
        phi_list=[]
        omega_list=[]
	bond_length=[]
	N_CA_length=[]
	CA_C_length=[]
	C_O_length=[]
	R_length=[]
	
	incompleteRes=[]
	occ_list=[]
	
	angle_list=[]
	CA_C_N2_angle=[]
	O_C_N2_angle=[]
	C_N2_CA2_angle=[]

	prev="0"
	
	X_angle=[]
	R_angle=[]

        N_prev="0"
        CA_prev="0"
        C_prev="0"
        O_prev="0"
        first="0"
        middle="0"
	for i in range(0,len(RESIDUES)):
		key = RESIDUES[i]
		if not key in resdict.values():
			print "unusual residue %s" % key
		else:
                        resNUM=int(RESNUMS[i])
                        res=chain[resNUM]
			res_list.append(resNUM)

                           
			acc = float( SOLACC[i] )
			acc_list.append( acc )
			
			rsa = float( SOLACC[i] ) / float( residue_max_acc[key] )
			rsa_list.append( rsa )
			
			seq_list.append(key)
			if key == '-':
				ss_list.append('n')
			else:
				ss_list.append(SS[i])

                        if(prev=="0"):
                                distance_list.append(0)
                                C_N2_CA2_angle.append(0)
                                phi_list.append(0)
                                omega_list.append(0)
                                
                                CA_prev=res['CA']
                                C_prev=res['C']
                                N_prev=res['N']
                                O_prev= res['O']

                                prev="1"

                        else:
                                CA=res['CA']
                                C=res['C']
                                N=res['N']
                                O=res['O']

                                CA_vector=CA.get_vector()
                                C_vector=C.get_vector()
                                N_vector= N.get_vector()
                                O_vector=O.get_vector()

                                CA_prev_vector= CA_prev.get_vector()
                                C_prev_vector= C_prev.get_vector()
                                N_prev_vector= N_prev.get_vector()
                                O_prev_vector= O_prev.get_vector()

                                CA_C_N2_angle.append(calc_angle(CA_prev_vector, C_prev_vector, N_vector))
                                O_C_N2_angle.append(calc_angle(O_prev_vector, C_prev_vector, N_vector))
                                C_N2_CA2_angle.append(calc_angle(C_prev_vector, N_vector, CA_vector))

                                distance= N-C_prev
                                distance_list.append(distance) 

                                psi= calc_dihedral(N_prev_vector, CA_prev_vector, C_prev_vector, N_vector)
                                phi= calc_dihedral(C_prev_vector, N_vector, CA_vector, C_vector)
                                omega= calc_dihedral(CA_prev_vector, C_prev_vector, N_vector, CA_vector)

                                psi_list.append(psi)
                                phi_list.append(phi)

                                omega_list.append(omega)
                                CA_prev=CA
                                C_prev=C
                                N_prev=N
                                O_prev=O

                        if(first=="0"):
                                first=res['CA'].get_vector()
                                angle_list.append(0)
                        elif(middle=="0"):
                                middle=res['CA'].get_vector()
                        else:
                                last=res['CA'].get_vector()
                                ang=calc_angle(first, middle, last)
                                angle_list.append(ang)
                                first=middle
                                middle=last
                                
			try:
                                r=[]
                                CA=res['CA']
                                C=res['C']
                                N=res['N']
                                O=res['O']

                                bond=(CA-C) + (CA-N) + (C- res['O'])

                                name=res.get_resname()

                                N_CA_length.append(N-CA)
                                CA_C_length.append(CA-C)
                                C_O_length.append(C- O)

                                ca=CA.get_vector()
                                c=C.get_vector()
                                n=N.get_vector()
                                o=O.get_vector()

                                X=[calc_dihedral(n, ca, c, o)]
                                R=[calc_angle(n, ca, c), calc_angle(ca, c, o)]
                                if(name=='ALA'):
                                        r.append(CA-res['CB'])
                                        
                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        
                                        X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
                                        incompleteRes.append(0)
                                elif(name=='ASN'):
                                        r.append(CA-res['CB'])
                                        
                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))

                                        r.append(res['CB']-res['CG'])
                                        r.append(res['CG']-res['OD1'])
                                        r.append(res['CG']-res['ND2'])

                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['OD1'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['ND2'].get_vector() ) )

                                        X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['OD1'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['ND2'].get_vector()))
                                        
                                        incompleteRes.append(0)
                                elif(name=='ASP'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['CG'])
                                        r.append(res['CG']-res['OD1'])
                                        r.append(res['CG']-res['OD2'])
                                        
                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['OD1'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['OD2'].get_vector() ) )

                                        X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['OD1'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['OD2'].get_vector()))
                                        
                                        incompleteRes.append(0)
                                elif(name=='ARG'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['CG'])
                                        r.append(res['CG']-res['CD'])
                                        r.append(res['CD']-res['NE'])
                                        r.append(res['NE']-res['CZ'])
                                        r.append(res['CZ']-res['NH1'])
                                        r.append(res['CZ']-res['NH2'])

                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector() ) )
                                        R.append(calc_angle(res['CG'].get_vector(), res['CD'].get_vector(), res['NE'].get_vector() ) )
                                        R.append(calc_angle(res['CD'].get_vector(), res['NE'].get_vector(), res['CZ'].get_vector() ) )
                                        R.append(calc_angle(res['NE'].get_vector(), res['CZ'].get_vector(), res['NH1'].get_vector() ) )
                                        R.append(calc_angle(res['NE'].get_vector(), res['CZ'].get_vector(), res['NH2'].get_vector() ) )

                                        X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector(), res['NE'].get_vector()))
                                        X.append(calc_dihedral(res['CG'].get_vector(), res['CD'].get_vector(), res['NE'].get_vector(), res['CZ'].get_vector()))
                                        X.append(calc_dihedral(res['CD'].get_vector(), res['NE'].get_vector(), res['CZ'].get_vector(), res['NH1'].get_vector()))
                                        X.append(calc_dihedral(res['CD'].get_vector(), res['NE'].get_vector(), res['CZ'].get_vector(), res['NH2'].get_vector()))

                                        incompleteRes.append(0)
                                elif(name=='CYS'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['SG'])

                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['SG'].get_vector()))

                                        X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['SG'].get_vector()))

                                        incompleteRes.append(0)
                                elif(name=='GLN'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['CG'])
                                        r.append(res['CG']-res['CD'])
                                        r.append(res['CD']-res['OE1'])
                                        r.append(res['CD']-res['NE2'])

                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector() ) )
                                        R.append(calc_angle(res['CG'].get_vector(), res['CD'].get_vector(), res['OE1'].get_vector() ) )
                                        R.append(calc_angle(res['CG'].get_vector(), res['CD'].get_vector(), res['NE2'].get_vector() ) )
                                        
                                        X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(),res['CD'].get_vector(), res['OE1'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(),res['CD'].get_vector(), res['NE2'].get_vector()))
                                        incompleteRes.append(0)
                                elif(name=='GLU'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['CG'])
                                        r.append(res['CG']-res['CD'])
                                        r.append(res['CD']-res['OE1'])
                                        r.append(res['CD']-res['OE2'])

                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector() ) )
                                        R.append(calc_angle(res['CG'].get_vector(), res['CD'].get_vector(), res['OE1'].get_vector() ) )
                                        R.append(calc_angle(res['CG'].get_vector(), res['CD'].get_vector(), res['OE2'].get_vector() ) )

					X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(),res['CD'].get_vector(), res['OE1'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(),res['CD'].get_vector(), res['OE2'].get_vector()))
                                        
                                        incompleteRes.append(0)
                                elif(name=='GLY'):
                                        r.append(0)
                                        incompleteRes.append(0)
                                elif(name=='HIS'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['CG'])
                                        r.append(res['CG']-res['ND1'])
                                        r.append(res['CG']-res['CD2'])
                                        r.append(res['ND1']-res['CE1'])
                                        r.append(res['CD2']-res['NE2'])
                                        r.append(res['NE2']-res['CE1'])

                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['ND1'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector() ) )
                                        R.append(calc_angle(res['CG'].get_vector(), res['ND1'].get_vector(), res['CE1'].get_vector() ) )
                                        R.append(calc_angle(res['CG'].get_vector(), res['CD2'].get_vector(), res['NE2'].get_vector() ) )
                                        R.append(calc_angle(res['ND1'].get_vector(), res['CE1'].get_vector(), res['NE2'].get_vector() ) )
                                        R.append(calc_angle(res['CD2'].get_vector(), res['NE2'].get_vector(), res['CE1'].get_vector() ) )
                                        
					X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
					X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['ND1'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['ND1'].get_vector(), res['CE1'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector(), res['NE2'].get_vector()))
                                        X.append(calc_dihedral(res['CG'].get_vector(), res['ND1'].get_vector(), res['CE1'].get_vector(), res['NE2'].get_vector()))

                                        incompleteRes.append(0)
                                elif(name=='ILE'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['CG1'])
                                        r.append(res['CB']-res['CG2'])
                                        r.append(res['CG1']-res['CD1'])

                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG1'].get_vector() ) )
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG2'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG1'].get_vector(), res['CD1'].get_vector() ) )

                                        X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG1'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG2'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG1'].get_vector(), res['CD1'].get_vector()))
                                        incompleteRes.append(0)
                                elif(name=='LEU'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['CG'])
                                        r.append(res['CG']-res['CD1'])
                                        r.append(res['CG']-res['CD2'])

                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD1'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector() ) )

                                        X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD1'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector()))
                                        
                                        incompleteRes.append(0)
                                elif(name=='LYS'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['CG'])
                                        r.append(res['CG']-res['CD'])
                                        r.append(res['CD']-res['CE'])
                                        r.append(res['CE']-res['NZ'])
                                                 
                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector() ) )
                                        R.append(calc_angle(res['CG'].get_vector(), res['CD'].get_vector(), res['CE'].get_vector() ) )
                                        R.append(calc_angle(res['CD'].get_vector(), res['CE'].get_vector(), res['NZ'].get_vector() ) )

					X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector(), res['CE'].get_vector()))
                                        X.append(calc_dihedral(res['CG'].get_vector(), res['CD'].get_vector(), res['CE'].get_vector(), res['NZ'].get_vector()))
                                        incompleteRes.append(0)
                                elif(name=='MET'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['CG'])
                                        r.append(res['CG']-res['SD'])
                                        r.append(res['SD']-res['CE'])

                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['SD'].get_vector() ) )
                                        R.append(calc_angle(res['CG'].get_vector(), res['SD'].get_vector(), res['CE'].get_vector() ) )

					X.append(calc_dihedral(n, c, ca, res['CB'].get_vector() ) )
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['SD'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['SD'].get_vector(), res['CE'].get_vector()))
                                        incompleteRes.append(0)
                                elif(name=='PHE'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['CG'])
                                        r.append(res['CG']-res['CD1'])
                                        r.append(res['CG']-res['CD2'])
                                        r.append(res['CD1']-res['CE1'])
                                        r.append(res['CD2']-res['CE2'])
                                        r.append(res['CE1']-res['CZ'])
                                        r.append(res['CE2']-res['CZ'])

                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD1'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector() ) )

                                        X.append( calc_dihedral(n, c, ca, res['CB'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD1'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['CD1'].get_vector(), res['CE1'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector(), res['CE2'].get_vector()))
                                        X.append(calc_dihedral(res['CG'].get_vector(), res['CD1'].get_vector(), res['CE1'].get_vector(), res['CZ'].get_vector()))

                                        incompleteRes.append(0)
                                elif(name=='PRO'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['CG'])
                                        r.append(res['CG']-res['CD'])
                                        r.append(res['CD']-N)
                                        
                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector() ) )
                                        R.append(calc_angle(res['CG'].get_vector(), res['CD'].get_vector(), N.get_vector() ) )
                                        
					X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
					X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['CD'].get_vector(), n))
					incompleteRes.append(0)
                                elif(name=='SER'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['OG'])
                                        
                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['OG'].get_vector() ) )

                                        X.append( calc_dihedral(n, c, ca, res['CB'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['OG'].get_vector() ))
                                        incompleteRes.append(0)
                                elif(name=='THR'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['OG1'])
                                        r.append(res['CB']-res['CG2'])
                                        
                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['OG1'].get_vector() ) )
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG2'].get_vector() ) )

                                        X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['OG1'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG2'].get_vector()))

                                        incompleteRes.append(0)
                                elif(name=='TRP'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['CG'])
                                        r.append(res['CG']-res['CD1'])
                                        r.append(res['CG']-res['CD2'])
                                        r.append(res['CD1']-res['NE1'])
                                        r.append(res['CD2']-res['CE2'])
                                        r.append(res['CD2']-res['CE3'])
                                        r.append(res['CE2']-res['NE1'])
                                        r.append(res['CE2']-res['CZ2'])
                                        r.append(res['CE3']-res['CZ3'])
                                        r.append(res['CZ3']-res['CH2'])
                                        r.append(res['CZ2']-res['CH2'])

                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD1'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector() ) )
                                        R.append(calc_angle(res['CG'].get_vector(), res['CD2'].get_vector(), res['CE3'].get_vector()  ) )
                                        R.append(calc_angle(res['NE1'].get_vector(), res['CE2'].get_vector(), res['CZ2'].get_vector()  ) )
                                        
					X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
					X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD1'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['CD1'].get_vector(), res['NE1'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector(), res['CE2'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector(), res['CE3'].get_vector()))
                                        X.append(calc_dihedral(res['CG'].get_vector(), res['CD2'].get_vector(), res['CE2'].get_vector(), res['CZ2'].get_vector()))
                                        X.append(calc_dihedral(res['CG'].get_vector(), res['CD2'].get_vector(), res['CE3'].get_vector(), res['CZ3'].get_vector()))
                                        X.append(calc_dihedral(res['CD2'].get_vector(), res['CE2'].get_vector(), res['CZ2'].get_vector(), res['CH2'].get_vector()))
					incompleteRes.append(0)
                                elif(name=='TYR'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['CG'])
                                        r.append(res['CG']-res['CD1'])
                                        r.append(res['CG']-res['CD2'])
                                        r.append(res['CD1']-res['CE1'])
                                        r.append(res['CD2']-res['CE2'])
                                        r.append(res['CE1']-res['CZ'])
                                        r.append(res['CE2']-res['CZ'])
                                        r.append(res['CZ']-res['OH'])

                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD1'].get_vector() ) )
                                        R.append(calc_angle(res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector() ) )
                                        R.append(calc_angle(res['CE1'].get_vector(), res['CZ'].get_vector(), res['OH'].get_vector() ) )
                                        R.append(calc_angle(res['CE2'].get_vector(), res['CZ'].get_vector(), res['OH'].get_vector() ) )

					X.append( calc_dihedral(n, c, ca, res['CB'].get_vector()))
					X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD1'].get_vector()))
                                        X.append(calc_dihedral(ca, res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['CD1'].get_vector(), res['CE1'].get_vector()))
                                        X.append(calc_dihedral(res['CB'].get_vector(), res['CG'].get_vector(), res['CD2'].get_vector(), res['CE2'].get_vector()))
                                        X.append(calc_dihedral(res['CG'].get_vector(), res['CD1'].get_vector(), res['CE1'].get_vector(), res['CZ'].get_vector()))
                                        X.append(calc_dihedral(res['CD1'].get_vector(), res['CE1'].get_vector(), res['CZ'].get_vector(), res['OH'].get_vector()))
                                        incompleteRes.append(0)
                                elif(name=='VAL'):
                                        r.append(CA-res['CB'])
                                        r.append(res['CB']-res['CG1'])
                                        r.append(res['CB']-res['CG2'])

                                        R.append(calc_angle(n, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(c, ca, res['CB'].get_vector()))
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG1'].get_vector() ) )
                                        R.append(calc_angle(ca, res['CB'].get_vector(), res['CG2'].get_vector() ) )
                                        
					X.append(calc_dihedral(n, c, ca, res['CB'].get_vector()))
					X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG1'].get_vector()))
                                        X.append(calc_dihedral(n, ca, res['CB'].get_vector(), res['CG2'].get_vector()))
					incompleteRes.append(0)
                                else:
                                        r.append([0])
                                        R.append([0])
                                        X.append([0])
                                        incompleteRes.append(1)
                                if(incompleteRes[-1]==0):
                                        num_atom=0.0
                                        occ_atom=0.0
                                        for atom in res:
                                                num_atom+=1
                                                occ_atom+=atom.get_occupancy()
                                        occ_list.append(occ_atom/num_atom)
                                else:
                                        occ_list.append(0)
                                for num in r:
                                        bond+=num
                                bond_length.append(bond)
                                R_length.append(r)
                                R_angle.append(R)
                                X_angle.append(X)
                        except KeyError:
                                incompleteRes.append(1)
                                R_length.append([0])
                                X_angle.append([0])
                                bond_length.append(0)
                                occ_list.append(0)
                                R_angle.append([0])
        angle_list.append(0)
        psi_list.append(0)


	os.unlink('tempfile.pdb')
	return seq_list, res_list, acc_list, occ_list, rsa_list, distance_list, angle_list, phi_list, psi_list, omega_list, ss_list, N_CA_length, CA_C_length, C_O_length, bond_length, R_length, incompleteRes, CA_C_N2_angle, O_C_N2_angle, C_N2_CA2_angle, R_angle, X_angle
	
	

def get_structure_header(pdb_id):
	structure = parsePDBStructure( pdb_id )
	if not structure:
		print "could not get structure"
		return False
	resolution = structure.header['resolution']
	if not resolution:
		resolution = 'NA'
	method = structure.header['structure_method']
	model_count = 0
	print "Resolution is ", resolution
	return (resolution, method, len(structure))
