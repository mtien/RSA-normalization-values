#!/usr/bin/python

import get_PDB
import sys


dssp_cmd = "/work/mtien/relative_solvent_acc/theoretical_max/Test/newStructure/GeoFolder/dssp_test/dssp"

count=0
c=0
fileP= open("cullpdb_pc20_res1.8_R0.25_d130517_chains3211")
fileP.readline()

##fileP= open("tester")
##fileP.readline()

fileError= open("Error_report_parse_alignment.txt", "wt")

fileALA= open("Ala_geo", 'w')
fileALA.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tN-Ca-Co-O\tN-Co-Ca-Cb\n")
fileASN=open("Asn_geo",'w')
fileASN.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg\tCg-Od1\tCg-Nd2\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg\tCb-Cg-Od1\tCb-Cg-Nd2\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Cg\tCa-Cb-Cg-Od1\tCa-Cb-Cg-Nd2\n")
fileASP=open("Asp_geo", 'w')
fileASP.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg\tCg-Od1\tCg-Od2\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg\tCb-Cg-Od1\tCb-Cg-Od2\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Cg\tCa-Cb-Cg-Od1\tCa-Cb-Cg-Od2\n")
fileARG= open("Arg_geo", 'w')
fileARG.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg\tCg-Cd\tCd-Ne\tNe-Cz\tCz-Nh1\tCz-Nh2\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg\tCb-Cg-Cd\tCg-Cd-Ne\tCd-Ne-Cz\tNe-Cz-Nh1\tNe-Cz-Nh2\tN-Ca-Co-O\tN-C-Ca-Cb\tN-Ca-Cb-Cg\tCa-Cb-Cg-Cd\tCb-Cg-Cd-Ne\tCg-Cd-Ne-Cz\tCd-Ne-Bz-Nh1\tCd-Ne-Cz-Nh2\n")
fileCYS=open("Cys_geo", 'w')
fileCYS.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Sg\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Sg\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Sg\n")
fileGLN=open("Gln_geo", 'w')
fileGLN.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg\tCg-Cd\tCd-Oe1\tCd-Ne2\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg\tCb-Cg-Cd\tCg-Cd-Oe1\tCg-Cd-Ne2\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Cg\tCa-Cb-Cg-Cd\tCb-Cg-Cd-Oe1\tCb-Cg-Cd-Ne2\n")
fileGLU=open("Glu_geo", 'w')
fileGLU.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg\tCg-Cd\tCd-Oe1\tCd-Oe2\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg\tCb-Cg-Cd\tCg-Cd-Oe1\tCg-Cd-Oe2\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Cg\tCa-Cb-Cg-Cd\tCb-Cg-Cd-Oe1\tCb-Cg-Cd-Oe2\n")
fileGLY=open("Gly_geo", 'w')
fileGLY.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tR-length\tN-Ca-Co\tCa-Co-O\tN-Ca-Co-O\n")
fileHIS=open("His_geo", 'w')
fileHIS.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg\tCg-Nd1\tCg-Cd2\tNd1-Ce1\tCd2-Ne2\tNe2-Ce1\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg\tCb-Cg-Nd1\tCb-Cg-Cd2\tCg-Nd1-Ce1\tCg-Cd2-Ne2\tNd1-Ce1-Ne2\tCd2-Ne2-Ce1\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Cg\tCa-Cb-Cg-Nd1\tCa-Cb-Cg-Cd2\tCb-Cg-Nd1-Ce1\tCb-Cg-Cd2-Ne2\tCg-Nd1-Ce1-Ne2\n")
fileILE=open("Ile_geo", 'w')
fileILE.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg1\tCb-Cg2\tCg1-Cd1\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg1\tCa-Cb-Cg2\tCb-Cg1-Cd1\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Cg1\tN-Ca-Cb-Cg2\tCa-Cb-Cg1-Cd1\n")
fileLEU=open("Leu_geo", 'w')
fileLEU.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg\tCg-Cd1\tCg-Cd2\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg\tCa-Cg-Cd1\tCb-Cg-Cd2\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Cg\tCa-Cb-Cg-Cd1\tCa-Cb-Cg-Cd2\n")
fileLYS=open("Lys_geo", 'w')
fileLYS.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg\tCg-Cd\tCd-Ce\tCe-Nz\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg\tCb-Cg-Cd\tCg-Cd-Ce\tCd-Ce-Nz\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Cg\tCa-Cb-Cg-Cd\tCb-Cg-Cd-Ce\tCg-Cd-Ce-Nz\n")
fileMET=open("Met_geo", 'w')
fileMET.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg\tCg-Sd\tSd-Ce\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg\tCb-Cg-Sd\tCg-Sd-Ce\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Cg\tCa-Cb-Cg-Sd\tCb-Cg-Sd-Ce\n")
filePHE=open("Phe_geo", 'w')
filePHE.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg\tCg-Cd1\tCg-Cd2\tCd1-Ce1\tCd2-Ce2\tCe1-Cz\tCe2-Cz\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg\tCb-Cg-Cd1\tCb-Cg-Cd2\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Cg\tCa-Cb-Cg-Cd1\tCa-Cb-Cg-Cd2\tCb-Cg-Cd1-Ce1\tCb-Cg-Cd2-Ce2\tCg-Cd1-Ce1-Cz\n")
filePRO=open("Pro_geo", 'w')
filePRO.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg\tCg-Cd\tCd-N\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg\tCb-Cg-Cd\tCg-Cd-N\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Cg\tCa-Cb-Cg-Cd\tCb-Cg-Cd-N\n")
fileSER=open("Ser_geo", 'w')
fileSER.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Og\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Og\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Og\n")
fileTHR=open("Thr_geo", 'w')
fileTHR.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Og1\tCb-CG2\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Og1\tCa-Cb-Cg2\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Og1\tN-Ca-Cb-Cg2\n")
fileTRP=open("Trp_geo", 'w')
fileTRP.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg\tCg-Cd1\tCg-Cd2\tCd1-Ne1\tCd2-Ce2\tCd2-Ce3\tCe2-Ne1\tCe2-Cz2\tCe3-Cz3\tCz3-Ch2\tCz2-Ch2\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg\tCb-Cg-Cd1\tCb-Cg-Cd2\tCg-Cd2-Ce3\tNe1-Ce2-Cz2\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Cg\tCa-Cb-Cg-Cd1\tCa-Cb-Cg-Cd2\tCb-Cg-Cd1-Ne1\tCb-Cg-Cd2-Ce2\tCb-Cg-Cd2-Ce3\tCg-Cd2-Ce2-Cz2\tCg-Cd2-Ce3-Cz3\tCd2-Ce2-Cz2-Ch2\n")
fileTYR=open("Tyr_geo", 'w')
fileTYR.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg\tCg-Cd1\tCg-Cd2\tCd1-Ce1\tCd2-Ce2\tCe1-Cz\tCe2-Cz\tCz-Oh\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg\tCb-Cg-Cd1\tCb-Cg-Cd2\tCe1-Cz-Oh\tCe2-Cz-Oh\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Cg\tCa-Cb-Cg-Cd1\tCa-Cb-Cg-Cd2\tCb-Cg-Cd1-Ce1\tCb-Cg-Cd2-Ce2\tCg-Cd1-Ce1-Cz\tCd1-Ce1-Cz-Oh\n")
fileVAL=open("Val_geo", 'w')
fileVAL.write("Seq\tSA\tRSA\tPDB_ID\tResNum\toccupancy\tpsi\tphi\tomega1\tomega2\tPep\ttotalbondLength\tN-Ca\tCa-Co\tC=O\tCa-Cb\tCb-Cg1\tCb-Cg2\tN-Ca-Co\tCa-Co-O\tN-Ca-Cb\tCo-Ca-Cb\tCa-Cb-Cg1\tCa-Cb-Cg2\tN-Ca-Co-O\tN-Co-Ca-Cb\tN-Ca-Cb-Cg1\tN-Ca-Cb-Cg2\n")

fileALA.close()
fileASN.close()
fileASP.close()
fileARG.close()
fileCYS.close()
fileHIS.close()
fileGLY.close()
fileGLU.close()
fileGLN.close()
fileILE.close()
fileLEU.close()
fileLYS.close()
fileTHR.close()
fileMET.close()
fileSER.close()
filePRO.close()
filePHE.close()
fileTYR.close()
fileVAL.close()
fileTRP.close()

def write_out(filename, bonds, tline):
    filename.write(tline)
    for z in bonds:
        filename.write("\t" + str(z))
    filename.write("\n")

for line in fileP:

    splits= line.split()
    IDChain= splits[0]
    ID=IDChain[0:4]
    chain=IDChain[4]

    try:
        seq_list, res_list, acc_list, occ_list, rsa_list, d_list, angle_list, phi_list, psi_list, omega_list, ss_list, N_CA_length, CA_C_length, C_O_length, bond_length, R_length, incompleteRes, CA_C_N2_angle, O_C_N2_angle, C_N2_CA2_angle, R_angle, X_angle=get_PDB.getStruInfo(ID, chain, dssp_cmd )

        if(len(seq_list)!=len(d_list)):
            count+=1
            fileError.write(IDChain + " lengths don't match")

        else:
            fileALA=open("Ala_geo", 'a')
            fileASN=open("Asn_geo", 'a')
            fileASP=open("Asp_geo", 'a')
            fileARG=open("Arg_geo", 'a')
            fileCYS=open("Cys_geo", 'a')
            fileGLN=open("Gln_geo", 'a')
            fileGLU=open("Glu_geo", 'a')
            fileGLY=open("Gly_geo", 'a')
            fileHIS=open("His_geo", 'a')
            fileILE=open("Ile_geo", 'a')
            fileLEU=open("Leu_geo", 'a')
            fileLYS=open("Lys_geo", 'a')
            fileMET=open("Met_geo", 'a')
            filePHE=open("Phe_geo", 'a')
            filePRO=open("Pro_geo", 'a')
            fileSER=open("Ser_geo", 'a')
            fileTHR=open("Thr_geo", 'a')
            fileTRP=open("Trp_geo", 'a')
            fileTYR=open("Tyr_geo", 'a')
            fileVAL=open("Val_geo", 'a')
            c+=1
            
            for i in range(1, len(seq_list)-1):
                am=seq_list[i]
                psi=psi_list[i]
                phi=phi_list[i]
                omega1=omega_list[i]
                omega2=omega_list[i+1]
                if(incompleteRes[i]==0 and incompleteRes[i-1]==0 and incompleteRes[i+1]==0 and occ_list[i]==1 and occ_list[i-1]==1 and occ_list[i+1]==1):
                    ##5 std. away
                    if(d_list[i]<1.37 and d_list[i]>1.29):
                        Nterm= seq_list[i-1]
                    else:
                        Nterm='-'

                    if(d_list[i+1]<1.37 and d_list[i+1]>1.29):
                        Cterm=seq_list[i+1]
                    else:
                        Cterm='-'
                        
                    tri=Nterm + am + Cterm

                    start=""
                    
                    if(tri.find('-')== -1 ):
                        start=tri+ "\t" + str(acc_list[i])+ "\t" + str(rsa_list[i]) + "\t" + IDChain + "\t"+ str(res_list[i]) + "\t" + str(occ_list[i]) +"\t"+ str(psi) +"\t"+ str(phi) +"\t"+ str(omega1) +"\t"+ str(omega2) +"\t"+ str(d_list[i]) +"\t"+ str(bond_length[i]) +"\t"+ str(N_CA_length[i]) +"\t"+ str(CA_C_length[i]) +"\t"+str(C_O_length[i])
                        info=R_length[i]+R_angle[i]+X_angle[i]
                        if(am=='A'):
                            write_out(fileALA, info, start)
                        elif(am=='N'):
                            write_out(fileASN, info, start)
                        elif(am=='D'):
                            write_out(fileASP, info, start)
                        elif(am=='R'):
                            write_out(fileARG, info, start)
                        elif(am=='C'):
                            write_out(fileCYS, info, start)
                        elif(am=='G'):
                            write_out(fileGLY, info, start)
                        elif(am=='Q'):
                            write_out(fileGLN, info, start)
                        elif(am=='E'):
                            write_out(fileGLU, info, start)
                        elif(am=='H'):
                            write_out(fileHIS, info, start)
                        elif(am=='I'):
                            write_out(fileILE, info, start)
                        elif(am=='L'):
                            write_out(fileLEU, info, start)
                        elif(am=='K'):
                            write_out(fileLYS, info, start)
                        elif(am=='M'):
                            write_out(fileMET, info, start)
                        elif(am=='P'):
                            write_out(filePRO, info, start)
                        elif(am=='F'):
                            write_out(filePHE, info, start)
                        elif(am=='S'):
                            write_out(fileSER, info, start)
                        elif(am=='T'):
                            write_out(fileTHR, info, start)
                        elif(am=='W'):
                            write_out(fileTRP, info, start)
                        elif(am=='Y'):
                            write_out(fileTYR, info, start)
                        elif(am=='V'):
                            write_out(fileVAL, info, start)
                        else:
                            fileError.write("Error in: " + IDChain +" incomplete Residue "+ am + str(i) +"\n")

            fileALA.close()
            fileASN.close()
            fileASP.close()
            fileARG.close()
            fileCYS.close()
            fileHIS.close()
            fileGLY.close()
            fileGLU.close()
            fileGLN.close()
            fileILE.close()
            fileLEU.close()
            fileLYS.close()
            fileTHR.close()
            fileMET.close()
            fileSER.close()
            filePRO.close()
            filePHE.close()
            fileTYR.close()
            fileVAL.close()
            fileTRP.close()
                    
    except:
        fileError.write("Error in: " + IDChain +"\n")
        fileError.write(sys.exc_info()[0])
    
fileError.close()
print str(count)
print str(c)
