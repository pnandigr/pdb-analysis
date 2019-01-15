#!/usr/bin/env python

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.PDB.Polypeptide import *
from collections import Counter
import Bio
from Bio.SeqUtils import seq1, seq3
from Bio.Seq import Seq
from Bio.Alphabet import ProteinAlphabet, ThreeLetterProtein
import pylab as pyl
import numpy as np
import sys
import os
import re
import argparse
import glob

#-------------------------print all pdb files in the current directory
#path = "/home/raj/clustalw-2.1-linux-x86_64-libcppstatic/calc_mutual_information_entropy/*.pdb"
#for f in glob.glob(path):
#    print(f)

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


#three_letter["S"] will now return "SER"
#three_letter = dict([[v,k] for k,v in one_letter.items()])

def parseArgs():
    """Parse command line arguments"""

    try:
        parser = argparse.ArgumentParser(
            description = 'Read and extract items from input PDB file and sequence file.')

        parser.add_argument('-i1',
                            '--inputPDB',
                            action='store',
                            required=True,
                            help='input PDB file in standard format')

        #parser.add_argument('-i2',
        #                    '--inputSeq',
        #                    action='store',
        #                    required=True,
        #                    help='input sequence file in standard format')

    except:
        print ("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()
        
    return parser.parse_args()


###############################################################################################
# Reads a PDB file and returns the residue name and coordinates for each C-alpha atom
# (the input argument for this routine is the pdb file name.)
###############################################################################################

def get_coordinates_PDB(PDB_File_In):
    try:
        fl = open(PDB_File_In,'r')
    except:
        print('Could not open input file {0}'.format(PDB_File_In))
        sys.exit()
    Res = []
    Resid = []
    Points = []
    #Getting from a PDB file
    for line in fl:
        if not(line.startswith('ATOM')):
            continue
        elif (line[13:15] != 'CA'):
            continue
        resname = line[17:20]
        resindex = line[23:26]
        xyz = re.findall('[-+]?\d+\.\d+', line)
        tmp = np.zeros(3)
        Res.append(resname)
        
        Resid.append(resindex)
        
        tmp[0] = float(xyz[0])
        tmp[1] = float(xyz[1])
        tmp[2] = float(xyz[2])
        Points.append(tmp)
        
    fl.close()
    return Res


#def shorten(x):
#    if len(x) % 3 != 0: 
#        raise ValueError('Input length should be a multiple of three')

#    y = ''
#    for i in range(len(int(x)/3)):
#            y += d[x[3*i:3*i+3]]
#    return y

#------------------------this function is not working yet
#def get_seq(seq_File_In):

#    try:
#        fl_seq = open(seq_File_In, 'r')
#        file_contents = fl_seq.read()
        #print(file_contents)
        #seq_record_c = list(SeqIO.parse(open(seq_File_In, "clustal"))) #all DprDIP cognate pair sequences
        #seq_record_c = SeqIO.parse(file_contents, "clustal") #all DprDIP cognate pair sequences
 #       for seq_record in SeqIO.parse(file_contents, "clustal"):
 #           print(repr(seq_record.seq))
 #           print(len(seq_record))
        
 #   except:
 #       print('Could not open input file {0}'.format(seq_File_In))
 #       sys.exit()

 #       msa_lst_c = []
 #       msa_seq_lst_c = []
 #       for i in seq_record_c:
 #           msa_lst_c.append(i.seq)

#        for i,j in enumerate(msa_lst_c):
#            msa_seq_lst_c.append(j)
#            
#    return msa_seq_lst_c
#---------------------------------------------------
            
        
def main():
    """Read and parse input PDB and Seq file."""

    #Parse arguments
    args = parseArgs()

    PDB_File_In = args.inputPDB
    #seq_File_In = args.inputSeq

    seq_record_c = list(SeqIO.parse("cognate.aln", "clustal")) #all DprDIP cognate pair sequences
    
    msa_lst_c = []
    msa_seq_lst_c = []
    
    for i in seq_record_c:
        msa_lst_c.append(i.seq)
        
    for i,j in enumerate(msa_lst_c):
        msa_seq_lst_c.append(j)
            
            
    with open("seq-out-35.txt", 'w') as fs:
        for seq,index in enumerate(msa_seq_lst_c[35]):
            fs.write('{0} \t {1} \n'.format(index,seq))

    seq_from_pdb = get_coordinates_PDB(PDB_File_In) #list object
    seq_from_msa = msa_seq_lst_c[0]  #string
    seq_from_pdb1 = ''.join(seq_from_pdb) #convert list to string
    #aa_s = shorten(seq_from_pdb1)

    #-------------------debug
    #print(seq_from_pdb1)
    #print(seq_from_msa[0])
    #print(len(seq_from_msa))
    #------------------------
    
    #one_lett_seq = Bio.SeqUtils.IUPACData.protein_letters_3to1[seq_from_pdb1]
    #print(one_lett_seq)
    seq_from_pdb_c = "".join([three_to_one(aa3) for aa3 in ["".join(g) for g in zip(*(iter(seq_from_pdb1),) * 3)]])

    #------------------debug
    #print(seq_from_pdb_c[0])
    #print(len(seq_from_pdb_c))
    #------------------------
    
    mod_indx = []
    rescaled = []

    for indx_p,res_p in enumerate(seq_from_msa):
        for indx_s,res_s in enumerate(seq_from_pdb_c):
            if res_p == res_s:
                print(res_p)
                mod_indx[indx_p] = indx_s
                break
            else:
                next
        rescaled.append(mod_indx)
    print(rescaled)


    #list3 = []

    #for i in range(len(seq_from_msa)):
    #    if seq_from_msa[i] not in seq_from_pdb_c:
    #        pass
    #    else:
    #        list3.append(seq_from_msa[i])
    #        print(list3)
                   
    

if __name__ == '__main__':
    main()
