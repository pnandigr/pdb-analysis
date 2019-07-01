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


###############################################################################################
# Reads a PDB file and returns the residue name and coordinates for each C-alpha atom
# (the input argument for this routine is the pdb file name.)
###############################################################################################

def get_coordinates_PDB(PDB_File_In):

    #list_of_files = glob.glob('./*.pdb')
    #for PDB_File_In in list_of_files:
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

    
def make_map(seq_from_msa):
    map_pdb_to_seq = {}
    map_seq_to_pdb = {}
    c = 0
    for i,rescode in enumerate(seq_from_msa):
        if rescode != '-':
            map_pdb_to_seq[c] = i
            map_seq_to_pdb[i] = c
            c += 1
    return map_pdb_to_seq, map_seq_to_pdb

        
def main():
    """Read and parse input PDB and Seq file."""

    #Parse arguments
    #args = parseArgs()

    #PDB_File_In = args.inputPDB
    #seq_File_In = args.inputSeq

    seq_record_c = list(SeqIO.parse("cognate.aln", "clustal")) #all DprDIP cognate pair sequences
    
    msa_lst_c = []
    msa_seq_lst_c = []
    
    for i in seq_record_c:
        msa_lst_c.append(i.seq)
        
    for i,j in enumerate(msa_lst_c):
        msa_seq_lst_c.append(j)
            
            
    with open("seq-out-mult.txt", 'w') as fs:
        for seq,index in enumerate(msa_seq_lst_c[0]):
            fs.write('{0} \t {1} \n'.format(index,seq))

    mult_pdb_list = []
    list_of_files = glob.glob('./cognate-pdbs/*.pdb')
    for file_name in list_of_files:
        seq_from_pdb = get_coordinates_PDB(file_name) #list object
        mult_pdb_list.append(seq_from_pdb)
        

    seq_from_msa = []
    for lst in msa_seq_lst_c:
        #for i in lst:
        seq_from_msa.append(lst)

    pdb_str = str(mult_pdb_list)
    seq_from_pdb1 = ''.join(pdb_str) #convert list to string
    print(seq_from_pdb1)
    
    seq_from_pdb_c = "".join([three_to_one(aa3) for aa3 in ["".join(g) for g in zip(*(iter(seq_from_pdb1),) * 3)]])  #convert three letter aa codes to one letter aa codes
    print(seq_from_pdb_c)
    
    msa_seq = list(seq_from_msa)
    pdb_seq = list(seq_from_pdb_c)
    
    #---------------call function to do mapping and inverse mapping
    #map,_ = make_map(seq_from_msa)
    #map,map_inv = make_map(seq_from_msa)

    map_pdb_to_seq, map_seq_to_pdb = make_map(seq_from_msa)

    #--------------------------write outout to file
    with open("map_pdb_to_seq-mult.txt", "w") as fl:
        for i,rescode in enumerate(pdb_seq):
            fl.write("%s \t %d \t %s \t %d \n"%(rescode, i, seq_from_msa[map_pdb_to_seq[i]], map_pdb_to_seq[i]))
            

    #--------------------------write outout to file
    with open("map_seq_to_pdb-mult.txt", "w") as fh:
        for i,rescode in enumerate(seq_from_msa):
            if i in map_seq_to_pdb:
                fh.write("%s \t %d \t %s \t %d \n"%(rescode, i, pdb_seq[map_seq_to_pdb[i]], map_seq_to_pdb[i]))

if __name__ == '__main__':
    main()
