#!/usr/bin/env python

from Bio import SeqIO
from Bio.PDB.Polypeptide import *
import numpy as np
import sys
import os
import re
import argparse
import glob
import traceback

# -------------------------print all pdb files in the current directory
# path = "/home/raj/clustalw-2.1-linux-x86_64-libcppstatic/calc_mutual_information_entropy/*.pdb"
# for f in glob.glob(path):
#     print(f)

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


# three_letter["S"] will now return "SER"
# three_letter = dict([[v,k] for k,v in one_letter.items()])

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

        # parser.add_argument('-i2',
        #                     '--inputSeq',
        #                     action='store',
        #                     required=True,
        #                     help='input sequence file in standard format')

    except:
        print("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()
        
    return parser.parse_args()


###############################################################################################
# Reads a PDB file and returns the residue name and coordinates for each C-alpha atom
# (the input argument for this routine is the pdb file name.)
###############################################################################################

def get_coordinates_PDB(PDB_File_In):
    try:
        fl = open(PDB_File_In, 'r')
    except:
        print('Could not open input file {0}'.format(PDB_File_In))
        sys.exit()
    Res = []
    Resid = []
    Points = []
    # Getting from a PDB file
    for line in fl:
        if not(line.startswith('ATOM')):
            continue
        elif line[13:15] != 'CA':
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
    for i, rescode in enumerate(seq_from_msa):
        if rescode != '-':
            map_pdb_to_seq[c] = i
            map_seq_to_pdb[i] = c
            c += 1
    return map_pdb_to_seq, map_seq_to_pdb

        
def main():
    """Read and parse input PDB and Seq file."""

    # Parse arguments
    args = parseArgs()

    PDB_File_In = args.inputPDB
    # seq_File_In = args.inputSeq

    # all DprDIP cognate pair sequences
    seq_record_c = list(SeqIO.parse("cognate.aln", "clustal")) 
    
    msa_lst_c = []
    msa_seq_lst_c = []
    
    for i in seq_record_c:
        msa_lst_c.append(i.seq)
        
    for i, j in enumerate(msa_lst_c):
        msa_seq_lst_c.append(j)
                
    # with open("seq-out-35.txt", 'w') as fs:
    #     for seq,index in enumerate(msa_seq_lst_c[35]):
    #         fs.write('{0} \t {1} \n'.format(index,seq))

    seq_from_pdb = get_coordinates_PDB(PDB_File_In)  # list object
    
    seq_from_pdb1 = ''.join(seq_from_pdb)  # convert list to string
    
    seq_from_pdb_c = "".join([three_to_one(aa3) for aa3 in ["".join(g) for g in zip(*(iter(seq_from_pdb1),) * 3)]])

    seq_from_msa = msa_seq_lst_c[20]  # string

    msa_seq = list(seq_from_msa)
    pdb_seq = list(seq_from_pdb_c)

    print(seq_from_msa)
    print(msa_seq)
    print(pdb_seq)

    # ---------------call function to do mapping and inverse mapping
    # map,_ = make_map(seq_from_msa)
    # map,map_inv = make_map(seq_from_msa)

    map_pdb_to_seq, map_seq_to_pdb = make_map(seq_from_msa)

    # ---------------------------write output to screen
    # for i,rescode in enumerate(pdb_seq):
    #     #print(rescode, "\t", i, "\t", seq_from_msa[map[i]], "\t", map[i])
    #     print(rescode, "\t", i, "\t", seq_from_msa[map_pdb_to_seq[i]], "\t", map_pdb_to_seq[i])

    # --------------------------write outout to file
    with open("map_pdb_to_seq-Dpr11DIPgamma.model.txt", "w") as fl:
        for i,rescode in enumerate(pdb_seq):
            # print(rescode, "\t", i, "\t", seq_from_msa[map[i]], "\t", map[i])
            # print(rescode, "\t", i, "\t", seq_from_msa[map_pdb_to_seq[i]], "\t", map_pdb_to_seq[i])
            fl.write("%s \t %d \t %s \t %d \n"%(rescode, i, seq_from_msa[map_pdb_to_seq[i]], map_pdb_to_seq[i]))
            
    # ---------------------------write output to screen
    # for i,rescode in enumerate(seq_from_msa):
    #     if i in map_inv:
    #         print(rescode, "\t", i, "\t", pdb_seq[map_inv[i]], "\t", map_inv[i])

    # --------------------------write outout to file
    with open("map_seq_to_pdb-Dpr11DIPgamma.model.txt", "w") as fh:
        for i,rescode in enumerate(seq_from_msa):
            if i in map_seq_to_pdb:
                # print(rescode, "\t", i, "\t", pdb_seq[map_seq_to_pdb[i]], "\t", map_seq_to_pdb[i])
                fh.write("%s \t %d \t %s \t %d \n"%(rescode, i, pdb_seq[map_seq_to_pdb[i]], map_seq_to_pdb[i]))


if __name__ == '__main__':
    main()

