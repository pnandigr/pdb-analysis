from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from collections import Counter
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
    return Res, Resid


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
#        return msa_seq_lst_c
            

seq_record_c = list(SeqIO.parse("cognate.aln", "clustal")) #all DprDIP cognate pair sequences

msa_lst_c = []
msa_seq_lst_c = []

for i in seq_record_c:
    msa_lst_c.append(i.seq)

for i,j in enumerate(msa_lst_c):
    msa_seq_lst_c.append(j)


#for i in range(len(seq_record_c)):
#for seq,i in enumerate(seq_record_c):
#    print(i,seq)

with open("seq-out-35.txt", 'w') as fs:
    for seq,index in enumerate(msa_seq_lst_c[35]):
        fs.write('{0} \t {1} \n'.format(index,seq))
        #print(index,seq)

def main():
    """Read and parse a provided PDB file."""


    #Parse arguments
    args = parseArgs()

    PDB_File_In = args.inputPDB
    #seq_File_In = args.inputSeq
    
    print(get_coordinates_PDB(PDB_File_In))
    #print(get_seq(seq_File_In))

if __name__ == '__main__':
    main()
