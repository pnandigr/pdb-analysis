from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from collections import Counter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import math
import random
import os
import sys
import warnings
import traceback
import argparse
from tqdm import tqdm
from itertools import product
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.interpolate import interp2d


def parseArgs():
    """Parse command line arguments"""

    try:
        parser = argparse.ArgumentParser(
            description='Calculate mutual/cross entropy between all possible pairs of residues in Dpr and DIP.')

        parser.add_argument('-a',
                            '--alignment',
                            action='store',
                            required=True,
                            help='The multiple sequence alignment (MSA) for all 189 DprDIP complexes or all cognate DprDIP complexes in any of the formats supported by Biopython\'s AlignIO.')
        parser.add_argument('-f',
                            '--alnformat',
                            action='store',
                            default='fasta',
                            help='Specify the format of the input MSA to be passed in to AlignIO.')
        parser.add_argument('-v',
                            '--verbose',
                            action='count',
                            default=0,
                            help='Verbose behaviour, printing parameters of the script.')
        parser.add_argument('-m',
                            '--makeplot',
                            action='store_true',
                            help='Plot the results via Matplotlib.')
        parser.add_argument('-o',
                            '--output',
                            action='store_true',
                            help='Save results to a text file.')
        
    except:
        print ("An exception occurred with argument parsing. Check your input arguments.")
        traceback.print_exc()

    return parser.parse_args()


def getEntropy(x):
    """Function to compute the informational entropy"""
    return -x * np.log(x)
    #return  x * (math.log(x,2))


class Prot_Seq_Entropy:
    """ This class will compute the mutual information, Shannon entropy, or joint entropy of 
        a family of protein sequences.
        Initialization requires a list of sequences and the length of the protein residues.
        Call the mutual information function to obtain the mutual information for a set of residues"""
    
    AA_resid = ['A', 'R', 'N', 'D', 'B', 'C', 
                'E', 'Q', 'Z', 'G', 'H', 'I', 
                'L', 'K', 'M', 'F', 'P', 'T', 
                'W', 'Y', 'V', 'S', 'X']
    # AA_zeros = [0] * len(AA_resid)
    ME_headers = AA_resid
    MI_headers = ["".join(map(str, prod)) for prod in product(AA_resid, repeat=2)]


    def __init__(self, seq_list, S_len):
        """Initialize with a list of sequences, and the length of a protein sequence."""
        self.seq_list = seq_list
        self.S_len = S_len
        self.MI_index = [np.arange(0, S_len).tolist(), np.arange(0, S_len).tolist()]
        self.midf_MI = pd.MultiIndex.from_product(self.MI_index, names=['resnum1', 'resnum2'])
        # Marginal entropy abbreviated ME
        # mutual information abbreviated MI
        self.df_MEC = pd.DataFrame(0.0, index=np.arange(0, S_len), columns=self.ME_headers)
        self.df_MEP  = pd.DataFrame(0.0, index=np.arange(0, S_len), columns=self.ME_headers)
        self.df_MIC = pd.DataFrame(0.0, index=self.midf_MI, columns=self.MI_headers)
        self.df_MIP  = pd.DataFrame(0.0, index=self.midf_MI, columns=self.MI_headers)
        self.df_MEC.index.rename('resnum', inplace=True)
        self.df_MEP.index.rename('resnum', inplace=True)
        self.df_JEC = pd.DataFrame(0.0, index=self.AA_resid, columns=self.ME_headers)

        
    def compute_p_resnum(self, resnum):
        """Computes the probability of an individual residue occuring in the sequence family."""
        self.df_MEC.loc[resnum, :] = 0.0
        for myseq in self.seq_list:
            self.df_MEC.loc[resnum, myseq[resnum]] = self.df_MEC.loc[resnum, myseq[resnum]] + 1.0
        self.df_MEP.loc[resnum] = (self.df_MEC.loc[resnum].multiply(1 / len(self.seq_list)))
        #self.df_MEP.loc[resnum] = (self.df_MEC.loc[resnum].multiply(1 / len(self.msa_lst_c)))
        return self.df_MEP.loc[resnum]


    def add_seq(self, new_seq_list):
        """Add additional sequences to your family."""
        self.seq_list.extend(new_seq_list)
        return new_seq_list
        
        
    def entropy_of_resnum(self, resnum):
        """Computes the Shannon entropy of a residue in a protein family. Residue number is a parameter."""
        self.compute_p_resnum(resnum)
        return self.df_MEP.loc[resnum].apply(getEntropy)
    
    
    def joint_entropy(self, resnum1, resnum2):
        """Computes the joint entropy of a residue pair, specified by parameters resnum1, and resnum2. 
            If an exception occurs here, the residue pair is not found. 
            Make sure all residue characters are in AA_resid"""
        try:
            for myseq in self.seq_list:
                pair = "".join(map(str, [myseq[resnum1], myseq[resnum2]]))
                if pair in self.df_MIC.columns:
                    self.df_MIC.loc[(resnum1, resnum2), pair] = self.df_MIC.loc[(resnum1, resnum2), pair] + 1.0
                else:
                    raise ValueError('amino acid pair '+ pair + ' not found')
        except Exception as error:
            print('Caught this error: ' + repr(error))
            
        # self.df_MIP.loc[(resnum1, resnum2)] = (self.df_MIC.loc[(resnum1, resnum2), :] 
        #                                        * (1 / len(self.seq_list)))
        self.df_MIP.loc[(resnum1, resnum2)] = self.df_MIC.loc[(resnum1, resnum2)].multiply(1/len(self.seq_list))
        return self.df_MIP.loc[(resnum1, resnum2)].apply(getEntropy)
 

    def mutual_information_pair(self, resA, resB):
        """Computes the mutual information entropy for a set of residues A, and set of residues B.
            Parameters are passed as residue indices. Returns a list object:
            [resnumA, resnumB, entropyA, entropyB, joint entropy, mutual information]"""
        try:
            self.MI_A = self.entropy_of_resnum(resA).fillna(0.0).sum()
            self.MI_B = self.entropy_of_resnum(resB).fillna(0.0).sum()
            self.MI_AB = self.joint_entropy(resA, resB).fillna(0.0).sum()
            #self.MI_tot = - ((self.MI_A + self.MI_B) - self.MI_AB)
            self.MI_tot = ((self.MI_A + self.MI_B) - self.MI_AB)
            if self.MI_tot > np.minimum(np.array(self.MI_A), np.array(self.MI_B)):
                raise ValueError('mutual information is out of bounds for ' + str(resA) + ' ' + str(resB))
            if (self.MI_A < 0.0) or (self.MI_B < 0.0) or (self.MI_AB < 0.0):
                raise ValueError('A marginal or joint entropy is out of bounds for' + str(resA) + ' or ' + str(resB))
        except Exception as error:
            print('Caught this error: ' + repr(error))
            pass
        return [resA, resB, self.MI_A, self.MI_B, self.MI_AB, self.MI_tot]
        #return [resA, resB, self.MI_tot]
        
    
    
    def mutual_information_list(self, reslistA, reslistB):
        """Computes the mutual information entropy for a set of residues A, and set of residues B.
            Parameters are passed as two lists of residues. Returns a list object."""
        self.df_MEC[:]=0.0
        self.df_MIC[:]=0.0
        self.df_MEP[:]=0.0
        self.df_MIP[:]=0.0
        self.MI_list = []
        for res_1 in reslistA:
            for res_2 in reslistB:
                self.MI_list.append(self.mutual_information_pair(res_1, res_2))
        return np.array(self.MI_list)    
    
    
    def JEC_matrix_element(self, resA, resB):
        #return self.df_MIC.loc[(resA, resB), self.df_MIC.loc[(resA, resB)]]
        return self.df_MIC.loc[(resA, resB), self.df_MIC.loc[(resA, resB)] != 0]
    
    
    def JEC_matrix(self, resnum1, resnum2):
        self.df_JEC[:] = 0.0
        for myseq in self.seq_list:
            Ai = self.AA_resid.index(myseq[resnum1])
            Bi = self.AA_resid.index(myseq[resnum2])
            self.df_JEC.iloc[Ai, Bi] = self.df_JEC.iloc[Ai, Bi] + 1.0
        return self.df_JEC

    def JEC_matrix_sparse(self, resnum1, resnum2):
        self.df_JEC_s = self.JEC_matrix(resnum1, resnum2)
        return self.df_JEC_s.loc[(self.df_JEC_s != 0).any(axis=1), (self.df_JEC_s != 0).any(axis=0)]
        #return self.df_JEC_s.loc[(self.df_JEC_s).any(axis=1), (self.df_JEC_s).any(axis=0)]

#-----------------------sequence alignment for cognate complexes
seq_record_c = list(SeqIO.parse("cognate-cured.aln", "clustal")) #all DprDIP cognate pair sequences

msa_lst_c = []

for i in seq_record_c:
    msa_lst_c.append(i.seq)    

#----------replace gaps "-" in alignment with "X"
msa_lst_c = [str(s).replace('-','X') for s in msa_lst_c] #all DprDIP cognate pairs

#------------------for all cognate pairs
msa_seq_lst_c = []
for i,j in enumerate(msa_lst_c):
    msa_seq_lst_c.append(j)
#--------------------------------------------------

myProtSeq = Prot_Seq_Entropy(msa_seq_lst_c,277)  #arguments are protein sequences and length of sequences (all alignment sequence must have same length)

lA = list(range(0,169))   #length of Dpr: 168 AA (including gaps)
lB = list(range(169,277)) #length of DIP: 109 AA (including gaps)

out_lst = myProtSeq.mutual_information_list(lA, lB)  #output DataFrame
np.save("data-cognates-positive-MI.npy",out_lst) #save output data for easy access later

#out_lst = np.load("dat-subset-cognates.npy")

with open("all-cognates-positive-MI.model.txt", 'w') as fh1:
    for x in out_lst:
        fh1.write("%d \t %d \t %f \n"%(int(x[0]),int(x[1]),x[5]))

lst1 = out_lst[:,0]
lst2 = out_lst[:,1]
order = np.argsort(out_lst[:,5])

sorted_lst = out_lst[order,:]

with open("all-cognates-sorted-positive-MI.model.txt", 'w') as fh2:
    for x in sorted_lst:
        fh2.write("%d \t %d \t %f \n"%(int(x[0]),int(x[1]),x[5]))        
    
#------------------------plot shannon entropy, mutual information, cross entropy
f1 = plt.figure()
#------------------------------create a colorbar
cmap1 = colors.ListedColormap(['red', 'blue','green', 'black', 'orange', 'yellow', 'purple'])
boundaries1 = [1.3, 1.1, 0.8, 0.5, 0.4, 0.3, 0.2, 0.0][::-1]
norm1 = colors.BoundaryNorm(boundaries1, cmap1.N, clip=True)
#----------------------------------------------
plt.scatter(out_lst[:, 0], out_lst[:, 1], c = out_lst[:,4], marker = 's', cmap = cmap1, norm = norm1)
plt.colorbar()
plt.xlabel('Dpr AA index', fontsize=12)
plt.ylabel('DIP AA index', fontsize=12)
plt.show()
f1.savefig("S_AB.cognates.pdf", bbox_inches='tight')

f2 = plt.figure()
#------------------------------create a colorbar
cmap2 = colors.ListedColormap(['red', 'blue','green', 'black', 'orange', 'yellow', 'purple'])
#boundaries2 = [-1.5, -1.2, -0.9, -0.6, -0.4, -0.3, -0.1, 0.0]
boundaries2 = [1.5, 1.2, 0.9, 0.6, 0.4, 0.3, 0.1, 0.0][::-1]
norm2 = colors.BoundaryNorm(boundaries2, cmap2.N, clip=True)
plt.scatter(out_lst[:,0],out_lst[0:,1], c = out_lst[0:,5], marker = 's', cmap = cmap2, norm = norm2)
plt.colorbar()
plt.xlabel('Dpr AA index', fontsize=12)
plt.ylabel('DIP AA index', fontsize=12)
plt.show()
f2.savefig("deltaS_AB.cognates.pdf", bbox_inches='tight')



