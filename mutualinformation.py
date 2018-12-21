# All the important Import statements

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from collections import Counter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import os
from tqdm import tqdm
from itertools import product

def getEntropy(x):
    """Function to compute the informational entropy"""
    return - x * np.log(x)


class Prot_Seq_Entropy:
    """ This class will compute the mutual information, Shannon entropy, or joint entropy of 
        a family of protein sequences.
        Initialization requires a list of the sequences and the length of the protein residues.
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
            Parameters are passed as residues indices. Returns a list object:
            [resnumA, resnumB, Mutual Information, joint entropy, entropyA, entropyB]"""
        try:
            self.MI_A = self.entropy_of_resnum(resA).fillna(0.0).sum()
            self.MI_B = self.entropy_of_resnum(resB).fillna(0.0).sum()
            self.MI_AB = self.joint_entropy(resA, resB).fillna(0.0).sum()
            self.MI_tot = - ((self.MI_A + self.MI_B) - self.MI_AB)
            if self.MI_tot > np.minimum(np.array(self.MI_A), np.array(self.MI_B)):
                raise ValueError('mutual information is out of bounds for ' + str(resA) + ' ' + str(resB))
            if (self.MI_A < 0.0) or (self.MI_B < 0.0) or (self.MI_AB < 0.0):
                raise ValueError('A marginal or joint entropy is out of bounds for' + str(resA) + ' or ' + str(resB))
        except Exception as error:
            print('Caught this error: ' + repr(error))
            pass
        return [resA, resB, self.MI_A, self.MI_B, self.MI_AB, self.MI_tot]
        
    
    
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
        
    
#---------------------------test with toy sequence
seq1 = list(Seq("GVGGGLGLLLG",IUPAC.protein))
seq2 = list(Seq("LLGLLGGGGGG",IUPAC.protein))
seq3 = list(Seq("VGGGGVGLLLG",IUPAC.protein))
seq4 = list(Seq("GGLGGLLLLLG",IUPAC.protein))
seq5 = list(Seq("GGLXGLLLLLG",IUPAC.protein))
seq_lst = [seq1, seq2, seq3, seq4, seq5]

#-----------------------sequence alignment for Dpr, DIP, and DprDIP
seq_record1 = list(SeqIO.parse("1.aln", "clustal"))
seq_record2 = list(SeqIO.parse("2.aln", "clustal"))
seq_record3 = list(SeqIO.parse("3.aln", "clustal"))
seq_record_c = list(SeqIO.parse("cognate.aln", "clustal"))

msa_lst1 = []
msa_lst2 = []
msa_lst3 = []
msa_lst_c = []

for i in seq_record1:
    msa_lst1.append(i.seq)
    
for i in seq_record2:
    msa_lst2.append(i.seq)

for i in seq_record3:
    msa_lst3.append(i.seq)

for i in seq_record_c:
    msa_lst_c.append(i.seq)    

#----------replace gaps "-" in alignment with "X"
msa_lst1 = [str(s).replace('-','X') for s in msa_lst1]
msa_lst2 = [str(s).replace('-','X') for s in msa_lst2]
msa_lst3 = [str(s).replace('-','X') for s in msa_lst3]
msa_lst_c = [str(s).replace('-','X') for s in msa_lst_c]

#sq_lst = [msa_lst3[0],msa_lst3[1],msa_lst3[2],msa_lst3[3],msa_lst3[4],msa_lst3[5]]

#------------------for all 189 complexes
msa_seq_lst = []
for i,j in enumerate(msa_lst3):
    msa_seq_lst.append(j)
#--------------------------------------------------

#------------------for all cognate pairs
msa_seq_lst_c = []
for i,j in enumerate(msa_lst_c):
    msa_seq_lst_c.append(j)
#--------------------------------------------------
    
#myProtSeq = Prot_Seq_Entropy(seq_lst,6)
myProtSeq = Prot_Seq_Entropy(msa_seq_lst_c, 277)  #arguments are protein sequences and length of sequences
# print(myProtSeq.midf_MI)
# print(myProtSeq.df_mutual_information_count)
# print(myProtSeq.df_mutual_information_count.loc[(0, 5), ['GL', 'LG']])

#---------------test with smaller sequence length
#lA = [0, 1, 2, 3, 4]
#lB = [5, 6, 7, 8, 9, 10]

# lA = [0]
# lB = [5]

#lA = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]
#lB = [27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55]

#lA = list(range(0,5))   
#lB = list(range(6,15)) 

lA = list(range(0,169))   #length of protein A
lB = list(range(169,277)) #length of protein B

#f1 = open('p_res-prot-A.txt','w')
#f2 = open('p_res-prot-B.txt','w')

print(myProtSeq.compute_p_resnum(169))  #argument is column number/position in a sequence in protein A
print(myProtSeq.compute_p_resnum(276))  #argument is column number/position in a sequence in protein B
#print(myProtSeq.mutual_information_list(lA_lst, lB_lst)) #arguments are all positions in sequence A and all positions in sequence B
print(myProtSeq.mutual_information_list(lA,lB)) #arguments are all positions in sequence A and all positions in sequence B
#print(myProtSeq.JEC_matrix_element(16,180)) #arguments are the positions in sequence A and sequence B
#mat36 = myProtSeq.JEC_matrix(2,8)  #arguments are the positions in sequence A and sequence B
print(myProtSeq.JEC_matrix_sparse(16,180))  #arguments are the positions in sequence A and sequence B
print(myProtSeq.JEC_matrix_sparse(100,245))
