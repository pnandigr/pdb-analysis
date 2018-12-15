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
                'W', 'Y', 'V', 'X']
    # AA_zeros = [0] * len(AA_resid)
    ME_headers = AA_resid
    MI_headers = ["".join(map(str, prod)) for prod in product(AA_resid, repeat=2)]

    def __init__(self, seq_list, A_len):
        """Initialize with a list of sequences, and the length of a protein sequence."""
        self.seq_list = seq_list
        self.A_len = A_len
        self.MI_index = [np.arange(0, 2*A_len).tolist(), np.arange(0, 2*A_len).tolist()]
        self.midf_MI = pd.MultiIndex.from_product(self.MI_index, names=['resnum1', 'resnum2'])
        # Marginal entropy abbreviated ME
        # mutual information abbreviated MI
        self.df_MEC = pd.DataFrame(0.0, index=np.arange(0, 2*A_len), columns=self.ME_headers)
        self.df_MEP  = pd.DataFrame(0.0, index=np.arange(0, 2*A_len), columns=self.ME_headers)
        self.df_MIC = pd.DataFrame(0.0, index=self.midf_MI, columns=self.MI_headers)
        self.df_MIP  = pd.DataFrame(0.0, index=self.midf_MI, columns=self.MI_headers)
        self.df_MEC.index.rename('resnum', inplace=True)
        self.df_MEP.index.rename('resnum', inplace=True)



        def compute_p_resnum(self, resnum):
            """Computes the probability of an individual reside occuring in the sequence family."""
            self.df_MEC.loc[resnum, :] = 0.0
            for myseq in self.seq_list:
                self.df_MEC.loc[resnum, myseq[resnum]] = self.df_MEC.loc[resnum, myseq[resnum]] + 1.0
            self.df_MEP.loc[resnum] = (self.df_MEC.loc[resnum].multiply(1 / len(self.seq_list)))
    
    
        def add_seq(self, new_seq_list):
            """Add additional sequences to your family."""
            self.seq_list.extend(new_seq_list)


        def entropy_of_resnum(self, resnum):
            """Computes the Shannon entropy of a residue in a protein family. Residue number is a parameter."""
            self.compute_p_resnum(resnum)
            return self.df_MEP.loc[resnum].apply(getEntropy)


        def joint_entropy(self, resnum1, resnum2):
            """Computes the joint entropy of a residue pair, specified by parameters resnum1, and resnum2. 
            If an exception occurs here, the residue pair is not found. 
            Make sure all residue characters are in AA_resid"""

            self.df_MIC.loc[resnum1,:] = 0.0
            self.df_MIC.loc[resnum2,:] = 0.0
            try:
                for myseq in self.seq_list:
                    pair = "".join(map(str, [myseq[resnum1], myseq[resnum2]]))
                    if pair in self.df_MIC.columns:
                        self.df_MIC.loc[(resnum1, resnum2), pair] = self.df_MIC.loc[(resnum1, resnum2), pair] +1.0
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
            Parameters are passed as residues indicies Returns a list object:
            [resnumA, resnumB, Mutual Information, joint entropy, entropyA, entropyB]"""
        try:
            self.MI_A = self.entropy_of_resnum(resA).fillna(0.0).sum()
            self.MI_B = self.entropy_of_resnum(resB).fillna(0.0).sum()
            self.MI_AB = self.joint_entropy(resA, resB).fillna(0.0).sum()
            self.MI_tot = (self.MI_B + self.MI_A) - self.MI_AB
            if self.MI_tot > np.minimum(np.array(self.MI_A), np.array(self.MI_B)):
                raise ValueError('mutual information is out of bounds for ' + str(resA) + ' ' + str(resB))
            if (self.MI_A < 0.0) or (self.MI_B < 0.0) or (self.MI_AB < 0.0):
                raise ValueError('A marginal or joint entropy is out of bounds for' + str(resA) + ' or ' + str(resB))
        except Exception as error:
            print('Caught this error: ' + repr(error))
            pass
        return [resA, resB, self.MI_tot, self.MI_AB, self.MI_A, self.MI_B]


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



#-------------test with toy sequence
list_AA = ['G','A','V','L','I','N','C','Q','H','M','F','P','S','T','W','Y','R','K','E','D','X']

seq1 = list(Seq("GGGGGLLLLL",IUPAC.protein))
seq2 = list(Seq("LLGLLGGGGG",IUPAC.protein))
seq3 = list(Seq("VGGGGVLLLL",IUPAC.protein))
seq4 = list(Seq("GGGGGLLLLL",IUPAC.protein))
seq_lst = [seq1, seq2, seq3, seq4]

N_t = len(seq_lst)

myProtSeq = Prot_Seq_Entropy(seq_lst, 5)
# print(myProtSeq.midf_MI)
# print(myProtSeq.df_mutual_information_count)
# print(myProtSeq.df_mutual_information_count.loc[(0, 5), ['GL', 'LG']])
# print(e15)
lA = [0, 1, 2, 3, 4]
lB = [5, 6, 7, 8, 9]
# lA = [0]
# lB = [5]
print(myProtSeq.mutual_information_list(lA, lB))


