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

seq_record1 = list(SeqIO.parse("1.aln", "clustal"))
seq_record2 = list(SeqIO.parse("2.aln", "clustal"))
seq_record3 = list(SeqIO.parse("3.aln", "clustal"))

#-------------test with toy sequence
list_AA = ['G','A','V','L','I','N','C','Q','H','M','F','P','S','T','W','Y','R','K','E','D','X']

seq1 = list(Seq("GGGGGLLLLL",IUPAC.protein))
seq2 = list(Seq("LLGLLGGGGG",IUPAC.protein))
seq3 = list(Seq("GGGGGLLLLL",IUPAC.protein))
seq4 = list(Seq("GGGGGLLLLL",IUPAC.protein))
seq_lst = [seq1, seq2, seq3, seq4]

N_t = 2
 
#-----------initialize empty probability list for each AA    
p_AA = np.zeros((len(seq1)))    
    

#columns = ['i','j','p_A','p_B','p_AB']
#index = range(0,len(seq1))

#df = pd.DataFrame(index=index,columns=columns)

#for i in range(len(seq1)):
#    for j in range(len(seq2)):
#        if seq1[i] == seq2[i]:

for i in range(len(seq1)):
    for myseq in seq_lst:
        if i < 5: 
            if myseq[i].eq.(myseq[i+5]):
                probs.
        
    

            
p1 = np.zeros((len(seq1)))
p2 = np.zeros((len(seq2)))
#s1 = np.zeros((len(seq1),len(seq2)))

#def entropy(i,j):
#    probs = []
#    for i,c1 in enumerate(seq1):
#        for j,c2 in enumerate(seq2):
                       
            #probs.append(np.mean(np.logical_and(X == c1, Y == c2)))

#    return np.sum(-p * np.log2(p) for p in probs)
    

cs = np.zeros((len(seq1),len(seq2)))

for i,z1 in enumerate(seq1):
    for j,z2 in enumerate(seq2):
        if z1 == z2:
            #if i == j:
                cs[i,j] += 1
                #print(i,j,z1,z2)
print(cs)
        #else:
        #        cs[i,j] = 0
            
plt.matshow(cs)
plt.colorbar()
#------------------------------------

msa_lst1 = []
msa_lst2 = []
msa_lst3 = []

for i in seq_record1:
    msa_lst1.append(i.seq)                                                                                                                                                           

for i in seq_record2:
    msa_lst2.append(i.seq)

for i in seq_record3:
    msa_lst3.append(i.seq)
    
msa_lst1 = [str(s).replace('-','X') for s in msa_lst1]
msa_lst2 = [str(s).replace('-','X') for s in msa_lst2]
msa_lst3 = [str(s).replace('-','X') for s in msa_lst3]

n_A = 169 #---number of residues in each aligned Dpr
n_B = 108 #---number of residues in each aligned DIP
n_tot = 277 #---number of residues in each aligned Dpr-DIP complex

print((msa_lst3[0][:n_A]))
print(len(msa_lst3[0][:n_A]))
print((msa_lst3[0][n_A+1:n_tot]))
print(len(msa_lst3[0][n_A:]))
print(msa_lst3[0])
print(len(msa_lst3[0]))


#print(len(msa_lst3[0]))
#print(len(msa_lst1[0]))
#print(len(msa_lst2[8]))

count = np.zeros((len(msa_lst1[0]),len(msa_lst2[0])))

#for l1 in msa_lst1: # 21
#    for l2 in msa_lst2: # 9
for i,s1 in enumerate(msa_lst1[0]):
    for j,s2 in enumerate(msa_lst2[0]):
        
        if s1 == s2:
            #if i == j:
                count[i,j] += 1
                #print(i,j,s1,s2)
        #else:
        #       count[i,j] = 0
                
#print(count)

plt.matshow(count)
plt.colorbar()
