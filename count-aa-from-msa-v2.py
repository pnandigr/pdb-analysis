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


list_AA = ['G','A','V','L','I','N','C','Q','H','M','F','P','S','T','W','Y','R','K','E','D','X']
list_AA_1 = ['G','A','V','L','I','N','C','Q','H','M','F','P','S','T','W','Y','R','K','E','D']

#for seq_record in SeqIO.parse("com.pir", "pir"):
    #if r"//" in seq_record:
    #    print(seq_record)
#    seq_record.split("//")
#    print(seq_record)
    #print(seq_record.id)
    #print(repr(seq_record.seq))
    #print(len(seq_record))
    
#align1 = AlignIO.read("Dpr-com.aln", "clustal")  #read Dpr MSA
#align2 = AlignIO.read("DIP-com.aln", "clustal")  #read DIP MSA
#align3 = AlignIO.read("DprDIP-com.aln", "clustal") #read DprDIP MSA


#with open("DprDIP-com.pir") as handle:
#    for record in SeqIO.parse(handle, "pir"):
        #print("%s length %i" % (record.id, len(record)))
#        print(handle)

seq_record1 = list(SeqIO.parse("Dpr.aln", "clustal"))
seq_record2 = list(SeqIO.parse("DIP.aln", "clustal"))
seq_record3 = list(SeqIO.parse("DprDIP.aln", "clustal"))

msa_lst1 = []
msa_lst2 = []
msa_lst3 = []

for i in seq_record1:
    msa_lst1.append(i.seq)

#print(msa_lst1[0])
    
for i in seq_record2:
    msa_lst2.append(i.seq)

for i in seq_record3:
    msa_lst3.append(i.seq)    

#----------replace gaps "-" in alignment with "X"
msa_lst1 = str(msa_lst1).replace('-','X')
msa_lst2 = str(msa_lst2).replace('-','X')
msa_lst3 = str(msa_lst3).replace('-','X')

#----------convert strings to lists
def convert(string):
    li = list(string.split( ))
    return li

msa_lst_1 = convert(msa_lst1)
msa_lst_2 = convert(msa_lst2)
msa_lst_3 = convert(msa_lst3)

print(msa_lst1[0])
for i in msa_lst1[0]:
    print(i)


#print(msa_lst_2[12])

#count = 0
#for i in msa_lst_1[0]:
#    for j in msa_lst_2[0]:
        
#aa_counts_1 = Counter(msa_lst_1[1])
#df = pd.DataFrame.from_dict(aa_counts_1, orient='index')
#df.plot(kind='bar')
#plt.show()
