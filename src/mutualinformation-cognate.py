from Bio import SeqIO
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd
# from itertools import product

# User coded modules for protein structures
# from protein_data_structure import parseArgs
# from protein_data_structure import getEntropy
from protein_data_structure import Prot_Seq_Entropy

from myconfig import *
# filename = "../data/cognate-cured.aln"
# lA_0 = 0
# lA_end = 169
# lB_0 = 169
# lB_end = 277
# outfile = "../output/data-cognates-positive-MI.npy"
# outfile_sorted = "../output/data-cognates-sorted-positive-MI.npy"
# figure1_file = "../output/S_AB.cognates.pdf"
# figure2_file = "../output/deltaS_AB.cognates.pdf"

# -----------------------sequence alignment for cognate complexes
seq_record_c = list(SeqIO.parse(filename, "clustal"))  # all DprDIP cognate pair sequences

msa_lst_c = []

for i in seq_record_c:
    msa_lst_c.append(i.seq)    

# ----------replace gaps "-" in alignment with "X"
msa_lst_c = [str(s).replace('-', 'X') for s in msa_lst_c]  # all DprDIP cognate pairs

# ------------------for all cognate pairs
msa_seq_lst_c = []
for i, j in enumerate(msa_lst_c):
    msa_seq_lst_c.append(j)
# --------------------------------------------------

# arguments are protein sequences and length of sequences (all alignment sequence must have same length)
myProtSeq = Prot_Seq_Entropy(msa_seq_lst_c, lB_end)

lA = list(range(lA_0, lA_end))     # length of Dpr: 168 AA (including gaps)
lB = list(range(lB_0, lB_end))   # length of DIP: 109 AA (including gaps)

out_lst = myProtSeq.mutual_information_list(lA, lB)  # output DataFrame
np.save(outfile, out_lst)    # save output data for easy access later

# out_lst = np.load("dat-subset-cognates.npy")

with open(outfile, 'w') as fh1:
    for x in out_lst:
        fh1.write("%d \t %d \t %f \n" % (int(x[0]), int(x[1]), x[5]))

lst1 = out_lst[:, 0]
lst2 = out_lst[:, 1]
order = np.argsort(out_lst[:, 5])

sorted_lst = out_lst[order, :]

with open(outfile-sorted, 'w') as fh2:
    for x in sorted_lst:
        fh2.write("%d \t %d \t %f \n" % (int(x[0]), int(x[1]), x[5]))
    
# ------------------------plot shannon entropy, mutual information, cross entropy
f1 = plt.figure()
# ------------------------------create a colorbar
cmap1 = colors.ListedColormap(['red', 'blue', 'green', 'black', 'orange', 'yellow', 'purple'])
boundaries1 = [1.3, 1.1, 0.8, 0.5, 0.4, 0.3, 0.2, 0.0][::-1]
norm1 = colors.BoundaryNorm(boundaries1, cmap1.N, clip=True)
# ----------------------------------------------
plt.scatter(out_lst[:, 0], out_lst[:, 1], c=out_lst[:, 4], marker='s', cmap=cmap1, norm=norm1)
plt.colorbar()
plt.xlabel('Dpr AA index', fontsize=12)
plt.ylabel('DIP AA index', fontsize=12)
plt.show()
f1.savefig(figure1_file, bbox_inches='tight')

f2 = plt.figure()
# ------------------------------create a colorbar
cmap2 = colors.ListedColormap(['red', 'blue', 'green', 'black', 'orange', 'yellow', 'purple'])
# boundaries2 = [-1.5, -1.2, -0.9, -0.6, -0.4, -0.3, -0.1, 0.0]
boundaries2 = [1.5, 1.2, 0.9, 0.6, 0.4, 0.3, 0.1, 0.0][::-1]
norm2 = colors.BoundaryNorm(boundaries2, cmap2.N, clip=True)
plt.scatter(out_lst[:, 0], out_lst[0:, 1], c=out_lst[0:, 5], marker='s', cmap=cmap2, norm=norm2)
plt.colorbar()
plt.xlabel('Dpr AA index', fontsize=12)
plt.ylabel('DIP AA index', fontsize=12)
plt.show()
f2.savefig(figure2_file, bbox_inches='tight')



