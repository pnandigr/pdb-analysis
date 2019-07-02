from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# User coded modules for protein structures
# then import the user configuration settings
from protein_data_structure import Prot_Seq_Entropy
from myconfig import *

# -----------------------sequence alignment for cognate complexes
seq_record_c = list(SeqIO.parse("cognate-cured.aln", "clustal"))  # all DprDIP cognate pair sequences

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
myProtSeq = Prot_Seq_Entropy(msa_seq_lst_c,
                             lB_end)


lA = list(range(lA_0, lA_end))  # length of Dpr: 168 AA (including gaps)
lB = list(range(lB_0, lB_end))  # length of DIP: 109 AA (including gaps)

# out_lst = myProtSeq.mutual_information_list(lA, lB)  #output DataFrame
# np.save("data-cognates-positive-MI.npy",out_lst) #save output data for easy access later

out_lst = np.load("data-cognates-positive-MI.npy")

sq_mat_joint1 = myProtSeq.JEC_matrix(targetA1, targetB1)

# sq_mat_A = myProtSeq.compute_p_resnum(71)
# print(sq_mat_joint)
# sns.set_style("whitegrid")

# ---------------attempt to discretize the plot
# df_q = pd.DataFrame()
# for col in sq_mat_joint:
#    df_q[col] = pd.to_numeric( pd.qcut(sq_mat_joint[col], 10, labels=list(range(10))) )

# sns.heatmap(df_q)
# ---------------------------------------
# sns.heatmap(sq_mat_joint, cmap='PiYG')
# plt.gca().invert_yaxis()  #put y-axis labels in reverse order
# plt.savefig("plotmat-MI-h10.pdf", bbox_inches='tight')

# -----------------------------------plot heatmap
# f = plt.figure()
# cmap1 = colors.ListedColormap(['red', 'blue','green', 'black', 'orange', 'yellow', 'purple'])
# boundaries1 = [100, 80, 60, 45, 35, 25, 10, 0.0][::-1]
# norm1 = colors.BoundaryNorm(boundaries1, cmap1.N, clip=True)
sns.heatmap(np.add(myProtSeq.JEC_matrix(targetA1, targetB1),
                   myProtSeq.JEC_matrix(targetA1, targetB2)), cmap='coolwarm', annot=True)
plt.gca().invert_yaxis()  # put y-axis labels in reverse order
# plt.savefig("plotmat-MI-h10.pdf", bbox_inches='tight')
plt.show()


