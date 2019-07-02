import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors

indx_A = []
indx_B = []
entropy = []
resA = []
resB = []

# with open ('interface_entropy_msa_indx.Dpr1DIPeta.dat') as fobj:
with open ('interface_entropy_pdb_indx.Dpr1DIPeta.dat') as fobj:
    for line in fobj:
        row = line.split()
        # print(row[0])
        indx_A.append(row[0])
        resA.append(row[1])
        indx_B.append(row[2])
        resB.append(row[3])
        entropy.append(row[4])

out_lst = zip(indx_A, resA, indx_B, resB, entropy)
        
# arr_resA = np.array(resA)
# arr_resB = np.array(resB)
# arr_entropy = np.array(entropy)

arr = np.array([(i1, r1, i2, r2, x) for i1, r1, i2, r2, x in out_lst])

lst1 = arr[:, 0]
indx1 = arr[:, 1]
lst2 = arr[:, 2]
indx2 = arr[:, 3]
order = np.argsort(arr[:, 4])

sorted_lst = arr[order, :]

# print(sorted_lst)

f1 = plt.figure()
cmap1 = colors.ListedColormap(['red', 'blue', 'green', 'black', 'orange', 'purple'])
# boundaries = [1.3, 1.1, 0.8, 0.5, 0.4, 0.3, 0.2, 0.0]
boundaries1 = [-0.77, -0.65, -0.55, -0.50, -0.45, -0.40, -0.35]
norm1 = colors.BoundaryNorm(boundaries1, cmap1.N, clip=True)
#----------------------------------------------
# plt.pcolormesh(out_lst[:,0],out_lst[0:,1],out_lst[0:,4], cmap=cmap, norm=norm)
# plt.scatter(out_lst[:, 0], out_lst[:, 1], c = out_lst[:,4], marker = 's')
# plt.scatter(sorted_resA, sorted_resB, c = sorted_entropy, marker = 's', cmap = cmap1, norm = norm1)
plt.scatter(arr[:, 1], arr[:, 3], c=arr[:, 4], marker='s', cmap=cmap1, norm=norm1)
plt.colorbar()
plt.xlabel('Dpr index in pdb file', fontsize=12)
plt.ylabel('DIP index in pdb file', fontsize=12)
plt.show()
f1.savefig("interface_pdb_indx.Dpr1DIPeta.model.pdf", bbox_inches='tight')
