import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import os

indx_A_pdb = []
indx_B_pdb = []
indx_A_map_pdb = []
indx_B_map_pdb = []
indx_A_map_seq = []
indx_B_map_seq = []
res_A = []
res_B = []
MI = []

fout = open("sorted-data.txt","w")

with open ('combined-data.txt') as fobj:
    for line in fobj:
        row = line.split()
        indx_A_pdb.append(row[0])       #actual index in pdb file for Dpr
        indx_B_pdb.append(row[1])       #actual index in pdb file for DIP
        indx_A_map_pdb.append(row[2])   #index of Dpr for pdb file in mapping
        indx_B_map_pdb.append(row[3])   #index of DIP for pdb file in mapping
        indx_A_map_seq.append(row[4])   #index of Dpr for seq alignment file in mapping
        indx_B_map_seq.append(row[5])   #index of DIP for seq alignment file in mapping
        res_A.append(row[6])            #AA for Dpr (same in pdb file and seq alignment)
        res_B.append(row[7])            #AA for DIP (same in pdb file and seq alignment)
        MI.append(row[8])               #calculated mutual information

out_lst = zip(indx_A_pdb,indx_B_pdb,indx_A_map_pdb,indx_B_map_pdb,indx_A_map_seq,indx_B_map_seq,res_A,res_B,MI)

arr = np.array([(i1_pdb,i2_pdb,i1_map_pdb,i2_map_pdb,i1_map_seq,i2_map_seq,r1,r2,entr) for i1_pdb,i2_pdb,i1_map_pdb,i2_map_pdb,i1_map_seq,i2_map_seq,r1,r2,entr in out_lst])

lstA_pdb = arr[:,0]
lstB_pdb = arr[:,1]
lstA_map_pdb = arr[:,2]
lstB_map_pdb = arr[:,3]
lstA_map_seq = arr[:,4]
lstB_map_seq = arr[:,5]
res_A = arr[:,6]
res_B = arr[:,7]
entropy = arr[:,8]

order = np.argsort(arr[:,8])

sorted_lst = arr[order,:]

for x in sorted_lst:
    #print(x[0],'\t',x[1],'\t',x[2],'\t',x[3],'\t',x[4],'\t',x[5],'\t',x[6])
    fout.write("%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n"%(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]))
