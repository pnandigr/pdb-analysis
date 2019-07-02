#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import os,sys
import numpy as np

seqList=[]
seqIdx=[]
pdbList=[]
pdbIdx=[]
res1_res2_map={}
seqtopdb = {}
seqtopdb_aa = {}

def readSeqTOPDBList(inF):
    fin=open(inF,"r")
    
    for line in fin:
        token=line.rstrip().split()
        try:
            seqList.append(token[0])
            seqIdx.append(int(token[1]))
            pdbList.append(token[2])
            pdbIdx.append(int(token[3]))
            #seqtopdb[int(token[1])] = int(token[3])
            #seqtopdb_aa[int(token[1])] = token[2]

            seqtopdb[int(token[3])] = int(token[1])
            seqtopdb_aa[int(token[3])] = token[2]  
            
        except IndexError:continue
        

def readInputSeq(inF):
    
    fin = open(inF,"r")
    fout = open('map-seq-to-pdb-1.txt',"w")
    
    for line in fin:
        token=line.rstrip().split()
        #try:
        if True:
            idx1 = int(token[0])-1
            idx2 = int(token[1])-1
            #idx1 = int(token[0])
            #idx2 = int(token[1])  
            entr = token[2]

            #print(idx2+1)
            #print(pdbIdx[idx1+1],pdbIdx[idx2+1])
            
            #res1_res2_map[str(seqIdx[idx1])+"_"+str(seqIdx[idx2])] = entr

            #print(res1_res2_map)
            
            #fout.write("%s \t %s \t %s \t %s \t %s \t %s \t %s \n"%(idx1+1,idx2+1,entr,seqtopdb[idx1]+1,seqtopdb[idx2]+1,seqtopdb_aa[idx1],seqtopdb_aa[idx2]))
            fout.write("%s \t %s \t %s \t %s \t %s \n"%(idx1+1,idx2+1,entr,seqtopdb[idx1]+1,seqtopdb[idx2]+1))
            #fout.write("%s \t %s \t %s \t %s \t %s \t %s \t %s \n"%(idx1,idx2,entr,seqtopdb[idx1],seqtopdb[idx2],seqtopdb_aa[idx1],seqtopdb_aa[idx2]))

            #print(seqtopdb[idx1]+1)
            #print(seqtopdb[idx2]+1)
            #print(seqtopdb_aa[idx1])
            #print(seqtopdb_aa[idx2])
            
    fin.close()
    fout.close()
                     
    
    """    
    ###may be create a new funciton to do the task below
    ####for reading
    ###read the "all_list.txt and make a map as follows"               
    newmap={}
    newMap[res1+"_"+res2] = float(dist)
    
    #### to query dist###
    newDist = newMap[str(seqIdx[idx1])+"_"+str(seqIdx[idx2])] ##answer you look for
    
    
    res1_res2_map[str(seqIdx[idx1])+"_"+str(seqIdx[idx2])] = newDist
    """
                  
if __name__=="__main__": 
    
    #readSeqTOPDBList("map_seq_to_pdb-Dpr2DIPeta.txt")
    readSeqTOPDBList("map_pdb_to_seq.txt")
    readInputSeq("top-18-entr.txt")
    #readMatrix("../all_lst.txt")
