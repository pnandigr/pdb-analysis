#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import os,sys
import numpy as np

pdbList = []
pdbIdx = []
seqList = []
seqIdx = []
resA = []
resB = []
res1_res2_map = {}


def readPDBTOSeqList(inF):
    fin = open(inF, "r")
    for line in fin:
        token = line.rstrip().split()
        try:
            pdbList.append(token[0])
            pdbIdx.append(int(token[1]))
            seqList.append(token[2])
            seqIdx.append(token[3])
        except IndexError: continue


def readInputSeq(inF):
    #    queryRes1=[]
    #    queryIdx1=[]
    #    queryRes2=[]
    #    queryIdx2=[]
    fin = open(inF, "r")
    fout = open('map-pdb-to-seq.txt',"w")
    
    for line in fin:
        token=line.rstrip().split()
        try:
            idx1 = int(token[1])-1
            # idx1 = int(token[1])
            res1 = str(token[2])     
            idx2 = int(token[4])-1
            # idx2 = int(token[4])
            res2 = str(token[5])   
            dist = token[7]
            res1_res2_map[str(seqIdx[idx1])+"_"+str(seqIdx[idx2])] = dist
            
            # print ("idx 1= "+str(idx1)+ " is "+str(seqIdx[idx1])+" res="+res1)
            # print ("idx 2= "+str(idx2)+ " is "+str(seqIdx[idx2])+" res="+res2)
            fout.write("%s \t %s \t %s \t %s \t %s \t %s \n" %
                       (idx1+1, idx2+1, seqIdx[idx1], seqIdx[idx2], res1, res2))
            
#            queryIdx1.append(int(token[1]))
#            queryRes1.append(token[2])
#            queryIdx2.append(int(token[4]))
#            queryRes2.append(token[5])
        except IndexError:continue
    fin.close()
    fout.close()
    

def readMatrix(infile1):
    in_f1 = open(infile1,"r")
    # in_f2 = open(infile2,"r")
    newMap = {}
    fout = open("entropy_msa_indx.txt","w")

    #for line in in_f2:
    #    token1=line.rstrip().split()
    #    try:
    #        resA = str(token1[2])     
    #        resB = str(token1[5])
    #    except IndexError:continue
    #print(resA,resB)
        
    for line in in_f1:
        token = line.rstrip().split()
        try:
            idx1 = token[0]
            idx2 = token[1]
            entr = token[2]
            newMap[idx1+"_"+idx2] = float(entr)
        except IndexError:continue

    for r in res1_res2_map:
        newDist = newMap[r]
        r1, r2 = r.split("_")[0], r.split("_")[1]


        fout.write("%s \t %s \t %s \n"%(r1,r2,newDist))
        # fout.write("%s \t %s \t %s \t %s \t %s \n"%(r1,r2,newDist,resA,resB))
        
    in_f1.close()
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


if __name__ =="__main__":
    
    readPDBTOSeqList("map_pdb_to_seq.txt")
    readInputSeq("contact-all_0.dat")
    readMatrix("../all_cognates-positive-MI.model.txt")

