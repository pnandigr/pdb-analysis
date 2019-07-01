#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import os,sys
import numpy as np

def readInputSeq(inF):
    fin = open(inF,"r")
    nM = {}
    fout = open('check-entr.txt',"w")

    for line in fin:
        token=line.rstrip().split()
        try:
            idx1 = (token[0])
            idx2 = (token[1])
            entr = token[2]
            nM[idx1+"_"+idx2] = entr
        #except IndexError:continue
        
            if(int(idx2) >= 273):
                #if(float(entr) > 0.75):
                fout.write("%s \t %s \t %s \n"%(idx1,idx2,entr))
                #print(idx1,idx2,entr)
        except IndexError:continue
        
    fin.close()
    fout.close()

if __name__=="__main__":

    readInputSeq("all-cognates-positive-MI.model.txt")
