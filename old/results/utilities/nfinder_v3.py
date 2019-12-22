#!/usr/bin/env python3

"""
Script to find region of non-N sequence of a FASTA file

Sina Majidian

Attention this version only works with on record per fasta file

"""



import sys

from Bio import SeqIO
from Bio.Seq import Seq

def nonNregion_extractor(seq):
    
    list_tuple_pos=[]
    for idx, char in enumerate(seq):
        if idx==0:       
            if char == 'N': 
                    start_pos=1 
        else:
            previous_char=seq[idx-1]
            if char == 'N': 
                if previous_char != 'N': # new N region started 
                    start_pos=idx+1 #  0-index to genomic pos
                    #else:  char =previous_char = 'N' 
            else: # char = 'A'
                if previous_char == 'N':
                    lastN_idx=idx-1
                    end_pos=lastN_idx+1
                    if end_pos-start_pos >13000 :
                        list_tuple_pos.append((start_pos,end_pos))
                    
   
    list_tuple_pos_nonN=[]
    if len(list_tuple_pos)>0:
        for idx_nregion, tuple_pos in enumerate(list_tuple_pos):
            previous_tuple_pos=list_tuple_pos[idx_nregion-1]
            list_tuple_pos_nonN.append((previous_tuple_pos[1],tuple_pos[0]))
        list_tuple_pos_nonN.append((tuple_pos[1],len(seq))) # last tuple
        list_tuple_pos_nonN[0]=(1,list_tuple_pos[0][0]) # first tuple should be edited  by this line
    else:
        list_tuple_pos_nonN=[(1,len(seq))]

    return list_tuple_pos_nonN

if __name__ == "__main__":

    file_fa= sys.argv[1]  # 'tst.fa'
    with open(file_fa, mode="r") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            seq = str(record.seq)
    
    list_tuple_pos_nonN=nonNregion_extractor(seq)
    print('The number of region around 50kbs is '+str(len(list_tuple_pos_nonN)))
    file_bed = open("nonNregion.bed","w") 
    for tuple_pos_nonN in list_tuple_pos_nonN:
        file_bed.write(str(tuple_pos_nonN[0])+'\t'+str(tuple_pos_nonN[1])+'\n')
    file_bed.close()

    
