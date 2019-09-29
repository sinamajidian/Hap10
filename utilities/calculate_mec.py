
# compute the MEC score of two haplotypes

import numpy as np
from sys import argv

read_num_all=[]
    

def read_hapblock(hapblock_file): # hapcut format
    hap_dic={}
    block_num=0

    with open(hapblock_file, 'r') as hbf:
        for line in hbf:
            elements = line.strip().split('\t')
            if 'BLOCK' in line: # new block started
                block_num+=1
                hap1 ={} # ['-']*num_snps
                hap2 ={} # ['-']*num_snps
                read_num_all.append(int(line.split(' ')[-1][:-1]))
                
            elif '*' in line: # block finished. 
                hap_dic[block_num]=[hap1,hap2]
            else: # 
                pos = int(elements[0])-1
                hap1[pos] = elements[1]
                hap2[pos] = elements[2]
            
    return hap_dic

def read_frag(frag_file):
    frags=[]
    with open(frag_file, 'r') as ff:
        for line in ff:
            line_strp = line.strip()       
            frags.append(line)
    return frags

    

def calculate_mec(hapblock_file,frags):

    MEC_all=[]
    block_num=len(hap_dic)
    for block_num in range(1,block_num+1):                               
        [hap1,hap2]=hap_dic[block_num]
        MEC = 0
        for line_strp in frags:
            el = line_strp.split()            
            num_blks = int(el[0]) # blks of each fragment 
            assert len(el) >= 3 + 2*num_blks     # technically should have qual score, but for our purposes do not enforce
            hap1_EC  = 0 # counter
            hap2_EC  = 0
            for i in range(0, num_blks):
                offset   = int(el[2+2*i]) - 1 # convert to 0-index
                seq      = el[3+2*i]
                l        = len(seq)
                #hap1_ref = hap1[offset:(offset+l)]   # relevant slice of hap1
                hap1_ref = ['-']*l
                hap2_ref = ['-']*l   # relevant slice of hap2
                for i_l in range(offset,offset+l):
                    if i_l in hap1.keys():
                        hap1_ref[i_l-offset]=hap1[i_l]
                        hap2_ref[i_l-offset]=hap2[i_l]
                    # hap1_EC  = 0 # counter
                    # hap2_EC  = 0
                for x,y in zip(seq, hap1_ref):
                    if (x == '1' and y == '0') or (x == '0' and y == '1'):
                        hap1_EC += 1
                for x,y in zip(seq, hap2_ref):
                    if (x == '1' and y == '0') or (x == '0' and y == '1'):
                        hap2_EC += 1
            MEC += min(hap1_EC, hap2_EC)
        MEC_all.append(MEC)
    return MEC_all

if __name__ == "__main__":


    hapblock_file=   argv[1]  #'tst_hap.txt'# 'out1.hap'  # 
    frag_file= argv[2]  #'tst_frag.txt'  #'frag.txt' #
    hap_dic=read_hapblock(hapblock_file)
    frags=read_frag(frag_file)
    #print('reading finished, calculating MEC started')
    MEC_all=calculate_mec(hapblock_file,frags)
    print('Mean MEC ',np.mean(MEC_all) )
    print('Sum MEC ',np.sum(MEC_all) )
    
    block_num=len(hap_dic)
    mec_r=[]
    for i in range(block_num):
        mec_r.append(MEC_all[i]/read_num_all[i])
    m_mec_r=np.mean(mec_r)
    if m_mec_r>.001:
        print('MEC over Nread',np.round(m_mec_r,3)) # per block
    else:
        print('MEC over Nread',np.round(m_mec_r,6)) # per block
        
    

