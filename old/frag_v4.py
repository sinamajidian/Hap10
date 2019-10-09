#!/usr/bin/env python3

import networkx as nx
from numpy import linalg as LA
import numpy as np
#import matplotlib.pyplot as plt
import ast
from sys import argv

import scipy
from scipy import sparse
from scipy import linalg
from scipy.sparse import diags
from scipy.sparse import lil_matrix



def parse_fragment_file(fragment_file_name):  
    
    """
    input: text file as the fragment file with the format of sdhap  (assumption: The lines are sorted based on the position of first allele)
    output: two python lists, 
                 First list:  each element is a list of snp pos (line number in the vcf file) correspond to each line of input fragment matrix  
                 Second list: each element is a line of input fragment matrix 
                 
                 Note that the alleles are not saved.
    """

    list_fragment_pos = []                              # pos is the line number of each variant in the VCF file (withouht counting #) as provided in fragment file
    list_fragment_lines = []  
    with open(fragment_file_name, 'r') as file_fragment: # credit https://github.com/EhsanMotazedi/Haplosim/blob/master/FragmentPoly.py    
            for line in file_fragment:                   #  each line of fragment file is a fragment 
                line_strip = line.strip()
                list_fragment_lines.append(line_strip)
                columns = line_strip.split()         # ' ' or "\t"     
                l_col = len(columns)-1
                if ((l_col >= 4) and (l_col%2 == 0)):
                    reads_start = columns[2:l_col:2]    # start pos for each segments in the line 
                    hetro_sites = columns[3:l_col:2]    # alleles for all segments in the line 
                    list_fragment_pos_line=[]           # list of all start pos for all segments
                    for i in range(0,len(reads_start)): 
                        start_pos = int(reads_start[i])
                        list_fragment_pos_line+=list(range(start_pos,start_pos+len(hetro_sites[i])))

                list_fragment_pos.append(list_fragment_pos_line)
                
                # if you need the alleles you may extract them as dic for each line
                    #hetro_dict = {}
#                     for i in range(0,len(reads_start)):
#                         start_pos = int(reads_start[i])
#                         for j in range(0,len(hetro_sites[i])):
#                             hetro_dict[start_pos-1] = int(hetro_sites[i][j])
#                             start_pos += 1
#                 list_fragment_dic.append(hetro_dict)
#     order_id=sorted(range(len(list_fragment_dic)),key=lambda k: min(list(list_fragment_dic[k].keys()))) # sorting list of dic based on the first pos of each line
#     list_fragment_dic=[list_fragment_dic[i] for i in order_id]

    # if not sorted
    order_id=sorted(range(len(list_fragment_pos)),key=lambda k: list_fragment_pos[k]) # sorting list of dic based on the first pos of each line
    list_fragment_pos=[list_fragment_pos[i] for i in order_id]
    list_fragment_lines=[list_fragment_lines[i] for i in order_id]
    
    
    return list_fragment_pos, list_fragment_lines



def extract_cc_id(list_fragment_pos):
    
    """ extract connected componenets 
    
    input:  A list of lists. Inner list contains the SNP indices of all alleles in the all segments in one fragment (line of fragment matrix)
    output: A list of lists. Inner list contains the line number of fragment in the fragment file for each cc  [[2], [6], [7, 5, 4, 3, 1, 0]]
    
    """
    
    list_cc_fragpos=[]
    list_cc_fragid=[]
    for fragid_k, list_fragment_pos_k in enumerate(list_fragment_pos): # the fragment are sorted in the list
        # fragid_k the 0-based index of this fragment in the fragment file 
        indices = set(list_fragment_pos_k)

        cc_id_overlapped_list=[]
        for cc_id, cc in enumerate(list_cc_fragpos):    # finding the index of cc that are overlapped with the new fragment
            intersction = indices & cc    
            if len(intersction):
                cc_id_overlapped_list.append(cc_id)

        cc_fragpos_new=indices                           # new_cc initialized with the indices of the new fragment
        cc_fragid_new=[fragid_k]
        for cc_id in cc_id_overlapped_list:             # adding those cc that are overlapped with the fragment to the new cc
            cc_fragpos_new |= list_cc_fragpos[cc_id] 
            cc_fragid_new  += list_cc_fragid[cc_id]

        for cc_id in sorted(cc_id_overlapped_list, reverse=True):      # removing those cc that are overlapped from the list_cc
            del list_cc_fragpos[cc_id]
            del list_cc_fragid[cc_id]

        list_cc_fragpos.append(cc_fragpos_new)
        list_cc_fragid.append(cc_fragid_new)

    return list_cc_fragid

        
def report_cc(list_fragment_lines, list_cc_fragid, fragment_file_name_prefx):
    """
    
    write several fragment files correspond to each connected components
    
    
    """
     
    print('Number of connected parts is '+str(len(list_cc_fragid))+'\n')

    for cc_id, cc_fragid in enumerate(list_cc_fragid):
        
        report_file.write("Number of fragments in the connected components "+str(cc_id)+" is "+str(len(cc_fragid))+"\n")

        cc_file_cc_id = open(fragment_file_name_prefx+'_'+str(cc_id)+'.txt', "w") 

        for fragid_k in cc_fragid:    
            line_fragid_k=list_fragment_lines[fragid_k]
            cc_file_cc_id.write(line_fragid_k+"\n")

        cc_file_cc_id.close()
        
    return 1
    



def calculate_fragment_weight_matrix(list_fragment_pos): 
    
    """
    input: a python list, each element is a dic {snp_indx:allele,..}
    output: weight matrix (scipy sparse matrix, lil_matrix), each node is the fragment. The weight is the number of shared SNPs between two fragments.
    
    """
    
    number_fragments = len(list_fragment_pos)
    W_f = lil_matrix((number_fragments,number_fragments), dtype=np.int16) 
    
    for i in range(number_fragments):      # if number_fragments=3 -> i=0,1,2
        
        fragment_i_pos = list_fragment_pos[i]  
        W_f[i,i] = len(fragment_i_pos)
        
        for j in range(i):                 # i=2, range(i)=0,1 
            fragment_j_pos = list_fragment_pos[j]
            num_shared_snp = len( set(fragment_i_pos) & set(fragment_j_pos) )
            
            #shared_snp_allel_num= len(set(fragment_i.items()) & set(fragment_j.items()))
            #if shared_snp_num>0: W[i,j]=(2*shared_snp_allel_num-num_shared_snp)/num_shared_snp else: W[i,j]=0
            #if W[i,j]<0:W[i,j]=0
            
            W_f[i,j] = num_shared_snp
            W_f[j,i] = W_f[i,j]
            
    return W_f


def bipartition_graph(W):
    """
    
    input: weight matrix of graph
    output: a dictionary contains two weight matrices (corresponding to two parts of graph after ncut)
    
    Based on paper "Normalized Cuts and Image Segmentation"
    
    """
    
    dim_W=np.shape(W)[0]
    diag_data=W.sum(axis=0)
    
    I = sparse.identity(dim_W, dtype=np.float32) # credit https://github.com/lmcinnes/umap/blob/master/umap/spectral.py
    D = sparse.spdiags(1.0 / np.sqrt(diag_data), 0, dim_W, dim_W)

    # Normalized Laplacian
    normalized_laplacian = I - D * W * D  # equals to D^(-0.5)*(D-W)*D^(-0.5)
    
    #print(' eigenvalues calculation is started')
    try:
        if dim_W < 8000:

            #num_lanczos_vectors = int(np.sqrt(dim_W) # #ncv=num_lanczos_vectors, 
            eig_vals, eig_vecs = sparse.linalg.eigsh(normalized_laplacian, k=2, which='SM',tol=1e-3, v0=np.ones(dim_W), maxiter=dim_W*2) 

        else:

            initial_rand=np.random.normal(size=(dim_W, 2)) # 2 is the number of eigen vlues  
            eig_vals, eig_vecs = sparse.linalg.lobpcg(normalized_laplacian, initial_rand, largest=False, tol=1e-3,maxiter=9000)
        #print('eigenvalues calculation is finished')
        idx = np.argsort(eig_vals)    
        second_eig=eig_vecs[:,idx[1]]

        #plt.plot(second_eig)
        #file_name_nct='second_eig_vector'+str(ii)+'.png'
        #plt.savefig(file_name_nct)
        #plt.close()
        #plt.show()  
        
        tresh_ncut = 0  # in future it can be optimzaized, maybe mean of vector 
        part1_indecis = list(np.where(second_eig > tresh_ncut)[0]) 
        part2_indecis = list(np.where(second_eig <= tresh_ncut)[0]) 

        
        W_part1_allrows = W.tocsc()[:,part1_indecis]
        W_part1 = W_part1_allrows.tocsr()[part1_indecis,:]

        W_part2_allrows = W.tocsc()[:,part2_indecis]
        W_part2 = W_part2_allrows.tocsr()[part2_indecis,:]

        dic_2parts={'W_p1':W_part1,'W_p2':W_part2,'ind_p1':part1_indecis,'ind_p2':part2_indecis }

        
        return dic_2parts
        
    except scipy.sparse.linalg.eigen.arpack.arpack.ArpackNoConvergence or scipy.sparse.linalg.ArpackError:
        warn( "The eigenvector solver failed. The output may not strangly connected\n ")
        
        return -1
        
  
    
    
def calcualte_ncut_value(ind_A,ind_B,W): 
    
    """
    
    input: two list of node indices  ind_A, ind_B and  weight matrix of graph W
    output: the value of ncut between two input sets of node graph (row/column of matrix W)
    
    Based on paper "Normalized Cuts and Image Segmentation"  page 889, equation 2
    
    """
    
    W_sliced_raw = W.tocsc()[:,ind_B]
    W_sliced = W_sliced_raw.tocsr()[ind_A,:]
    
    cut = W_sliced.sum()
    W_restA = W.tocsr()[ind_A,:] # extracting those rows
    W_restB = W.tocsr()[ind_B,:]
    assoc_A = W_restA.sum()
    assoc_B = W_restB.sum()
    if assoc_A and assoc_B:
        ncut_value=cut/assoc_A+cut/assoc_B 
    else:
        ncut_value=1000      # cut is not needed. The output will be the same as input
        
    return ncut_value







if __name__ == "__main__":
    
    """
    
    input: fragment file with the format of haptree dic as text file. It can be generated by
            python bamprocess.py file.bam file.vcf
            mv Hap...txt frag.txt
            python2 FragmentPoly.py -f frag.txt  -o frag_dic.txt -x HapTree
    
    output:    
        
        
    """
    
    
    ii=0 # iteration for ouput parts
    fragment_file_name = 'data/frag_sp_0_weight_cc0.npz'#'data/frag_sp_0.txt'  #argv[1]  # 'data/frag_weight.npz'#
#     list_frag_dic = []
    

    report_file = open("report_ncut_values.txt", "a") 
    
    if fragment_file_name[-3:] == 'txt':
        
        list_fragment_pos, list_fragment_lines = parse_fragment_file(fragment_file_name)        # list_fragment_dic[2][5] shows the allele of snp_idx=5 of 2nd fragment
        report_file.write("Number of fragments in the txt file is "+str(len(list_fragment_lines))+"\n")

 
    
    out_mode = "strong"
    
    if out_mode == "cc":
        list_cc_fragid= extract_cc_id(list_fragment_pos)
        if len(list_cc_fragid)>1:
            report_cc(list_fragment_lines, list_cc_fragid, fragment_file_name[:-4])
        else:
            print("The input file contains only one connected components. No output file is generated.")

        
   
    if out_mode == "strong":
        
        if fragment_file_name[-3:] == "npz":  # one connected components npz
            W_f = sparse.load_npz(fragment_file_name)
            report_file.write("Number of fragments in the npz file is "+str(np.shape(W_f)[0])+"\n")
            
        if fragment_file_name[-3:]=="txt":
            
            list_cc_fragid= extract_cc_id(list_fragment_pos)
            
            print("There are "+str(len(list_cc_fragid))+" connected components.\n") 
            for cc_id, cc_fragid in enumerate(list_cc_fragid):  
                print("Working on connected components "+str(cc_id)+"\n") 
                list_fragment_pos_cc_id=[]
                for fragid_k in cc_fragid: 
                    list_fragment_pos_cc_id.append(list_fragment_pos[fragid_k])
                
                report_file.write("Number of fragments in the connected components "+str(cc_id)+" is "+str(len(list_fragment_pos_cc_id))+"\n")
                W_f = calculate_fragment_weight_matrix(list_fragment_pos_cc_id) 
                sparse.save_npz(fragment_file_name[:-4]+"_weight_cc"+str(cc_id)+".npz", sparse.coo_matrix(W_f,dtype=np.int16))
        
    report_file.close()
                                                   
                                   

            
                
        