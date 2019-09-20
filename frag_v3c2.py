#!/usr/bin/env python3



# warning the ouput may not be connected. use connected extractor another round

# the input is connected!

# in this version if within weakly spliting, become disconnected it will be printed
import networkx as nx
from numpy import linalg as LA
import numpy as np
#import matplotlib.pyplot as plt
import ast
from sys import argv
from scipy import sparse
from scipy import linalg
from scipy.sparse import diags


thresh=.031
    
def check_is_connected(W):
    read_graph=nx.convert_matrix.from_numpy_matrix(W, parallel_edges=False)
    check=nx.is_connected(read_graph)
    return check


def bipartition_graph(W):
    dic_2parts={}
#     print(W.todense())
    d1=W.sum(axis=0)
    d_half=np.power(d1,-.5)
    D_half=diags(np.array(d_half), [0])

    D=diags(np.array(d1), [0])
    mat_=D_half.dot(D-W)  
    mat=mat_.dot(D_half) # normalized laplacian
    dim_w=np.shape(W)[0]
    #print(' eigenvalues calculation is started')
    if dim_w < 10000:
        eig_vals, eig_vecs = sparse.linalg.eigsh(mat,k=2,which='SM')
    else:
        initial_rand=np.random.normal(size=(dim_w, 2)) # 2 is the number of eigen vlues  check the https://github.com/lmcinnes/umap/blob/master/umap/spectral.py
        eig_vals, eig_vecs = sparse.linalg.lobpcg(mat,initial_rand, largest=False, tol=1e-8)#sparse.linalg.eigsh(mat,k=2,which='SM', tol=1e-2,maxiter=10000)
    #print(' eigenvalues are calculated ')
    idx = np.argsort(eig_vals)
    second_eig=eig_vecs[:,idx[1]]
    #print(second_eig)
    #plt.plot(second_eig)
    #file_name_nct='second_eig_vector'+str(ii)+'.png'
    #plt.savefig(file_name_nct)
    #plt.close()
    #plt.show()    
    
    part1_indecis=np.where(second_eig>0)
    part2_indecis=np.where(second_eig<=0)
    part1_list_ind=list(part1_indecis[0])
    part2_list_ind=list(part2_indecis[0])
    #print('part 1  ',part1_list_ind)
    #print('part 2  ',part2_list_ind)
    W_part1_raw = W.tocsc()[:,part1_list_ind]
    W_part1 = W_part1_raw.tocsr()[part1_list_ind,:]

    W_part2_raw = W.tocsc()[:,part2_list_ind]
    W_part2 = W_part2_raw.tocsr()[part2_list_ind,:]
#     W_part1=W[part1_list_ind][:,part1_list_ind] #A[[0,1]] extracts first and secnond row of matrix A
#     W_part2=W[part2_list_ind][:,part2_list_ind]
    dic_2parts={'W_p1':W_part1,'W_p2':W_part2,'ind_p1':part1_list_ind,'ind_p2':part2_list_ind }
    return dic_2parts

def Ncut_graph(ind_A,ind_B,W): # eq 2 page 889 shi malik paper
    # W 2d numpy array
    W_sliced_raw = W.tocsc()[:,ind_B]
    W_sliced = W_sliced_raw.tocsr()[ind_A,:]
    
    #W_sliced=W[ind_A][:, ind_B]
    cut=W_sliced.sum()
    W_restA=W.tocsr()[ind_A,:] # extracting those rows
    W_restB=W.tocsr()[ind_B,:]
    assoc_A=W_restA.sum()
    assoc_B=W_restB.sum()
    Ncut_graph=cut/assoc_A+cut/assoc_B
    return Ncut_graph

def partitioning(fragments_connected_matrix,fragments_connected_indecis): #partitioning(W):
    """
    a
    """
    global ii
    list_branch_working=[fragments_connected_matrix]
    #list_branch_glob_indices_working=[list(range(np.shape(W)[0]))]
    list_branch_glob_indices_working=[fragments_connected_indecis]

    list_branch_finished=[]    
    list_branch_glob_indices_finished=[]
    
    while len(list_branch_working): 
        
        report_file.write('*******'+str(ii)+'*******\n')
        report_file.write('The working list (number of fragments) is'+str([np.shape(W1)[0] for W1 in list_branch_working])+"\n")
        report_file.write('The finsihed list (number of fragments) is '+str([np.shape(W1)[0] for W1 in list_branch_finished])+"\n") 
        report_file.flush()

        W=list_branch_working[0]
        glob_indices=list_branch_glob_indices_working[0]

        
        ###??if check_is_connected(W) & (np.shape(W)[0]>2) :
        if (np.shape(W)[0]>2) :
            dic_2parts= bipartition_graph(W) # dic_2parts['W_p1'] and dic_2parts['ind_p1']
            ncut=Ncut_graph(dic_2parts['ind_p1'],dic_2parts['ind_p2'],W)
            report_file.write('Ncut value between first (l='+str(np.shape(dic_2parts['W_p1'])[0])+') and second part (l='+str(np.shape(dic_2parts['W_p2'])[0])+') is '+str(np.round(ncut,6))+"\n")
            report_file.flush()
            if ((ncut>thresh) | (np.shape(dic_2parts['W_p1'])[0] ==0) | (np.shape(dic_2parts['W_p2'])[0]==0)):
                list_branch_finished.append(W)
                list_branch_glob_indices_finished.append(glob_indices)
                ii=ii+1
                filename_out= 'part'+str(ii)+'.txt' #'part'+str(len(list_branch_finished))+'.txt'
                with open(filename_out, 'w+') as file_out: 
                    for line_num in glob_indices:
                        file_out.write(str(line_num+1)+'\n')
                #print('it goes to finsihed')
            else:
                list_branch_working.append(dic_2parts['W_p1'])
                glob_indices_p1=[glob_indices[i] for i in dic_2parts['ind_p1']] 
                #print('glob ind 1 is',glob_indices_p1)
                list_branch_glob_indices_working.append(glob_indices_p1)
                list_branch_working.append(dic_2parts['W_p2'])
                glob_indices_p2=[glob_indices[i] for i in dic_2parts['ind_p2']] 
                list_branch_glob_indices_working.append(glob_indices_p2)
                #print('glob ind 2 is',glob_indices_p2)
                #print(list_branch_glob_indices_working)
        else:
            list_branch_finished.append(W)
            list_branch_glob_indices_finished.append(glob_indices)
            ii=ii+1
            #print(glob_indices)
            filename_out='part'+str(ii)+'.txt'#'part'+str(len(list_branch_finished))+'.txt'
            #print('sina \n \ n')
            with open(filename_out, 'w+') as file_out: 
                for line_num in glob_indices:
                    file_out.write(str(line_num+1)+'\n')
            file_out.close()
            #print('it is not connected and goes to finished')
        
        del list_branch_glob_indices_working[0]
        del list_branch_working[0]
              
    return list_branch_finished


def connected_components_extractor(fragment_graph):
    """
    input: networx graph
    output: a tuple of two lists: list of weight matrices of  connected comppentes
                                  list of indecis of connected comppentes,  
    """
    fragments_connected_indecis_list_set = sorted(nx.connected_components(fragment_graph), key=lambda x: min(x))
    # list of set of nodes [{0, 1, 2}, {3, 4}]
    fragments_connected_indecis_list = [] # list of list [[0, 1, 2], [3, 4]]
    fragments_connected_matrix_list = []
    for index, fragments_connected_indecis_set in enumerate(fragments_connected_indecis_list_set):
        fragments_connected_indecis_sorted = sorted(fragments_connected_indecis_set)
        fragments_connected_indecis_list.append(fragments_connected_indecis_sorted)
        fragments_connected_matrix_list.append(W_f[fragments_connected_indecis_sorted][:,fragments_connected_indecis_sorted]) #A[[0,1]] extracts first and secnond row of matrix A

    return (fragments_connected_matrix_list, fragments_connected_indecis_list)


#def snp_weight_matrix(list_fragment_dic):  # is not used ??

def fragment_weight_matrix(list_fragment_dic): 
    """
    input: a python list, each element is a dic {snp_indx:allele,..}
    output: weight matrix numpy, each node is the fragment.
    """
    number_fragments = len(list_fragment_dic)
    W_f=np.zeros((number_fragments,number_fragments))
    for i in range(number_fragments):  # if number_fragments=3 -> i=0,1,2
        fragment_i = list_fragment_dic[i]  # fragment_i is dictinary
        W_f[i,i]=len(fragment_i)
        for j in range(i):  # i=2, range(i)=0,1 is not 
            fragment_j = list_fragment_dic[j]
            num_shared_snp = len(set(fragment_i.keys())&set(fragment_j.keys()))
            #shared_snp_allel_num= len(set(fragment_i.items()) & set(fragment_j.items()))
            #if shared_snp_num>0: W[i,j]=(2*shared_snp_allel_num-num_shared_snp)/num_shared_snp else: W[i,j]=0
            #if W[i,j]<0:W[i,j]=0
            W_f[i,j]=num_shared_snp
            W_f[j,i]=W_f[i,j]
    return W_f


def parse_fragment_file(fragment_file_name):  
    """
    input: fragment file with the format of haptree dic as text file. 
    output: a python list, each element is a dic {snp_indx:allele,..}
    """
    list_fragment_dic = [] # output list
    with open(fragment_file_name, 'r') as file_fragment: # each line of fragment file is a fragment
            for fragment in file_fragment: 
                fragment_dic = ast.literal_eval(fragment)
                list_fragment_dic.append(fragment_dic)
    return list_fragment_dic


if __name__ == "__main__":
    """
    input: fragment file with the format of haptree dic as text file. It can be generated by
        python bamprocess.py file.bam file.vcf
        mv Hap...txt frag.txt
        python2 FragmentPoly.py -f frag.txt  -o frag_dic.txt -x HapTree
    """
    
    ii=0 # iteration for ouput parts
    fragment_file_name =   argv[1]  # 'small_exmp2.txt'  'old_data/small_exmp2_weight_f.npy'
    list_frag_dic = []
    report_file_name="report_.txt"#"+fragment_file_name[:-4]+"
    report_file = open(report_file_name, "a") 
    if fragment_file_name[-3:]=='txt':
        list_fragment_dic = parse_fragment_file(fragment_file_name) # list_fragment_dic[2][5] shows the allele of snp_idx=5 of 2nd fragment
        report_file.write('Number of fragments in the txt file is'+str(len(list_fragment_dic))+"\n")
        W_f_np = fragment_weight_matrix(list_fragment_dic)
        W_f=sparse.coo_matrix(W_f_np,dtype=np.int8)
        filename_numpy_matirx = fragment_file_name[:-4]+'_weight_f'
        #np.save(filename_numpy_matirx, W_f)
        sparse.save_npz('frag_weight.npz', W_f)
    if fragment_file_name[-3:]=='npy':
        W_f_np=np.load(fragment_file_name)
        W_f=sparse.coo_matrix(W_f_np,dtype=np.int8)
        sparse.save_npz('frag_weight.npz', W_f)
    if fragment_file_name[-3:]=='npz':
        W_f = sparse.load_npz(fragment_file_name)
        
    
    report_file.write('Number of fragments in the npy file is'+str(np.shape(W_f)[0])+"\n")
    
     
    
    
    #fragment_graph = nx.convert_matrix.from_numpy_matrix(W_f, parallel_edges=False)
    #connected = nx.is_connected(fragment_graph)
        
#     if not connected:
#         fragments_connected_matrix_list, fragments_connected_indecis_list =connected_components_extractor(fragment_graph)

#     if connected:
#         fragments_connected_matrix_list = [W_f]
#         fragments_connected_indecis_list = [list(range(np.shape(W_f)[0]))]
     
#     #print(fragments_connected_matrix_list)
#     #print(fragments_connected_indecis_list)    
#     report_file.write('Number of connected componnets in the file is '+str(len(fragments_connected_indecis_list))+"\n")
#     report_file.flush()
    
  
     
#     for connected_i in range(len(fragments_connected_matrix_list)):
        
#         report_file.write('==='+str(connected_i)+'== \n')
#         report_file.flush()
#         fragments_connected_matrix=fragments_connected_matrix_list[connected_i]
#         fragments_connected_indecis=fragments_connected_indecis_list[connected_i]
#         #print(fragments_connected_indecis)
#         
#         #print('number of partitions is ',len(a))
    
#     print(fragment_file_name+" is finished")
    #return check 
    #a= partitioning(W)
    #print('>>><<<')
    
    fragments_connected_matrix = W_f
    fragments_connected_indecis = list(range(np.shape(W_f)[0]))
    a= partitioning(fragments_connected_matrix,fragments_connected_indecis)

    
