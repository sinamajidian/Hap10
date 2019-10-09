#!/usr/bin/env python3


#import matplotlib.pyplot as plt

from sys import argv
import networkx as nx
import numpy as np
from scipy import sparse
from scipy import linalg




def parse_fragment_file(fragment_file_name):

    """
    input: text file as the fragment file with the format of hapcut
    output: two python lists,
                 First list:  each element is a list of snp pos (line number in the vcf file) correspond to each line of input fragment matrix
                 Second list: each element is a string which is the line of input fragment matrix

                 Note that the alleles are not saved.
    """

    list_varids_infragments = []                        # varid is the line number of each variant in the VCF file (withouht counting #) as provided in fragment file
    list_lines = []                                     # all lines of the input fragment matrix

    with open(fragment_file_name, 'r') as fragment_file: # credit https://github.com/EhsanMotazedi/Haplosim/blob/master/FragmentPoly.py
                for line in fragment_file:                   
                    line_strip = line.strip()            #  each line of fragment file is a fragment
                    list_lines.append(line_strip)
                    columns = line_strip.split()         # ' ' or "\t"
                    col_num = len(columns)-1
                    if ((col_num >= 4) and (col_num%2 == 0)):
                        list_start_varids = columns[2:col_num:2]        # start var_id for all segments in the fragment
                        list_alleles_segments = columns[3:col_num:2]    # alleles for all segments in the fragment

                        # 33 01 42 000
                        #list_start_var_id=['33','42']

                        list_varids_infragment=[]           # list of all var_id for all segments in the fragment
                        for i in range(0,len(list_start_varids)):

                            start_varids = int(list_start_varids[i])
                            alleles_segment = list_alleles_segments[i]
                            list_varids_infragment += list(range(start_varids,start_varids+len(alleles_segment)))

                        list_varids_infragments.append(list_varids_infragment)
                    else:
                        print('Error in format of input fragment file ')
                        exit(1)


    order_line_number=sorted(range(len(list_varids_infragments)),key=lambda k: list_varids_infragments[k]) # sorting list of dic based on the first pos of each line
    list_varids_infragments_sorted=[list_varids_infragments[i] for i in order_line_number]
    list_lines_sorted=[list_lines[i] for i in order_line_number]


    return list_varids_infragments, list_lines_sorted



def extract_ccs(list_varids_infragments):

    """ extracting connected componenets

    input:  A list of lists. Inner list contains the variant indices (varids) of all alleles in the all segments in one fragment (line of fragment matrix)
    output: A list of lists. Inner list contains the line number of fragment in the fragment file for each cc  [[2], [6], [7, 5, 4, 3, 1, 0]]

    """

    list_ccs_varids=[]
    list_ccs_fragids=[] # list_fragment_pos list_varid_fragment
    
    
    for fragid, list_varids_infragment in enumerate(list_varids_infragments):   # each iteration is for one fragment in list_varids_infragments 
         
        # fragid the 0-based index of this fragment in the list_varid_fragments
        # list_varid_fragment is a list of numbers (varids) 
        set_varids_infragment = set(list_varids_infragment)

        # obtaining overlap of current fragment with existing (cc)s and reporting the id of those cc. (the cc_id is chanign in each outer for loop)
            
        cc_id_overlapped_list=[]    
        for cc_id, cc in enumerate(list_ccs_varids):                          # finding the index of cc that are overlapped with the new fragment
            intersction = set_varids_infragment & cc
            if len(intersction):
                cc_id_overlapped_list.append(cc_id)

        cc_varids_new=set_varids_infragment                                  # new_cc initialized with the indices of the new fragment
        cc_fragids_new=[fragid]
        for cc_id in cc_id_overlapped_list:                                   # adding those cc that are overlapped with the fragment to the new cc
            cc_varids_new |= list_ccs_varids[cc_id]
            cc_fragids_new  += list_ccs_fragids[cc_id]

        for cc_id in sorted(cc_id_overlapped_list, reverse=True):             # removing those cc that are overlapped from the list_cc
            del list_ccs_varids[cc_id]
            del list_ccs_fragids[cc_id]

        list_ccs_varids.append(cc_varids_new)
        list_ccs_fragids.append(cc_fragids_new)

    return list_ccs_fragids






def report_frag_file(list_fragment_lines, list_parts_fragid, fragment_file_name_prefx,tresh_size_report):
    """

    write several fragment files correspond to each (strongly) connected components.


    """

    #print('Number of connected parts is '+str(len(list_cc_fragid))+'\n')

    for part_id, part_fragid in enumerate(list_parts_fragid):

        
        if len(part_fragid) > tresh_size_report:
            report_file.write("Number of fragments in part with id "+str(part_id)+" is "+str(len(part_fragid))+"\n")

            file_part_id = open(fragment_file_name_prefx+'_'+str(part_id)+'.txt', "w")

            for fragid_k in part_fragid:
                line_fragid_k=list_fragment_lines[fragid_k]
                file_part_id.write(line_fragid_k+"\n")

            file_part_id.close()
            
        else:
            report_file.write("Number of fragments in part with id "+str(part_id)+" is "+str(len(part_fragid))+". So, it's not reported as output file.\n")

    return 1



#W_f_cc = calculate_fragment_weight_matrix(list_varids_infragments_cc)

def calculate_fragment_weight_matrix(list_varids_infragments_cc):

    """
    input: a python list, each element is a dic {snp_indx:allele,..}
    output: weight matrix (scipy sparse matrix, lil_matrix), each node is the fragment. The weight is the number of shared SNPs between two fragments.

    """

    number_fragments = len(list_varids_infragments_cc)
    report_file.write("Calculating weight matrix of "+str(number_fragments)+" fragments.\n")
    W_f = sparse.lil_matrix((number_fragments,number_fragments), dtype=np.int16)

    for i in range(number_fragments):      # if number_fragments=3 -> i=0,1,2

        varids_infragment_i = list_varids_infragments_cc[i]
        W_f[i,i] = len(varids_infragment_i)

        for j in range(i):                 # i=2, range(i)=0,1
            varids_infragment_j = list_varids_infragments_cc[j]
            num_shared_vars = len( set(varids_infragment_i) & set(varids_infragment_j) )

            W_f[i,j] = num_shared_vars
            W_f[j,i] = W_f[i,j]

    return W_f




def extract_strongly_cc(W_f_cc):

    # To implement the spectral clustering, consider a tree. at each step it grows and creat two brances. A simplay way to traverse the tree.

    list_node_matrix_working = [W_f_cc]
    list_node_indices_working = [list(range(np.shape(W_f_cc)[0]))]  # here the values of indices is based on the input matrix
    #list_node_matrix_finished = []     # we don't need this
    list_node_indices_finished = []



    while len(list_node_matrix_working):

        cut_is_not_needed = False

        W_f = list_node_matrix_working[0]
        indices = list_node_indices_working[0]



        dic_2parts = bipartition_graph(W_f)  # dic_2parts['W_p1'] and dic_2parts['ind_p1']

        if dic_2parts == -1:                # the eiganvlue decompostion in bipartition_graph function faces error, we can only say cut is not needed.
            cut_is_not_needed =True

        else:
            ind_part1 = dic_2parts['ind_p1']
            ind_part2 = dic_2parts['ind_p2']
            ncut_value = calcualte_ncut_value(ind_part1,ind_part2,W_f)

            if (np.shape(dic_2parts['W_p1'])[0] < tresh_size_part) | (np.shape(dic_2parts['W_p2'])[0] < tresh_size_part):
                cut_is_not_needed = True

            if ncut_value > thresh_ncut:
                cut_is_not_needed = True

            report_file.write('Ncut value between first (l='+str(np.shape(dic_2parts['W_p1'])[0])+') and second part (l='+str(np.shape(dic_2parts['W_p2'])[0])+') is '+
                              str(np.round(ncut_value,6))+", therefore, ") #(cut is needed)="+str(not cut_is_not_needed)+"\n
            report_file.flush()


        # remove the current from working and then

        del list_node_matrix_working[0]
        del list_node_indices_working[0]



        if cut_is_not_needed:            # move the current to the finished list  (the node is the ending node of the tree)
            report_file.write('Cut is not needed.\n')
            #list_node_matrix_finished.append(W)
            list_node_indices_finished.append(indices)

        else:                            # put the new two smaller parts to working list (as two new nodes)
            report_file.write('Cut is needed.\n')
            W_part1 = dic_2parts['W_p1']
            W_part2 = dic_2parts['W_p2']

            list_node_matrix_working.append(W_part1)
            list_node_matrix_working.append(W_part2)


            ind_part1_local = dic_2parts['ind_p1']
            ind_part2_local = dic_2parts['ind_p2']

            ind_part1 = [indices[i] for i in ind_part1_local]
            ind_part2 = [indices[i] for i in ind_part2_local]

            list_node_indices_working.append(ind_part1)
            list_node_indices_working.append(ind_part2)
            
    return list_node_indices_finished




def bipartition_graph(W_f):
    """

    input: weight matrix of graph
    output: a dictionary contains two weight matrices (corresponding to two parts of graph after ncut)

    Based on paper "Normalized Cuts and Image Segmentation"

    """


    dim_W=np.shape(W_f)[0]
    report_file.write("Currently on bipartition graph with the size of "+str(dim_W)+"\n")
    diag_data=W_f.sum(axis=0)

    I = sparse.identity(dim_W, dtype=np.float32) # credit https://github.com/lmcinnes/umap/blob/master/umap/spectral.py
    D = sparse.spdiags(1.0 / np.sqrt(diag_data), 0, dim_W, dim_W)

    # Normalized Laplacian
    normalized_laplacian = I - D * W_f * D  # equals to D^(-0.5)*(D-W)*D^(-0.5)

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

        thresh_bi_ncut = 0  # in future it can be optimzaized, maybe mean of vector
        part1_indecis = list(np.where(second_eig > thresh_bi_ncut)[0])
        part2_indecis = list(np.where(second_eig <= thresh_bi_ncut)[0])


        W_part1_allrows = W_f.tocsc()[:,part1_indecis]
        W_part1 = W_part1_allrows.tocsr()[part1_indecis,:]

        W_part2_allrows = W_f.tocsc()[:,part2_indecis]
        W_part2 = W_part2_allrows.tocsr()[part2_indecis,:]

        dic_2parts={'W_p1':W_part1,'W_p2':W_part2,'ind_p1':part1_indecis,'ind_p2':part2_indecis }


        return dic_2parts

    except sparse.linalg.eigen.arpack.arpack.ArpackNoConvergence or sparse.linalg.ArpackError:
        warn( "The eigenvector solver failed. The output may not strangly connected\n ")

        return -1




def calcualte_ncut_value(ind_A,ind_B,W_f):

    """

    input: two list of node indices  ind_A, ind_B and  weight matrix of graph W_f
    output: the value of ncut between two input sets of node graph (row/column of matrix W_f)

    Based on paper "Normalized Cuts and Image Segmentation"  page 889, equation 2

    """

    W_sliced_raw = W_f.tocsc()[:,ind_B]
    W_sliced = W_sliced_raw.tocsr()[ind_A,:]

    cut = W_sliced.sum()
    W_restA = W_f.tocsr()[ind_A,:] # extracting those rows
    W_restB = W_f.tocsr()[ind_B,:]
    assoc_A = W_restA.sum()
    assoc_B = W_restB.sum()
    if assoc_A and assoc_B:
        ncut_value=cut/assoc_A+cut/assoc_B
    else:
        ncut_value=1000      # cut is not needed. The output will be the same as input

    return ncut_value







if __name__ == "__main__":

    """

    input:  fragment file with the format of hapcut. 
            Here, each line is called a fragment. Each fragment consists of several segments.  
    
            2 aa 33 01 42 000
           
            Number of segments=2
            fragment id= aa
            segment 1: start_var_id=33, alleles=01  -> 33:0 34:1
            segment 2: start_var_id=42, alleles=000 - > 42:0 43:0 44:0
           
            var_id is the index of variant in the vcf file (line number in the vcf file after removing #)

    output:  several fragment files which are (strongly) connected components  
    

    run:    python3  extract_scc.py fragment.txt cc  ./out # connected componnets
            python3  extract_scc.py fragment.txt scc ./out # strongly connected componnets

    """
 
    print('You can check the report file. It is updating when the code is running.')

   
    tresh_size_part=100       # If a part id smaller than that, we won't split it.
    tresh_size_report=10     # If a part id smaller than that, we won't report it in outputs.  **
     
    
    thresh_ncut=.031
    #ii=0 # iteration for ouput parts

    out_mode = argv[2]
    input_fragment_file_name = argv[1] #'data_1m01_2_c10/frag_sp_weight_cc1.npz' #'data/frag_sp_0.txt'  #'data/frag_sp_0_weight_cc0.npz'# #argv[1]  # 'data/frag_weight.npz'#

    output_files_prefx = argv[3] #'data_1m01_2_c10/out' 
    report_file_name = output_files_prefx+"_report.txt"
    report_file = open(report_file_name, "a")

    list_varids_infragments, list_lines = parse_fragment_file(input_fragment_file_name)               # list_fragment_dic[2][5] shows the allele of snp_idx=5 of 2nd fragment

    #list_fragment_pos, list_fragment_lines = parse_fragment_file(fragment_file_name)        # list_fragment_dic[2][5] shows the allele of snp_idx=5 of 2nd fragment
    report_file.write("Number of fragments in the txt file is "+str(len(list_lines))+"\n")


    if out_mode == "cc": # All lines of output fragment files are overlapped. Ncut is not used.
        
        list_ccs_fragids =  extract_ccs(list_varids_infragments)

        report_file.write("The input file contains "+str(len(list_ccs_fragids))+" connected component(s).\n")
        report_frag_file(list_lines, list_ccs_fragids, output_files_prefx,tresh_size_report)                   # it may contain only one cc
    
    
    
    if out_mode == "scc": # Extracting strongly connected components

        list_ccs_fragids = extract_ccs(list_varids_infragments) 
        report_file.write("The input file contains "+str(len(list_ccs_fragids))+" connected component(s).\n")
    
         # extracting varids of all fragments in each cc from list_varids_infragments (one list for input file) and export it as list_varids_infragments_cc.
        list_varids_infragments_ccs=[]
        for cc_id, list_cc_fragids in enumerate(list_ccs_fragids):   # list_cc_fragids is the list of line number which are connected to each other.
            
            list_varids_infragments_cc=[]
            for fragid in list_cc_fragids:
                list_varids_infragments_cc.append(list_varids_infragments[fragid])
                
            list_varids_infragments_ccs.append(list_varids_infragments_cc)
            
            
        for cc_id, list_varids_infragments_cc in enumerate(list_varids_infragments_ccs): 
            
            report_file.write("*** Working on connected component index "+str(cc_id)+"***\n")
            list_cc_fragids=list_ccs_fragids[cc_id]   # it is used for exporting output fragment file

            if len(list_varids_infragments_cc) > tresh_size_part:    
                
                #report_file.write("Number of fragments in the connected components with index "+str(cc_id)+" is "+str(len(list_fragment_pos_cc_id))+"\n")
                W_f_cc = calculate_fragment_weight_matrix(list_varids_infragments_cc)
                #sparse.save_npz(fragment_file_name[:-4]+"_weight_cc"+str(cc_id)+".npz", sparse.coo_matrix(W_f_cc,dtype=np.int16))

                list_strong_parts=extract_strongly_cc(W_f_cc) # list of list
                #print("cc with id"+str(cc_id)+"contains "+str(len(list_strong_parts))+"stronlgy connected components.")
                # in a rare case, the output of this step was not connected components. Check for next version. 
                
                # extract their global line number
                list_parts_fragids=[]
                for part_id, strong_part in enumerate(list_strong_parts): # strng_part contains the indices of cc_fragid
                    strong_fragids=[list_cc_fragids[i] for i in strong_part]
                    list_parts_fragids.append(strong_fragids)
                   
            else:
                report_file.write("The cc is very small, so it is considered as strongly cc.\n")
                list_parts_fragids=[list_cc_fragids]
            
        
            # then report them as several files
            report_file.write("The connected component with index "+str(cc_id)+", contains "+str(len(list_parts_fragids))+" strongly connected component(s)\n.")
            fragment_file_out_prefx=output_files_prefx+str(cc_id)
            report_frag_file(list_lines, list_parts_fragids, fragment_file_out_prefx,tresh_size_report)
            

    report_file.close()
 
