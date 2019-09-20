#!/usr/bin/env python3

"""
This code is for reading a fragment file of 10x data end split those rows of the sme barcode but came from different molecules

usage: python3

Sina Majdiain
Start: 20 Jan

"""
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
from sys import argv
import numpy as np
from sklearn.cluster import MeanShift


def parse_pos_file(pos_free):
    """ parsing a position file

        input: a file containing positions  of variants from freebayes
        output: a list of varinats positions

        """
    #  Parsing the file containing the  positions  of variants
    position_list = []  # initializing the list of variants positions
    with open(pos_free, 'r') as infile:
        for line in infile:  # each row of the file is a variant position
            position_list.append(int(line.strip()))

    return position_list


def parse_frag_file(frag_file):
    """ parsing a fragment file

    input: a fragment file
    output: a dictionary, key: barcode  values: [list of parts, quality]], each part is (65,'1100') which 65 is the starting point


    consider  22  1100 -> it means taht 22:1 23:1 24:0 25:0

    The fragment file is the output of bamprocess.py of E. Motazedi.
    The starting point is the index of variant not the real position  of reference genome
    """

    #  Parsing the file containing fragments
    dic_frag = {}
    with open(frag_file, 'r') as infile:
        for line in infile:
            ls = line.strip().split('\t')  # ls: separate all parts of the row based on tab.
            barcode = ls[1][4:-2]  # the barcode of the row
            start_point_list = [int(i) for i in ls[::2][1:-1]]  # list of all starting points of all rows. ls[::2][0] is number of parts, ls[::2][end] is quality
            allele_part_list = ls[3::2]  # list of all alleles of all parts of the row. Each element of list is string '100'

            list_parts = []  # list of
            for idx_part in range(len(start_point_list)):
                part = (start_point_list[idx_part], allele_part_list[idx_part])  # tuple( )
                list_parts.append(part)
            quality = ls[-1]
            dic_frag[barcode] = [list_parts, quality]

    return dic_frag


def convert_indx_position(list_parts, position_list):
    """ finding the positions of each part

    input: a fragment file and list of position
    output:  genomic positions of each part,  [[100,110,112],[..], ...]

    """
    list_parts = list_parts_quality[0]  # a list of tuples, each tuple is a part like (65,'1100') start index and alleles
    pos_list = []
    for idx_prt, tuple_part in enumerate(list_parts):  # tuple_part = (64, '1100')
        start_idx = tuple_part[0]-1  # -1 is to make as index starting form zero.
        alleles_part = tuple_part[1]  # '1100'
        pos_part = list()
        for i_allele in range(len(alleles_part)):
            allele_idx = start_idx + i_allele
            allele_pos = position_list[allele_idx]
            pos_part.append(allele_pos)
        pos_list.append(pos_part) # save the list of each part for each barcode. each part contains list of positions
    return pos_list


def cluster_parts(pos_list, mol_length):
    """ clustering each barcode based on the list of start_end

    input: a  list of list genomic positions [[100,110,112],[..], ...]
    output:  list of list of indices of parts for each cluster [[0,1,2],[3,4,5]]

    split those rows with the same barcode but came from different molecules

    """
    # for each barcode, we have a few clusters.
    #clusters_parts = {}  # each cluster corresponds to a molecule  (10x)
    # Here, the index of parts are clustered based on the same molecules

    posistions = []
    for part in pos_list:
        posistions += part
    positions_np = np.array(posistions)
    positions_np_reshaped = positions_np.reshape(-1, 1)
    
    
    
    bandwidth = mol_length/2
    #clustering = MeanShift(bandwidth=bandwidth).fit(positions_np_reshaped) 
    
    clustering = MeanShift(bandwidth=bandwidth,bin_seeding=True).fit(positions_np_reshaped) 
    
    #if enough memory and run on one dataset, -1 means using all processors. 
    #clustering = MeanShift(bandwidth=bandwidth,bin_seeding=True n_jobs= -1).fit(positions_np_reshaped) #for better time use ,

    # cluster labels are not in order, the first one may be 1, next 0 next 2!
    cluster_labels = clustering.labels_ 
    #cluster_num = cluster_labels[-1]+1  # number of clusters

    cluster_list = []  # for each cluster find the start/end genomic position

    clst_unq, clst_indx = np.unique(cluster_labels, return_index=True)
    clst_unq_order_preserved = clst_unq[np.argsort(clst_indx)]


    for cluster_indx  in clst_unq_order_preserved:
        idx_pos_cluster = np.where(cluster_labels == cluster_indx)
        pos_this_cluster = positions_np[idx_pos_cluster]
        start_pos_cluster = pos_this_cluster[0]
        end_pos_cluster = pos_this_cluster[-1]
        cluster_list.append([start_pos_cluster, end_pos_cluster])

    # Now, split the position list to clusters based on the start/end in following list
    cluster_idx = 0
    [start_pos_clust, end_pos_clust] = cluster_list[cluster_idx]
    clusters_parts = [[]] # [[0,1,2,3],[4,5,6]] for all clusters of a barcode
    for idx_part, part in enumerate(pos_list):

        score_part = 0  # number of allele in a position that belong to this cluster
        for each_pos in part:
            if each_pos >= start_pos_clust and each_pos <= end_pos_clust:
                score_part += 1

        if score_part >= len(part)/2:  # if most allele in a part belong to that cluster, add the index of this part to list
            clusters_parts[cluster_idx].append(idx_part)
        else:  # else consider the next cluster
            cluster_idx += 1
            [start_pos_clust, end_pos_clust] = cluster_list[cluster_idx]
            clusters_parts.append([idx_part])  # adding new inner list (correspond to new cluster) to list

    return clusters_parts


def split_frg(idx_parts, list_parts_quality):
    """ splitting a fragment dictionary

    input: a fragment dictionary and a dic_idx_parts={brcd:{0:{1,2,3}, 1:{4,5,6},...  },..}
    output: a dictionary, key: barcode  values:

    split those rows with the same barcode but came from different molecules

    """
    dic_splited={}
    clusters_parts = idx_parts# dic_idx_parts[barcode]  # for each barcode, parts are clustered. the list of indecis of parts,
    quality_all = list_parts_quality[1]  # quality of each barcode
    parts = list_parts_quality[0]
    accumulated_length_old = 0
    for idx_cluster, cluster_indx_part in enumerate(clusters_parts):
        parts_cluster = [parts[i] for i in cluster_indx_part]
        length_this_cluster = sum([len(part[1]) for part in parts_cluster])
        quality_cluster = quality_all[accumulated_length_old:accumulated_length_old+length_this_cluster]
        accumulated_length_old += length_this_cluster
        dic_splited[idx_cluster] = [parts_cluster, quality_cluster]

    return dic_splited


def write_new_frag_file(barcode_cluster,list_parts_quality, new_fragment_file):
    """ writing a file

    input: a  dictionary of splited
    output: a dictionary, key: barcode  values:

    """
    
    alleles_parts = list_parts_quality[0]
    quality = list_parts_quality[1]
    if len(quality)>1:  # it can be improved by adding the single allele to previous part ??
        new_fragment_file.write("{} BX::{}".format(len(alleles_parts), barcode_cluster))
        for idx, part in enumerate(alleles_parts):
            new_fragment_file.write(" {} {}".format(part[0], part[1]))
        new_fragment_file.write(" {}\n".format(quality))
    return 1


if __name__ == "__main__":

    #address = argv[1]
    #filename_fragment = address+argv[2]
    #pos_freebys = address+argv[3]
    filename_fragment = argv[1] #'frag1a.txt'
    pos_freebys = argv[2] #'pos_freebayes.txt'
    molecule_length = int(argv[3])* 1000   # the molecule length


    pos_list = parse_pos_file(pos_freebys)
    dic_frg = parse_frag_file(filename_fragment)
    
    out_file_name='frag_sp.txt' #filename_fragment[:-4]+
    new_fragment_file=open(out_file_name, "w")

    for barcode, list_parts_quality in dic_frg.items():
        list_parts=list_parts_quality[0]
        pos_parts = convert_indx_position(list_parts, pos_list)
        idx_parts = cluster_parts(pos_parts, molecule_length)
        dic_splited = split_frg(idx_parts, list_parts_quality)
        for idx_cluster, list_parts_quality in dic_splited.items():  
            barcode_cluster=barcode + '_' + str(idx_cluster)
            write_new_frag_file(barcode_cluster,list_parts_quality, new_fragment_file)
