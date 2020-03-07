#!/usr/bin/env python3

"""
This code is for reading a fragment file of 10x data end split those rows of the same barcode but came from different molecules.

Step 3 of Hap10 paper. Extracting strongly connected components of fragments



usage: python3 splitter.py frag_linked.txt var_het.vcf 50

Sina Majidian
Start: 20 Jan 2019

"""


from sys import argv
import numpy as np
from sklearn.cluster import MeanShift



def parse_frag_file(frag_file):
    """ parsing a fragment file

    input: a fragment file
    output: a dictionary, key: barcode  values: [list of segments, qualities]], each segment is (65,'1100') which 65 is the starting point


    22  1100 -> it means that 22:1 23:1 24:0 25:0

    The fragment file is the output ofextractHairs.
    The starting point is the index of variant not the real position  of reference genome.
    """

    dic_frag = {}
    with open(frag_file, 'r') as infile:
        for line in infile:

            ls = line.strip().split()                            # ls: separate all parts of the row based on tab.
            barcode_raw = ls[1]

            if barcode_raw.startswith('B'):
                barcode = barcode_raw[4:-2]                      # the barcode of the row
            else:
                barcode = barcode_raw[:-2]

            start_point_list = [int(i) for i in ls[::2][1:-1]]   # list of all starting points of all rows. ls[::2][0] is number of parts, ls[::2][end] is qualities
            allele_segment_list = ls[3::2]                       # list of all alleles of all parts of the row. Each element of list is string '100'

            list_segments = []  # list of
            for idx_segment in range(len(start_point_list)):

                segment = (start_point_list[idx_segment], allele_segment_list[idx_segment])  # tuple( )
                list_segments.append(segment)
            qualities = ls[-1]

            fragment=[list_segments, qualities]
            dic_frag[barcode] = fragment

    return dic_frag




def parse_vcf_file(filename_vcf):
    """ parsing a vcf file

        input: a vcf file containing positions of variants output of variant caller in it's second column.
        output: a list of variants positions

        """


    vcf_file = open(filename_vcf,'r')
    list_positions=[]

    for line in vcf_file:
        line_strip=line.strip()

        if line_strip.startswith('#'):
            pass
        else:

            line_parts=line_strip.split('\t')
            list_positions.append(int(line_parts[1]))

    return list_positions




def convert_varid_position(list_segments, position_list):
    """ finding the positions of each segment

    input:   a list of segments [(6, '00')]  and list of variants positions
    output:  genomic positions of each segment,  [[100,110,112],[..], ...]

    """

    #list_segments = list_segments_quality[0]  # a list of tuples, each tuple is a segment like (65,'1100') start index and alleles
    list_pos_segments = []
    for tuple_segment in list_segments:  # tuple_segment = (64, '1100')

        start_id = tuple_segment[0]-1  # -1 is to make 0-based indexing.
        alleles_segment = tuple_segment[1]  # '1100'

        pos_segment = list()
        for i_allele in range(len(alleles_segment)):

            var_id = start_id + i_allele
            var_pos = position_list[var_id]
            pos_segment.append(var_pos)
        list_pos_segments.append(pos_segment) # save the list of each segment for each barcode. each segment contains list of positions
    return list_pos_segments





def cluster_segments(list_pos_segments, molecule_length):
    """
    split each fragment with the same barcode but came from different molecules
    clustering each list of positions of barcode-based fragment to molecule-based list of positions
    for each barcode, we will have a few clusters.
    the indices of segments are clustered based on the same molecules

    input: a  list of list genomic positions [[100,110,112],[..], ...]
    output:  list of list of indices of segments for each cluster [[0,1,2],[3,4,5]]


    """

    posistions = []
    for segment in list_pos_segments:   # combine all segments to one list
        posistions += segment
    positions_np = np.array(posistions)
    positions_np_reshaped = positions_np.reshape(-1, 1)


    bandwidth = molecule_length/2
    try:
        clustering = MeanShift(bandwidth=bandwidth,bin_seeding=True).fit(positions_np_reshaped)  # n_jobs= -1 using all cores.

    except MemoryError:  # how to handle this?
        log_excepttion(MemoryError)
        exit(1)

    clusters_labels=clustering.labels_
    # cluster labels are not in order, the first one may be 1, next 0 next 2!  # array([1, 1, 0, 0, 0, 2, 2])
    cluster_unique, cluster_id = np.unique(clusters_labels, return_index=True)
    cluster_unique_order_preserved = cluster_unique[np.argsort(cluster_id)]

    # for each cluster, find the boundries (start/end genomic position)
    list_boundries_clusters = []
    for cluster_id  in cluster_unique_order_preserved:
        indices_this_cluster = np.where(clusters_labels == cluster_id)  # indices of alleles in the combined list which are estimated to be in this cluster
        pos_this_cluster = positions_np[indices_this_cluster]
        start_pos_this_cluster = pos_this_cluster[0]
        end_pos_this_cluster = pos_this_cluster[-1]
        list_boundries_clusters.append([start_pos_this_cluster, end_pos_this_cluster])

    # Now, split the position list to clusters based on the start/end in list_boundries_clusters
    cluster_id = 0
    [start_pos_this_cluster, end_pos_this_cluster] = list_boundries_clusters[cluster_id]

    indices_clusters_segments = [[]] # [[0,1,2,3],[4,5,6]] for all clusters of a barcode
    for idx_segment, segment in enumerate(list_pos_segments):

        score_segment = 0  # number of alleles that are in within the boundry belonging to this cluster
        for each_pos in segment:
            if each_pos >= start_pos_this_cluster and each_pos <= end_pos_this_cluster:
                score_segment += 1

        if score_segment >= len(segment)/2:  # if most allele in a segment belong to that cluster, add the idx of this segment to this cluster
            indices_clusters_segments[cluster_id].append(idx_segment)

        # Since, values in the list are ordered, we do not need to check all clusters
        else:  #  consider the next cluster.
            cluster_id += 1
            [start_pos_this_cluster, end_pos_this_cluster] = list_boundries_clusters[cluster_id]
            indices_clusters_segments.append([idx_segment])  # adding new inner list (correspond to new cluster) to list

    return indices_clusters_segments



def split_frg(indices_clusters_segments, list_segments, qualities):
    """ splitting a fragment based on the clustered indecis

    input: a fragment dictionary and a list of indices for each cluster={[1,2,3], [4,5,6],...  ]
    output: a list of list, inner list for each cluster  [[[(1, '00')], 'BB'], [[(32, '010')], 'AAA']]

    """
    list_fragment_splitted=[]

    accumulated_length = 0
    for cluster_id, indices_cluster_segments in enumerate(indices_clusters_segments):

        segments_cluster = [list_segments[i] for i in indices_cluster_segments]

        length_this_cluster = sum([len(segment[1]) for segment in segments_cluster])

        quality_cluster = qualities[accumulated_length:accumulated_length+length_this_cluster]
        accumulated_length += length_this_cluster

        list_fragment_splitted.append([segments_cluster, quality_cluster])

    return list_fragment_splitted

def write_new_frag_file(molecule_name,fragment_splitted , out_file):
    """ writing a new fragment file

    input:
    output: a file

    """

    alleles_segments = fragment_splitted[0]
    quality = fragment_splitted[1]
    if len(quality)>1:  # it can be improved by adding the single allele to previous segment

        out_file.write("{} {}".format(len(alleles_segments), molecule_name))

        for segment in alleles_segments:
            out_file.write(" {} {}".format(segment[0], segment[1]))
        out_file.write(" {}\n".format(quality))
    return 1



if __name__ == "__main__":


    filename_fragment = argv[1] #'data/frag_out.txt'#argv[1]      # Barcode-specific fragment file
    filename_vcf= argv[2] #'data/a.vcf' #argv[2]                  #  VCf file, we need second column, genomic position of variant for clustering
    molecule_length = int(argv[3])* 1000 # 50000           # length of 10x DNA molecule

    dic_frg = parse_frag_file(filename_fragment)
    list_positions=parse_vcf_file(filename_vcf)


    filename_out= filename_fragment[:-4]+'_sp.txt'
    out_file=open(filename_out, "w")




    for barcode, fragment in dic_frg.items():
        [list_segments, qualities]=fragment

        # example: list_segments=[(1, '00'),(30, '010')]; qualities='BBAAA'

        list_pos_segments = convert_varid_position(list_segments, list_positions)

        indices_clusters_segments = cluster_segments(list_pos_segments, molecule_length)

        list_fragment_splitted = split_frg(indices_clusters_segments, list_segments,qualities)

        for idx_cluster, fragment_splitted in enumerate(list_fragment_splitted):
            molecule_name=barcode + '_' + str(idx_cluster)
            write_new_frag_file(molecule_name, fragment_splitted, out_file)

    out_file.close()
