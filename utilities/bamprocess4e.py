# Make the SNP-fragment list from a multi-sample barcoded 10x-Genomics bam file and its corresponding VCF file with input options: mmq, mbq and qoffset.
# Written by Ehsan Motazedi, Wageningen UR, 04-11-2016.
# Last mofified: 18-07-2018.

import copy
import difflib
import os
import pysam
import re
import subprocess
import sys
from sys import argv
from genotypes import getAllelesPop, getGenotypesPop
from math import log
from reads import Read
# from memory_profiler import profile
import numba as nb
import numpy as np




@nb.njit
def numba_search(arr, element):
    for i in range(arr.shape[0]):
        if arr[i] == element:
            return i
    return -1



# @profile
def adjust_seq(read):
	""" Adjust the read sequence according to the mapped positions, i.e. get rid of the insertions and clipped bases."""
	read_cigarstring=read.cigarstring
	read_seq=read.seq
	read_qual=read.qual
	if read_cigarstring=='151M' or read_cigarstring=='128M' :
		adj_seq=read_seq
		adj_qual=read_qual
	else:
		cig = list(_x for _x in re.split('[0-9]{1,}', read_cigarstring) if _x)
		cign = list(_x for _x in re.split('[^0-9]', read_cigarstring) if _x)
		cig = list(cig[_n]*int(cign[_n]) for _n in range(0, len(cig)))
		cig = ''.join(_x for _x in cig if 'D' not in _x) # deleted nucleotides from the reference are not present in the read sequence 
		adj_seq = []
		adj_qual = []
		for _n, _x in enumerate(read_seq):
			if cig[_n]=='M': # Allow only match/mismatch of the read nucleotides, i.e. no clipping, no insertion.
				adj_qual.append(read_qual[_n])
				adj_seq.append(_x)
		adj_seq=''.join(adj_seq)
	return adj_seq, adj_qual

#@profile
def frag_gen(varpos, allelelst, genolst, coordinates, nucleotides, qscores):
	""" Generate SNP-fragments from (linked) reads.""" 
	if (not coordinates) or (not nucleotides):
		return Read(), {}
	var_codes = []
	var_num = []
	var_q = []
	#print('length coordinates is',len(coordinates))
	#print('length varpos is ',len(varpos))
	#print('10th element of varpos is ',varpos[10])
	coordinates_set=set(coordinates)
	coordinates_np=np.array(coordinates, dtype=np.int64)
	for _n, _x in enumerate(varpos):
		if _x < coordinates[0]:
			continue
		if _x > coordinates[-1]:
			break
		if _x in coordinates_set:  # _x in coordinates:
			if set(genolst[_n].GetGenes()).intersection(set(['.','-'])): # Throw away missing genotypes or genotypes with one or more missing allele(s)
				continue
			if len(set(genolst[_n].GetGenes()))<2: # Do not include in the SNP-fragments belonging to a population member its homozygous alleles
				continue
			try:
				coordinates_index_x=numba_search(coordinates_np, _x)
				#coordinates_index_x=coordinates.index(_x)
				var_codes.append(allelelst[_n][2][nucleotides[coordinates_index_x]])
				var_num.append(allelelst[_n][0])
				var_q.append(qscores[coordinates_index_x])
			except (ValueError,KeyError): # if the called nucleotide is wrong, i.e. does not exist in VCF alleles or _x not in coordinates
				pass
	try: # return the reads {SNP number: allele} and the associated quality scores {SNP number: Qscore}
		return Read({_x:str(_y) for _x, _y in zip([var_num[0], var_num[1]]+var_num[2:], var_codes)}), {_x:str(_y) for _x, _y in zip(var_num, var_q)}
	except IndexError: # throw away reads with less than 2 SNPs
		return Read(), {}

class InputError(Exception):
	""" Handle invalid input specifications.""" 
	def __init__(self, msg, *args):
		super(InputError, self).__init__(args)
		self.msg = msg
	def __str__(self):
		return "InputError: {}\n".format(self.msg)
	def __repr__(self):
		return (self.msg,)+self.args+('\n',)

#@profile
def get_frags(bamfile, vcffile, mbq=13, mmq=20, qoffset=33):
	""" 
mmq : minimum read mapping quality to consider a read for phasing, default 20\n
qvoffset <33/64> : quality value offset, 33/64 depending on how quality values were encoded, default is 33\n
mbq  : minimum base quality to consider a base for haplotype fragment, default 13\n
	"""
	try:
		all_reads = pysam.Samfile(bamfile, 'rb')
	except IOError:
		raise InputError('The input BAM file was not found!')
	ReadHeader = subprocess.Popen(["samtools","view","-H", bamfile], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
	header, err_header = ReadHeader.communicate()
	if ReadHeader.returncode!=0:
		raise InputError('Failed to read the header from the bam file! Original error message:\n'+err_header)
	if isinstance(header, bytes):
		header=bytes.decode(header)
	else:
		pass
	RGIDs, SMnames = [], []
	for _headerline in header.splitlines(): # pasre the header of the bam file to extract the ID and SM fields of each Read Group
		if _headerline[0:3]=='@RG':
			RGID_added, SM_added = False, False
			for _n, _field in enumerate(_headerline.split()):
				if 'ID' in _field:
					if not RGID_added:
						RGIDs.append(''.join(_headerline.split()[_n].split(':')[1:])) # add read group ID
						RGID_added = True
					else:
						raise InputError('Double ID fields detected in @RG header line!')
				elif 'SM' in _field:
					if not SM_added:
						SMnames.append(''.join(_headerline.split()[_n].split(':')[1:])) # add the sample name associated with the read group ID
						SM_added = True
					else:
						raise InputError('Double SM fields detected in @RG header line!')
			if SM_added and RGID_added:
				pass
			elif SM_added:
				raise InputError('ID field missing in @RG header line!')
			elif RGID_added:
				raise InputError('SM field missing in @RG header line!')
			else:
				raise InputError('ID and SM fields missing in @RG header line!')
	if len(RGIDs)!=len(set(RGIDs)):
		raise InputError('Duplicate read group IDs detected in the bam header!')
	GroupedReadsWithID = [[] for _id in RGIDs] # separate reads belonging to each Read Group
	for _read in all_reads:
		try:
			GroupedReadsWithID[RGIDs.index(dict(_read.get_tags())['RG'])].append(_read)
		except ValueError: # this error occurs if a read group in the read tags is not present in the header. This is NOT expected to happen though!
			SMnames.append(SMnames[RGIDs.index(difflib.get_close_matches(dict(_read.get_tags())['RG'],RGIDs)[0])]) # Choose a sample name for the new read group from the header read groups the most close to it 
			RGIDs.append(dict(_read.get_tags())['RG'])
			GroupedReadsWithID.append([_read])
	GroupedReads, GroupedSM = [], [] # combine the reads with different RGID but the same SM as they are assumed to belong to the same sample
	for _SMname, _ReadGroup in zip(SMnames, GroupedReadsWithID):
		if _SMname not in GroupedSM:
			GroupedReads.append(_ReadGroup)
			GroupedSM.append(_SMname)
		else:
			GroupedReads[GroupedSM.index(_SMname)]+=_ReadGroup
	del GroupedReadsWithID
	try:
		genolst = getGenotypesPop(vcffile, GroupedSM)
		allelelst = getAllelesPop(vcffile, GroupedSM)
	except IOError:
		raise InputError('The VCF file was not found!')
	except:
		raise
	varpos = list(_x[1] for _x in allelelst)
	for _group, reads in enumerate(GroupedReads):
		out=open('HapCUT_frags_{}.matrix'.format(GroupedSM[_group]),'w')
		BX_reads = []
		Chromium_barcode = []
		for _r in reads:
			try:
				Chromium_barcode.append(_r.get_tag('BX'))
				BX_reads.append(_r)
			except KeyError:
				pass
		#garbage = sys.stderr.write('WARNING: {0:d} reads were not considered as they lacked an adjusted chromium barcode!\n'.format(len(reads)-len(BX_reads)))
		reads = BX_reads
		linked_reads = {_barcode:[] for _barcode in Chromium_barcode}
		frag_lst = {_barcode:[] for _barcode in Chromium_barcode} 
		q_lst = {_barcode:[] for _barcode in Chromium_barcode}
		HapCUT =  {_barcode:[] for _barcode in Chromium_barcode}
		for _rnum, _barcode in enumerate(Chromium_barcode): # scan through the alignments to find the linked reads)
			read = copy.deepcopy(reads[_rnum])
			if read.is_unmapped or read.is_duplicate or (read.flag & 0x900): # throw away unmapped reads, duplicates and secondary/supplementary alignments
				continue
			if read.mapping_quality >= mmq:
				try:
					adj_seq, adj_qual =  adjust_seq(read)
					coordinates, nucleotides, quals = list(zip(*[(int(_x), _y, _z) for _x, _y, _z in zip(read.positions, adj_seq.upper(), list(ord(_x)-qoffset for _x in adj_qual)) if _z>=mbq]))
				except ValueError as e:
					if e.args[0][0:len("need more than 0 values to unpack")]=="need more than 0 values to unpack" or e.args[0][0:len("not enough values to unpack")]=="not enough values to unpack":
						coordinates, nucleotides, quals= [(), (), ()] 
					else:
						raise
				linked_reads[_barcode].append((coordinates, nucleotides, quals))
			#break # for testing only one read 
		
		
		
		
		#print('len ', len(linked_reads))
		#print('first ', linked_reads[0])
		#print('len  set', len(set(linked_reads)))
		barcode_indx=0; number_barcodes=len(linked_reads);print('number of barcodes is',number_barcodes)
		for _barcode, _rec in linked_reads.items(): 
			try: 
				barcode_indx+=1;
				#if barcode_indx%10 == 0:	print(''+str(barcode_indx)+' barcodes out of '+str(number_barcodes),'are considered')
				_rec = [tuple(zip(*sorted(zip(*_r), key = lambda _x: int(_x[0])))) for _r in _rec if len(_r[0])>0]
				_rec = sorted(_rec, key = lambda x: x[0][0])
			except IndexError as e:
				if "tuple index out of range" in e.args[0]:
					_rec = [(), (), ()] 
				else:
					raise
			else: # Take care of overlaps between paired-end reads or between reads in barcoded molecules
				_rn = 0
				while _rn < len(_rec)-1:
					if _rec[_rn][0][-1]<_rec[_rn+1][0][0]:
						_rn+=1
					else:
						coordinates = [_c for _c in _rec[_rn][0]]+[_c for _c in _rec[_rn+1][0]]
						nucleotides = [_n for _n in _rec[_rn][1]]+[_n for _n in _rec[_rn+1][1]]
						quals = [_q for _q in _rec[_rn][2]]+[_q for _q in _rec[_rn+1][2]]
						coordinates, nucleotides, quals = [_t for _t in zip(*sorted(zip(coordinates,nucleotides,quals), key=lambda x: int(x[0])))]
						unique_q = []
						unique_c = []
						unique_n = []
						for _n, _c in enumerate(coordinates): # remove the duplicates from overlapping positions 
							try:
								if unique_c[-1]!=_c:
									#print("\nN:\n",nucleotides, "\nQ:\n", quals)
									unique_c.append(_c)
									unique_n.append(nucleotides[_n])
									unique_q.append(quals[_n])
								elif unique_n[-1]==nucleotides[_n]:
									unique_q[-1] = min(126-qoffset, unique_q[-1]+quals[_n])
								else: # if the called nucleotides differ at overlapping sites, use the one with the highest Phred score and adjust the Phred score.
									if quals[_n]>unique_q[-1]:
										_new_q_score = round(-10*log(1-10**(-unique_q[-1]/10)*(1-10**(-quals[_n]/10)), 10), 5) # Q=-10log(p,10)
										if _new_q_score >= mbq:
											unique_n[-1] = nucleotides[_n]
											unique_q[-1] = _new_q_score
										else:
											del(unique_c[-1], unique_n[-1], unique_q[-1])
									else:
										_new_q_score = round(-10*log(1-(1-10**(-unique_q[-1]/10))*10**(-quals[_n]/10), 10), 5)
										if _new_q_score >= mbq:
											unique_q[-1] = _new_q_score
										else:
											del(unique_c[-1], unique_n[-1], unique_q[-1])
							except IndexError:
								unique_c.append(_c)
								unique_n.append(nucleotides[_n])
								unique_q.append(quals[_n])
						_rec[_rn+1] = (tuple(unique_c), tuple(unique_n), tuple(unique_q))
						_rec[_rn] = ((),(),())
						_rn+=1
				_rec = [tuple(zip(*sorted(zip(*_r), key = lambda _x: int(_x[0])))) for _r in _rec if len(_r[0])>0]
				coordinates = [int(_c)+1 for _r in _rec for _c in _r[0]] # Convert the zero-based BAM coordinates to 1-based, as the coordinates are 1-based in the VCF (like the SAM format).
				nucleotides = [_n for _r in _rec for _n in _r[1]]
				quals = [_q for _r in _rec for _q in _r[2]]
				frag_lst_barcode, q_lst_barcode = frag_gen(varpos, allelelst, genolst[_group], coordinates, nucleotides, quals)
				_HapCUT_alleles = [[]]
				_HapCUT_starts = []
				_snp_alleles = frag_lst_barcode.GetalleleSet(frag_lst_barcode.GetBegin(), frag_lst_barcode.GetEnd())
				if len([_x for _x in _snp_alleles if _x!='-'])>1: #keep only SNP-fragments that have a length of at least 2
					for _n, _allele in enumerate(_snp_alleles): # prepare records from the fragment list and the quality list to be written in HapCUT format
						if _allele == '-' and _HapCUT_alleles[-1]!=[]:
							_HapCUT_alleles.append([])
						elif _allele != '-' and _HapCUT_alleles[-1]==[]:
							_HapCUT_alleles[-1].append(_allele)
							_HapCUT_starts.append(str(_n+frag_lst_barcode.GetBegin()+1))
						elif _allele != '-':
							_HapCUT_alleles[-1].append(_allele)
					_HapCUT_alleles = _HapCUT_alleles[0:len(_HapCUT_starts)] # eliminate empty allele sets
					for _n in range(0, len(_HapCUT_starts)):
						_HapCUT_alleles[_n] = ''.join(_HapCUT_alleles[_n])
					HapCUT_barcode = str(len(_HapCUT_starts))+'\t'+'BX::'+_barcode+'\t'+'\t'.join('\t'.join([_start, _alleles]) for _start, _alleles in zip(_HapCUT_starts, _HapCUT_alleles))+'\t'+''.join(chr(int(_q)+qoffset) for _q in q_lst_barcode.values())
					out.write(HapCUT_barcode+os.linesep)
		#break # for testing only one barcode
		out.close()

if __name__ == "__main__":
	
	bamfile=argv[1]
	vcffile=argv[2]
	#get_frags('../20180327/cultivar1_C1R5_BXsorted.bam','../20180327/cultivar1_C1R5_recoded.vcf', mmq=9, mbq=12)
	get_frags(bamfile,vcffile)#, mmq=9, mbq=12)


