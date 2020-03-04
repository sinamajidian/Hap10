#!/usr/bin/env python3


"""
    This code has to modes
    
    1) haplogen:
    converting output of haplogenerator (input.hap) to a fake VCF
    
    python hap2vcf.py input.hap haplogen
    
    2) sdhap_genomic:    
    Combine outpout of convertalleSDHap (input.hap) with a real VCF
    
    python hap2vcf.py input.hap input.vcf
    
""" 
    
    



from sys import argv

def read_haplo_gen(haplos_file_address):
    
    """ reading the output of haplogenerator
    
    input name of file 
    output list
    """
    haplos_file = open(haplos_file_address,'r'); 
    var_list=[]
    for line in haplos_file:
        line_strip=line.strip() 
        if line_strip.startswith('B'):
            header_line=line_strip
        else:
            line_split=line_strip.split('\t')
            chrom= line_split[1]
            genomic_position=int(line_split[2])
            ref_allele=line_split[3]
            alt_allele=line_split[4]
            hap_values=line_split[5:]

            var_list.append([chrom,genomic_position,ref_allele,alt_allele,hap_values])

    return var_list


def writ_vcf_gen(vcf_file_address,out_unphased):

    
    """ writing a vcf file 
    
    input name of output file 
    output 1
    """
    
    out_unphased=True # if you need unphased,
    vcf_file = open(vcf_file_address,'w'); 

    vcf_file.write('##fileformat=VCFv4.2\n##source=Haplogenerator\n')
    vcf_file.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="test">\n')
    vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')


    for var in var_list:
        [chrom,genomic_position,ref_allele,alt_allele,hap_values]=var
        ID='.'
        quality=str(100)
        filtr='PASS'
        formt='GT'
        info='NS=1'
        if out_unphased:
            hap='/'.join(hap_values)
        else:
            hap='|'.join(hap_values)

        var_out=[chrom,str(genomic_position),ID,ref_allele,alt_allele,quality,filtr,info,formt,hap]
        vcf_file.write('\t'.join(var_out)+'\n')
    vcf_file.close()
    return 1




def read_vcf_file(vcf_file_address):    

    header_lines=[]
    var_lines=[]
    var_pos_list=[]

    vcf_file = open(vcf_file_address,'r')

    for line in vcf_file:
        
        
        line_strip = line.strip() 
              
        if line_strip.startswith('#'):
            header_lines.append(line_strip)
        else:
            var_lines.append(line_strip)
            
            line_parts=line_strip.split('\t') 
            var_pos_list.append(int(line_parts[1]))
            
    vcf_file.close()
    
    return header_lines, var_lines, var_pos_list 
            
    
def read_haplo_sdhap(haplos_file_address):
    
    """ reading the output of sdhap  with genomic positon (output of convertallele sdhap )
    
    input name of file 
    output list
    """
    haplos_file = open(haplos_file_address,'r'); 
    phasings={}
    for line in haplos_file:
        line_strip=line.strip() 
        if line_strip.startswith('B'): # first line of each haplotype block 
            # new block started            
            first_var_block=True
        else:

            line_split=line_strip.split('\t')
            genomic_position=int(line_split[1])

            if first_var_block:
                block_id=genomic_position
                
            first_var_block=False
            
            alleles=line_split[2:]
            phase= "|".join(alleles)
            
            phasings[genomic_position]=[phase,block_id]
        
    return phasings
    
    
def convert_vcf(var_lines, var_pos, phasings):
    var_lines_phasing=[]
    for var_i, var_pos in enumerate(var_pos_list):

        var_line=var_lines[var_i]

        if var_pos in phasings:

            phasing=phasings[var_pos]
            phase_value=phasing[0]+':'+str(phasing[1])

            var_line_split=var_line.split("\t")
            var_line_new_list=var_line_split[:8]+["GT|PS"]+[phase_value]
            var_line_new="\t".join(var_line_new_list)
        else:

            var_line_new=var_line

        var_lines_phasing.append(var_line_new)

    return var_lines_phasing

    
    
def writ_vcf_sdhap(vcf_file_out,header_lines,var_lines_phasing):
    
    vcf_file = open(vcf_file_out,'w'); 
    for header_line in header_lines:
        vcf_file.write(header_line+"\n")
    for var_line_phasing in var_lines_phasing:
        vcf_file.write(var_line_phasing+"\n")
    vcf_file.close()
    
    return 1
    
    

if __name__ == "__main__":

    
    

    
    haplos_file_address= argv[1]  # input name 'sim_varianthaplos.txt'#

    argv2= argv[2]
    
    input_type="sdhap_genomic" # default
    
    
    if argv2=="haplogen":
        input_type = "haplogen"
        vcf_file_out=haplos_file_address+".vcf"  # output name 'out2.vcf'#
    else:
        vcf_file_address=argv2
            
    if input_type == "haplogen":
        var_list=read_haplo_gen(haplos_file_address)
        out_unphased=True # if you need unphased,
        result=writ_vcf_gen(vcf_file_out,out_unphased)
        
        print("A vcf file \""+vcf_file_out+"\" is generated based on the provided haplotype.")
        
    elif input_type=="sdhap_genomic": # convertic 
        
        
        header_lines, var_lines, var_pos_list = read_vcf_file(vcf_file_address)
        
        phasings=read_haplo_sdhap(haplos_file_address)
        
        var_lines_phasing=convert_vcf(var_lines, var_pos, phasings)

        vcf_file_out= haplos_file_address+".vcf"
        
        result=writ_vcf_sdhap(vcf_file_out,header_lines,var_lines_phasing)
 
        print("The new vcf \""+vcf_file_out+"\" contains phase information.")
    
    

