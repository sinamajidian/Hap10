#!/bin/bash

cc_extr='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly6'
frag_weak='/mnt/LTR_userdata/majid001/software/code_py/frag_weak/frag_v3c2.py'
bamprocess='/mnt/LTR_userdata/majid001/software/bamprocess_chaning/bamprocess4e.py'
fragpoly='/mnt/LTR_userdata/majid001/software/code_py/FragmentPoly.py' # python2
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'
split_molec='/mnt/LTR_userdata/majid001/software/code_py/split_v3b.py'



# input 

#frag_path='../frb/frag.txt' # $1  
#vcf_path= './frb/var_het.vcf' # $2 
#mol_len=50   #$3 




#mkdir pre; cd pre

#cp  frag_path .  &&
#grep -v "#"  $vcf_path | cut -f 2  > pos_freebayes.txt


## clustering barcode-specific fragment to molecule-specific fragment
#python3 $split_molec frag.txt pos_freebayes.txt ${mol_len}   &&


## extracting  connected componnets (cc) of (part of sdhap code )
## output of this step is the line number of each cc
#python2 $fragpoly -f frag_sp.txt  -o frag_sd.txt -x SDhaP  &&
#$cc_extr frag_sd.txt        &&  

#csplit --quiet -z connected_dic.txt -f cc -n 1 /block/ {*}  


csplit --quiet -z ../connected_dic_small.txt  -f cc -n 1 /block/ {*}   #connected_dic_good.txt



cc_num=$(ls | grep "cc" | wc -l)

for (( i=0; i<$cc_num; i++)); do # each cc 

    # For whole genome case, extract good cc 
    
    # working
    # if  cat cc${i} | wc -l    ??   gt  100
    
    mkdir ${i}; cd ${i} &&
    mkdir other
    # extracting fragment of cc based on the line number found beforehand
    mv ../cc${i} . &&
    grep -v "block" cc${i} >  cc${i}_pure &&
    awk 'NR == FNR{a[$0]; next};FNR in a' cc${i}_pure ../frag_sp.txt > frag${i}.txt 

    #extracting strongly connected components 
    
    #  note that the ouput may not connected
    # for sdhap it's fine. for matlab code you need to uses another version of this code. 
    python2 $fragpoly -f frag${i}.txt -o frag${i}_dic.txt -x HapTree &&    
    python3 $frag_weak frag${i}_dic.txt  &&
    part_num_i=$(ls | grep "part" | wc -l)

    for ((j=1; j<=$part_num_i; j++)); do    
        awk 'NR == FNR{a[$0]; next};FNR in a' part${j}.txt frag${i}.txt > frag${i}_${j}.txt  &&
        python2 $fragpoly -f frag${i}_${j}.txt -o frag${i}_${j}_sd.txt -x SDhaP &&
        mv frag${i}_${j}.txt other
    done
    
    mv cc* other
    mv part* other
    mv frag_weight.npz other
    mv report_.txt other
    #mv  notcc* other
    
    cd ..

done

wait
pwd
