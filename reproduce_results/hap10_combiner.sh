#!/bin/bash


ConvertAllelesSDhaP='/mnt/LTR_userdata/majid001/software/code_py/ConvertAllelesSDhaP.py' # python2

cc_num=$(($(find ./ -maxdepth 1 -type d| wc -l)-1))

echo -n '' > raw.hap
for (( i=0; i<$cc_num; i++)); do  
    #ls
    cd ${i}
    #part_num_i=$(ls | grep ".hap" | grep -v "raw" | wc -l)
    #for ((j=1; j<=$part_num_i; j++)); do
    part_num_2cc=$(($(ls | grep "pure" | wc -l)-1)) # number of parts after second round of extracting connected
    part_num=$(ls | grep "part" | wc -l) # number of output of strongly cc
    part_num_i=$(($part_num_2cc + $part_num ))  # the max possible index of fragment files
    for ((j=1; j<=$part_num_i; j++)); do
        #if [ -f frag${i}_${j}.hap ]; then
            cat  frag${i}_${j}.hap | sed  "s/Total MEC/frag${i}_${j}.hapTotal MEC/g" >> ../raw.hap
            

            #cat  raw${i}_${j}.hap >> ../raw.hap # for wcc
        #fi
    done
    cd ..
done

python2 $ConvertAllelesSDhaP -p raw.hap -o haplotype.hap -v ../frb/var_het.vcf 
pwd
