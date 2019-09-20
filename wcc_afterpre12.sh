#!/bin/bash

ConvertAllelesSDhaP='/mnt/LTR_userdata/majid001/software/code_py/ConvertAllelesSDhaP.py' # python2
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'

k=$2 # ploidy level
cd $1  &&
#cd wcc &&

cc_num=$(($(find ./ -maxdepth 1 -type d| wc -l)-1))
#for (( i=(($cc_num-1)); i>=0; i--)); do  
for (( i=0; i<$cc_num; i++)); do  
    cd ${i} &&
    #part_num_i=$(ls | grep "_sd.txt" | wc -l)
    part_num_2cc=$(($(ls | grep "pure" | wc -l)-1)) # number of parts after second round of extracting connected
    part_num=$(ls | grep "part" | wc -l) # number of output of strongly cc
    part_num_i=$(($part_num_2cc + $part_num ))  # the max possible index of fragment files
    (for ((j=1; j<=$part_num_i; j++)); do
        if [ -f frag${i}_${j}_sd.txt ]; then
            $sdhap frag${i}_${j}_sd.txt  raw${i}_${j}.hap $k  
        fi
    done) &
    cd ..
done



wait # after all finished. combining estimated haplotype from each folder

echo 'finished'

echo -n '' > hap.hap
for (( i=0; i<$cc_num; i++)); do  
    cd ${i}
    part_num_2cc=$(($(ls | grep "pure" | wc -l)-1)) # number of parts after second round of extracting connected
    part_num=$(ls | grep "part" | wc -l) # number of output of strongly cc
    part_num_i=$(($part_num_2cc + $part_num ))  # the max possible index of fragment files
    for ((j=1; j<=$part_num_i; j++)); do
        if [ -f raw${i}_${j}.hap ]; then
            cat  raw${i}_${j}.hap >> ../hap.hap
        fi
    done
    cd ..
done

python2 $ConvertAllelesSDhaP -p hap.hap -o haplotype.hap -v ../frb/var_het.vcf 
pwd
