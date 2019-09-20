#!/bin/bash


# a pralell version of hapcompare

hapcompare='/mnt/LTR_userdata/majid001/software/code_py/compare/hapcompare_v1.py' # python2

#input name is  haplotype.hap

haplotype_file=$1

mkdir parts;
cp $haplotype_file parts/input_hap   &&
cd parts 
csplit -z input_hap -f hap -n 1 /Block/ {*}  &&
p_num=$(ls | grep "hap" | wc -l)  &&

#for (( i=0; i<$p_num; i++)); do   
for (( i=0; i<$(($p_num-1)); i++)); do

    len_i=$(cat hap${i} | wc -l )
    if [ $len_i -gt 2 ]; then 
    python2 $hapcompare  ../../haplogen/hap.variants_het.txt hap${i} -t -v > res${i}.txt  &
    fi
    
    
done
wait

echo -n '' > ../results_rr.txt
echo -n ' ' > ../results_ver.txt

p_num=$(ls | grep "^hap" | wc -l) # start with hap (not including input_hap)
#for ((j=0; j<$p_num; j++)); do
for (( j=0; j<$p_num; j++)); do
    grep "Allelic correlation"  res${j}.txt | head -n 1 | sed 's/Allelic correlation between the corresponding blocks (including variants with discordant dosages):/    /g'  >> ../results_rr.txt
    grep "Vector Error"  res${j}.txt | head -n 1 | sed 's/Vector Error Rates for common variants of each block pair: /    /g'  >> ../results_ver.txt

done


cd ..
cat results_rr.txt| cut -f 2  | sed 's/    /\t/g' | cut -f 2 |sed 's/ ,/\n/g' | tr '\n' ,  > results_rr_pure.txt
cat results_ver.txt| cut -f 2  | sed 's/    /\t/g' | cut -f 2 |sed 's/ ,/\n/g' | tr '\n' ,  > results_ver_pure.txt


echo -n '' > length
grep "lock" $haplotype_file | cut -f 2 | sed 's/ Length of haplotype block / /g' | tr '\n' , >> length


#grep -n "lock" $haplotype_file | cut -f1 -d: | tr '\n' ,  >line_number


#cat $haplotype_file | wc -l >>line_number


wait
pwd

