#!/bin/bash

cc_extr='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly4'
frag_weak='/mnt/LTR_userdata/majid001/software/code_py/frag_weak/frag_v3c2.py'
bamprocess='/mnt/LTR_userdata/majid001/software/bamprocess_chaning/bamprocess4e.py'
fragpoly='/mnt/LTR_userdata/majid001/software/code_py/FragmentPoly.py' # python2
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'
split_molec='/mnt/LTR_userdata/majid001/software/code_py/split_v3b.py'


cd $1  &&
mkdir pre12; cd pre12 &&

mol_len=$2


cp ../frb/frag.txt . &&
cp ../frb/pos_freebayes.txt . &&
python3 $split_molec frag.txt pos_freebayes.txt ${mol_len}   &&
python2 $fragpoly -f frag_sp.txt  -o frag_sd.txt -x SDhaP  &&

$cc_extr frag_sd.txt        &&  # extracting connected componnets (part of sdhap code )
csplit --quiet -z connected_dic.txt -f cc -n 1 /block/ {*}  &&
cc_num=$(ls | grep "cc" | wc -l)

#for (( i=(($cc_num-1)); i>=0; i--)); do  
for (( i=0; i<$cc_num; i++)); do
    mkdir ${i}; cd ${i} &&
    mv ../cc${i} . &&
    grep -v "block" cc${i} >  cc${i}_pure &&
    awk 'NR == FNR{a[$0]; next};FNR in a' cc${i}_pure ../frag_sp.txt > frag${i}.txt 
    python2 $fragpoly -f frag${i}.txt -o frag${i}_dic.txt -x HapTree &&
    
    #extracting weak parts  # the ouput may not connected
    python3 $frag_weak frag${i}_dic.txt  &&
    part_num_i=$(ls | grep "part" | wc -l)
    counter_newcc=$part_num_i
    for ((j=1; j<=$part_num_i; j++)); do
    
        awk 'NR == FNR{a[$0]; next};FNR in a' part${j}.txt frag${i}.txt > frag${i}_${j}.txt  &&
        python2 $fragpoly -f frag${i}_${j}.txt -o frag${i}_${j}_sd.txt -x SDhaP &&
        $cc_extr frag${i}_${j}_sd.txt  &&
        mv connected_dic.txt  connected_dic${i}_${j}.txt &&
        cc_ij=$(grep "block" connected_dic${i}_${j}.txt  | wc -l)
        echo $cc_ij
        if [ $cc_ij -gt 1 ]; then  # if a fragment (output of frag_weak.py) is not connected, the file is renamed, and split to some. the index starts from the  
            mv frag${i}_${j}.txt notcc${i}_${j}.txt &&
            rm frag${i}_${j}_sd.txt
            csplit -z connected_dic${i}_${j}.txt  -f cc${i}_${j}_ -n 1 /block/ {*}  &&
            cc_num=$(ls | grep "cc${i}_${j}_" | wc -l)
            for (( ii=$(($cc_num-1)); ii>=0; ii--)); do  
                counter_newcc=$(($counter_newcc+1)) &&
                #grep -v "block" cc${i}_${j}_${ii} >  cc${i}_${j}_${ii}_pure &&
                #awk 'NR == FNR{a[$0]; next};FNR in a' cc${i}_${j}_${ii}_pure notcc${i}_${j}.txt > frag${i}_${counter_newcc}.txt
                grep -v "block" cc${i}_${j}_${ii} | awk 'NR == FNR{a[$0]; next};FNR in a' notcc${i}_${j}.txt > frag${i}_${counter_newcc}.txt
                python2 $fragpoly -f frag${i}_${counter_newcc}.txt -o frag${i}_${counter_newcc}_sd.txt -x SDhaP
                 
            done   
        fi
        
    done
    cd ..
done

wait
pwd
