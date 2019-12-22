#!/bin/bash

#10x pipline
#triploid

hap_generator_path='/mnt/LTR_userdata/majid001/software/haplogenerator/haplogenerator.py'


sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'
bamprocess='/mnt/LTR_userdata/majid001/software/bamprocess_chaning/bamprocess4e.py'
fragpoly='/mnt/LTR_userdata/majid001/software/code_py/FragmentPoly.py' # python2
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'
split_molec='/mnt/LTR_userdata/majid001/software/code_py/split_v3b.py'
ConvertAllelesSDhaP='/mnt/LTR_userdata/majid001/software/code_py/ConvertAllelesSDhaP.py' # python2
hapcompare='/mnt/LTR_userdata/majid001/software/code_py/compare/hapcompare_v1.py' # python2
lrsim='/home/majid001/Install/LRSIM/simulateLinkedReads.pl' # simulateLinkedReads_notRandom-m.pl'



genomeGenerating()
{    
    mkdir ref
    # soltub_200k.fasta  soltub_1M_2_up.fasta solanum.fa solanum_5m_R1.fa solanum_5m_a.fa solanum_5m_noN.fa
    cp /mnt/LTR_userdata/majid001/genomes/solanum_tuberosum/solanum_1me.fa ref/ref.fasta
    #head -n 1000 ref/ref1.fasta > ref/ref.fasta
    samtools faidx ref/ref.fasta
    longranger mkref  ref/ref.fasta
    
    #printf '=%.0s' {1..100}; printf "\n"; printf "haplogenerator started\n"; printf '=%.0s' {1..100}; printf "\n"
    mkdir haplogen
    python2 $hap_generator_path -f ref/ref.fasta -o haplogen/sim -m "{'A':'C','C':'A','G':'T','T':'G'}"  -p 3 -v --model poisson -s [.01,0,0] 
    #--dosage "[0.5,0.23,0.27]"
    cd haplogen
    #grep -P "\tN\t" sim_varianthaplos.txt | wc -l;
    #grep -P "\t2\t" sim_varianthaplos.txt  | wc -l
    grep -vP "1\t1\t1" sim_varianthaplos.txt | grep -vP "0\t0\t0" > hap.variants_het.txt # for triploid !!
    cut -f 3  hap.variants_het.txt | grep -v "contig" > pos_truth.txt
    mkdir edited; cp *.fa edited; cd edited
    for idx in {1..3}; do sed -i "s/^>/>hap${idx}_/g" sim_hap$idx.fa; done # there is just one part in ref fasta
    for idx in {1..3}; do cat sim_hap${idx}.fa >> genome.fa; done
    cd ../..
}
    
ReadGenerating()
{
    printf '=%.0s' {1..100}; printf "\n"; printf "LRSIM started\n"; printf '=%.0s' {1..100}; printf "\n"
    mkdir out_lrsim
    #simulateLinkedReads_notRandom-m.pl  # number of molecules is not random    -f change , also change the spliter
    
    
    $lrsim -g haplogen/edited/genome.fa -p out_lrsim/sim \
     -o -x 0.193 -f 50 -m 10 -t 0.54 # for coverage 60 triploid 1mbp # it is sensitive to space ! do not put excessive space in the middle
     #coverage 50.3 for 4.8mb     -o -x .78 -f 50 -m 10 -t 2.2 #-x 9.3 -f 20 -m 10 -t 66.5 #for 88.6mb
    
    
    
    cd out_lrsim; mkdir edited
    cp sim_S1_L00* edited/
    cd edited
    gunzip *
    sed 's/\/1/ /g' sim_S1_L001_R1_001.fastq > sim_ed_S1_L001_R1_001.fastq;
    sed 's/\/2/ /g' sim_S1_L001_R2_001.fastq > sim_ed_S1_L001_R2_001.fastq;
    gzip -k sim_ed_S1_L001_R1_001.fastq;
    gzip -k sim_ed_S1_L001_R2_001.fastq
    mkdir filegz
    mv *.gz filegz
    cd ../../
    printf '=%.0s' {1..100}; printf "\n"; printf "LRSIM finished\n"; printf '=%.0s' {1..100}; printf "\n"    

}

Aligning()
{
    printf '=%.0s' {1..100}; printf "\n"; printf "LongRanger started\n"; printf '=%.0s' {1..100}; printf "\n"    
    cd ./$2
    longranger align --id=out_longranger --fastqs=out_lrsim/edited/filegz/ --reference=refdata-ref
    printf '=%.0s' {1..100}; printf "\n"; printf "LongRanger finished\n"; printf '=%.0s' {1..100}; printf "\n"    

}

VariantCalling()
{
    printf '=%.0s' {1..100}; printf "\n"; printf "Variant calling  started\n"; printf '=%.0s' {1..100}; printf "\n"
    mkdir frb
    cd frb
    cp ../out_longranger/outs/possorted_bam.bam .
    samtools index possorted_bam.bam
    cp ../haplogen/pos_truth.txt .
    freebayes -f ../ref/ref.fasta -p 3  possorted_bam.bam  > var.vcf
    /mnt/LTR_userdata/majid001/software/pipelines/break_vcf.sh var.vcf var_break.vcf
    cat var_break.vcf  | grep -v "1/1/1" | grep -v "0/0/0"  > var_het.vcf
    grep -v "#" var_het.vcf | cut -f 2  > pos_freebayes.txt
    cd ..
}


HaplotypingSdhap()
{
    cd frb
    python $bamprocess possorted_bam.bam var_het.vcf
    mv HapCUT_frags_out_longranger.matrix frag.txt
    #time(
    #python3 $split_molec frag.txt pos_freebayes.txt 50 
    #cd ..;
    #mkdir sdhap;  cp frb/frag_sp.txt sdhap;  cd sdhap;
    #python2 $fragpoly -f frag_sp.txt  -o frag_sd.txt -x SDhaP
    #$sdhap frag_sd.txt  out_sd_raw.hap  3)

    #python2 $ConvertAllelesSDhaP -p out_sd_raw.hap -o out_sdhap.hap -v ../frb/var_het.vcf
    #cut -f 2 out_sdhap.hap > pos_sdhap.txt 
    #sed -i '1d' pos_sdhap.txt 
    #cd ../..
    #python2 $hapcompare haplogen/hap.variants_het.txt sdhap/out_sdhap.hap   -t -v   > dhap/res_sdhap.txt
}





######  main
#usage:  ./code.sh 31

mkdir $1
cd $1

#mkdir 65c_m10b
#cd  65c_m10b
#cp -r ../65c_1/ref .
#cp -r ../65c_1/haplogen .

#longranger mkref  ref/ref.fasta


genomeGenerating
ReadGenerating
Aligning
VariantCalling

HaplotypingSdhap  





pwd


#python3 /mnt/scratch/majid001/software/code_py/frag_weak/frag_v3c2.py  

#HaplotypingSdhapSplited



