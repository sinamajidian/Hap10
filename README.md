Hap10
======




The goal is to reconstruct accurate and long haplotype from polyploid genome using linked reads (10x genomics reads).


# Workflow of Hap10

Step 0. Preparation procedure

Step 1. Extracting haplotype information

Step 2. Extracting molecule-specific fragments

Step 3. Extracting strongly connected components of fragments

Step 4. Haplotyping



# Usage

## Step 0. Preparation procedure

Consider that you have a reference genome (`ref/ref.fa`) and two FASTQ files (`reads/R1.fastq.gz` and `reads/R2.fastq.gz` corresponding to Illumina paired-end read). As mentioned in the paper you need to align them to the reference genome using [Longrander]((https://support.10xgenomics.com/genome-exome/software/pipelines/latest/installation)). Then, variants can be called using [freebayes](https://github.com/ekg/freebayes).

```
longranger mkref ref/ref.fasta
mkdir longranger_output
longranger align --id=longranger_output/ --fastqs=reads/ --reference=refdata-ref/
freebayes -f ref/ref.fasta -p $k longranger_output/outs/possorted_bam.bam  > var.vcf
```

in which k is the ploidy level.

## Step 1. Extracting haplotype information

Here, a fragment file (as a text file) is generated using BAM and filtered VCF file. Firstly, you need to build extract_poly, an edited version (under test) of [extracthairs](https://github.com/vibansal/HapCUT2) in which polyploid genomes are also allowed. The makefile will attempt to build samtools and htslib as git submodules. The output of this step is a binary file `extractHAIRS` in folder `build`.

```
git clone https://github.com/smajidian/extract_poly
cd extract_poly
make
```

The input of this step is two files:

- BAM file for an individual containing reads aligned to a reference genome
- VCF file containing only heterzygous SNVs. Complex SNVs and indels are not allowed for this version. The VCF file should only contain genotyped SNPs.


```
utilities/break_vcf.sh  var.vcf var_break.vcf

cat var_break.vcf  | grep -v "1/1/1" | grep -v "0/0/0" grep -v "/2" > var_het.vcf   #This is for triploid.

./extract_poly/build/extractHAIRS --10X 1 --bam out_longranger/outs/possorted_bam.bam --VCF var_het.vcf --out unlinked_fragment_file
```

Here, each line corresponds to one read. Then we link those read with the same barcode and generate barcode-specific fragment file.
```
grep -v "NULL"  unlinked_fragment_file > unlinked_fragment_file_filtered  # removing those reads without barcode

python3 utilities/LinkFragments_brcd_based.py  unlinked_fragment_file_filtered frag.txt
```




## Step 2.  Extracting molecule-specific fragments


```
mean_10x_molecule_length=50      # kilobase

python3 utilities/splitter.py frag.txt $mean_10x_molecule_length frag_sp.txt
```


## Step 3.  Extracting strongly connected components of fragments

```
python3 utilities/extract_scc.py frag_sp.txt scc ./out
```
If you want to have larger haplotype block you can use `cc` instead of `scc`. There are some parameters in the code, you can change them based on your preference.




## Step 4.  Haplotyping

-Fast mode (Hap++):

 You need to install [sdhap](https://sourceforge.net/projects/sdhap/) using this [instruction](https://github.com/smajidian/sdhapc).

```
# triploid

python2 utilities/FragmentPoly.py -f frag_sp.txt  -o frag_sd.txt -x SDhaP  &&

./sdhap/hap_poly frag_sd.txt  out.hap 3

python2 utilities/ConvertAllelesSDhaP.py -p out.hap -o out_with_genomic_position.hap -v var_het.vcf  
```



Since, we know that all variants are hetro, we may remove homo variants in estimated haplotype.

```
cat out_with_genomic_position.hap | grep -v "0\t0\t0" |grep -v "1\t1\t1" > out_filtered.hap
```


-Accurate mode (Hap10):

You may refer to [this folder](https://github.com/smajidian/Hap10/tree/master/accurate_mode).










## Citation:

[Extracthairs](https://github.com/vibansal/HapCUT2)

[Haplosim](https://github.com/EhsanMotazedi/Haplosim)

[Hap10] Hap10: reconstructing accurate and long polyploid haplotypes using linked reads
