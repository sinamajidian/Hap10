Hap10
======


The goal of Hap10 is to reconstruct accurate and long haplotypes polyploid genome using linked reads.

Here, we suppose that you prepared the aligned reads (as BAM files per chrs) and the called variants (as a VCF file). Additionally, you have installed [extract_poly](https://github.com/smajidian/extract_poly) and [SDhaP](https://github.com/smajidian/sdhapc) beforehand. For a complete description see the [Wiki page](https://github.com/smajidian/Hap10/wiki/Hap10-Wiki-page). We provide a bash script summarising the steps for a test case in [test_data folder](https://github.com/smajidian/Hap10/tree/master/test_data).


Sina Majidian, Mohammad Hossein Kahaei, and Dick de Ridder. ["Hap10: reconstructing accurate and long polyploid haplotypes using linked reads."](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03584-5)


## Usage


Firstly, you need to clone the Hap10 package.

```
git clone https://github.com/smajidian/Hap10.git
cd Hap10
```

You can run the following steps using a single bash script `utilities/hap.sh`.


* important* 
Please run the pipeline for each chromosome separately.

As the first step, we generate the fragment file using BAM and VCF files.
```
./extract_poly/build/extractHAIRS --10X 1 --bam out_longranger/outs/possorted_bam.bam --VCF var_het.vcf --out unlinked_fragments
python3 utilities/LinkFragments_brcd_based.py  unlinked_fragments frag.txt
```

Then, we extract the molecule-specific fragments in which `m` is the mean 10X molecule length (in Kb) which can be set as 50.

```
python3 utilities/splitter.py frag.txt var_het.vcf $m 
```

Now, extracting strongly connected components of fragments using the following.

```
python3 utilities/extract_scc.py frag_sp.txt scc ./out
```

We are ready to use the assembly core:

```
python2 utilities/FragmentPoly.py -f frag_sp.txt  -o frag_sd.txt -x SDhaP
./sdhap/hap_poly frag_sd.txt  out.hap 3
python2 utilities/ConvertAllelesSDhaP.py -p out.hap -o out_with_genomic_position.hap -v var_het.vcf  
python utilities/hap2vcf.py haplotype_with_genomic_position.hap $vcf
```
The output of the provided bash scripts is a VCF file in which the haplotypes are provided in GT field.

The above mentioned is the pipeline for Hap++. For accurate mode (Hap10), you may refer to this [page](https://github.com/smajidian/Hap10/tree/master/accurate_mode).





## Copyright

This package is distributed under the Creative Commons Attribution-ShareAlike 4.0 International Public License.

